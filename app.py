#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bundle CFD app : the local server behind mesh_explorer.html
============================================================
``mesh_explorer.html`` on its own is a mesh tool: open it as a file and the
Geometry tab works with no server at all.  Running it through this server adds
the other three tabs, because a browser page cannot launch Fluent.

    python3 app.py                 # serve on 127.0.0.1:8737, open a browser
    python3 app.py --port 9000     # somewhere else
    python3 app.py --backend mock  # force the offline mock, no Fluent needed
    python3 app.py --no-open       # do not touch the browser

Stdlib only, plus ``ansys-fluent-core`` when the real backend is used.  It
binds to the loopback address and nothing else: this drives a local solver and
writes local files, so it has no business being reachable from the network.

One job at a time.  That is not a limitation worth engineering around for a
tool driving one Fluent session on one workstation, and it keeps the state the
UI has to show down to something a person can hold in their head.

The mesh is built HERE, by mesh_explorer.py, from the same parameters the
Geometry tab is showing.  The browser never uploads a mesh - they run to
hundreds of megabytes - and the two implementations are checked against each
other node for node, so building it on this side changes nothing.
"""

from __future__ import annotations

import argparse
import json
import os
import socket
import sys
import threading
import time
import traceback
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import fluent_case as FC            # noqa: E402
import mesh_explorer as ME          # noqa: E402

PAGE = os.path.join(HERE, "mesh_explorer.html")

#  browser parameter name -> Case keyword.  The Geometry tab speaks JavaScript
#  camelCase; the mesh generator speaks Python.  One table, both directions.
PARAM_MAP = {
    "D": "D", "ST": "ST", "SL": "SL", "nRows": "n_rows", "nCols": "n_cols",
    "lIn": "l_in", "lOut": "l_out", "H": "H", "helix": "helix", "rMid": "r_mid",
    "sector": "sector", "nAz": "n_az", "nRad": "n_rad", "firstLayer": "first_layer",
    "nZ": "n_z", "nxIn": "n_x_in", "nxOut": "n_x_out", "halfRods": "half_rods",
    "firstOffset": "first_offset", "wrap": "wrap", "unit": "unit",
}
INT_PARAMS = {"n_rows", "n_cols", "n_az", "n_rad", "n_z", "n_x_in", "n_x_out"}
BOOL_PARAMS = {"half_rods", "first_offset", "wrap"}


def case_from(geometry, params):
    """Build a mesh_explorer.Case from what the Geometry tab is showing."""
    kw = {}
    for js, py in PARAM_MAP.items():
        if js not in (params or {}):
            continue
        v = params[js]
        if py in BOOL_PARAMS:
            kw[py] = bool(v)
        elif py == "unit":
            kw[py] = str(v)
        elif py in INT_PARAMS:
            kw[py] = int(round(float(v)))
        else:
            kw[py] = float(v)
    return ME.Case(geometry, **kw)


# =============================================================================
#  THE JOB
# =============================================================================
class Job(object):
    """One end-to-end run: mesh -> launch -> set up -> initialise -> iterate.

    Every stage reports itself, so a run that dies has a stage attached to it
    rather than just a stack trace.
    """

    STAGES = ["queued", "meshing", "launching", "setup", "initialising",
              "iterating", "finished"]

    def __init__(self, jid, geometry, params, settings, backend, out_dir):
        self.id = jid
        self.geometry = geometry
        self.params = params
        self.settings = settings
        self.backend = backend
        self.out_dir = out_dir
        self.stage = "queued"
        self.error = None
        self.log_lines = []
        self.started = time.time()
        self.finished_at = None
        self.mesh_stats = None
        self.mesh_path = None
        self.driver = None
        self.case = None
        self.lock = threading.Lock()
        self.thread = None

    # -- reporting --------------------------------------------------------
    def log(self, msg):
        with self.lock:
            self.log_lines.append("%7.1fs  %s" % (time.time() - self.started, msg))
            if len(self.log_lines) > 4000:
                del self.log_lines[:1000]

    def status(self, log_from=0):
        with self.lock:
            res = list(self.driver.residuals) if self.driver else []
            lines = self.log_lines[log_from:]
            n_lines = len(self.log_lines)
        return {
            "id": self.id, "stage": self.stage, "error": self.error,
            "geometry": self.geometry, "backend": self.backend,
            "mock": bool(self.driver and self.driver.mock),
            "elapsed": round(time.time() - self.started, 2),
            "finished": self.finished_at is not None,
            "mesh": self.mesh_stats, "mesh_path": self.mesh_path,
            "residuals": res, "log": lines, "log_next": n_lines,
            "surfaces": self.driver.surfaces() if self._can_report() else [],
        }

    def _can_report(self):
        return self.driver is not None and self.stage in ("iterating", "finished")

    # -- the run ----------------------------------------------------------
    def start(self):
        self.thread = threading.Thread(target=self._run, name="job-" + self.id)
        self.thread.daemon = True
        self.thread.start()

    def _run(self):
        try:
            self.stage = "meshing"
            self.case = case_from(self.geometry, self.params)
            self.case.validate()
            self.log("geometry %s, %s" % (self.geometry, self.case.name()))
            mesh = ME.Mesh(self.case)
            mesh.build()
            rep = ME.run_checks(mesh, quiet=True)
            self.mesh_stats = {k: rep[k] for k in
                               ("cells", "nodes2d", "quads", "points", "cracks",
                                "maxdev", "vol_err", "dev_max", "ar_max", "patches")}
            self.log("mesh: %d cells, %d cracks, volume error %.1e"
                     % (rep["cells"], rep["cracks"], rep["vol_err"]))
            name = self.case.name(mesh.az_full[1])
            self.mesh_path = os.path.join(self.out_dir, name + ".msh")
            if self.backend == "mock":
                self.log("MOCK backend: skipping the .msh write")
            else:
                ME.write_fluent(mesh, self.mesh_path)
                self.log("wrote %s (%.1f MB)"
                         % (os.path.basename(self.mesh_path), ME.mb(self.mesh_path)))

            problems = FC.validate(self.settings)
            for p in problems:
                self.log("settings: " + p)

            self.driver = FC.make_driver(self.backend, self.case, self.settings,
                                         self.mesh_path, self.log)
            self.stage = "launching"
            self.driver.launch()
            self.stage = "setup"
            self.driver.setup()
            self.stage = "initialising"
            self.driver.initialize()
            self.stage = "iterating"
            self.driver.iterate(int(self.settings["run"]["iterations"]))
            self.stage = "finished"
            self.log("done")
        except Exception as exc:                        # noqa: BLE001
            self.error = "%s: %s" % (type(exc).__name__, exc)
            self.log("FAILED in stage %r - %s" % (self.stage, self.error))
            for line in traceback.format_exc().splitlines()[-12:]:
                self.log("    " + line)
            self.stage = "finished"
        finally:
            self.finished_at = time.time()

    def stop(self):
        if self.driver is not None:
            self.log("stop requested")
            try:
                self.driver.interrupt()
            except Exception as exc:                    # noqa: BLE001
                self.log("interrupt failed: %s" % exc)

    def close(self):
        if self.driver is not None:
            self.driver.close()


# =============================================================================
#  THE SERVER
# =============================================================================
class App(object):
    def __init__(self, backend, out_dir):
        self.backend = backend
        self.out_dir = out_dir
        self.job = None
        self.counter = 0
        self.lock = threading.Lock()

    def info(self):
        return {
            "ok": True,
            "fluent_available": FC.fluent_available(),
            "backend": self.backend,
            "out_dir": self.out_dir,
            "settings": FC.SETTINGS,
            "defaults": FC.default_settings(),
            "variables": FC.VARIABLES,
            "patches": FC.PATCHES,
            "geometries": sorted(ME.GEOMETRIES),
            "reports": ["area-weighted-avg", "facet-min", "facet-max",
                        "mass-flow-rate", "area"],
        }

    def start_job(self, body):
        with self.lock:
            if self.job is not None and not self.job.finished_at:
                raise ValueError("a run is already going; stop it first")
            if self.job is not None:
                self.job.close()
            self.counter += 1
            settings = FC.merge_settings(body.get("settings"))
            backend = body.get("backend") or self.backend
            self.job = Job("j%d" % self.counter, body.get("geometry", "rod-inline"),
                           body.get("params") or {}, settings, backend, self.out_dir)
            self.job.start()
            return self.job.status()

    def need_job(self):
        if self.job is None:
            raise ValueError("nothing has been run yet")
        return self.job

    def need_results(self):
        job = self.need_job()
        if job.driver is None or job.stage not in ("iterating", "finished"):
            raise ValueError("no results yet - the run is at stage %r" % job.stage)
        if job.error:
            raise ValueError("the run failed: " + job.error)
        return job


class Handler(BaseHTTPRequestHandler):
    server_version = "BundleCFD/1.0"
    app = None                       # set on the server instance

    def log_message(self, fmt, *args):       # quieter than the default
        if "--verbose" in sys.argv:
            sys.stderr.write("  %s\n" % (fmt % args))

    # -- plumbing ---------------------------------------------------------
    def _send(self, code, body, ctype="application/json; charset=utf-8"):
        if not isinstance(body, (bytes, bytearray)):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            pass                                        # the tab went away

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj, ensure_ascii=False, allow_nan=False,
                                    default=float))

    def _fail(self, exc, code=400):
        self._json({"ok": False, "error": "%s: %s" % (type(exc).__name__, exc)}, code)

    def _body(self):
        n = int(self.headers.get("Content-Length") or 0)
        if not n:
            return {}
        return json.loads(self.rfile.read(n).decode("utf-8"))

    # -- routes -----------------------------------------------------------
    def do_GET(self):                                   # noqa: N802
        path = self.path.split("?")[0]
        try:
            if path in ("/", "/index.html", "/mesh_explorer.html"):
                with open(PAGE, "rb") as f:
                    return self._send(200, f.read(), "text/html; charset=utf-8")
            if path == "/api/info":
                return self._json(self.app.info())
            if path == "/api/job":
                q = self.path.split("?", 1)
                frm = 0
                if len(q) > 1:
                    for kv in q[1].split("&"):
                        if kv.startswith("log_from="):
                            frm = int(kv.split("=")[1] or 0)
                job = self.app.job
                if job is None:
                    return self._json({"ok": True, "job": None})
                return self._json({"ok": True, "job": job.status(frm)})
            return self._send(404, "not found", "text/plain; charset=utf-8")
        except Exception as exc:                        # noqa: BLE001
            return self._fail(exc, 500)

    def do_POST(self):                                  # noqa: N802
        path = self.path.split("?")[0]
        try:
            body = self._body()
            if path == "/api/run":
                return self._json({"ok": True, "job": self.app.start_job(body)})
            if path == "/api/stop":
                self.app.need_job().stop()
                return self._json({"ok": True})
            if path == "/api/mesh":
                return self._json({"ok": True, "mesh": self._mesh_only(body)})
            if path == "/api/field":
                job = self.app.need_results()
                f = job.driver.field(body["surface"], body["variable"])
                return self._json({"ok": True, "field": f})
            if path == "/api/report":
                job = self.app.need_results()
                v = job.driver.report(body.get("kind", "area-weighted-avg"),
                                      body.get("surfaces") or ["inlet"],
                                      body.get("variable", "pressure"))
                return self._json({"ok": True, "value": v,
                                   "mock": bool(job.driver.mock)})
            if path == "/api/journal":
                settings = FC.merge_settings(body.get("settings"))
                case = case_from(body.get("geometry", "rod-inline"),
                                 body.get("params") or {})
                #  name it exactly as the real run would, azimuthal count and all
                mesh_path = os.path.join(
                    self.app.out_dir, case.name(ME.quick_az(case)[1]) + ".msh")
                return self._json({"ok": True,
                                   "script": FC.journal(body.get("geometry"),
                                                        body.get("params"), settings,
                                                        mesh_path)})
            return self._send(404, "not found", "text/plain; charset=utf-8")
        except Exception as exc:                        # noqa: BLE001
            return self._fail(exc)

    def _mesh_only(self, body):
        """Write the .msh without running anything, for people who only want
        the mesh on this machine."""
        case = case_from(body.get("geometry", "rod-inline"), body.get("params") or {})
        case.validate()
        mesh = ME.Mesh(case)
        mesh.build()
        rep = ME.run_checks(mesh, quiet=True)
        name = case.name(mesh.az_full[1])
        path = os.path.join(self.app.out_dir, name + ".msh")
        ME.write_fluent(mesh, path)
        return {"path": path, "name": name, "mb": round(ME.mb(path), 2),
                "cells": rep["cells"], "cracks": rep["cracks"],
                "vol_err": rep["vol_err"]}


def free_port(preferred):
    s = socket.socket()
    try:
        s.bind(("127.0.0.1", preferred))
        s.close()
        return preferred
    except OSError:
        s.close()
        s = socket.socket()
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
        s.close()
        return port


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8737)
    ap.add_argument("--backend", choices=["auto", "fluent", "mock"], default="auto",
                    help="auto uses Fluent when ansys-fluent-core is importable")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "runs"),
                    help="where meshes, cases and data are written")
    ap.add_argument("--no-open", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args(argv)

    out_dir = os.path.abspath(a.out_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(PAGE):
        print("[ERROR] %s is missing" % PAGE)
        return 2

    port = free_port(a.port)
    Handler.app = App(a.backend, out_dir)
    httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = "http://127.0.0.1:%d/" % port

    have = FC.fluent_available()
    print("=" * 66)
    print(" bundle CFD app   %s" % url)
    print(" backend        : %s%s" % (a.backend,
                                      "" if have else "   (ansys-fluent-core NOT installed"
                                                      " - the mock will be used)"))
    print(" writing to     : %s" % out_dir)
    if port != a.port:
        print(" note           : port %d was busy" % a.port)
    print(" loopback only; Ctrl-C to stop")
    print("=" * 66)
    if not a.no_open:
        threading.Timer(0.4, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping")
    finally:
        if Handler.app.job is not None:
            Handler.app.job.close()
        httpd.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
