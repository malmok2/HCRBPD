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
import re
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

#  a plane name becomes a Fluent surface name and a key in the UI, so keep it
#  to something both can hold without quoting
PLANE_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_\-]{0,30}$")

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


def _jsonable(x):
    """json.dumps' last resort, for objects it does not know.

    This used to be plain `float`, which is wrong twice over: it silently
    turns anything numeric-looking into a number, and on a numpy ARRAY it
    fails with "only 0-dimensional arrays can be converted to Python scalars"
    - a message about numpy, three layers below the field that actually went
    unconverted.  Convert what can be converted; name what cannot.
    """
    tolist = getattr(x, "tolist", None)         # ndarray and numpy scalars
    if tolist is not None:
        return tolist()
    item = getattr(x, "item", None)
    if item is not None:
        return item()
    raise TypeError("cannot send a %s over the wire: %r"
                    % (type(x).__name__, x))


#  what each extension in the output folder is, so the panel can say rather
#  than leaving the user to recognise them
KINDS = [
    (".fields.json", "fields", "saved field - opens with no Fluent"),
    (".cas.h5", "case", "Fluent case (mesh + set-up)"),
    (".dat.h5", "data", "Fluent data (the solution)"),
    (".msh", "mesh", "Fluent mesh"),
    (".vtu", "vtu", "ParaView"),
    (".stl", "stl", "STL"),
    (".py", "journal", "PyFluent script"),
]


def kind_of(name):
    low = name.lower()
    for ext, key, label in KINDS:
        if low.endswith(ext):
            return key, label
    return "other", ""


def list_outputs(out_dir, limit=200):
    """Everything the app has written, newest first."""
    rows = []
    try:
        names = os.listdir(out_dir)
    except OSError:
        return rows
    for n in names:
        full = os.path.join(out_dir, n)
        try:
            st = os.stat(full)
        except OSError:
            continue
        if not os.path.isfile(full):
            continue
        key, label = kind_of(n)
        rows.append({"name": n, "path": full, "bytes": st.st_size,
                     "mtime": st.st_mtime, "kind": key, "label": label})
    rows.sort(key=lambda r: r["mtime"], reverse=True)
    return rows[:limit]


def open_in_file_manager(path):
    """Show a folder in the desktop's own file manager.

    The server is the only side that can do this, and it only ever opens its
    own output directory - the path is not taken from the request.
    """
    import subprocess
    if sys.platform.startswith("win"):
        os.startfile(path)                      # noqa: S606  (Windows only)
        return "explorer"
    cmd = "open" if sys.platform == "darwin" else "xdg-open"
    subprocess.Popen([cmd, path], stdout=subprocess.DEVNULL,
                     stderr=subprocess.DEVNULL)
    return cmd


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
        self.notes = []                 # settings the server had to correct
        self.load_path = None           # set when reopening a saved case
        self.loaded = False
        self.snapshot = False           # ... and that case was a saved field
        self.planes = {}                # name -> {axis, value}
        self._bbox = None

    # -- reporting --------------------------------------------------------
    def log(self, msg):
        with self.lock:
            self.log_lines.append("%7.1fs  %s" % (time.time() - self.started, msg))
            if len(self.log_lines) > 4000:
                del self.log_lines[:1000]

    def status(self, log_from=0):
        with self.lock:
            res = list(self.driver.residuals) if self.driver else []
            mon = list(getattr(self.driver, "monitors", []) or []) if self.driver else []
            lines = self.log_lines[log_from:]
            n_lines = len(self.log_lines)
        return {
            "id": self.id, "stage": self.stage, "error": self.error,
            "geometry": self.geometry, "backend": self.backend,
            "mock": bool(self.driver and self.driver.mock),
            "elapsed": round(time.time() - self.started, 2),
            "finished": self.finished_at is not None,
            "mesh": self.mesh_stats, "mesh_path": self.mesh_path,
            "residuals": res, "monitors": mon,
            "log": lines, "log_next": n_lines,
            "surfaces": self.surface_list(),
            "planes": dict(self.planes), "loaded": self.loaded,
            "snapshot": bool(self.driver is not None
                             and getattr(self.driver, "snapshot", False)),
            "variables": (self.driver.variables() if self._can_report() else None),
            "out_dir": self.out_dir,
        }

    def _can_report(self):
        return self.driver is not None and self.stage in ("iterating", "finished")

    def surface_list(self):
        """Every surface the Results tab may show.

        The planes are added from this job's own registry rather than taken on
        trust from the driver: they were created here, so losing them because
        a release will not enumerate surfaces would be the app forgetting its
        own work.
        """
        if not self._can_report():
            return []
        names = list(self.driver.surfaces())
        for n in self.planes:
            if n not in names:
                names.append(n)
        return names

    def bbox(self):
        """The domain box the plane sliders span, measured once."""
        if self._bbox is None:
            self._bbox = self.driver.bbox()
        return self._bbox

    # -- the run ----------------------------------------------------------
    def start(self):
        target = self._run_load if self.load_path else self._run
        self.thread = threading.Thread(target=target, name="job-" + self.id)
        self.thread.daemon = True
        self.thread.start()

    def _run(self):
        try:
            #  first line of every log: paste one back and the code that made
            #  it is not in doubt
            self.log("code %s  |  server up since %s" % (version_line(VERSION), STARTED))
            for n in self.notes:
                self.log("settings corrected: " + n)
            self.stage = "meshing"
            self.case = case_from(self.geometry, self.params)
            self.case.validate()
            mesh = ME.Mesh(self.case)
            mesh.build()
            rep = ME.run_checks(mesh, quiet=True)
            self.mesh_stats = {k: rep[k] for k in
                               ("cells", "nodes2d", "quads", "points", "cracks",
                                "maxdev", "vol_err", "dev_max", "ar_max", "patches")}
            self.log("mesh: %d cells, %d cracks, volume error %.1e"
                     % (rep["cells"], rep["cracks"], rep["vol_err"]))
            name = self.case.name(mesh.az_full[1])
            self.log("geometry %s, %s" % (self.geometry, name))
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

    def _run_load(self):
        """Reopen a case+data that was written earlier.

        Same Job, same stages, same everything downstream - it just skips
        meshing and iterating, which is the whole point: a solution that took
        an hour should not have to be produced twice to be looked at twice.
        """
        try:
            self.log("code %s  |  server up since %s" % (version_line(VERSION), STARTED))
            self.log("reopening %s" % os.path.basename(self.load_path))
            self.stage = "launching"
            if self.snapshot:
                #  a saved field needs no solver at all - that is the point
                self.driver = FC.SnapshotDriver(
                    FC.read_snapshot(self.load_path), self.log)
                self.driver.launch()
            else:
                self.driver = FC.make_driver(self.backend, self.case,
                                             self.settings, self.load_path,
                                             self.log)
                self.driver.launch()
                self.stage = "setup"
                self.driver.load_case(self.load_path)
            self.mesh_path = self.load_path
            self.stage = "finished"
            self.loaded = True
            self.log("done - this is a reopened solution, nothing was re-solved")
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
        installed, detail = FC.fluent_installed()
        return {
            "ok": True,
            "version": VERSION,
            "started": STARTED,
            "fluent_available": FC.fluent_available(),
            "fluent_installed": installed,
            "fluent_detail": detail,
            "effective_backend": ("mock" if self.backend == "mock"
                                  else ("fluent" if installed else "mock")),
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
            #  the browser keeps the last case in local storage, so a choice
            #  the schema has since corrected can still arrive here
            notes = []
            settings = FC.merge_settings(body.get("settings"), notes)
            backend = body.get("backend") or self.backend
            self.job = Job("j%d" % self.counter, body.get("geometry", "rod-inline"),
                           body.get("params") or {}, settings, backend, self.out_dir)
            self.job.notes = notes
            self.job.start()
            return self.job.status()

    def need_job(self):
        if self.job is None:
            raise ValueError("nothing has been run yet")
        return self.job

    def resolve_output(self, body, key="file"):
        """A file IN the output folder, named by its basename.

        The request names a file, never a path: nothing outside the folder the
        app writes to is reachable through it.
        """
        name = str(body.get(key) or "")
        if not name or os.path.basename(name) != name:
            raise ValueError("name a file in the output folder, not a path")
        full = os.path.join(self.out_dir, name)
        if not os.path.isfile(full):
            raise ValueError("%s is not in %s" % (name, self.out_dir))
        return name, full

    def start_load(self, body):
        """Open a saved case+data, or a saved field.

        A .cas.h5 needs Fluent; a snapshot does not need anything, which is
        the point of it.
        """
        name, full = self.resolve_output(body)
        with self.lock:
            if self.job is not None and not self.job.finished_at:
                raise ValueError("a run is already going; stop it first")
            if self.job is not None:
                self.job.close()
            self.counter += 1
            settings = FC.merge_settings(body.get("settings"))
            geometry = body.get("geometry", "rod-inline")
            self.job = Job("j%d" % self.counter, geometry,
                           body.get("params") or {}, settings,
                           body.get("backend") or self.backend, self.out_dir)
            #  the mock needs a Case to cut planes through; the real driver
            #  takes its geometry from the file it is about to read
            self.job.case = case_from(geometry, body.get("params") or {})
            self.job.load_path = full
            self.job.snapshot = name.endswith(FC.SNAPSHOT_EXT)
            self.job.start()
            return self.job.status()

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
                                    default=_jsonable))

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
            if path == "/api/files":
                return self._json({"ok": True, "out_dir": self.app.out_dir,
                                   "files": list_outputs(self.app.out_dir)})
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
            if path == "/api/load":
                return self._json({"ok": True, "job": self.app.start_load(body)})
            if path == "/api/snapshot":
                job = self.app.need_results()
                want = [str(x) for x in (body.get("surfaces") or [])]
                live = set(job.surface_list())
                bad = [x for x in want if x not in live]
                if bad:
                    raise ValueError("no such surface: %s" % ", ".join(bad))
                if not want:
                    raise ValueError("pick at least one surface to save")
                vars_ = [str(x) for x in (body.get("variables") or [])]
                known = set(job.driver.variables())
                vars_ = [v for v in vars_ if v in known] or list(job.driver.variables())
                base = os.path.splitext(os.path.splitext(
                    os.path.basename(job.mesh_path or "fields"))[0])[0]
                out = os.path.join(self.app.out_dir, base + FC.SNAPSHOT_EXT)
                data = FC.snapshot(job.driver, want, vars_, meta={
                    "geometry": job.geometry, "case": base,
                    "params": job.params, "settings": job.settings})
                FC.write_snapshot(out, data)
                job.log("saved %s (%d surface(s), %d variable(s), %.1f MB)"
                        % (os.path.basename(out), len(want), len(vars_),
                           os.path.getsize(out) / 1e6))
                return self._json({"ok": True, "file": os.path.basename(out),
                                   "bytes": os.path.getsize(out),
                                   "surfaces": want, "variables": vars_})
            if path == "/api/open_dir":
                how = open_in_file_manager(self.app.out_dir)
                return self._json({"ok": True, "opened": self.app.out_dir,
                                   "with": how})
            if path == "/api/outline":
                job = self.app.need_results()
                if job.case is None:
                    raise ValueError("this result has no geometry to outline")
                return self._json({"ok": True, "lines": FC.outline(job.case)})
            if path == "/api/bbox":
                job = self.app.need_results()
                return self._json({"ok": True, "bbox": job.bbox()})
            if path == "/api/plane":
                job = self.app.need_results()
                name = str(body.get("name") or "").strip()
                axis = str(body.get("axis") or "x")
                if not name or not PLANE_NAME.match(name):
                    raise ValueError("a plane name may use letters, digits, "
                                     "'-' and '_' only")
                if name in FC.PATCHES:
                    raise ValueError("%s is a boundary; pick another name" % name)
                value = float(body.get("value"))
                job.driver.make_plane(name, axis, value)
                job.planes[name] = {"axis": axis, "value": value}
                return self._json({"ok": True, "planes": job.planes,
                                   "surfaces": job.surface_list()})
            if path == "/api/plane_delete":
                job = self.app.need_results()
                name = str(body.get("name") or "")
                job.driver.drop_plane(name)
                job.planes.pop(name, None)
                return self._json({"ok": True, "planes": job.planes,
                                   "surfaces": job.surface_list()})
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


def code_version():
    """What this server process is actually running.

    Python imports a module once.  Refreshing the browser re-fetches the page
    and nothing else, so a change to fluent_case.py or mesh_explorer.py only
    takes effect when this process is restarted - and a run that keeps failing
    on a bug that was fixed is, nearly always, a server that was never
    restarted.  Stamping the revision on the console, on /api/info and on the
    first line of every run log makes that visible instead of a guess.
    """
    import subprocess
    out = {"commit": "", "subject": "", "date": "", "dirty": None, "how": ""}
    try:
        def git(*a):
            return subprocess.check_output(("git", "-C", HERE) + a,
                                           stderr=subprocess.DEVNULL,
                                           timeout=10).decode("utf-8", "replace").strip()
        out["commit"] = git("rev-parse", "--short", "HEAD")
        out["date"], out["subject"] = git("log", "-1", "--format=%cs%n%s").split("\n", 1)
        out["dirty"] = bool(git("status", "--porcelain", "--", "*.py", "*.html"))
        out["how"] = "git"
    except Exception:                                   # noqa: BLE001
        #  no git, or not a checkout: the newest source file still says
        #  whether the running process could possibly be current
        newest = 0.0
        for f in ("app.py", "fluent_case.py", "mesh_explorer.py", "mesh_explorer.html"):
            try:
                newest = max(newest, os.path.getmtime(os.path.join(HERE, f)))
            except OSError:
                pass
        out["date"] = time.strftime("%Y-%m-%d %H:%M", time.localtime(newest))
        out["subject"] = "newest source file"
        out["how"] = "mtime"
    return out


def version_line(v):
    """One line naming the code, for the console and the run log."""
    if v["how"] == "git":
        return "%s%s (%s) %s" % (v["commit"], "+local edits" if v["dirty"] else "",
                                 v["date"], v["subject"])
    return "no git here; %s %s" % (v["subject"], v["date"])


VERSION = code_version()
STARTED = time.strftime("%Y-%m-%d %H:%M:%S")


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

    installed, detail = FC.fluent_installed()
    print("=" * 66)
    print(" bundle CFD app   %s" % url)
    print(" code           : %s" % version_line(VERSION))
    if a.backend == "mock":
        print(" backend        : mock (forced) - nothing it produces is a result")
    elif installed:
        print(" backend        : %s -> Fluent %s" % (a.backend, detail))
    else:
        print(" backend        : %s -> no Fluent found, the MOCK will be used" % a.backend)
        print("                  %s" % detail)
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
