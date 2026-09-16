# -*- coding: utf-8 -*-
"""The parametric campaign - stages 2 and 3.

A study is a DEFINITION plus a set of RESULTS, kept apart on purpose.

  study.json    what is to be run.  Small, exact, reviewable, and written
                before anything is launched - so the matrix is a decision
                somebody made rather than whatever happened to get run.
  results.json  one record per case as it comes back, appended as it goes.
                Resumable: a study that stops after nine of fifty cases picks
                up at the tenth.

Nothing here invents a number.  A case that fails is recorded as failed, a
case run on the mock backend is recorded as mock and is excluded from every
fit and every comparison - a MOCK value is plumbing evidence, never a result.

TWO STUDY KINDS
---------------
mesh   one geometry, one flow condition, the mesh refined step by step.  The
       answer is a grid-convergence index and the coarsest mesh that is
       already inside tolerance, because the point of a mesh study is to find
       the CHEAPEST adequate mesh, not the finest one affordable.

sweep  one mesh resolution, the arrangement and the flow condition varied.
       The answer is Eu per row over the matrix, against the correlations,
       and the coefficients of a fit of our own.

WHY THE SIDE AND END WALLS ARE SYMMETRY PLANES
----------------------------------------------
Every correlation in correlations.py is for a bank that is infinitely wide and
made of infinitely long tubes.  A box with four real walls is not that: the
side walls add a passage of their own and the top and bottom add a boundary
layer the correlation never saw, and both biases grow as the domain shrinks.
Running them as symmetry planes - with half rods on the sides, so the mirror
lands on a rod centreline - gives back exactly the bank the correlation
describes.  It also makes the problem effectively two-dimensional, which is
why these cases are cheap enough to run dozens of.

That is a deliberate choice about what is being compared, and it is recorded
in the study definition rather than left in somebody's head.
"""

import json
import math
import os
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import correlations as CR              # noqa: E402
import mesh_explorer as ME             # noqa: E402

STUDIES = os.path.join(HERE, "studies")
KINDS = ("mesh", "sweep")


# =============================================================================
#  THE MESHING RULES THE STUDY APPLIES PER CASE
# =============================================================================
#  A parametric study must not change mesh QUALITY while it changes the
#  geometry, or the trend it measures is partly the mesh's.  Two rules keep it
#  fixed, and both are applied here rather than left to whoever fills in the
#  panel:
#
#    * the first cell is placed at a fixed y+, which means its height changes
#      with the velocity and with the pitch (both move u_max), and
#    * enough radial layers are used to get from that first cell out to the
#      gap at a bounded growth ratio, which means the layer COUNT changes too.
#
#  Holding n_rad fixed instead would let the growth ratio run from 1.2 at the
#  slowest case to 3 at the fastest, and the fastest cases would come back
#  wrong for a reason that has nothing to do with the physics being studied.

def first_layer_for(st, y_plus, export_scale):
    """First CELL height in model units, for a target y+.

    Same rule the Geometry tab's y+ button uses: Blasius skin friction on the
    gap velocity, friction velocity from it, and the cell centre at the target
    y+ - so the cell is twice that.
    """
    Re = max(st["Re"], 1.0)
    cf = 0.079 * Re ** -0.25                    # Fanning, Blasius
    u_tau = st["umax"] * math.sqrt(max(cf, 1e-9) / 2.0)
    y1 = 2.0 * y_plus * st["nu"] / max(u_tau, 1e-12)      # metres
    return y1 / export_scale


def _layer_span(first, q, n):
    """How far n layers reach, starting at `first` and growing by q."""
    return first * n if abs(q - 1.0) < 1e-12 else first * (q ** n - 1.0) / (q - 1.0)


def n_rad_for(clearance, first_layer, max_growth, floor=4, ceiling=80):
    """The fewest radial layers that get from `first_layer` out to `clearance`
    without any step growing by more than `max_growth`.

    The reach has to be checked, not just the ratio.  `growth_ratio` returns
    1.0 both when a uniform spacing happens to fit AND when no ratio in its
    bracket can span the distance - and the second case is a mesh that stops
    short of the gap, which would sail through a test on the ratio alone.
    That is exactly what happened at Re 1e5, where the first cell is 4 um and
    four layers cannot reach a 1.25 mm gap however fast they grow.
    """
    for n in range(floor, ceiling + 1):
        q = ME.growth_ratio(clearance, first_layer, n)
        if q <= max_growth and _layer_span(first_layer, q, n) >= clearance * 0.999:
            return n
    return ceiling


def case_state(geometry, params, settings):
    """The flow state a case will be in, before it is run.

    This is what decides the first layer and what the correlations are asked
    about, so it is computed from the same Case the mesher will build rather
    than from the numbers in the panel.
    """
    import app as APP                            # local: app imports us back
    case = APP.case_from(geometry, params)
    k = case.export_scale
    mat, inl = settings["material"], settings["inlet"]
    nu = mat["viscosity"] / mat["density"]
    return case, CR.flow_state(
        D=case.D * k, ST=case.ST * k, SL=case.SL * k, n_rows=case.n_rows,
        u_in=inl["velocity"], rho=mat["density"], nu=nu,
        staggered=case.stagger, dT=2.0 * case.bE * k, dL=2.0 * case.aE * k,
        helix=(case.helix if case.family == "helical" else 0.0),
        coil_alternating=False)


def apply_mesh_rules(geometry, params, settings, y_plus=1.0, max_growth=1.2):
    """Fill in firstLayer and nRad for one case.  Returns what it decided."""
    case, st = case_state(geometry, params, settings)
    fl = first_layer_for(st, y_plus, case.export_scale)
    n_rad = n_rad_for(case.clearance, fl, max_growth)
    params = dict(params)
    params["firstLayer"] = float("%.4g" % fl)
    params["nRad"] = int(n_rad)
    growth = ME.growth_ratio(case.clearance, params["firstLayer"], n_rad)
    return params, {"y_plus": y_plus, "first_layer": params["firstLayer"],
                    "n_rad": n_rad, "growth": round(growth, 4),
                    "clearance": round(case.clearance, 5),
                    "umax": st["umax"], "Re": st["Re"],
                    "XT": st["XT"], "XL": st["XL"],
                    "diagonal": st["diagonal"]}


def velocity_for_Re(geometry, params, settings, Re_target):
    """The inlet velocity that puts Re_max on target.

    Re is built on u_max, and u_max is u_in times S_T/gap, so the velocity
    that lands on a wanted Re depends on the pitch.  A sweep stated in Re -
    which is the only way to state one, since Re is what the correlations are
    written on - has to solve this for every case.
    """
    probe = dict(settings)
    probe["inlet"] = dict(settings["inlet"])
    probe["inlet"]["velocity"] = 1.0
    case, st = case_state(geometry, params, probe)
    return Re_target * st["nu"] / (case.D * case.export_scale) / st["umax"]


# =============================================================================
#  THE STUDY OBJECT
# =============================================================================
class Study(object):
    """A definition on disk, and the results beside it."""

    def __init__(self, d):
        self.d = d
        for key in ("name", "kind", "geometry", "base", "cases"):
            if key not in d:
                raise ValueError("a study needs %r" % key)
        if d["kind"] not in KINDS:
            raise ValueError("kind must be one of %s" % ", ".join(KINDS))
        ids = [c["id"] for c in d["cases"]]
        if len(set(ids)) != len(ids):
            raise ValueError("two cases share an id")

    # -- where it lives ---------------------------------------------------
    @property
    def dir(self):
        return os.path.join(STUDIES, self.d["name"])

    @property
    def name(self):
        return self.d["name"]

    @property
    def kind(self):
        return self.d["kind"]

    @property
    def cases(self):
        return self.d["cases"]

    def case(self, cid):
        for c in self.cases:
            if c["id"] == cid:
                return c
        raise KeyError("no case %r in %s" % (cid, self.name))

    #  a case is the base with its own overrides on top, never a fresh dict:
    #  a matrix where each case restated the whole set-up would drift entry by
    #  entry, and the difference between two cases is the whole point
    def params_for(self, cid):
        p = dict(self.d["base"].get("params") or {})
        p.update(self.case(cid).get("params") or {})
        return p

    def settings_for(self, cid):
        import fluent_case as FC
        s = {}
        for src in (self.d["base"].get("settings") or {},
                    self.case(cid).get("settings") or {}):
            for gid, grp in src.items():
                s.setdefault(gid, {}).update(grp or {})
        return FC.merge_settings(s)

    # -- on disk ----------------------------------------------------------
    def save(self):
        os.makedirs(self.dir, exist_ok=True)
        path = os.path.join(self.dir, "study.json")
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.d, fh, ensure_ascii=False, indent=1, sort_keys=False)
            fh.write("\n")
        return path

    @staticmethod
    def load(name):
        path = os.path.join(STUDIES, name, "study.json")
        if not os.path.isfile(path):
            raise ValueError("no study %r in %s" % (name, STUDIES))
        with open(path, encoding="utf-8") as fh:
            return Study(json.load(fh))

    @staticmethod
    def list_all():
        out = []
        if not os.path.isdir(STUDIES):
            return out
        for name in sorted(os.listdir(STUDIES)):
            if os.path.isfile(os.path.join(STUDIES, name, "study.json")):
                try:
                    s = Study.load(name)
                except Exception:                       # noqa: BLE001
                    continue
                rows = s.results()
                out.append({"name": name, "kind": s.kind,
                            "title": s.d.get("title", name),
                            "geometry": s.d["geometry"],
                            "cases": len(s.cases),
                            "done": len([r for r in rows if r.get("ok")]),
                            "failed": len([r for r in rows if not r.get("ok")]),
                            "mock": len([r for r in rows if r.get("mock")])})
        return out

    # -- results ----------------------------------------------------------
    @property
    def results_path(self):
        return os.path.join(self.dir, "results.json")

    def results(self):
        try:
            with open(self.results_path, encoding="utf-8") as fh:
                return json.load(fh).get("rows") or []
        except (IOError, ValueError):
            return []

    def record(self, row):
        """Append or replace one case's result, on disk, immediately.

        Immediately, because a campaign that only wrote its results at the end
        would lose a night's runs to one crash on the last case.
        """
        rows = [r for r in self.results() if r.get("id") != row["id"]]
        rows.append(row)
        order = {c["id"]: i for i, c in enumerate(self.cases)}
        rows.sort(key=lambda r: order.get(r.get("id"), 1e9))
        os.makedirs(self.dir, exist_ok=True)
        with open(self.results_path, "w", encoding="utf-8") as fh:
            json.dump({"study": self.name, "rows": rows}, fh,
                      ensure_ascii=False, indent=1)
            fh.write("\n")

    def clear_results(self):
        if os.path.isfile(self.results_path):
            os.unlink(self.results_path)

    def usable(self):
        """The rows a conclusion may be drawn from: finished, not mock."""
        return [r for r in self.results() if r.get("ok") and not r.get("mock")]


# =============================================================================
#  GRID CONVERGENCE
# =============================================================================
def _thin(rows, keep):
    """Every nth row, ends included, so a history fits in a result file."""
    n = len(rows)
    if n <= keep:
        return list(rows)
    step = (n - 1) / float(keep - 1)
    out = [rows[int(round(i * step))] for i in range(keep)]
    out[-1] = rows[-1]
    return out


def gci(h, phi, safety=1.25):
    """Roache's grid-convergence index, by the procedure in

        I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman and
        P. E. Raad, "Procedure for estimation and reporting of uncertainty due
        to discretization in CFD applications", J. Fluids Eng. 130 (2008)
        078001,

    which is also what ASME V&V 20 asks for.  `h` and `phi` are three meshes
    fine-to-coarse: h1 < h2 < h3.

    The observed order p solves

        p = |ln|eps32/eps21| + q(p)| / ln(r21),
        q(p) = ln((r21^p - s)/(r32^p - s)),   s = sign(eps32/eps21)

    which is implicit whenever the two refinement ratios differ - and here they
    always do, because a block mesh refines by integer division counts and
    never lands on the ratio asked for.  Fixed-point iteration, as the paper
    prescribes.

    Oscillatory convergence (s < 0) is reported rather than hidden: it means
    the three points do not sit on a monotone curve and the extrapolation is
    not meaningful, which a reader has to be told.
    """
    h1, h2, h3 = [float(x) for x in h]
    p1, p2, p3 = [float(x) for x in phi]
    if not (h1 < h2 < h3):
        raise ValueError("gci wants three meshes fine to coarse (h1 < h2 < h3)")
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = p2 - p1, p3 - p2
    out = {"r21": r21, "r32": r32, "eps21": e21, "eps32": e32,
           "phi1": p1, "phi2": p2, "phi3": p3, "h1": h1, "h2": h2, "h3": h3}
    if abs(e21) < 1e-14:
        out["note"] = "the two finest meshes give the same answer to 14 digits"
        out.update(p=None, phi_ext=p1, gci_fine=0.0, gci_coarse=0.0,
                   e_a21=0.0, e_ext21=0.0, oscillatory=False)
        return out
    ratio = e32 / e21
    s = 1.0 if ratio > 0 else -1.0
    out["oscillatory"] = ratio < 0
    p = 2.0
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(ratio)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-12:
            p = p_new
            break
        p = 0.5 * (p + p_new)                    # damped: it can oscillate
    out["p"] = p
    rp = r21 ** p
    out["phi_ext"] = (rp * p1 - p2) / (rp - 1.0) if abs(rp - 1.0) > 1e-12 else p1
    out["e_a21"] = abs((p1 - p2) / p1) if p1 else None
    out["e_ext21"] = (abs((out["phi_ext"] - p1) / out["phi_ext"])
                      if out["phi_ext"] else None)
    out["gci_fine"] = (safety * out["e_a21"] / (rp - 1.0)
                       if out["e_a21"] is not None and abs(rp - 1.0) > 1e-12
                       else None)
    rp32 = r32 ** p
    e_a32 = abs((p2 - p3) / p2) if p2 else None
    out["gci_coarse"] = (safety * e_a32 / (rp32 - 1.0)
                         if e_a32 is not None and abs(rp32 - 1.0) > 1e-12
                         else None)
    return out


def mesh_analysis(study, tol=0.01, key="dp_bundle"):
    """Every consecutive triplet of a ladder, and which mesh is good enough.

    "Good enough" is the COARSEST mesh whose distance from the Richardson
    extrapolation is inside `tol` - the point of a mesh study is the cheapest
    adequate mesh, not the finest one that fits in the night.
    """
    rows = [r for r in study.usable() if r.get(key) is not None and r.get("h")]
    out = {"ladders": [], "tol": tol, "key": key}
    #  A grid-convergence index measures the difference between CONVERGED
    #  solutions on different meshes.  Fed a case whose residuals stalled, it
    #  measures the difference between two half-solved ones and returns a
    #  number that looks exactly like an answer.  A case that did not meet its
    #  own residual criterion is therefore kept in the table - so it is
    #  visible - and kept out of every extrapolation.
    out["unconverged"] = [r["id"] for r in rows if r.get("converged") is False]
    ladders = {}
    for r in rows:
        ladders.setdefault(r.get("ladder", "main"), []).append(r)
    #  in the order the study defines them, not alphabetically: the
    #  representative ladder is meant to be read before the corner one
    order = []
    for c in study.cases:
        lad = c.get("ladder", "main")
        if lad not in order:
            order.append(lad)
    for name in [x for x in order if x in ladders] + \
                [x for x in sorted(ladders) if x not in order]:
        rs = sorted(ladders[name], key=lambda r: r["h"])       # fine first
        good = [r for r in rs if r.get("converged") is not False]
        lad = {"ladder": name, "levels": rs, "converged": [r["id"] for r in good],
               "triplets": [], "chosen": None}
        if len(good) < 3:
            lad["blocked"] = ("fewer than three converged meshes: no grid "
                              "convergence index can be formed")
            lad["blocked_ko"] = ("수렴한 격자가 세 개 미만입니다. 격자 수렴 "
                                 "지수를 만들 수 없습니다.")
            out["ladders"].append(lad)
            continue
        rs = good
        for i in range(len(rs) - 2):
            trio = rs[i:i + 3]
            try:
                g = gci([t["h"] for t in trio], [t[key] for t in trio])
            except Exception as exc:                            # noqa: BLE001
                g = {"error": str(exc)}
            g["ids"] = [t["id"] for t in trio]
            lad["triplets"].append(g)
        finest = lad["triplets"][0] if lad["triplets"] else None
        if finest and finest.get("phi_ext"):
            ref = finest["phi_ext"]
            lad["phi_ext"] = ref
            #  coarsest first, so the first one inside tolerance is the answer
            for r in sorted(rs, key=lambda r: -r["h"]):
                r["dev_from_ext"] = abs(r[key] - ref) / abs(ref)
                if lad["chosen"] is None and r["dev_from_ext"] <= tol:
                    lad["chosen"] = r["id"]
            if lad["chosen"] is None:
                lad["chosen"] = rs[0]["id"]
                lad["chosen_note"] = ("no mesh is inside the tolerance; the "
                                      "finest is named so that something is")
        out["ladders"].append(lad)
    return out


def sweep_analysis(study, keys=None):
    """Every usable case against every evaluable correlation, and our own fit."""
    rows = [r for r in study.usable() if r.get("eu_row")]
    keys = keys or CR.ENCODED
    out = {"rows": [], "keys": keys, "fit": None, "ratios": {}}
    for r in rows:
        st = CR.flow_state(D=r["D"], ST=r["ST"], SL=r["SL"],
                           n_rows=r["n_rows"], u_in=r["u_in"], rho=r["rho"],
                           nu=r["nu"], staggered=r["staggered"])
        rec = {"id": r["id"], "XT": st["XT"], "XL": st["XL"], "Re": st["Re"],
               "XT_d": st["XT_d"], "XL_d": st["XL_d"],
               "cos_eps": st["cos_eps"], "coil_K": st["coil_K"],
               "helix": st["helix"],
               "eu_row": r["eu_row"], "dp_bundle": r["dp_bundle"],
               "staggered": st["staggered"], "corr": {}}
        for k in keys:
            e = CR.evaluate(k, st)
            rec["corr"][k] = {"eu_row": e["eu_row"],
                              "ratio": (r["eu_row"] / e["eu_row"]
                                        if e["eu_row"] else None),
                              "outside": e["outside"]}
        #  The BAND: the span of the correlations that are in range here.
        #  Stage 1 measured that the published correlations differ from one
        #  another by 1.0x to 2.2x, so "within 20 % of Jakob" is not a test of
        #  anything - the literature is not that precise.  What can be tested
        #  is whether the CFD falls inside the band they span, and that is the
        #  number this study reports.
        span = [c["eu_row"] for c in rec["corr"].values()
                if c["eu_row"] and not c["outside"]]
        if len(span) >= 2:
            rec["band"] = [min(span), max(span)]
            rec["band_n"] = len(span)
            rec["in_band"] = rec["band"][0] <= r["eu_row"] <= rec["band"][1]
            rec["band_dev"] = (0.0 if rec["in_band"] else
                               (r["eu_row"] / rec["band"][1] - 1.0
                                if r["eu_row"] > rec["band"][1]
                                else r["eu_row"] / rec["band"][0] - 1.0))
        out["rows"].append(rec)
    for k in keys:
        rs = [x["corr"][k]["ratio"] for x in out["rows"]
              if x["corr"][k]["ratio"]]
        if rs:
            out["ratios"][k] = {"n": len(rs), "min": min(rs), "max": max(rs),
                                "mean": sum(rs) / len(rs),
                                "rms": math.sqrt(sum((v - 1.0) ** 2 for v in rs)
                                                 / len(rs))}
    banded = [r for r in out["rows"] if "in_band" in r]
    if banded:
        out["band"] = {
            "n": len(banded),
            "inside": len([r for r in banded if r["in_band"]]),
            "worst": max(abs(r["band_dev"]) for r in banded),
            "width_mean": (sum(r["band"][1] / r["band"][0] for r in banded)
                           / len(banded)),
        }
    if len(out["rows"]) >= 6:
        helical = any(r.get("cos_eps", 1.0) < 1.0 for r in out["rows"])
        try:
            out["fit"] = CR.fit_form(out["rows"], with_helix=helical)
            #  and the same fit with Shen's own restriction, p = q, so the
            #  question "does splitting the pitch exponents earn its keep?"
            #  is answered on this data rather than argued about
            out["fit_tied"] = CR.fit_form(out["rows"], with_helix=helical,
                                          tie_pitch=True)
        except Exception as exc:                                # noqa: BLE001
            out["fit_error"] = str(exc)
    return out


# =============================================================================
#  BUILDING THE CAMPAIGN
# =============================================================================
#  The base case.  Water at 20 C, a ten-row bank, and every wall that is not a
#  rod made a symmetry plane - see the module header for why that is the
#  configuration a correlation can be compared against at all.
WATER = {"name": "water-liquid", "density": 998.2, "viscosity": 1.003e-3}

BASE_SETTINGS = {
    "general": {"steady": True, "energy": False},
    "turbulence": {"viscous": "k-omega", "k_omega_variant": "sst"},
    "material": dict(WATER),
    "inlet": {"spec": "components", "velocity": 1.0,
              "turb_spec": "Intensity and Viscosity Ratio",
              "intensity": 5.0, "visc_ratio": 10.0},
    "outlet": {"gauge_pressure": 0.0, "prevent_reverse_flow": True},
    "zones": {"wall_rods": "wall",
              "wall_side_y0": "symmetry", "wall_side_ymax": "symmetry",
              "wall_bottom": "symmetry", "wall_top": "symmetry"},
    "methods": {"flow_scheme": "Coupled", "gradient": "least-square-cell-based",
                "pressure": "second-order", "momentum": "second-order-upwind",
                "turb": "second-order-upwind", "pseudo_time": True},
    "run": {"init_method": "hybrid", "iterations": 800,
            "residual_criterion": 1e-5, "monitor_every": 20,
            "write_case": False},
}

BASE_PARAMS = {
    "D": 10.0, "nRows": 10, "nCols": 4,
    "lIn": 60.0, "lOut": 120.0,
    #  a thin slab between two symmetry planes: the case is effectively 2-D,
    #  so z carries two cells and no more
    "H": 2.0, "nZ": 2,
    "halfRods": True, "firstOffset": False, "unit": "mm",
    "nAz": 29, "nRad": 12, "firstLayer": 0.03, "nxIn": 18, "nxOut": 29,
}

#  The refinement ladder.  In-plane counts rise by about 1.35 a step; the
#  radial direction refines by TIGHTENING the growth ratio rather than by
#  adding layers directly, because the first cell is pinned at y+ = 1 and must
#  not move - refining it would change the near-wall treatment, not the mesh.
#
#  Celik recommends a refinement ratio of at least 1.3 between levels.  A block
#  mesh cannot hit a ratio exactly - the counts are integers and the radial
#  layer count is decided by the growth rule, not chosen - so the ratios that
#  come out of this ladder run between about 1.25 and 1.40.  That is why the
#  GCI procedure carries an implicit p at all, and why the report prints the
#  ratio each triplet actually had and says so when one falls short.
LADDER = [
    {"id": "L1", "nAz": 14, "nxIn": 8, "nxOut": 14, "growth": 1.45},
    {"id": "L2", "nAz": 20, "nxIn": 12, "nxOut": 20, "growth": 1.35},
    {"id": "L3", "nAz": 29, "nxIn": 17, "nxOut": 29, "growth": 1.26},
    {"id": "L4", "nAz": 42, "nxIn": 25, "nxOut": 42, "growth": 1.18},
    {"id": "L5", "nAz": 60, "nxIn": 36, "nxOut": 60, "growth": 1.12},
]

Y_PLUS = 1.0                    # k-omega SST integrates to the wall


def _case(geometry, cid, label, params, settings, meta=None, ladder=None):
    params = dict(params)
    growth = params.pop("_growth", 1.25)
    p, rule = apply_mesh_rules(geometry, params, settings, Y_PLUS, growth)
    c = {"id": cid, "label": label, "params": p, "mesh_rule": rule}
    if settings.get("inlet", {}).get("velocity") is not None:
        c["settings"] = {"inlet": {"velocity": settings["inlet"]["velocity"]}}
    if meta:
        c["meta"] = meta
    if ladder:
        c["ladder"] = ladder
    return c


def _settings_at(velocity):
    import fluent_case as FC
    s = FC.merge_settings(BASE_SETTINGS)
    s["inlet"]["velocity"] = velocity
    return s


def make_mesh_study(geometry, name=None):
    """A refinement ladder at a representative condition, and a shorter one at
    the demanding corner of the stage-3 matrix.

    Two ladders, because a mesh chosen at X = 1.5 and Re = 1e4 is not
    automatically adequate at X = 1.25 and Re = 1e5, where the gap is half as
    wide and the first cell an order of magnitude thinner.  The second ladder
    is the check that the first one's answer travels.
    """
    name = name or (geometry + "-mesh")
    conditions = [
        ("main", 1.5, 1.0e4, "대표 조건 · representative", LADDER),
        ("corner", 1.25, 1.0e5, "가장 까다로운 조건 · the demanding corner",
         LADDER[2:]),
    ]
    cases = []
    for lad_name, X, Re, label, levels in conditions:
        base = dict(BASE_PARAMS)
        base.update(ST=10.0 * X, SL=10.0 * X)
        vel = velocity_for_Re(geometry, base, _settings_at(1.0), Re)
        st = _settings_at(vel)
        for lv in levels:
            p = dict(base)
            p.update(nAz=lv["nAz"], nxIn=lv["nxIn"], nxOut=lv["nxOut"])
            p["_growth"] = lv["growth"]
            cases.append(_case(
                geometry, "%s-%s" % (lad_name, lv["id"]),
                "%s · %s" % (label, lv["id"]), p, st,
                meta={"X": X, "Re_target": Re, "velocity": vel},
                ladder=lad_name))
    d = {
        "name": name, "kind": "mesh", "geometry": geometry,
        "title": "2단계 · %s 격자 민감도" % geometry,
        "created": time.strftime("%Y-%m-%d"),
        "phi": "dp_bundle",
        "notes": [
            "물 20 C, k-omega SST, y+ = 1, 정상상태.",
            "rod 외 모든 벽은 대칭면 - 상관식이 기술하는 무한 다발과 같게 하기 위함.",
            "격자는 면내 방향으로 단계당 약 1.35배 세밀해지고, 반경 방향은 "
            "첫 셀을 y+ = 1에 고정한 채 성장비를 조여서 세밀해집니다.",
            "판정은 Celik 등(2008) / ASME V&V 20 의 GCI 절차로 합니다.",
        ],
        "base": {"params": dict(BASE_PARAMS), "settings": dict(BASE_SETTINGS)},
        "cases": cases,
    }
    return Study(d)


#  Stage 3.  X_T and X_L are varied INDEPENDENTLY on purpose: on the square
#  diagonal (X_T = X_L) the two pitch terms of any correlation of the form
#  (X_T-1)^-p X_L^q are collinear and neither coefficient can be identified.
#  A sweep that only walked the diagonal would produce a fit that could not be
#  fitted, which is a mistake that is invisible until the solve is singular.
SWEEP_XT = (1.25, 1.5, 2.0)
SWEEP_XL = (1.25, 1.5, 2.0)
SWEEP_RE = (2.0e3, 1.0e4, 4.0e4)
SWEEP_RE_EXTRA = 1.0e5          # on the diagonal only - outside Jakob's range


def make_sweep_study(geometry, level="L3", name=None):
    """The arrangement/condition matrix, on one mesh resolution."""
    name = name or (geometry + "-sweep")
    lv = [x for x in LADDER if x["id"] == level]
    if not lv:
        raise ValueError("no ladder level %r" % level)
    lv = lv[0]
    cases = []
    grid = [(xt, xl, re) for xt in SWEEP_XT for xl in SWEEP_XL
            for re in SWEEP_RE]
    grid += [(x, x, SWEEP_RE_EXTRA) for x in SWEEP_XT]
    for xt, xl, re in grid:
        p = dict(BASE_PARAMS)
        p.update(ST=10.0 * xt, SL=10.0 * xl,
                 nAz=lv["nAz"], nxIn=lv["nxIn"], nxOut=lv["nxOut"])
        p["_growth"] = lv["growth"]
        vel = velocity_for_Re(geometry, p, _settings_at(1.0), re)
        cid = "T%s_L%s_Re%s" % (("%g" % xt).replace(".", ""),
                                ("%g" % xl).replace(".", ""),
                                ("%.0e" % re).replace("+0", "").replace("+", ""))
        cases.append(_case(
            geometry, cid,
            "X_T %.2f · X_L %.2f · Re %.0e" % (xt, xl, re),
            p, _settings_at(vel),
            meta={"XT": xt, "XL": xl, "Re_target": re, "velocity": vel}))
    d = {
        "name": name, "kind": "sweep", "geometry": geometry,
        "title": "3단계 · %s 압력강하 스윕" % geometry,
        "created": time.strftime("%Y-%m-%d"),
        "phi": "eu_row",
        "notes": [
            "물 20 C, k-omega SST, y+ = 1, 정상상태, rod 외 벽은 모두 대칭면.",
            "격자 해상도는 2단계에서 고른 수준(%s)을 씁니다." % level,
            "X_T 와 X_L 을 독립적으로 바꿉니다. 정사각 피치만 훑으면 피치 항 두 개가 "
            "서로 구분되지 않아 상관식 적합 자체가 불가능해집니다.",
            "Re 1e5 점은 Jakob 적용 범위 밖이며, 1단계에서 두 상관식이 크게 "
            "갈라지는 것으로 확인된 영역을 일부러 지나갑니다.",
        ],
        "base": {"params": dict(BASE_PARAMS), "settings": dict(BASE_SETTINGS)},
        "cases": cases,
    }
    return Study(d)


# =============================================================================
#  RUNNING IT
# =============================================================================
class Runner(object):
    """One study, one case after another, through the app's own Job.

    It drives the same Job the Run tab drives and hands it to the app as the
    current job, so a campaign is watchable: the residuals, the pressure-drop
    monitor, the log and the Results tab all follow whichever case is running.
    Nothing about a study case is special except that something else pressed
    the button.
    """

    def __init__(self, app, study, only=None, redo=False, force=False):
        self.app = app
        self.study = study
        self.redo = bool(redo)
        #  A campaign is a bet that the case SET-UP is right, repeated N times.
        #  This one cost 55 minutes to discover that every case stalled at a
        #  residual of 7.5e-2, because nothing looked at the first result
        #  before starting the second.  Unless overridden, the first case has
        #  to converge before the rest are allowed to run.
        self.force = bool(force)
        self.stopped_reason = None
        done = {r["id"] for r in study.results() if r.get("ok")} \
            if not redo else set()
        wanted = set(only) if only else None
        self.queue = [c["id"] for c in study.cases
                      if (wanted is None or c["id"] in wanted)
                      and c["id"] not in done]
        self.skipped = [c["id"] for c in study.cases
                        if (wanted is None or c["id"] in wanted)
                        and c["id"] in done]
        self.index = 0
        self.current = None
        self.stopping = False
        self.finished = False
        self.error = None
        self.started = time.time()
        self.log_lines = []
        self.thread = None

    # -- reporting --------------------------------------------------------
    def log(self, msg):
        self.log_lines.append("%7.1fs  %s" % (time.time() - self.started, msg))
        del self.log_lines[:-500]

    def status(self):
        return {
            "study": self.study.name, "kind": self.study.kind,
            "total": len(self.queue), "index": self.index,
            "current": self.current, "queue": list(self.queue),
            "skipped": list(self.skipped),
            "stopping": self.stopping, "finished": self.finished,
            "error": self.error, "stopped_reason": self.stopped_reason,
            "elapsed": round(time.time() - self.started, 1),
            "log": self.log_lines[-60:],
            "rows": self.study.results(),
        }

    # -- the loop ---------------------------------------------------------
    def start(self):
        self.thread = threading.Thread(target=self._run,
                                       name="study-" + self.study.name)
        self.thread.daemon = True
        self.thread.start()

    def stop(self):
        self.stopping = True
        job = self.app.job
        if job is not None and not job.finished_at:
            job.stop()

    def _run(self):
        try:
            if self.skipped:
                self.log("%d case(s) already done, skipping them: %s"
                         % (len(self.skipped), ", ".join(self.skipped)))
            for cid in self.queue:
                if self.stopping:
                    self.log("stopped before " + cid)
                    break
                self.current = cid
                self.index += 1
                self.log("[%d/%d] %s" % (self.index, len(self.queue), cid))
                try:
                    row = self._one(cid)
                except Exception as exc:                        # noqa: BLE001
                    row = {"id": cid, "ok": False,
                           "error": "%s: %s" % (type(exc).__name__, exc)}
                    self.log("  %s FAILED: %s" % (cid, row["error"]))
                self.study.record(row)
                if row.get("ok"):
                    self.log("  %s: dp_bundle %.4g Pa, Eu_row %.4f, %d cells%s"
                             % (cid, row["dp_bundle"] or float("nan"),
                                row["eu_row"] or float("nan"), row["cells"],
                                "   (MOCK)" if row.get("mock") else ""))
                #  the smoke test: one case is enough to tell whether the
                #  set-up produces a converged solution at all
                if self.index == 1 and not self.force and not self.stopping:
                    why = self._first_case_verdict(row)
                    if why:
                        self.stopped_reason = why
                        self.log("STOPPING after the first case: " + why)
                        self.log("  nothing is wrong with the study definition "
                                 "- the SET-UP does not produce a usable "
                                 "solution, and 7 more of the same would not "
                                 "change that. Fix it, or re-run with force.")
                        break
            self.current = None
        except Exception as exc:                                # noqa: BLE001
            self.error = "%s: %s" % (type(exc).__name__, exc)
            self.log("the runner itself failed: " + self.error)
        finally:
            self.finished = True
            self.log("done - %d of %d case(s) attempted"
                     % (self.index, len(self.queue)))

    @staticmethod
    def _first_case_verdict(row):
        """Why the campaign should not continue past case one, or None."""
        if not row.get("ok"):
            return "it failed: %s" % (row.get("error") or "?")
        if row.get("mock"):
            return None                     # the mock is for testing the plumbing
        if row.get("converged") is False:
            return ("it did not converge - final residual %.2e%s against a "
                    "criterion of %.0e"
                    % (row.get("residual_worst") or 0.0,
                       " in " + row["residual_worst_eq"]
                       if row.get("residual_worst_eq") else "",
                       row.get("criterion") or 0.0))
        if row.get("dp_bundle") is None:
            return "the bundle pressure drop could not be measured"
        if (row.get("dp_drift") or 0.0) > 0.05:
            return ("the pressure drop was still moving by %.1f %% over the "
                    "final stretch" % (100 * row["dp_drift"]))
        return None

    def _one(self, cid):
        """Run one case and measure it."""
        import app as APP
        import fluent_case as FC
        s = self.study
        c = s.case(cid)
        params = s.params_for(cid)
        settings = s.settings_for(cid)
        t0 = time.time()

        job = self.app.new_job(s.d["geometry"], params, settings)
        job.log("study %s, case %s" % (s.name, cid))
        job.start()
        while job.finished_at is None:
            time.sleep(0.4)
            if self.stopping:
                job.stop()
        row = {"id": cid, "ok": False, "label": c.get("label"),
               "ladder": c.get("ladder", "main"), "meta": c.get("meta"),
               "seconds": round(time.time() - t0, 1),
               "when": time.strftime("%Y-%m-%d %H:%M:%S"),
               "code": APP.version_line(APP.VERSION),
               "mock": bool(job.driver and job.driver.mock),
               "mesh_path": (os.path.basename(job.mesh_path)
                             if job.mesh_path else None),
               "error": job.error}
        if job.error:
            return row

        case = job.case
        k = case.export_scale
        mat, inl = settings["material"], settings["inlet"]
        st = CR.flow_state(D=case.D * k, ST=case.ST * k, SL=case.SL * k,
                           n_rows=case.n_rows, u_in=inl["velocity"],
                           rho=mat["density"], nu=mat["viscosity"] / mat["density"],
                           staggered=case.stagger,
                           dT=2.0 * case.bE * k, dL=2.0 * case.aE * k,
                           helix=(case.helix if case.family == "helical"
                                  else 0.0))
        row.update(D=st["D"], ST=st["ST"], SL=st["SL"], n_rows=st["n_rows"],
                   u_in=st["u_in"], rho=st["rho"], nu=st["nu"],
                   staggered=st["staggered"], umax=st["umax"], Re=st["Re"],
                   XT=st["XT"], XL=st["XL"], XT_d=st["XT_d"], XL_d=st["XL_d"],
                   helix=st["helix"], diagonal=st["diagonal"])

        ms = job.mesh_stats or {}
        row["cells"] = ms.get("cells")
        row["nodes"] = ms.get("nodes2d")
        row["skew_max"] = ms.get("dev_max")
        row["ar_max"] = ms.get("ar_max")
        row["volume"] = ms.get("vtot")
        #  representative cell size.  The top and bottom are symmetry planes
        #  and z carries a fixed number of cells, so this is a two-dimensional
        #  refinement and h is an AREA per cell, not a volume: taking the cube
        #  root of V/N would make every level look like it refined by less
        #  than it did, and the observed order would come out wrong.
        n_z = int(params.get("nZ") or 1)
        if row["volume"] and row["cells"]:
            area = row["volume"] / max(case.H * k, 1e-12)
            row["h"] = math.sqrt(area / max(row["cells"] / n_z, 1.0))
            row["h_kind"] = "2d"
        row["dp_bundle"] = None
        row["dp_total"] = None

        #  the pressure drop across the BUNDLE, the way the correlation defines
        #  it: two planes on the bundle faces, half a pitch outside the first
        #  and last rod centres.  inlet - outlet would span the boxes too.
        names = ["study_bundle_in", "study_bundle_out"]
        xs = [case.x_b0 * k, case.x_b1 * k]
        made = []
        try:
            for nm, x in zip(names, xs):
                job.driver.make_plane(nm, "x", x)
                made.append(nm)
            p = [job.driver.report("area-weighted-avg", [nm], "pressure")
                 for nm in names]
            row["p_bundle_in"], row["p_bundle_out"] = p[0], p[1]
            row["dp_bundle"] = p[0] - p[1]
            row["eu_row"] = row["dp_bundle"] / (st["n_rows"] * st["q"])
            row["dp_per_row"] = row["dp_bundle"] / max(1, st["n_rows"])
        finally:
            for nm in made:
                try:
                    job.driver.drop_plane(nm)
                except Exception:                               # noqa: BLE001
                    pass
        try:
            pi = job.driver.report("area-weighted-avg", ["inlet"], "pressure")
            po = job.driver.report("area-weighted-avg", ["outlet"], "pressure")
            row["dp_total"] = pi - po
        except Exception:                                       # noqa: BLE001
            pass

        res = list(getattr(job.driver, "residuals", []) or [])
        if res:
            last = res[-1]
            #  which equation stalled is the whole diagnosis: continuity says
            #  the pressure-velocity coupling never closed, turbulence says the
            #  model is struggling, and the two are fixed differently.  Keeping
            #  only the maximum threw that away.
            eqs = [(k, v) for k, v in last.items() if k != "iter" and v is not None]
            worst_eq, worst = (max(eqs, key=lambda kv: kv[1]) if eqs
                               else (None, 0.0))
            row["iterations"] = len(res)
            row["residual_worst"] = worst
            row["residual_worst_eq"] = worst_eq
            row["residuals_last"] = {k: v for k, v in eqs}
            row["criterion"] = float(settings["run"]["residual_criterion"])
            row["converged"] = worst < row["criterion"]
            #  The SHAPE of the residual curve is the diagnosis - flat means
            #  the solver has nothing to converge to, still-descending means it
            #  just wanted more iterations - and it lived only in the server
            #  process's memory.  Restart the server and the one thing needed
            #  to tell those apart was gone.  It goes in the record now.
            row["residual_history"] = _thin(res, 150)
        mon = list(getattr(job.driver, "monitors", []) or [])
        row["dp_history"] = _thin(mon, 150)
        if len(mon) >= 2:
            #  how much the answer was still moving over the last fifth of the
            #  run.  Residuals settling is not the answer settling, and a study
            #  that recorded only the final number could not tell them apart.
            tail = mon[max(0, len(mon) - max(2, len(mon) // 5)):]
            vals = [m["dp"] for m in tail if m.get("dp") is not None]
            if len(vals) >= 2 and vals[-1]:
                row["dp_drift"] = (max(vals) - min(vals)) / abs(vals[-1])
        row["ok"] = True
        #  deliberately NOT closed here.  The next case's new_job() closes the
        #  previous session before it launches its own, so a licence is still
        #  held by one case at a time - but the LAST case's session is left
        #  alive, and that is the one somebody wants to look at in the Results
        #  tab when the campaign stops.  Closing it here also meant every
        #  session was closed twice, once by this line and once by new_job.
        return row


# =============================================================================
#  PLOTS
# =============================================================================
#  SVG, written out here rather than drawn on a canvas in the browser: a study
#  report is read on paper as often as on screen, an SVG prints at the
#  printer's resolution instead of the screen's, and the same function serves
#  the tab and the command line.  It is about eighty lines; a plotting
#  dependency for that would be a poor trade in a project that ships a stdlib
#  server on purpose.
PALETTE = ["#2563eb", "#dc2626", "#059669", "#d97706", "#7c3aed", "#0891b2"]


def _ticks(lo, hi, log):
    if log:
        out = []
        d = int(math.floor(math.log10(lo)))
        while 10.0 ** d <= hi * 1.0001:
            for m in (1, 2, 5):
                v = m * 10.0 ** d
                if lo * 0.9999 <= v <= hi * 1.0001:
                    out.append(v)
            d += 1
        return out or [lo, hi]
    span = hi - lo
    if span <= 0:
        return [lo]
    step = 10.0 ** math.floor(math.log10(span))
    for m in (1, 2, 5, 10):
        if span / (m * step) <= 6:
            step *= m
            break
    out, v = [], math.ceil(lo / step) * step
    while v <= hi * 1.0001:
        out.append(v)
        v += step
    return out


def _tlabel(v, log):
    if log:
        e = math.log10(v)
        if abs(e - round(e)) < 1e-9:
            return "1e%d" % round(e)
        return ("%g" % v) if v < 1e4 else "%.0e" % v
    return "%g" % float("%.4g" % v)


def svg_plot(series, xlabel="", ylabel="", xlog=False, ylog=False,
             width=620, height=330, title=""):
    """One chart.  `series` are {label, points:[(x,y)], dash, marker}."""
    pts = [p for s in series for p in s["points"]]
    pts = [(x, y) for x, y in pts if x is not None and y is not None
           and (not xlog or x > 0) and (not ylog or y > 0)]
    if not pts:
        return ""
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    if xlog:
        x0, x1 = x0 / 1.25, x1 * 1.25
    else:
        pad = (x1 - x0) * 0.06 or max(abs(x0), 1.0) * 0.06
        x0, x1 = x0 - pad, x1 + pad
    if ylog:
        y0, y1 = y0 / 1.3, y1 * 1.3
    else:
        pad = (y1 - y0) * 0.10 or max(abs(y0), 1.0) * 0.10
        y0, y1 = y0 - pad, y1 + pad
    #  the legend sits ABOVE the frame, not inside it: in a log-log Eu-Re plot
    #  the curves run straight through the top-left corner, which is the only
    #  place a boxed legend fits
    n_leg = len([s for s in series if s.get("label")])
    L, R, T, B = 66, 14, 16 + (16 if title else 0) + (14 if n_leg else 0), 44
    iw, ih = width - L - R, height - T - B

    def px(v):
        f = ((math.log10(v) - math.log10(x0)) / (math.log10(x1) - math.log10(x0))
             if xlog else (v - x0) / (x1 - x0))
        return L + f * iw

    def py(v):
        f = ((math.log10(v) - math.log10(y0)) / (math.log10(y1) - math.log10(y0))
             if ylog else (v - y0) / (y1 - y0))
        return T + ih - f * ih

    o = ['<svg class="plot" viewBox="0 0 %d %d" width="100%%" '
         'xmlns="http://www.w3.org/2000/svg" font-family="inherit">' % (width, height)]
    if title:
        o.append('<text x="%d" y="14" font-size="11.5" font-weight="640" '
                 'fill="#33415a">%s</text>' % (L, title))
    for v in _ticks(x0, x1, xlog):
        x = px(v)
        o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e8edf4"/>'
                 % (x, T, x, T + ih))
        o.append('<text x="%.1f" y="%d" font-size="10" fill="#5a6678" '
                 'text-anchor="middle">%s</text>'
                 % (x, T + ih + 15, _tlabel(v, xlog)))
    for v in _ticks(y0, y1, ylog):
        y = py(v)
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e8edf4"/>'
                 % (L, y, L + iw, y))
        o.append('<text x="%d" y="%.1f" font-size="10" fill="#5a6678" '
                 'text-anchor="end">%s</text>' % (L - 6, y + 3, _tlabel(v, ylog)))
    o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" '
             'stroke="#c8d2df"/>' % (L, T, iw, ih))
    o.append('<text x="%.1f" y="%d" font-size="10.5" fill="#33415a" '
             'text-anchor="middle">%s</text>'
             % (L + iw / 2.0, height - 6, ylabel and xlabel or xlabel))
    o.append('<text x="12" y="%.1f" font-size="10.5" fill="#33415a" '
             'text-anchor="middle" transform="rotate(-90 12 %.1f)">%s</text>'
             % (T + ih / 2.0, T + ih / 2.0, ylabel))
    for i, s in enumerate(series):
        col = s.get("color") or PALETTE[i % len(PALETTE)]
        ps = [(x, y) for x, y in s["points"] if x is not None and y is not None
              and (not xlog or x > 0) and (not ylog or y > 0)]
        if not ps:
            continue
        if s.get("line", True) and len(ps) > 1:
            d = " ".join(("%s%.1f,%.1f" % ("M" if j == 0 else "L", px(x), py(y)))
                         for j, (x, y) in enumerate(ps))
            o.append('<path d="%s" fill="none" stroke="%s" stroke-width="1.7"%s/>'
                     % (d, col, ' stroke-dasharray="5 3"' if s.get("dash") else ""))
        if s.get("marker", True):
            for x, y in ps:
                o.append('<circle cx="%.1f" cy="%.1f" r="2.9" fill="%s"/>'
                         % (px(x), py(y), col))
    #  legend: one row, above the frame
    lx, ly = L, T - 6
    for i, s in enumerate(series):
        if not s.get("label"):
            continue
        col = s.get("color") or PALETTE[i % len(PALETTE)]
        o.append('<rect x="%.1f" y="%.1f" width="15" height="2.6" fill="%s"/>'
                 % (lx, ly - 3, col))
        o.append('<text x="%.1f" y="%.1f" font-size="10" fill="#33415a">%s</text>'
                 % (lx + 20, ly + 1, s["label"]))
        lx += 26 + 6.0 * len(s["label"])
    o.append("</svg>")
    return "".join(o)


# =============================================================================
#  THE STAGE REPORTS
# =============================================================================
def _esc(v):
    return (str(v).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _n(v, fmt="%.4g"):
    return "—" if v is None else fmt % v


def _warnings(study, rows, ko):
    """What a reader must be told before reading any number below."""
    out = []
    all_rows = study.results()
    mock = [r for r in all_rows if r.get("mock")]
    bad = [r for r in all_rows if not r.get("ok")]
    missing = [c["id"] for c in study.cases
               if c["id"] not in {r["id"] for r in all_rows}]
    if mock:
        out.append(("bad", "MOCK 백엔드로 실행된 케이스가 %d개 있습니다. 이 값들은 "
                           "결과가 아니며 아래 표·그림·적합에서 모두 제외되었습니다."
                    % len(mock) if ko else
                    "%d case(s) ran on the MOCK backend. Those are not results "
                    "and are excluded from every table, plot and fit below."
                    % len(mock)))
    if bad:
        out.append(("bad", "실패한 케이스 %d개: %s" % (len(bad),
                    ", ".join(r["id"] for r in bad[:8])) if ko else
                    "%d case(s) failed: %s"
                    % (len(bad), ", ".join(r["id"] for r in bad[:8]))))
    if missing:
        out.append(("warn", "아직 실행되지 않은 케이스 %d개." % len(missing) if ko
                    else "%d case(s) have not been run yet." % len(missing)))
    unconv = [r for r in rows if r.get("converged") is False]
    if unconv:
        eqs = sorted({r.get("residual_worst_eq") for r in unconv
                      if r.get("residual_worst_eq")})
        worst = max(r.get("residual_worst") or 0.0 for r in unconv)
        out.append(("bad",
                    "잔차 기준에 도달하지 못한 케이스 %d개 (%s). 가장 큰 잔차 "
                    "%.2e, 주로 %s 방정식. 수렴하지 않은 해로는 격자 수렴 지수를 "
                    "만들 수 없으므로 아래 외삽과 판정에서 제외했습니다 — 표에는 "
                    "보이도록 남겨 두었습니다."
                    % (len(unconv), ", ".join(r["id"] for r in unconv[:8]),
                       worst, ", ".join(eqs) or "?")
                    if ko else
                    "%d case(s) did not reach the residual criterion (%s). "
                    "Worst residual %.2e, mostly in %s. A grid convergence "
                    "index cannot be built on unconverged solutions, so they "
                    "are excluded from every extrapolation below - and left in "
                    "the table so that they are visible."
                    % (len(unconv), ", ".join(r["id"] for r in unconv[:8]),
                       worst, ", ".join(eqs) or "?")))
    drift = [r for r in rows if (r.get("dp_drift") or 0) > 0.01]
    if drift:
        out.append(("warn", "마지막 구간에서 Δp가 1 %% 이상 움직인 케이스 %d개 - "
                            "잔차가 내려가도 답은 아직 움직이고 있었습니다."
                    % len(drift) if ko else
                    "%d case(s) had Dp still moving by more than 1 %% over the "
                    "final stretch - residuals settling is not the answer "
                    "settling." % len(drift)))
    return out


def _head(study, rows, ko, subtitle=""):
    t = (lambda a, b: a if ko else b)
    h = ['<h1>%s</h1>' % _esc(study.d.get("title") or study.name)]
    h.append('<div class="sub"><code>%s</code> · %s · %s</div>'
             % (_esc(study.name), _esc(study.d["geometry"]),
                _esc(subtitle or time.strftime("%Y-%m-%d %H:%M"))))
    for level, msg in _warnings(study, rows, ko):
        h.append('<p class="%s">%s</p>'
                 % ("bad" if level == "bad" else "warn", _esc(msg)))
    if study.d.get("notes"):
        h.append("<h2>%s</h2>" % t("이 연구의 설정", "How it was set up"))
        h.append("<ul>")
        for nline in study.d["notes"]:
            h.append("<li>%s</li>" % _esc(nline))
        h.append("</ul>")
    return h


def mesh_report_body(study, lang="ko", tol=0.01):
    ko = (lang == "ko")
    t = (lambda a, b: a if ko else b)
    rows = study.usable()
    an = mesh_analysis(study, tol=tol)
    h = _head(study, rows, ko)
    if not rows:
        h.append("<p>%s</p>" % t("아직 사용할 수 있는 결과가 없습니다.",
                                 "No usable results yet."))
        return "\n".join(h)

    for lad in an["ladders"]:
        levels = lad["levels"]
        meta = (levels[0].get("meta") or {})
        h.append("<h2>%s · %s</h2>" % (
            t("격자 사다리", "ladder"), _esc(lad["ladder"])))
        h.append("<p class=\"muted\">%s X = %s, Re<sub>max</sub> = %s</p>"
                 % (t("조건:", "condition:"), _n(meta.get("X")),
                    _n(meta.get("Re_target"), "%.0e")))
        if lad.get("blocked"):
            h.append('<p class="bad">%s</p>'
                     % _esc(lad["blocked_ko"] if ko else lad["blocked"]))
        h.append('<table class="grid"><tr><th>%s</th><th>%s</th><th>h [m]</th>'
                 '<th>Δp<sub>bundle</sub> [Pa]</th><th>Eu<sub>row</sub></th>'
                 '<th>%s</th><th>%s</th><th>%s</th><th>%s</th></tr>'
                 % (t("레벨", "level"), t("셀 수", "cells"),
                    t("외삽값 대비", "vs extrapolated"),
                    t("최종 잔차", "final residual"),
                    t("수렴", "converged"), t("시간", "time")))
        for r in sorted(levels, key=lambda r: -r["h"]):
            dev = r.get("dev_from_ext")
            bad = r.get("converged") is False
            resid = ("—" if r.get("residual_worst") is None else
                     "%.2e%s" % (r["residual_worst"],
                                 " " + r["residual_worst_eq"]
                                 if r.get("residual_worst_eq") else ""))
            h.append('<tr><td>%s%s</td><td class="n">%s</td><td class="n">%s</td>'
                     '<td class="n">%s</td><td class="n">%s</td>'
                     '<td class="n">%s</td><td class="n">%s</td>'
                     '<td class="n"%s>%s</td><td class="n">%s s</td></tr>'
                     % (_esc(r["id"]),
                        ' <b>&larr;</b>' if r["id"] == lad.get("chosen") else "",
                        "{:,}".format(r["cells"]), _n(r.get("h"), "%.3e"),
                        _n(r.get("dp_bundle")), _n(r.get("eu_row"), "%.4f"),
                        "—" if dev is None else "%.2f %%" % (100 * dev),
                        resid,
                        ' style="color:#9b1c1c;font-weight:640"' if bad else "",
                        t("아니오", "no") if bad else t("예", "yes"),
                        _n(r.get("seconds"), "%.0f")))
        h.append("</table>")

        if lad.get("triplets"):
            h.append("<h3>%s</h3>" % t("GCI (Celik 등 2008 / ASME V&amp;V 20)",
                                       "GCI (Celik et al. 2008 / ASME V&amp;V 20)"))
            h.append('<table class="grid"><tr><th>%s</th><th>r<sub>21</sub></th>'
                     '<th>r<sub>32</sub></th><th>p</th>'
                     '<th>&phi;<sub>ext</sub></th><th>e<sub>a</sub><sup>21</sup></th>'
                     '<th>GCI<sub>fine</sub></th></tr>'
                     % t("격자 3개 (조밀→성김)", "triplet (fine to coarse)"))
            for g in lad["triplets"]:
                if g.get("error"):
                    h.append('<tr><td>%s</td><td colspan="6">%s</td></tr>'
                             % (_esc(", ".join(g.get("ids") or [])),
                                _esc(g["error"])))
                    continue
                h.append('<tr><td>%s</td><td class="n">%.3f</td>'
                         '<td class="n">%.3f</td><td class="n">%s%s</td>'
                         '<td class="n">%s</td><td class="n">%s</td>'
                         '<td class="n">%s</td></tr>'
                         % (_esc(" → ".join(g["ids"])), g["r21"], g["r32"],
                            _n(g.get("p"), "%.2f"),
                            " ⚠" if g.get("oscillatory") else "",
                            _n(g.get("phi_ext")),
                            "—" if g.get("e_a21") is None
                            else "%.2f %%" % (100 * g["e_a21"]),
                            "—" if g.get("gci_fine") is None
                            else "%.2f %%" % (100 * g["gci_fine"])))
            h.append("</table>")
            thin = [g for g in lad["triplets"] if (g.get("r21") or 9) < 1.3
                    or (g.get("r32") or 9) < 1.3]
            if thin:
                h.append('<p class="warn">%s</p>' % t(
                    "세밀화 비가 1.3 미만인 삼중항이 있습니다 (Celik 등의 권고값). "
                    "블록 격자는 분할 수가 정수라 비를 정확히 맞출 수 없습니다. "
                    "GCI 절차 자체는 비가 서로 달라도 성립하지만, 비가 작을수록 "
                    "관측 차수 p 가 불안정해집니다.",
                    "Some triplets refine by less than the 1.3 Celik et al. "
                    "recommend. A block mesh cannot hit a ratio exactly - the "
                    "division counts are integers. The procedure itself is "
                    "valid for unequal ratios, but the smaller the ratio the "
                    "less stable the observed order p."))
            if any(g.get("oscillatory") for g in lad["triplets"]):
                h.append('<p class="warn">%s</p>' % t(
                    "⚠ 로 표시된 삼중항은 진동 수렴입니다: 세 점이 단조 곡선 위에 "
                    "있지 않으므로 그 외삽은 의미가 없습니다.",
                    "The triplets marked ⚠ converge non-monotonically: the "
                    "three points are not on a monotone curve, so that "
                    "extrapolation does not mean anything."))

        if lad.get("phi_ext"):
            pts = [(r["h"], r["dp_bundle"]) for r in
                   sorted(levels, key=lambda r: r["h"])]
            hs = [p[0] for p in pts]
            h.append(svg_plot(
                [{"label": t("CFD", "CFD"), "points": pts},
                 {"label": t("Richardson 외삽 (h→0)", "Richardson h -> 0"),
                  "points": [(min(hs) * 0.7, lad["phi_ext"]),
                             (max(hs) * 1.1, lad["phi_ext"])],
                  "dash": True, "marker": False, "color": "#64748b"}],
                xlabel="h [m]", ylabel="Δp_bundle [Pa]",
                title=t("격자 수렴", "grid convergence")))
            h.append('<p class="muted">%s</p>' % t(
                "가로축은 대표 셀 크기입니다. 위·아래가 대칭면이고 z 방향 셀 수가 "
                "고정이므로 이 세밀화는 2차원이고, h 는 셀 하나의 면적의 제곱근입니다. "
                "(V/N)<sup>1/3</sup> 을 쓰면 각 단계가 실제보다 덜 세밀해진 것처럼 "
                "보여 관측 차수가 틀리게 나옵니다.",
                "The abscissa is the representative cell size. Top and bottom "
                "are symmetry planes and the z count is fixed, so this "
                "refinement is two-dimensional and h is the square root of a "
                "cell's area. Using (V/N)^(1/3) would make every level look "
                "less refined than it is and the observed order would come "
                "out wrong."))

        #  the residual curve of the finest case, because its SHAPE is what
        #  says whether a stalled run wanted more iterations or had nothing to
        #  converge to.  Drawn for the finest mesh: if that one is flat, no
        #  coarser one is going to be better.
        finest = min(levels, key=lambda r: r["h"])
        hist = finest.get("residual_history") or []
        if hist:
            eqs = sorted({k for row in hist for k in row
                          if k != "iter" and row.get(k)})
            series = [{"label": eq,
                       "points": [(r["iter"], r[eq]) for r in hist if r.get(eq)]}
                      for eq in eqs]
            crit = finest.get("criterion")
            if crit:
                xs = [r["iter"] for r in hist]
                series.append({"label": t("판정 기준", "criterion"),
                               "points": [(min(xs), crit), (max(xs), crit)],
                               "dash": True, "marker": False, "color": "#94a3b8"})
            h.append(svg_plot(series, xlabel=t("반복", "iteration"),
                              ylabel=t("잔차", "residual"), ylog=True,
                              title="%s · %s" % (_esc(finest["id"]),
                                                 t("잔차", "residuals"))))
            h.append('<p class="muted">%s</p>' % t(
                "곡선이 <b>평평</b>하면 정상상태 솔버가 수렴할 해가 없다는 뜻이고 "
                "(뭉툭한 물체 다발에서는 보통 유동이 비정상이라는 뜻입니다), "
                "<b>끝까지 내려가고</b> 있으면 반복 횟수가 모자랐다는 뜻입니다. "
                "둘은 다른 문제이고 고치는 방법도 다릅니다.",
                "A <b>flat</b> curve means the steady solver has nothing to "
                "converge to - on a bluff-body bank that usually means the "
                "flow is unsteady. A curve <b>still descending</b> at the end "
                "means it simply wanted more iterations. They are different "
                "problems with different fixes."))
        dph = finest.get("dp_history") or []
        if len(dph) >= 2:
            h.append(svg_plot(
                [{"label": "\u0394p", "points": [(r["iter"], r["dp"])
                                                  for r in dph if r.get("dp")]}],
                xlabel=t("반복", "iteration"), ylabel="\u0394p [Pa]",
                title="%s · %s" % (_esc(finest["id"]),
                                   t("압력강하 이력", "pressure drop history"))))

        chosen = lad.get("chosen")
        if chosen:
            r = [x for x in levels if x["id"] == chosen][0]
            h.append('<div class="tiles">')
            h.append('<div class="tile"><div class="k">%s</div>'
                     '<div class="v">%s</div></div>'
                     % (t("선택된 격자", "the mesh to use"), _esc(chosen)))
            h.append('<div class="tile"><div class="k">%s</div>'
                     '<div class="v">%s</div></div>'
                     % (t("셀 수", "cells"), "{:,}".format(r["cells"])))
            h.append('<div class="tile"><div class="k">%s</div>'
                     '<div class="v">%s</div></div>'
                     % (t("외삽값 대비", "from extrapolated"),
                        "—" if r.get("dev_from_ext") is None
                        else "%.2f %%" % (100 * r["dev_from_ext"])))
            h.append("</div>")
            h.append("<p>%s</p>" % t(
                "격자 연구의 목적은 가장 세밀한 격자가 아니라 <b>충분히 정확한 가장 "
                "싼 격자</b>를 찾는 것입니다. 위는 Richardson 외삽값에서 %.1f %% 이내에 "
                "드는 가장 성긴 격자입니다." % (100 * tol),
                "The point of a mesh study is not the finest mesh but the "
                "<b>cheapest adequate</b> one. That is the coarsest mesh "
                "within %.1f %% of the Richardson extrapolation."
                % (100 * tol)))
            if lad.get("chosen_note"):
                h.append('<p class="warn">%s</p>' % _esc(lad["chosen_note"]))
    h.append('<div class="foot">study.py · %s</div>'
             % _esc(time.strftime("%Y-%m-%d %H:%M:%S")))
    return "\n".join(h)


def sweep_report_body(study, lang="ko"):
    ko = (lang == "ko")
    t = (lambda a, b: a if ko else b)
    rows = study.usable()
    an = sweep_analysis(study)
    h = _head(study, rows, ko)
    if not an["rows"]:
        h.append("<p>%s</p>" % t("아직 사용할 수 있는 결과가 없습니다.",
                                 "No usable results yet."))
        return "\n".join(h)
    keys = an["keys"]
    names = {k: CR.BY_KEY[k]["name"] for k in keys}

    #  headline: how far the CFD is from each correlation
    h.append("<h2>%s</h2>" % t("상관식과의 비교", "Against the correlations"))
    band = an.get("band")
    if band:
        h.append("<p>%s</p>" % t(
            "1단계에서 측정한 대로 기존 상관식들은 서로 1.0–2.2배 차이가 납니다. "
            "따라서 '어느 한 상관식의 몇 % 이내'는 판정 기준이 될 수 없습니다. "
            "여기서 쓰는 기준은 <b>그 조건에서 적용 범위 안에 있는 상관식들이 "
            "만드는 띠 안에 CFD가 들어오는가</b> 입니다.",
            "As stage 1 measured, the published correlations differ from one "
            "another by 1.0x to 2.2x. 'Within N % of one of them' is therefore "
            "not a test. The test used here is whether the CFD falls "
            "<b>inside the band spanned by the correlations that are in range "
            "at that condition</b>."))
    h.append('<div class="tiles">')
    h.append('<div class="tile"><div class="k">%s</div><div class="v">%d</div></div>'
             % (t("사용된 케이스", "cases used"), len(an["rows"])))
    if band:
        h.append('<div class="tile"><div class="k">%s</div>'
                 '<div class="v">%d<span class="u">/ %d</span></div></div>'
                 % (t("띠 안에 든 케이스", "inside the band"),
                    band["inside"], band["n"]))
        h.append('<div class="tile"><div class="k">%s</div>'
                 '<div class="v">%.0f<span class="u">%%</span></div></div>'
                 % (t("띠 밖 최대 이탈", "worst excursion"),
                    100 * band["worst"]))
        h.append('<div class="tile"><div class="k">%s</div>'
                 '<div class="v">%.2f<span class="u">×</span></div></div>'
                 % (t("띠의 평균 폭", "mean band width"), band["width_mean"]))
    for k in keys:
        s = an["ratios"].get(k)
        if not s:
            continue
        h.append('<div class="tile"><div class="k">%s</div>'
                 '<div class="v">%.2f<span class="u">×</span></div></div>'
                 % (_esc(names[k]) + " " + t("평균비", "mean ratio"), s["mean"]))
    h.append("</div>")
    h.append('<table class="grid"><tr><th>%s</th><th>n</th><th>%s</th><th>%s</th>'
             '<th>%s</th><th>%s</th></tr>'
             % (t("상관식", "correlation"), t("최소", "min"), t("최대", "max"),
                t("평균", "mean"), t("RMS 편차", "rms deviation from 1")))
    for k in keys:
        s = an["ratios"].get(k)
        if not s:
            continue
        h.append('<tr><th>%s</th><td class="n">%d</td><td class="n">%.2f</td>'
                 '<td class="n">%.2f</td><td class="n">%.2f</td>'
                 '<td class="n">%.1f %%</td></tr>'
                 % (_esc(names[k]), s["n"], s["min"], s["max"], s["mean"],
                    100 * s["rms"]))
    h.append("</table>")
    h.append('<p class="muted">%s</p>' % t(
        "비 = CFD Eu<sub>row</sub> / 상관식 Eu<sub>row</sub>. 1보다 크면 CFD가 "
        "더 큰 압력강하를 냅니다.",
        "ratio = CFD Eu_row / correlation Eu_row. Above 1 means the CFD drops "
        "more pressure than the correlation predicts."))

    #  Eu against Re, one curve per pitch pair
    h.append("<h2>%s</h2>" % t("Eu – Re 곡선", "Eu against Re"))
    by_pitch = {}
    for r in an["rows"]:
        by_pitch.setdefault((round(r["XT"], 3), round(r["XL"], 3)), []).append(r)
    for i, key in enumerate(sorted(by_pitch)):
        rs = sorted(by_pitch[key], key=lambda r: r["Re"])
        if len(rs) < 2:
            continue
        series = [{"label": "CFD", "points": [(r["Re"], r["eu_row"]) for r in rs],
                   "color": PALETTE[0]}]
        for j, k in enumerate(keys):
            pts = [(r["Re"], r["corr"][k]["eu_row"]) for r in rs
                   if r["corr"][k]["eu_row"]]
            if pts:
                series.append({"label": names[k], "points": pts, "dash": True,
                               "color": PALETTE[(j + 1) % len(PALETTE)]})
        h.append(svg_plot(series, xlabel="Re_max", ylabel="Eu per row",
                          xlog=True, ylog=True,
                          title="X_T = %g, X_L = %g" % key))

    #  parity
    h.append("<h2>%s</h2>" % t("일대일 비교", "Parity"))
    series = []
    for j, k in enumerate(keys):
        pts = [(r["corr"][k]["eu_row"], r["eu_row"]) for r in an["rows"]
               if r["corr"][k]["eu_row"]]
        if pts:
            series.append({"label": names[k], "points": pts, "line": False,
                           "color": PALETTE[(j + 1) % len(PALETTE)]})
    if series:
        allv = [v for s in series for p in s["points"] for v in p]
        lo, hi = min(allv), max(allv)
        series.insert(0, {"label": "1:1", "points": [(lo, lo), (hi, hi)],
                          "marker": False, "color": "#94a3b8"})
        h.append(svg_plot(series, xlabel=t("상관식 Eu", "correlation Eu"),
                          ylabel=t("CFD Eu", "CFD Eu"), xlog=True, ylog=True))

    #  the matrix itself
    h.append("<h2>%s</h2>" % t("전체 결과", "Every case"))
    h.append('<table class="grid"><tr><th>%s</th><th>X<sub>T</sub></th>'
             '<th>X<sub>L</sub></th>'
             '<th>Re<sub>max</sub></th><th>Eu CFD</th><th>%s</th>'
             % (t("케이스", "case"), t("띠", "band"))
             + "".join("<th>%s</th><th>%s</th>" % (_esc(names[k]), t("비", "ratio"))
                       for k in keys) + "</tr>")
    for r in sorted(an["rows"], key=lambda r: (r["XT"], r["XL"], r["Re"])):
        cells = ""
        for k in keys:
            c = r["corr"][k]
            #  a dagger where the correlation was asked outside its own limits:
            #  that number is an extrapolation and its ratio is not evidence
            mark = " †" if c["outside"] else ""
            cells += ('<td class="n">%s%s</td><td class="n">%s</td>'
                      % (_n(c["eu_row"], "%.4f"), mark, _n(c["ratio"], "%.2f")))
        verdict = ("—" if "in_band" not in r else
                   (t("안", "in") if r["in_band"]
                    else "%+.0f %%" % (100 * r["band_dev"])))
        h.append('<tr><td>%s</td><td class="n">%.2f</td><td class="n">%.2f</td>'
                 '<td class="n">%.3g</td><td class="n">%.4f</td>'
                 '<td class="n">%s</td>%s</tr>'
                 % (_esc(r["id"]), r["XT"], r["XL"], r["Re"], r["eu_row"],
                    verdict, cells))
    h.append("</table>")
    h.append('<p class="muted">%s</p>' % t(
        "† 는 그 상관식이 자기 적용 범위 밖에서 계산되었다는 뜻이고, 그 값과 비는 "
        "외삽입니다. 띠 계산에서는 제외됩니다.",
        "A dagger means that correlation was evaluated outside its own quoted "
        "limits; that value and its ratio are an extrapolation, and it is left "
        "out of the band."))

    #  our own fit
    h.append("<h2>%s</h2>" % t("우리 상관식 (4단계 준비)",
                               "A correlation of our own - toward stage 4"))
    fit = an.get("fit")
    if not fit:
        h.append("<p>%s %s</p>" % (
            t("아직 적합할 수 없습니다.", "Not fittable yet."),
            _esc(an.get("fit_error", ""))))
    else:
        h.append("<p><code>%s</code></p>" % CR.FIT_FORM["expr"])
        h.append("<p class=\"muted\">%s</p>" % t(
            "형태는 새로 만든 것이 아니라 Shen 등(2024) 식 (12) 의 골격이고, "
            "우리 데이터로 계수만 다시 맞춘 것입니다. 오른쪽 열이 그 논문의 값입니다.",
            "The shape is not new: it is the skeleton of Shen et al. (2024) "
            "eq. (12), refitted on our own data. The right-hand column is that "
            "paper's own value."))
        h.append('<table class="grid"><tr><th>%s</th><th>%s</th><th>%s</th>'
                 "<th>%s</th></tr>"
                 % (t("계수", "coefficient"), t("우리 적합", "our fit"),
                    "Shen 2024", t("비", "ratio")))
        for cname in CR.FIT_FORM["coefficients"]:
            ref = CR.FIT_FORM["shen_values"].get(cname)
            held = cname in (fit.get("fixed") or [])
            h.append('<tr><th>%s</th><td class="n">%.5g%s</td>'
                     '<td class="n">%s</td><td class="n">%s</td></tr>'
                     % (cname, fit[cname],
                        (" <i>(%s)</i>" % t("고정", "held")) if held else "",
                        _n(ref, "%.5g"),
                        "—" if not ref else "%.2f" % (fit[cname] / ref)))
        h.append('<tr><th>%s</th><td colspan="3" class="n">%d</td></tr>'
                 % (t("사용 점 수", "points"), fit["n_points"]))
        h.append('<tr><th>%s</th><td colspan="3" class="n">%.2f %%</td></tr>'
                 % (t("평균 절대 오차", "mean absolute error"),
                    100 * fit["mean_abs_error"]))
        h.append('<tr><th>%s</th><td colspan="3" class="n">%.2f %%</td></tr>'
                 % (t("최대 절대 오차", "worst absolute error"),
                    100 * fit["max_abs_error"]))
        h.append("</table>")
        local = (lambda k: fit.get(k + "_ko", fit.get(k)) if ko
                 else fit.get(k))
        if fit.get("fixed_note"):
            h.append('<p class="muted">%s</p>' % _esc(local("fixed_note")))
        if fit.get("laminar_dropped"):
            h.append('<p class="warn">%s</p>' % _esc(local("laminar_dropped")))

        tied = an.get("fit_tied")
        if tied:
            h.append("<h3>%s</h3>" % t("피치 지수를 나눌 값이 있는가",
                                       "Is splitting the pitch exponents worth it?"))
            h.append("<p>%s</p>" % t(
                "Shen 식은 두 피치비를 곱으로만 씁니다 (X_T X_L)^-p. 그 논문의 "
                "해석은 전부 S_T = S_L 이었으므로 두 지수를 나눌 데이터가 "
                "없었습니다 — 우리 행렬은 일부러 그 공선성을 피했으므로, 나눌 "
                "값이 있는지 여기서 답할 수 있습니다.",
                "Shen's form uses the two pitch ratios only as a product, "
                "(X_T X_L)^-p. Every case in that paper had S_T = S_L, so no "
                "data of theirs could split the exponents. Ours varies them "
                "independently on purpose, so the question can be answered "
                "here."))
            h.append('<table class="grid"><tr><th>%s</th><th>p</th><th>q</th>'
                     "<th>%s</th><th>%s</th></tr>"
                     % (t("적합", "fit"), t("평균 오차", "mean error"),
                        t("최대 오차", "worst")))
            for label, fc in ((t("p, q 각각", "p and q free"), fit),
                              (t("p = q (Shen)", "p = q (Shen's)"), tied)):
                h.append('<tr><th>%s</th><td class="n">%.3f</td>'
                         '<td class="n">%.3f</td><td class="n">%.2f %%</td>'
                         '<td class="n">%.2f %%</td></tr>'
                         % (label, fc["p"], fc["q"],
                            100 * fc["mean_abs_error"],
                            100 * fc["max_abs_error"]))
            h.append("</table>")
            better = tied["mean_abs_error"] / max(fit["mean_abs_error"], 1e-12)
            h.append("<p>%s</p>" % t(
                "나누면 평균 오차가 %.1f 배 줄어듭니다. p 가 q 보다 크다는 것은 "
                "<b>횡방향 피치가 압력강하를 지배하고 종방향 피치는 거의 영향이 "
                "없다</b>는 뜻이며, u_max 가 횡방향 간극으로 정해진다는 사실과 "
                "맞습니다." % better,
                "Splitting cuts the mean error by a factor of %.1f. That p "
                "comes out larger than q says the <b>transverse pitch governs "
                "the pressure drop and the longitudinal one barely enters</b>, "
                "which is consistent with u_max being set by the transverse "
                "gap." % better))
        h.append('<p class="muted">%s</p>' % t(
            "이것은 <b>직선 rod</b> 데이터만으로 맞춘 것이며, 나선 코일 항 r 은 아직 "
            "포함되어 있지 않습니다. 4단계에서 코일 데이터가 더해지면 같은 형태에 "
            "항 하나가 붙고, 경사각 0에서 이 계수들로 정확히 돌아가야 합니다 — 그것이 "
            "코일 적합이 맞는지 확인하는 방법입니다.",
            "This is fitted on the STRAIGHT-ROD data alone; the coil term r is "
            "not in it. Stage 4 adds one factor to the same form, and at zero "
            "lean it must collapse back to exactly these coefficients - which "
            "is how the coil fit gets checked."))
        #  parity, not Eu against Re: with nine pitch pairs on one Re axis the
        #  points stack in vertical columns and nothing about the FIT is
        #  visible.  Fitted against measured, with the 1:1 line, shows it.
        pts = [(CR.fit_eu(fit, {"Re": r["Re"], "XT": r["XT"], "XL": r["XL"]}),
                r["eu_row"]) for r in an["rows"]]
        allv = [v for p in pts for v in p]
        lo, hi = min(allv), max(allv)
        h.append(svg_plot(
            [{"label": "1:1", "points": [(lo, lo), (hi, hi)], "marker": False,
              "color": "#94a3b8"},
             {"label": t("적합 대 CFD", "fit against CFD"), "points": pts,
              "line": False, "color": PALETTE[1]}],
            xlabel=t("적합 Eu", "fitted Eu"), ylabel=t("CFD Eu", "CFD Eu"),
            xlog=True, ylog=True,
            title=t("우리 상관식의 적합도", "how well our form fits")))
    h.append('<div class="foot">study.py · %s</div>'
             % _esc(time.strftime("%Y-%m-%d %H:%M:%S")))
    return "\n".join(h)


def report_body(study, lang="ko"):
    return (mesh_report_body(study, lang) if study.kind == "mesh"
            else sweep_report_body(study, lang))


def write_report(study, lang="ko", path=None):
    body = report_body(study, lang)
    path = path or os.path.join(study.dir, "report.html")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(CR.standalone(body, study.d.get("title") or study.name, lang))
    return path


# =============================================================================
#  COMMAND LINE
# =============================================================================
STRAIGHT = ("rod-inline", "rod-staggered")


def make_campaign(level="L3"):
    """The straight-rod campaign: a mesh study and a sweep for each
    arrangement.  Helical coils are stage 4 and are deliberately not here."""
    out = []
    for geom in STRAIGHT:
        out.append(make_mesh_study(geom))
    for geom in STRAIGHT:
        out.append(make_sweep_study(geom, level=level))
    return out


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)

    def opt(flag, default=None):
        if flag in argv:
            i = argv.index(flag)
            if len(argv) > i + 1 and not argv[i + 1].startswith("--"):
                return argv[i + 1]
        return default

    if "--make" in argv:
        level = opt("--level", "L3")
        for s in make_campaign(level):
            path = s.save()
            print("wrote %s  (%d cases)" % (path, len(s.cases)))
        return 0
    if "--list" in argv:
        rows = Study.list_all()
        if not rows:
            print("no studies in %s - run --make" % STUDIES)
            return 0
        print("%-24s %-6s %-16s %5s %5s %6s %5s"
              % ("name", "kind", "geometry", "cases", "done", "failed", "mock"))
        for r in rows:
            print("%-24s %-6s %-16s %5d %5d %6d %5d"
                  % (r["name"], r["kind"], r["geometry"], r["cases"],
                     r["done"], r["failed"], r["mock"]))
        return 0
    name = opt("--report") or opt("--show")
    if "--report" in argv:
        if not name:
            print("--report needs a study name")
            return 2
        s = Study.load(name)
        lang = "en" if "--en" in argv else "ko"
        path = write_report(s, lang)
        print("wrote %s (%.1f kB)" % (path, os.path.getsize(path) / 1e3))
        return 0
    if "--show" in argv:
        s = Study.load(name)
        print("%s  (%s, %s, %d cases)"
              % (s.name, s.kind, s.d["geometry"], len(s.cases)))
        for c in s.cases:
            r = c.get("mesh_rule") or {}
            print("  %-20s %-34s az %-3s rad %-3s first %-9s Re %.3g"
                  % (c["id"], (c.get("label") or "")[:34],
                     c["params"].get("nAz"), r.get("n_rad"),
                     r.get("first_layer"), r.get("Re") or 0))
        return 0
    print(__doc__.strip().splitlines()[0])
    print("\n  python3 study.py --make [--level L3]   build the campaign")
    print("  python3 study.py --list                what is defined and how far it got")
    print("  python3 study.py --show NAME           the cases of one study")
    print("  python3 study.py --report NAME [--en]  write its report")
    print("\nRunning a study is done from the app's 파라메트릭 tab, because it "
          "needs Fluent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
