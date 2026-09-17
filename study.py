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
import re
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


#  The time step is DERIVED, not typed.  Same rule as the first cell height,
#  which follows from a y+ target rather than from somebody's judgement: a
#  transient bundle case has one time scale that matters, the shedding period,
#  and every other choice is a count of it.
#
#      f = St u_max / D,   T = 1/f,   dt = T / steps_per_period
#
#  St = 0.2 is the standard value for a circular cylinder over the whole
#  sub-critical range, and in a bank the passing frequency is close enough to
#  it for setting a time step.  Getting it wrong by 30 % costs 30 % of the run
#  time; getting it wrong by a factor of ten loses the oscillation entirely,
#  which is what typing a number would eventually do.
SHEDDING = {
    "strouhal": 0.2,
    "steps_per_period": 25,     # 25 steps resolves the cycle for 2nd-order
    "periods": 20,              # the shedding floor: never fewer than this
    "average_periods": 15,      # of which this many are averaged
    "flush_flowthroughs": 2.0,  # and the bundle is swept this many times first
    "max_iter_per_step": 12,
}


def transient_controls(st, cfg=None, bundle_length=None):
    """Time step, step count and averaging window for one flow state.

    TWO clocks, and the first version of this only had one.

    The time step comes from the shedding period, D / (St u_max): twenty-five
    steps resolve the cycle.  The RUN LENGTH was taken from the same clock -
    twenty periods, five of them nominally to "flush the start-up" - and that
    is the wrong clock for FLUSHING.  A shedding period is set by one rod and
    the gap velocity; the time the bundle takes to establish is set by its
    whole length and the bulk velocity, and the ratio between the two moves
    with the pitch and the Reynolds number rather than staying put.  Measured
    over the 76 cases this campaign defines, five shedding periods came to
    between 0.25 and 1.00 sweeps of the case's own bundle - every one of the
    76 began averaging before the flow had crossed the bundle once, so the
    last rows were still carrying the initial condition of the first, by
    differing amounts across a matrix whose whole purpose is comparison.

    The AVERAGING window stays on the shedding clock, and should: fifteen
    periods is a statistical window over a periodic signal, and how many
    bundle sweeps that happens to be does not change how well it averages.
    It is the flush, not the average, that has to be measured against the
    bundle.

    So the run is long enough for BOTH: sweep the bundle `flush_flowthroughs`
    times, then average over `average_periods` shedding periods, and never
    fewer steps than the shedding floor.  The flushing length is the BUNDLE,
    not the whole domain - the inlet box carries uniform flow and has nothing
    to develop, and the outlet box is downstream of both measuring planes.

    Costs a factor of about two in steps at the demanding end and nothing at
    the easy end, and makes the amount of physics behind each number the same
    across the matrix, which is what a sweep compares.
    """
    c = dict(SHEDDING)
    c.update(cfg or {})
    f = c["strouhal"] * st["umax"] / st["D"]
    period = 1.0 / f if f > 0 else 1.0
    dt = period / c["steps_per_period"]
    dt = float("%.4g" % dt)
    keep = int(round(c["average_periods"] * c["steps_per_period"]))
    floor = int(round(c["periods"] * c["steps_per_period"]))
    t_flow = None
    steps = floor
    if bundle_length and st.get("u_in"):
        t_flow = bundle_length / float(st["u_in"])
        steps = max(floor, int(math.ceil(
            (c["flush_flowthroughs"] * t_flow + keep * dt) / dt)))
    out = {
        "time_step": dt,
        "time_steps": steps,
        "average_last": min(keep, steps),
        "max_iter_per_step": int(c["max_iter_per_step"]),
        "shedding_hz": f,
        "period_s": period,
        "flow_time_s": steps * dt,
    }
    if t_flow:
        out["bundle_flowthrough_s"] = t_flow
        out["flowthroughs"] = steps * dt / t_flow
        #  what the ANSWER is averaged over, which is the number that matters
        out["flowthroughs_averaged"] = out["average_last"] * dt / t_flow
    return out


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
        #  several workers append to one result file, and record() is a
        #  read-modify-write of the whole thing
        self._lock = threading.Lock()
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
    def save(self, force=False):
        """Write the definition.  Refuses to overwrite one that has results.

        A definition and the results beside it are one object: every case in a
        study is meant to have been run the same way, and that is the only
        reason the cases can be compared with each other.  Rewriting the
        definition under existing results - which `--make` and the tab's
        rebuild both did, silently - leaves a folder whose study.json
        describes a run that did not happen, and whose remaining cases would
        be run differently from the ones already in it.  That is not a stale
        file; it is a study that has quietly become two studies.
        """
        if not force and os.path.isfile(self.results_path) and self.results():
            raise ValueError(
                "%s already has %d recorded case(s); rewriting its definition "
                "would describe a run that did not happen and would run the "
                "rest differently. Clear the results first, or save with force."
                % (self.name, len(self.results())))
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
        with self._lock:
            rows = [r for r in self.results() if r.get("id") != row["id"]]
            rows.append(row)
            order = {c["id"]: i for i, c in enumerate(self.cases)}
            rows.sort(key=lambda r: order.get(r.get("id"), 1e9))
            os.makedirs(self.dir, exist_ok=True)
            tmp = self.results_path + ".part"
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump({"study": self.name, "rows": rows}, fh,
                          ensure_ascii=False, indent=1)
                fh.write("\n")
            #  written aside and moved into place: a worker crashing mid-write
            #  would otherwise leave a truncated file where the campaign's
            #  whole record used to be
            os.replace(tmp, self.results_path)

    def clear_results(self):
        if os.path.isfile(self.results_path):
            os.unlink(self.results_path)

    def usable(self):
        """The rows a conclusion may be drawn from: finished, not mock."""
        return [r for r in self.results() if r.get("ok") and not r.get("mock")]


# =============================================================================
#  GRID CONVERGENCE
# =============================================================================
#  How much of each history the record keeps.
#
#  These two are not the same kind of data and must not share a budget.
#
#  The residual history is a DIAGNOSIS: what is wanted from it is the shape -
#  falling, flat, or falling and then flat - and a few hundred points draw
#  that as well as six thousand do.  A transient run reports every inner
#  iteration, so 500 time steps x 12 inner iterations is 6000 rows for one
#  case; thinning them costs nothing.
#
#  The Dp history is the RESULT.  A transient run samples it once per time
#  step, and the trace is a shedding oscillation resolved at 25 steps per
#  period.  Thin 500 samples to 150 and the same signal is left at 7.5
#  points per period - too few to read a frequency off honestly, and thinning
#  at a non-integer stride aliases the oscillation rather than merely
#  coarsening it.  So the cap is set above any run this campaign defines
#  (20 periods x 25 steps = 500), and nothing is dropped in practice.  The
#  cost is about 60 kB per case.
RES_KEEP = 400
DP_KEEP = 2000


def mean_drift(row):
    """Has the TIME AVERAGE stopped moving?  None if it cannot be told.

    Not the same question as how far the instantaneous Dp swings.  A settled
    vortex street swings ten per cent of its own mean every period and will do
    so for ever; a criterion built on that swing rejects converged answers,
    which is what the first version of this did.  What has to stop moving is
    the mean, so the averaging window is halved and the two half-means are
    compared.

    The measurement itself is fluent_case.mean_drift - one implementation, so
    what the driver records and what a report recomputes cannot disagree.
    Recomputed from dp_history when the recorded field is absent, so a case
    that ran before this existed is judged by the same rule as one that ran
    after rather than being re-run for a number already in its record.
    """
    import fluent_case as FC            # lazy: FC does not import this module
    if row.get("dp_mean_drift") is not None:
        return row["dp_mean_drift"]
    hist = row.get("dp_history") or []
    key = ("dp_bundle" if any(m.get("dp_bundle") is not None for m in hist)
           else "dp")
    vals = [m[key] for m in hist if m.get(key) is not None]
    keep = row.get("dp_averaged_over") or 0
    if keep and keep <= len(vals):
        vals = vals[-keep:]
    return FC.mean_drift(vals)


def _thin(rows, keep):
    """Every nth row, ends included, so a history fits in a result file."""
    n = len(rows)
    if n <= keep:
        return list(rows)
    step = (n - 1) / float(keep - 1)
    out = [rows[int(round(i * step))] for i in range(keep)]
    out[-1] = rows[-1]
    return out


#  Directions a mesh ladder can refine in, and where to read each one off a
#  case's parameters.  The wall-normal one is the first cell height, which is
#  a LENGTH and so refines when it gets smaller - hence the flag.
REFINEMENT_DIRS = [
    ("azimuthal", "nAz", False),
    ("radial", "nRad", False),
    ("streamwise", "nxIn", False),
    ("wall-normal", "firstLayer", True),
]


def refinement_factors(study, ids):
    """Per-direction refinement factor from each level to the next.

    Returns {direction: [factor, ...]} over the ids given, coarsest first.
    """
    out = {}
    levels = []
    for cid in ids:
        try:
            levels.append(study.params_for(cid))
        except Exception:                                   # noqa: BLE001
            return out
    for name, key, inverse in REFINEMENT_DIRS:
        fs = []
        for a, b in zip(levels, levels[1:]):
            x, y = a.get(key), b.get(key)
            if not x or not y:
                fs = []
                break
            fs.append(float(x) / float(y) if inverse else float(y) / float(x))
        if fs:
            out[name] = fs
    return out


def anisotropy(study, ids, tol=0.15):
    """Directions whose refinement factor differs from the fastest one.

    The grid-convergence index rests on ONE refinement ratio: the error is
    assumed to go as h^p for a single h, which is only meaningful if every
    direction was refined by the same factor.  Refine the bulk and leave the
    wall spacing alone and the reported h - built from cells per unit area -
    follows the direction that moved, while the error is set partly by the one
    that did not.  The extrapolation is then reading a mixture, and nothing in
    the number says so.

    Returns [(direction, mean factor)] for the directions that lag, plus the
    leader, or [] when the ladder is uniform enough.
    """
    fac = refinement_factors(study, ids)
    if len(fac) < 2:
        return []
    mean = {k: sum(v) / len(v) for k, v in fac.items()}
    lead = max(mean.values())
    if lead <= 1.0:
        return []
    lag = [(k, m) for k, m in sorted(mean.items(), key=lambda kv: kv[1])
           if (lead - m) / lead > tol]
    return [(k, m) for k, m in lag] + [("__lead__", lead)] if lag else []


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
    out["unconverged"] = [r["id"] for r in rows
                          if r.get("converged") is False and not r.get("transient")]
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
        good = [r for r in rs
                if r.get("transient") or r.get("converged") is not False]
        lad = {"ladder": name, "levels": rs, "converged": [r["id"] for r in good],
               "triplets": [], "chosen": None}
        #  coarsest first for the refinement factors: a ladder is built by
        #  refining, and the factors only read as ">1 means finer" that way
        order = sorted(rs, key=lambda r: -r["h"])
        lad["refinement"] = refinement_factors(study, [r["id"] for r in order])
        lad["anisotropy"] = anisotropy(study, [r["id"] for r in order])
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

#  TRANSIENT, not steady.  The first attempt at stage 2 ran every case steady
#  and all eight stalled at a residual of 7.5e-2 in continuity.  That is not a
#  mesh being too coarse: a strictly two-dimensional bank at Re_max of 1e4 and
#  above has no steady solution to converge to - a 2-D bluff-body wake is
#  time-periodic above Re of order 200 and these are fifty to five hundred
#  times that.  The domain's own symmetry planes, which are there so the
#  geometry matches what a correlation describes, remove the spanwise
#  decorrelation that would otherwise break the vortices up, so this is the
#  most shedding-prone version of the problem rather than the least.  Shen et
#  al. (2024) solved the same problem with URANS and reported time averages.
BASE_SETTINGS = {
    "general": {"steady": False, "energy": False},
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
    #  time_step, time_steps and average_last are filled in per case from the
    #  shedding period - see transient_controls - so what is here is only the
    #  part that does not depend on the flow
    "run": {"init_method": "hybrid", "iterations": 800,
            "residual_criterion": 1e-5, "monitor_every": 20,
            "max_iter_per_step": SHEDDING["max_iter_per_step"],
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
    over = {}
    if settings.get("inlet", {}).get("velocity") is not None:
        over["inlet"] = {"velocity": settings["inlet"]["velocity"]}
    if not settings["general"]["steady"]:
        case, st = case_state(geometry, p, settings)
        tc = transient_controls(st, bundle_length=case.l_bund * case.export_scale)
        over["run"] = {k: tc[k] for k in ("time_step", "time_steps",
                                          "average_last", "max_iter_per_step")}
        #  the flushing figures are not settings Fluent takes; they are the
        #  reason the step count is what it is, and they ride in c["transient"]
        #  beside it rather than being recomputed by every reader
        c["transient"] = tc
    if over:
        c["settings"] = over
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
    #  A mesh study writes its case files.  Eight of them is affordable, and
    #  when a ladder does something unexpected the thing you want is to open
    #  the offending mesh and look at the field - which a result record cannot
    #  give you, because it holds histories and numbers and no solution.  The
    #  sweep leaves them off: thirty cases of this is a different question.
    base_settings = {k: dict(v) for k, v in BASE_SETTINGS.items()}
    base_settings["run"]["write_case"] = True
    d = {
        "name": name, "kind": "mesh", "geometry": geometry,
        "title": "2단계 · %s 격자 민감도" % geometry,
        "created": time.strftime("%Y-%m-%d"),
        "phi": "dp_bundle",
        "notes": [
            "물 20 C, k-omega SST, y+ = 1, 비정상 (URANS).",
            "정상상태로 풀 수 없기 때문입니다. 엄밀 2차원 다발의 Re_max 1e4 이상 "
            "에서는 정상해가 존재하지 않고, 처음 시도에서 8개 케이스 전부 "
            "continuity 잔차 7.5e-2 에서 멈췄습니다.",
            "시간 간격은 케이스마다 와류 방출 주기에서 유도합니다 "
            "(St = 0.2, 주기당 25 스텝, 20 주기, 마지막 15 주기를 평균).",
            "결과는 마지막 순간값이 아니라 Δp 의 시간 평균입니다.",
            "rod 외 모든 벽은 대칭면 - 상관식이 기술하는 무한 다발과 같게 하기 위함.",
            "격자는 면내 방향으로 단계당 약 1.35배 세밀해지고, 반경 방향은 "
            "첫 셀을 y+ = 1에 고정한 채 성장비를 조여서 세밀해집니다.",
            "판정은 Celik 등(2008) / ASME V&V 20 의 GCI 절차로 합니다.",
            "case 파일을 함께 씁니다 - 사다리가 이상하게 나올 때 그 격자를 열어 "
            "유동장을 봐야 하기 때문입니다. 결과 기록만으로는 유동장을 되살릴 수 "
            "없습니다.",
        ],
        "base": {"params": dict(BASE_PARAMS), "settings": base_settings},
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
            "물 20 C, k-omega SST, y+ = 1, 비정상 (URANS), rod 외 벽은 모두 대칭면.",
            "시간 간격은 케이스마다 와류 방출 주기에서 유도하고, 결과는 Δp 의 "
            "시간 평균입니다.",
            "격자 해상도는 2단계에서 고른 수준(%s)을 씁니다." % level,
            "케이스 수가 많아 case 파일은 쓰지 않습니다. 결과 기록(잔차·Δp 이력, "
            "격자 수, 측정된 압력강하)은 남으므로 파라메트릭 탭에서 다시 열 수 "
            "있지만, 유동장은 남지 않습니다.",
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
def run_controls(study, cid):
    """The solver controls one case would be run with, as recorded on a row."""
    st = study.settings_for(cid)
    out = {k: st["run"].get(k)
           for k in ("time_step", "time_steps", "average_last",
                     "max_iter_per_step", "iterations", "residual_criterion")}
    out["steady"] = bool(st["general"]["steady"])
    return out


def definition_drift(study):
    """Recorded cases whose run no longer matches the definition on disk.

    A study is a study only because every case in it was run the same way;
    that is the whole basis for putting them in one table.  A definition can
    change under existing results - the rule that sets the step count was
    wrong once and had to be corrected - and the next run would then quietly
    append cases advanced differently from the ones already there.  Neither
    half is wrong on its own.  The mixture is, and nothing in the table shows
    it.

    Returns [(id, field, recorded, defined)], empty when they agree or when
    the recorded rows predate this check.
    """
    out = []
    for r in study.results():
        was = r.get("run_controls")
        if not was or not r.get("ok"):
            continue
        try:
            now = run_controls(study, r["id"])
        except Exception:                                   # noqa: BLE001
            continue
        for k in sorted(now):
            a, b = was.get(k), now.get(k)
            if a is None or b is None:
                continue
            if isinstance(a, float) or isinstance(b, float):
                if abs(float(a) - float(b)) <= 1e-9 * max(1.0, abs(float(b))):
                    continue
            elif a == b:
                continue
            out.append((r["id"], k, a, b))
    return out


class Runner(object):
    """One study, one case after another, through the app's own Job.

    It drives the same Job the Run tab drives and hands it to the app as the
    current job, so a campaign is watchable: the residuals, the pressure-drop
    monitor, the log and the Results tab all follow whichever case is running.
    Nothing about a study case is special except that something else pressed
    the button.
    """

    #  a launch failure that says any of these is the licence server, not the
    #  case: the ramp stops adding workers rather than failing every remaining
    #  case one at a time
    LICENCE_WORDS = ("licen", "flexlm", "ansyslmd", "lmgrd", "-15", "1055")

    def __init__(self, app, study, only=None, redo=False, force=False,
                 workers=1, cores=None):
        self.app = app
        self.study = study
        self.redo = bool(redo)
        #  How many Fluent sessions to keep going at once.  One case cannot
        #  fill a workstation - the solve is a fraction of each case and the
        #  meshes are small enough that four ranks is already generous - so
        #  the way to use the machine is more CASES, not more cores per case.
        #  What limits it is solver TASKS in the licence, which this cannot
        #  know, so it finds out: workers are started one at a time and each
        #  has to prove a session will launch before the next is added.
        self.workers = max(1, min(int(workers or 1), 16))
        #  and the other half of that sentence, which used to be unreachable.
        #  A campaign took whatever the settings schema defaulted to - four
        #  ranks - for every case, and no panel could change it, so the advice
        #  above could not be followed even by someone who agreed with it.
        #  Four ranks on a 13 504-cell mesh is 3 376 cells a rank; below
        #  roughly 50 000 the partition boundaries cost more than the cells
        #  inside them save, and this whole study is under that at four.
        self.cores = int(cores) if cores else None
        self.capped = None              # why the ramp stopped, if it did
        self.ready = threading.Event()
        self.qlock = threading.Lock()
        self.live = {}                  # worker -> the case it is on
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
        self.cursor = 0
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

    def _cores(self):
        """Ranks per case: what this run asked for, or the schema's default."""
        if self.cores:
            return self.cores
        try:
            return int(self.study.settings_for(self.queue[0])["launch"]
                       ["processors"])
        except Exception:                               # noqa: BLE001
            return 1

    def cells_per_core(self):
        """{case id: cells per rank} for the cases that have been meshed.

        Estimated from whatever has already run; a case that has not been
        meshed yet has no cell count to divide, so it is simply absent.
        """
        n = self._cores()
        out = {}
        for r in self.study.results():
            if r.get("cells") and r["id"] in self.queue:
                out[r["id"]] = int(r["cells"] / max(n, 1))
        return out

    def status(self):
        return {
            "study": self.study.name, "kind": self.study.kind,
            "total": len(self.queue), "index": self.index,
            "current": self.current, "queue": list(self.queue),
            "skipped": list(self.skipped),
            "stopping": self.stopping, "finished": self.finished,
            "error": self.error, "stopped_reason": self.stopped_reason,
            "workers": self.workers, "cores": self.cores,
            "capped": self.capped,
            "live": dict(self.live),
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
        """One case alone, then as many at once as the licence will bear.

        The order is not an accident.  The first case is the smoke test: if the
        set-up does not produce a usable solution, running seven more of it in
        parallel only wastes the machine faster.  Only once it has passed does
        the ramp start, and it starts one worker at a time because nothing here
        can know how many solver tasks the licence has - so each new worker has
        to prove a session will launch before the next is added.
        """
        try:
            if self.skipped:
                self.log("%d case(s) already done, skipping them: %s"
                         % (len(self.skipped), ", ".join(self.skipped)))
            if not self.queue:
                return
            self.log("%d case(s) to run, up to %d at once, %d core(s) each "
                     "(%d cores in total)"
                     % (len(self.queue), self.workers, self._cores(),
                        self.workers * self._cores()))
            #  the arithmetic that decides whether more ranks per case helps.
            #  Said once, at the top, from the meshes this study actually
            #  defines - not as a rule of thumb the reader has to apply.
            per = self.cells_per_core()
            if per:
                lo, hi = min(per.values()), max(per.values())
                if hi < CELLS_PER_CORE:
                    self.log("  %d-%d cells per core: below about %d the "
                             "partition boundaries cost more than the cells "
                             "inside them save, so fewer cores per case and "
                             "more cases at once is the faster arrangement"
                             % (lo, hi, CELLS_PER_CORE))

            #  --- case one, alone ---
            first = self._next()
            row = self._run_case(first, worker=0)
            if not self.stopping:
                why = None if self.force else self._first_case_verdict(row)
                if why:
                    self.stopped_reason = why
                    self.log("STOPPING after the first case: " + why)
                    self.log("  nothing is wrong with the study definition "
                             "- the SET-UP does not produce a usable solution, "
                             "and more of the same would not change that. Fix "
                             "it, or re-run with force.")
                    return

            #  --- and now the rest, ramped ---
            threads = []
            for w in range(self.workers):
                if self.stopping or self._empty():
                    break
                t = threading.Thread(target=self._worker, args=(w,),
                                     name="study-w%d" % w, daemon=True)
                t.start()
                threads.append(t)
                if w + 1 >= self.workers:
                    break
                #  wait for this one to get a session up before adding another
                if not self.ready.wait(timeout=900):
                    self.capped = "worker %d never got a session up" % w
                    break
                self.ready.clear()
                if self.capped:
                    break
            if self.capped:
                self.log("running %d at a time: %s" % (len(threads), self.capped))
            for t in threads:
                t.join()
        except Exception as exc:                                # noqa: BLE001
            self.error = "%s: %s" % (type(exc).__name__, exc)
            self.log("the runner itself failed: " + self.error)
        finally:
            self.current = None
            self.live = {}
            self.finished = True
            self.log("done - %d of %d case(s) attempted"
                     % (self.index, len(self.queue)))

    # -- the queue, shared -------------------------------------------------
    def _next(self):
        with self.qlock:
            if self.cursor >= len(self.queue):
                return None
            cid = self.queue[self.cursor]
            self.cursor += 1
            self.index = self.cursor
            self.current = cid
            return cid

    def _empty(self):
        with self.qlock:
            return self.cursor >= len(self.queue)

    def _worker(self, w):
        while not self.stopping:
            cid = self._next()
            if cid is None:
                break
            self._run_case(cid, worker=w)
        self.live.pop(w, None)

    def _run_case(self, cid, worker):
        """One case, recorded whatever happens to it."""
        self.live[worker] = cid
        self.log("[%d/%d] %s%s" % (self.index, len(self.queue), cid,
                                   "  (w%d)" % worker if self.workers > 1 else ""))
        try:
            row = self._one(cid, worker=worker)
        except Exception as exc:                                # noqa: BLE001
            row = {"id": cid, "ok": False,
                   "error": "%s: %s" % (type(exc).__name__, exc)}
            self.log("  %s FAILED: %s" % (cid, row["error"]))
        self.study.record(row)
        if row.get("ok"):
            self.log("  %s: dp_bundle %.4g Pa, Eu_row %.4f, %d cells%s"
                     % (cid, row["dp_bundle"] or float("nan"),
                        row["eu_row"] or float("nan"), row["cells"] or 0,
                        "   (MOCK)" if row.get("mock") else ""))
        self.live.pop(worker, None)
        #  a finished case is also proof that a session could be had, and the
        #  ramp must not sit on `ready` waiting for a signal from a case that
        #  has already come and gone
        self.ready.set()
        return row

    @staticmethod
    def _first_case_verdict(row):
        """Why the campaign should not continue past case one, or None."""
        if not row.get("ok"):
            return "it failed: %s" % (row.get("error") or "?")
        if row.get("mock"):
            return None                     # the mock is for testing the plumbing
        #  before either branch: a run that measured nothing is no use however
        #  well its residuals behaved
        if row.get("dp_bundle") is None:
            return "the bundle pressure drop could not be measured"
        if row.get("transient"):
            #  A transient run's inner iterations are not supposed to reach a
            #  steady criterion, so convergence is judged on the ANSWER: did
            #  the time average settle?  A trace still spanning half its own
            #  mean has not been run long enough.
            spread = row.get("dp_spread")
            if spread is not None and spread > 0.5:
                return ("the pressure drop is still swinging %.0f %% of its "
                        "own mean over the averaging window - the run is too "
                        "short, or it has not settled" % (100 * spread))
            if row.get("dp_averaged_over", 0) < 20:
                return ("only %d samples in the averaging window"
                        % row.get("dp_averaged_over", 0))
            #  and NOT the residual criterion.  A transient run's inner
            #  iterations are not meant to reach it - continuity in particular
            #  floors out at the level the pressure-velocity coupling can hold
            #  within one time step - so the answer is judged on the answer.
            md = mean_drift(row)
            if md is not None and md > 0.01:
                return ("the time average is still moving by %.1f %% across "
                        "the averaging window" % (100 * md))
            return None
        if row.get("converged") is False:
            return ("it did not converge - final residual %.2e%s against a "
                    "criterion of %.0e"
                    % (row.get("residual_worst") or 0.0,
                       " in " + row["residual_worst_eq"]
                       if row.get("residual_worst_eq") else "",
                       row.get("criterion") or 0.0))
        if (row.get("dp_drift") or 0.0) > 0.05:
            return ("the pressure drop was still moving by %.1f %% over the "
                    "final stretch" % (100 * row["dp_drift"]))
        return None

    def _licence_trouble(self, msg):
        m = (msg or "").lower()
        return any(w in m for w in self.LICENCE_WORDS)

    def _one(self, cid, worker=0):
        """Run one case and measure it."""
        import app as APP
        s = self.study
        c = s.case(cid)
        params = s.params_for(cid)
        settings = s.settings_for(cid)
        if self.cores:
            settings["launch"] = dict(settings.get("launch") or {},
                                      processors=self.cores)
        t0 = time.time()

        #  Worker 0's case is handed to the app as the current job, so the Run
        #  tab's residuals, monitor and log follow it.  The others are not:
        #  there is one live view and it cannot show four cases at once, and
        #  quietly flipping it between them would be worse than not having it.
        if worker == 0:
            job = self.app.new_job(s.d["geometry"], params, settings)
        else:
            job = APP.Job("%s-w%d" % (cid, worker), s.d["geometry"], params,
                          settings, self.app.backend, self.app.out_dir)
        job.log("study %s, case %s" % (s.name, cid))
        job.start()
        signalled = False
        while job.finished_at is None:
            time.sleep(0.4)
            #  a session is up: the next worker may start
            if not signalled and job.stage not in ("queued", "meshing",
                                                   "launching"):
                signalled = True
                self.ready.set()
            if self.stopping:
                job.stop()
        if not signalled:
            #  it never got past launching.  If the licence said no, stop
            #  adding workers rather than failing every remaining case on it.
            if job.error and self._licence_trouble(job.error):
                self.capped = ("the licence would not give worker %d a session "
                               "(%s)" % (worker, job.error.split(":")[-1].strip()[:80]))
            self.ready.set()
        row = {"id": cid, "ok": False, "label": c.get("label"),
               "ladder": c.get("ladder", "main"), "meta": c.get("meta"),
               "seconds": round(time.time() - t0, 1),
               "when": time.strftime("%Y-%m-%d %H:%M:%S"),
               "code": APP.version_line(APP.VERSION),
               "mock": bool(job.driver and job.driver.mock),
               "mesh_path": (os.path.basename(job.mesh_path)
                             if job.mesh_path else None),
               "error": job.error, "worker": worker}
        if job.error:
            if worker != 0:
                job.close()
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

        #  The pressure drop across the BUNDLE, the way the correlation defines
        #  it: two planes on the bundle faces, half a pitch outside the first
        #  and last rod centres.  inlet - outlet would span the boxes too.
        #  The driver stands those planes up during setup now, so the monitor
        #  can follow the bundle drop while the case runs; the study uses them
        #  if they are there and makes its own if they are not.
        drv = job.driver
        transient = bool(getattr(drv, "transient", lambda: False)())
        row["transient"] = transient
        #  which spelling of the time mode this Fluent took.  Different
        #  installations of the same release allow different ones, and
        #  1st-order implicit damps the very oscillation being measured, so
        #  the scheme belongs beside the number it produced.
        row["time_scheme"] = getattr(drv, "time_scheme", None)
        #  what this case was actually advanced with, so a later run can tell
        #  whether it is about to add cases to a study that no longer matches
        row["run_controls"] = {k: settings["run"].get(k)
                               for k in ("time_step", "time_steps",
                                         "average_last", "max_iter_per_step",
                                         "iterations", "residual_criterion")}
        row["run_controls"]["steady"] = bool(settings["general"]["steady"])
        #  how much of the bundle's own flow-through the run covered before it
        #  began averaging.  Under one sweep the last rows were still carrying
        #  the first rows' initial condition when the average opened.
        tc = c.get("transient") or {}
        for k in ("bundle_flowthrough_s", "flowthroughs",
                  "flowthroughs_averaged"):
            if tc.get(k) is not None:
                row[k] = tc[k]
        if tc.get("bundle_flowthrough_s"):
            _dt = float(settings["run"]["time_step"])
            row["flush_flowthroughs"] = (
                (int(settings["run"]["time_steps"])
                 - int(settings["run"]["average_last"])) * _dt
                / tc["bundle_flowthrough_s"])
        names = list(drv.bundle_surfaces or [])
        made = []
        try:
            if not names:
                names = ["study_bundle_in", "study_bundle_out"]
                for nm, x in zip(names, [case.x_b0 * k, case.x_b1 * k]):
                    drv.make_plane(nm, "x", x)
                    made.append(nm)
            p = [drv.report("area-weighted-avg", [nm], "pressure")
                 for nm in names]
            row["p_bundle_in"], row["p_bundle_out"] = p[0], p[1]
            row["dp_bundle_last"] = p[0] - p[1]
            row["dp_bundle"] = row["dp_bundle_last"]
        finally:
            for nm in made:
                try:
                    drv.drop_plane(nm)
                except Exception:                               # noqa: BLE001
                    pass

        #  On a transient run the ANSWER is the time average, not whatever the
        #  last instant happened to be: the last instant is one sample of an
        #  oscillation.  The instantaneous value is kept beside it so the two
        #  can be compared, which is also how you see whether the averaging
        #  window was long enough.
        avg = getattr(drv, "dp_time_average", lambda: None)()
        if transient and avg:
            row["dp_bundle"] = avg["mean"]
            row["dp_averaged_over"] = avg["n"]
            row["dp_samples"] = avg["of"]
            row["dp_spread"] = avg["spread"]
            row["dp_mean_drift"] = avg.get("mean_drift")
            row["dp_average_of"] = avg["key"]
        if row.get("dp_bundle") is not None:
            row["eu_row"] = row["dp_bundle"] / (st["n_rows"] * st["q"])
            row["dp_per_row"] = row["dp_bundle"] / max(1, st["n_rows"])
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
            row["residual_history"] = _thin(res, RES_KEEP)
        mon = list(getattr(job.driver, "monitors", []) or [])
        row["dp_history"] = _thin(mon, DP_KEEP)
        if len(mon) >= 2:
            #  how much the answer was still moving over the last fifth of the
            #  run.  Residuals settling is not the answer settling, and a study
            #  that recorded only the final number could not tell them apart.
            tail = mon[max(0, len(mon) - max(2, len(mon) // 5)):]
            vals = [m["dp"] for m in tail if m.get("dp") is not None]
            if len(vals) >= 2 and vals[-1]:
                row["dp_drift"] = (max(vals) - min(vals)) / abs(vals[-1])
        row["ok"] = True
        #  Worker 0's session is deliberately NOT closed: the next case's
        #  new_job() closes it before launching its own, so one licence is
        #  held at a time and the LAST case stays open for the Results tab.
        #  Every other worker owns its session outright and nothing else will
        #  ever close it, so it closes its own.
        if worker != 0:
            job.close()
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
#  Roughly where a partitioned solve stops paying for itself.  Below this
#  many cells on a rank the halo exchange each iteration costs more than the
#  interior cells it saves, so adding ranks to one case makes it slower while
#  running another case alongside it makes the machine faster.  A round number
#  from practice, not a measurement of this machine - which is why it is used
#  to say something to the reader rather than to decide anything on its own.
CELLS_PER_CORE = 50000

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
    #  the legend wraps.  Eight cases on one row ran off the right-hand edge
    #  and the last labels were simply cut, which is the one failure mode a
    #  legend must not have - a chart you cannot read the key of is a chart
    #  that says nothing.
    L, R, B = 66, 14, 44
    legend, row, used = [], [], 0.0
    for i, sr in enumerate(series):
        if not sr.get("label"):
            continue
        w = 26 + 6.0 * len(sr["label"])
        if row and used + w > width - L - R:
            legend.append(row)
            row, used = [], 0.0
        row.append((i, sr))
        used += w
    if row:
        legend.append(row)
    T = 16 + (16 if title else 0) + 14 * len(legend)
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
    #  legend: as many rows as it takes, above the frame
    for ri, row in enumerate(legend):
        lx = L
        ly = T - 6 - 14 * (len(legend) - 1 - ri)
        for i, sr in row:
            col = sr.get("color") or PALETTE[i % len(PALETTE)]
            o.append('<rect x="%.1f" y="%.1f" width="15" height="2.6" '
                     'fill="%s"/>' % (lx, ly - 3, col))
            o.append('<text x="%.1f" y="%.1f" font-size="10" '
                     'fill="#33415a">%s</text>' % (lx + 20, ly + 1, sr["label"]))
            lx += 26 + 6.0 * len(sr["label"])
    o.append("</svg>")
    return "".join(o)


#  The two charts that ARE the result of a stage, built in one place because
#  the stage report and the at-a-glance panel both draw them and a second copy
#  is a second thing to keep in step.
def _grid_series(lad, levels, t):
    """Δp against cell size, with the Richardson extrapolation as a line."""
    pts = [(r["h"], r["dp_bundle"]) for r in sorted(levels, key=lambda r: r["h"])]
    hs = [p[0] for p in pts]
    return [{"label": t("CFD", "CFD"), "points": pts},
            {"label": t("Richardson 외삽 (h→0)", "Richardson h -> 0"),
             "points": [(min(hs) * 0.7, lad["phi_ext"]),
                        (max(hs) * 1.1, lad["phi_ext"])],
             "dash": True, "marker": False, "color": "#64748b"}]


def _eu_re_series(rs, keys, names):
    """Our Eu against Re, and every correlation that is in range, dashed."""
    series = [{"label": "CFD", "points": [(r["Re"], r["eu_row"]) for r in rs],
               "color": PALETTE[0]}]
    for j, k in enumerate(keys):
        pts = [(r["Re"], r["corr"][k]["eu_row"]) for r in rs
               if r["corr"][k]["eu_row"]]
        if pts:
            series.append({"label": names[k], "points": pts, "dash": True,
                           "color": PALETTE[(j + 1) % len(PALETTE)]})
    return series


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
    #  A TRANSIENT run is not judged on the residual criterion and never was:
    #  mesh_analysis keeps it in the extrapolation, and the criterion belongs
    #  to a steady solver.  Continuity in particular floors out at whatever
    #  the pressure-velocity coupling can hold inside one time step, so
    #  "did not reach 1e-5" is a statement about the time step, not about the
    #  answer.  This banner said those cases had been thrown out of the GCI.
    #  They had not been, and they should not be; it was telling the reader
    #  something untrue about the analysis directly below it.
    unconv = [r for r in rows if r.get("converged") is False
              and not r.get("transient")]
    tr_floor = [r for r in rows if r.get("transient")
                and r.get("converged") is False]
    if tr_floor:
        out.append(("warn",
                    "비정상 해석 케이스 %d개는 정상상태 잔차 기준(%.0e)에 닿지 "
                    "않았습니다. 시간 스텝 안에서 continuity 잔차는 더 내려가지 "
                    "않는 바닥이 있으므로 이는 정상입니다 — 판정은 잔차가 아니라 "
                    "Δp 시간평균이 멈췄는지로 하고, 이 케이스들은 아래 외삽에 "
                    "그대로 들어갑니다."
                    % (len(tr_floor), tr_floor[0].get("criterion") or 0.0)
                    if ko else
                    "%d transient case(s) did not reach the steady residual "
                    "criterion (%.0e). That is expected - continuity floors "
                    "out at what the pressure-velocity coupling can hold "
                    "inside one time step - so they are judged on whether the "
                    "Δp time average stopped moving, and they are included in "
                    "the extrapolation below."
                    % (len(tr_floor), tr_floor[0].get("criterion") or 0.0)))
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
    #  the same distinction again: on a transient run the thing that has to
    #  have stopped moving is the MEAN, not the instantaneous value.  A
    #  settled vortex street swings ten per cent every period and always will.
    drift = []
    for r in rows:
        if r.get("transient"):
            md = mean_drift(r)
            if md is not None and md > 0.01:
                drift.append(r)
        elif (r.get("dp_drift") or 0) > 0.01:
            drift.append(r)
    if drift:
        out.append(("warn", "Δp 가 아직 1 %% 이상 움직인 케이스 %d개 (%s) - "
                            "비정상 해석은 시간평균 기준, 정상 해석은 마지막 "
                            "구간 기준입니다."
                    % (len(drift), ", ".join(r["id"] for r in drift[:8]))
                    if ko else
                    "%d case(s) still had Δp moving by more than 1 %% (%s) - "
                    "measured on the time average for a transient run and "
                    "over the final stretch for a steady one."
                    % (len(drift), ", ".join(r["id"] for r in drift[:8]))))
    #  Flushing.  A statistical average over a flow that has not crossed the
    #  bundle yet is an average over the initial condition, and how far short
    #  it falls differs from case to case, so it is not even a consistent bias.
    short = [r for r in rows
             if r.get("flush_flowthroughs") is not None
             and r["flush_flowthroughs"] < 1.0]
    if short:
        worst = min(r["flush_flowthroughs"] for r in short)
        out.append(("bad",
                    "\ud3c9\uade0\uc744 \uc2dc\uc791\ud558\uae30 \uc804\uc5d0 \ub2e4\ubc1c\uc744 \ud55c \ubc88\ub3c4 \ud1b5\uacfc\ud558\uc9c0 \ubabb\ud55c \ucf00\uc774\uc2a4 "
                    "%d\uac1c (\ucd5c\uc18c %.2f\ud68c). \ub4b7\uc904\uc774 \uc544\uc9c1 \uc55e\uc904\uc758 \ucd08\uae30\uc870\uac74\uc744 \uc9c0\ub098\uac00\ub294 "
                    "\uc911\uc5d0 \ud3c9\uade0\uc744 \ub0c8\ub2e4\ub294 \ub73b\uc774\uace0, \ucf00\uc774\uc2a4\ub9c8\ub2e4 \uadf8 \uc815\ub3c4\uac00 \ub2ec\ub77c "
                    "\uc77c\uad00\ub41c \ud3b8\ud5a5\ub3c4 \uc544\ub2d9\ub2c8\ub2e4."
                    % (len(short), worst) if ko else
                    "%d case(s) began averaging before the flow had swept the "
                    "bundle once (worst %.2f sweeps). The last rows were "
                    "averaging over the first rows' initial condition, by "
                    "different amounts in different cases, so it is not even "
                    "a consistent bias." % (len(short), worst)))

    #  Which spelling of the time mode a given Fluent accepted is not
    #  bookkeeping.  First-order implicit damps the shedding oscillation that
    #  the time step was sized to resolve and that the time average is taken
    #  over, so a run that fell back to it measured something else - and two
    #  cases on different schemes are not comparable with each other at all.
    schemes = sorted({r.get("time_scheme") for r in rows
                      if r.get("transient") and r.get("time_scheme")})
    if len(schemes) > 1:
        out.append(("bad", "케이스마다 시간 이산화가 다릅니다 (%s). 같은 연구 안의 "
                           "케이스끼리 비교할 수 없습니다."
                    % ", ".join(schemes) if ko else
                    "the cases did not all advance time the same way (%s); "
                    "cases in one study on different schemes are not "
                    "comparable with each other." % ", ".join(schemes)))
    elif schemes and schemes[0] == "unsteady-1st-order":
        out.append(("warn", "시간 이산화가 1차 음해법입니다. 시간 간격은 와류 이탈 "
                            "주기를 분해하도록 잡혀 있는데 1차는 그 진동을 감쇠시키므로 "
                            "Δp 진폭이 실제보다 작게 나옵니다."
                    if ko else
                    "time was advanced with first-order implicit. The time "
                    "step is sized to resolve the shedding period and "
                    "first-order damps that oscillation, so the Δp swing is "
                    "smaller here than it should be."))
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
        #  Said ABOVE the table, because it decides whether the table below it
        #  means anything.  A grid-convergence index rests on one refinement
        #  ratio; a ladder that refines the bulk and leaves the wall spacing
        #  alone has several, the reported h follows whichever moved, and the
        #  extrapolation reads a mixture without saying so.
        aniso = lad.get("anisotropy") or []
        if aniso:
            lead = dict(aniso).get("__lead__")
            lag = [(k, v) for k, v in aniso if k != "__lead__"]
            names_ko = {"azimuthal": "원주", "radial": "반경",
                        "streamwise": "유동", "wall-normal": "벽면 수직"}
            h.append('<p class="bad">%s</p>' % _esc(
                "이 사다리는 방향마다 다른 비율로 세밀해집니다: 가장 빠른 방향이 "
                "단계당 %.2f배인데 %s. 격자 수렴 지수는 단일 세밀화 비를 전제로 "
                "하므로 (오차 ~ h^p, h 하나), 아래 외삽값은 여러 비가 섞인 값을 "
                "읽고 있습니다. 특히 벽면 수직 방향이 1.00배면 벽 근처 격자는 "
                "전혀 세밀해지지 않은 것입니다."
                % (lead, ", ".join("%s %.2f배" % (names_ko.get(k, k), v)
                                   for k, v in lag))
                if ko else
                "this ladder refines by a different factor in each direction: "
                "the fastest is %.2f per level while %s. A grid-convergence "
                "index assumes ONE refinement ratio (error ~ h^p for a single "
                "h), so the extrapolation below is reading a mixture. A "
                "wall-normal factor of 1.00 in particular means the near-wall "
                "mesh was not refined at all."
                % (lead, ", ".join("%s is %.2f" % (k, v) for k, v in lag))))
            h.append('<table class="grid"><tr><th>%s</th>%s</tr>'
                     % (t("방향", "direction"),
                        "".join("<th>%d→%d</th>" % (i + 1, i + 2)
                                for i in range(max(
                                    (len(v) for v in
                                     (lad.get("refinement") or {}).values()),
                                    default=0)))))
            for k, v in (lad.get("refinement") or {}).items():
                h.append("<tr><th>%s</th>%s</tr>"
                         % (_esc(names_ko.get(k, k) if ko else k),
                            "".join('<td class="n">%.2f</td>' % x for x in v)))
            h.append("</table>")
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
            h.append(svg_plot(_grid_series(lad, levels, t),
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
        series = _eu_re_series(rs, keys, names)
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


# =============================================================================
#  AT A GLANCE
# =============================================================================
#  A stage report is written to be read end to end, and it draws the residual
#  curve of ONE case - the finest - because that is the one its argument turns
#  on.  That is the wrong shape for the question actually asked most often,
#  which is "how did the whole batch go?"  Answering it by opening eight or
#  thirty cases one at a time is how a bad set-up survives a night.
#
#  So: small multiples.  Every case on one screen, every panel on the SAME
#  axes, because the comparison is the entire point - a flat residual curve
#  next to a descending one is obvious, and a flat curve on its own is not.
#  Nothing here is computed that the report does not already compute; this is
#  the same record, arranged for the eye instead of for the argument.
OV_STATUS = {
    #  key:      (colour,    ko,            en)
    "queued":    ("#94a3b8", "대기",        "queued"),
    "failed":    ("#9b1c1c", "실패",        "failed"),
    "mock":      ("#7c3aed", "MOCK",        "MOCK"),
    "good":      ("#15803d", "정상",        "settled"),
    "moving":    ("#d97706", "아직 이동 중", "still moving"),
    "stalled":   ("#b91c1c", "미수렴",      "not converged"),
}


def case_status(row):
    """One word for how a case went, by the same rules the runner judges it.

    Deliberately delegates to Runner._first_case_verdict rather than
    re-deriving the test: two places deciding what "converged" means is two
    places to disagree, and the one that decides whether a campaign keeps
    running has to be the one that wins.
    """
    if row is None:
        return "queued"
    if not row.get("ok"):
        return "failed"
    if row.get("mock"):
        return "mock"
    if Runner._first_case_verdict(row) is None:
        return "good"
    #  the verdict says there is a problem; these three say WHICH, in the order
    #  that matters.  A steady run that never met its criterion is stalled; a
    #  run that finished but measured nothing is no more use than one that
    #  crashed; everything else is an answer that had not stopped moving.
    if row.get("converged") is False and not row.get("transient"):
        return "stalled"
    if row.get("dp_bundle") is None:
        return "failed"
    return "moving"


def _worst_trace(hist):
    """The largest residual at each recorded point.

    Six equations on one panel at thumbnail size is a smudge.  The envelope is
    what the convergence test looks at anyway, so the panel shows exactly the
    number the verdict was made on.
    """
    out = []
    for r in hist:
        vals = [v for k, v in r.items() if k != "iter"
                and isinstance(v, (int, float)) and v > 0]
        if vals:
            out.append((r["iter"], max(vals)))
    return out


def _dp_deviation(hist):
    """The Dp trace as a percentage of its own settled mean.

    Cases at different Reynolds numbers have pressure drops orders of
    magnitude apart, so a shared Pa axis would show one curve and five flat
    lines.  What is being compared is not the value but whether it stopped
    moving, and that is scale-free.  The reference is the mean of the last
    fifth, which is the window the run itself averages over.
    """
    vals = [(r["iter"], r["dp"]) for r in hist if r.get("dp") is not None]
    if len(vals) < 2:
        return [], None
    tail = vals[max(0, len(vals) - max(2, len(vals) // 5)):]
    mean = sum(v for _, v in tail) / len(tail)
    if not mean:
        return [], None
    return [(i, 100.0 * (v / mean - 1.0)) for i, v in vals], mean


def svg_spark(series, width=236, height=86, ylog=False, yrange=None,
              xrange=None, hline=None, band=None, note=""):
    """A panel of one small multiple: the curve, a frame, and nothing else.

    The ranges are passed IN rather than fitted here - that is what makes a
    wall of these comparable.  A panel whose data leaves the given range is
    clipped at the frame rather than rescaled, because a panel that quietly
    rescaled itself would break the only promise this layout makes.
    """
    pts = [p for s in series for p in s["points"]]
    if not pts:
        return ('<svg class="spark" viewBox="0 0 %d %d" width="100%%" '
                'xmlns="http://www.w3.org/2000/svg"><text x="%d" y="%d" '
                'font-size="10.5" fill="#94a3b8" text-anchor="middle">%s'
                '</text></svg>' % (width, height, width // 2, height // 2 + 3,
                                   _esc(note or "-")))
    x0, x1 = xrange or (min(p[0] for p in pts), max(p[0] for p in pts))
    y0, y1 = yrange or (min(p[1] for p in pts), max(p[1] for p in pts))
    if x1 <= x0:
        x1 = x0 + 1
    if y1 <= y0:
        y1 = y0 + (abs(y0) or 1.0) * 0.1
    L, R, T, B = 4, 4, 4, 4

    def px(v):
        return L + (v - x0) / float(x1 - x0) * (width - L - R)

    def py(v):
        f = ((math.log10(max(v, 1e-30)) - math.log10(y0))
             / (math.log10(y1) - math.log10(y0)) if ylog
             else (v - y0) / float(y1 - y0))
        return T + (1.0 - min(max(f, 0.0), 1.0)) * (height - T - B)

    o = ['<svg class="spark" viewBox="0 0 %d %d" width="100%%" '
         'xmlns="http://www.w3.org/2000/svg" font-family="inherit">'
         % (width, height)]
    if band:
        lo, hi = band
        o.append('<rect x="%d" y="%.1f" width="%.1f" height="%.1f" '
                 'fill="#dcfce7"/>' % (L, py(hi), width - L - R,
                                       max(py(lo) - py(hi), 1.0)))
    if hline is not None and y0 < hline < y1:
        o.append('<line x1="%d" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#94a3b8" '
                 'stroke-width="1" stroke-dasharray="4 3"/>'
                 % (L, py(hline), width - R, py(hline)))
    for s in series:
        ps = s["points"]
        if len(ps) < 2:
            continue
        d = " ".join("%s%.1f,%.1f" % ("M" if j == 0 else "L", px(x), py(y))
                     for j, (x, y) in enumerate(ps))
        o.append('<path d="%s" fill="none" stroke="%s" stroke-width="%s"/>'
                 % (d, s.get("color") or "#2563eb", s.get("width") or "1.5"))
    o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" '
             'stroke="#dde3ea"/>' % (L, T, width - L - R, height - T - B))
    o.append("</svg>")
    return "".join(o)


def _ov_cards(study, rows, ko, t):
    """One panel per case, all on shared axes."""
    by = {r["id"]: r for r in rows}
    cases = [(c["id"], by.get(c["id"])) for c in study.cases]
    #  the shared ranges, over every case that has anything to show
    rv = [v for _, r in cases if r for _, v in
          _worst_trace(r.get("residual_history") or [])]
    crits = [r.get("criterion") for _, r in cases
             if r and r.get("criterion")]
    if crits:
        rv.append(min(crits))
    ry = ((min(rv) / 2.0, max(rv) * 2.0) if len(rv) > 1 else None)
    devs = [abs(v) for _, r in cases if r for _, v in
            _dp_deviation(r.get("dp_history") or [])[0]]
    #  at least +-2 % so a settled batch does not get magnified into noise,
    #  at most +-50 % so one wild case does not flatten all the others - and
    #  snapped to a round number, because the limit is written in the panel
    #  label and "+-34.2455 %" is a number nobody asked for
    want = min(max(max(devs) * 1.15 if devs else 2.0, 2.0), 50.0)
    lim = next(v for v in (2, 5, 10, 20, 50) if v >= want - 1e-9)
    #  the x axis is the FRACTION of each run, not the iteration count: cases
    #  that ran for different numbers of iterations are still comparable in
    #  shape, which is what a wall of these is read for
    h = ['<div class="ovgrid">']
    for cid, r in cases:
        key = case_status(r)
        col, kko, ken = OV_STATUS[key]
        h.append('<div class="ovcard" data-case="%s" style="border-left-color:%s">'
                 % (_esc(cid), col))
        h.append('<div class="ovhead"><b>%s</b>'
                 '<span class="ovst" style="color:%s">%s</span></div>'
                 % (_esc(cid), col, _esc(kko if ko else ken)))
        if r is None:
            h.append('<div class="ovnums">%s</div></div>'
                     % t("아직 실행되지 않았습니다", "not run yet"))
            continue
        bits = []
        if r.get("eu_row") is not None:
            bits.append("Eu %s" % _n(r["eu_row"], "%.3f"))
        if r.get("cells"):
            bits.append("{:,}".format(r["cells"]).replace(",", " ") +
                        t(" 셀", " cells"))
        if r.get("Re"):
            bits.append("Re %s" % _n(r["Re"], "%.3g"))
        if r.get("residual_worst") is not None:
            bits.append("%s %s" % (t("최종 잔차", "final res"),
                                   _n(r["residual_worst"], "%.1e")))
        #  the number the verdict was actually made on, so the colour of the
        #  border is never something the reader has to take on trust
        if r.get("transient"):
            #  the swing is context; the MEAN DRIFT is the criterion.  Showing
            #  only the swing invited the reading that a ten-per-cent
            #  oscillation was a ten-per-cent uncertainty, which it is not.
            if r.get("dp_spread") is not None:
                bits.append("%s %s" % (t("Δp 진폭", "Δp swing"),
                                       "%.0f %%" % (100 * r["dp_spread"])))
            md = mean_drift(r)
            if md is not None:
                bits.append("%s %s" % (t("평균 이동", "mean drift"),
                                       "%.2f %%" % (100 * md)))
        elif r.get("dp_drift") is not None:
            bits.append("%s %s" % (t("Δp 이동", "Δp drift"),
                                   "%.1f %%" % (100 * r["dp_drift"])))
        h.append('<div class="ovnums">%s</div>' % _esc(" · ".join(bits) or "—"))
        res = _thin(_worst_trace(r.get("residual_history") or []), 140)
        h.append('<div class="ovplot"><span>%s</span>%s</div>'
                 % (t("잔차(최대)", "residual (worst)"),
                    svg_spark([{"points": _fracx(res), "color": col}],
                              ylog=True, yrange=ry, xrange=(0.0, 1.0),
                              #  a steady criterion drawn across a transient
                              #  panel would be a line the run was never asked
                              #  to cross
                              hline=(None if r.get("transient")
                                     else r.get("criterion")),
                              note=t("이력 없음", "no history"))))
        dev, mean = _dp_deviation(r.get("dp_history") or [])
        h.append('<div class="ovplot"><span>%s</span>%s</div>'
                 % (t("Δp 정착 ±%g%%" % lim, "Δp settling ±%g%%" % lim),
                    svg_spark([{"points": _fracx(_thin(dev, 140)),
                                "color": col}],
                              yrange=(-lim, lim), xrange=(0.0, 1.0),
                              band=(-1.0, 1.0), hline=0.0,
                              note=t("이력 없음", "no history"))))
        h.append("</div>")
    h.append("</div>")
    return "\n".join(h)


def _fracx(pts):
    """Re-index a trace onto 0..1 so runs of different length line up."""
    if len(pts) < 2:
        return list(pts)
    x0, x1 = pts[0][0], pts[-1][0]
    if x1 <= x0:
        return [(0.0, y) for _, y in pts]
    return [((x - x0) / float(x1 - x0), y) for x, y in pts]


def overview_body(study, lang="ko"):
    """Every case of one study on a single screen."""
    ko = lang == "ko"
    t = lambda a, b: a if ko else b                          # noqa: E731
    rows = study.results()
    h = ['<h1>%s · %s</h1>' % (_esc(study.d.get("title") or study.name),
                               t("한눈에 보기", "at a glance"))]
    counts = {}
    for c in study.cases:
        r = next((x for x in rows if x["id"] == c["id"]), None)
        k = case_status(r)
        counts[k] = counts.get(k, 0) + 1
    h.append('<p class="sub">%s</p>' % _esc(" · ".join(
        "%s %d" % ((OV_STATUS[k][1] if ko else OV_STATUS[k][2]), n)
        for k, n in sorted(counts.items(), key=lambda kv: -kv[1]))))
    for level, msg in _warnings(study, [r for r in rows if r.get("ok")], ko):
        h.append('<p class="%s">%s</p>'
                 % ("bad" if level == "bad" else "warn", _esc(msg)))
    h.append('<p class="muted">%s</p>' % t(
        "패널 하나가 케이스 하나입니다. 모든 패널의 축은 <b>동일</b>하므로 "
        "옆칸과 바로 비교할 수 있습니다. 가로축은 반복 횟수가 아니라 그 런의 "
        "<b>진행률(0→1)</b>이라 길이가 다른 런도 모양으로 비교됩니다. "
        "Δp 패널의 초록 띠는 ±1 % 입니다. 정상 해석 패널의 점선은 수렴 판정선이고, "
        "비정상 해석에는 그 선이 없습니다 — 시간 스텝 안의 잔차는 정상상태 기준에 "
        "닿도록 만든 것이 아니라서, 판정은 Δp 시간평균이 멈췄는지로 합니다. "
        "패널을 클릭하면 그 케이스가 해석 탭에 복원됩니다.",
        "One panel per case. Every panel is on the <b>same</b> axes, so a "
        "panel can be read against the one beside it. The abscissa is each "
        "run's <b>progress, 0 to 1</b>, not its iteration count, so runs of "
        "different length still compare by shape. The dashed line in a "
        "steady panel is its convergence criterion; a transient panel has no "
        "such line, because residuals inside a time step were never meant to "
        "reach a steady criterion - what is judged there is whether the Δp "
        "time average stopped moving. The green band in a Δp panel is ±1 %. "
        "Clicking a panel puts that case back on the Run tab."))
    h.append(_ov_cards(study, rows, ko, t))

    #  the same traces again, overlaid.  The wall says which case; the overlay
    #  says how far apart they are, which a wall of separate frames cannot.
    good = [r for r in rows if r.get("ok") and (r.get("residual_history")
                                                or r.get("dp_history"))]
    if len(good) > 1:
        h.append("<h2>%s</h2>" % t("겹쳐 보기", "overlaid"))
        wide = lambda svg: '<div class="ovwide">%s</div>' % svg   # noqa: E731
        #  Coloured by STATUS, not by case, and with a status key rather than
        #  thirty names.  A thirty-case sweep has five times more curves than
        #  the palette has colours, so a per-case legend would put the same
        #  blue against six different names - and the question a stack of
        #  curves is read for is not which one is main-L3, it is which ones
        #  went wrong.  The colours are the same ones the panels above are
        #  bordered with, so the two read as one picture.
        h.append('<p class="muted">%s</p>' % t(
            "곡선 색은 위 패널의 테두리 색과 같은 상태 색입니다 (케이스별 색이 "
            "아닙니다). 어느 케이스인지는 위 패널에서 보세요.",
            "Curve colour is the status colour the panels above are bordered "
            "with, not a per-case colour. Which case is which is in the panels."))

        def _stack(trace_of, **kw):
            out, seen = [], []
            for r in good:
                tr = trace_of(r)
                if not tr:
                    continue
                key = case_status(r)
                col = OV_STATUS[key][0]
                out.append({"points": tr, "marker": False, "color": col})
                if key not in seen:
                    seen.append(key)
            for key in seen:
                col, kko, ken = OV_STATUS[key]
                out.append({"label": kko if ko else ken, "points": [],
                            "color": col})
            return out

        series = _stack(lambda r: _fracx(_thin(
            _worst_trace(r.get("residual_history") or []), 140)))
        if series:
            crit = next((r.get("criterion") for r in good if r.get("criterion")),
                        None)
            if crit:
                series.append({"label": t("판정 기준", "criterion"),
                               "points": [(0.0, crit), (1.0, crit)],
                               "dash": True, "marker": False, "color": "#94a3b8"})
            h.append(wide(svg_plot(
                series, xlabel=t("진행률", "progress"),
                ylabel=t("최대 잔차", "worst residual"), ylog=True, width=880,
                height=340, title=t("케이스별 잔차", "residual, every case"))))
        series = _stack(lambda r: _fracx(_thin(
            _dp_deviation(r.get("dp_history") or [])[0], 140)))
        if series:
            h.append(wide(svg_plot(
                series, xlabel=t("진행률", "progress"),
                ylabel=t("자기 평균 대비 [%]", "from own mean [%]"),
                width=880, height=340, title=t("Δp 정착", "Δp settling"))))
            h.append('<p class="muted">%s</p>' % t(
                "0 에 가까이 붙어 평평해진 곡선은 정착한 것이고, 끝에서 "
                "기울어져 있으면 아직 이동 중, 규칙적으로 진동하면 와류 "
                "이탈입니다 — 비정상 해석에서는 정상입니다.",
                "A curve that flattens onto 0 has settled; one still sloping "
                "at the right-hand edge has not; a regular oscillation is "
                "vortex shedding, which on a transient run is what should "
                "happen."))
    h.extend(_ov_result(study, ko, t))
    h.append('<div class="foot">study.py · %s</div>'
             % _esc(time.strftime("%Y-%m-%d %H:%M:%S")))
    return "\n".join(h)


def _ov_result(study, ko, t):
    """The one chart the stage exists to produce, at the bottom of the wall.

    The panels above say whether the runs are trustworthy; this says what they
    measured.  Both on one screen, because the second is worth nothing without
    the first and reading them in different places is how that gets forgotten.
    """
    h = []
    try:
        an = (mesh_analysis(study) if study.kind == "mesh"
              else sweep_analysis(study))
    except Exception:                                       # noqa: BLE001
        return h
    if study.kind == "mesh":
        for lad in an.get("ladders", []):
            levels = lad.get("levels") or []
            if not lad.get("phi_ext") or len(levels) < 2:
                continue
            h.append("<h2>%s · %s</h2>"
                     % (t("결과", "the result"), _esc(lad["ladder"])))
            h.append('<div class="ovwide">%s</div>' % svg_plot(
                _grid_series(lad, levels, t), width=880, height=340,
                xlabel="h [m]", ylabel="Δp_bundle [Pa]",
                title=t("격자 수렴", "grid convergence")))
            if lad.get("chosen"):
                h.append('<p class="muted">%s</p>'
                         % t("선택된 격자: <b>%s</b>" % _esc(lad["chosen"]),
                             "the mesh to use: <b>%s</b>" % _esc(lad["chosen"])))
        return h
    rows = an.get("rows") or []
    if not rows:
        return h
    keys = an.get("keys") or []
    names = {k: CR.BY_KEY[k]["name"] for k in keys}
    by_pitch = {}
    for r in rows:
        by_pitch.setdefault((r["XT"], r["XL"]), []).append(r)
    drawn = []
    for key in sorted(by_pitch):
        rs = sorted(by_pitch[key], key=lambda r: r["Re"])
        if len(rs) >= 2:
            drawn.append('<div class="ovwide">%s</div>' % svg_plot(
                _eu_re_series(rs, keys, names), width=880, height=340,
                xlabel="Re_max", ylabel="Eu per row", xlog=True, ylog=True,
                title="X_T = %g, X_L = %g" % key))
    if drawn:
        h.append("<h2>%s</h2>" % t("결과", "the result"))
        h.extend(drawn)
    return h


def report_body(study, lang="ko"):
    return (mesh_report_body(study, lang) if study.kind == "mesh"
            else sweep_report_body(study, lang))


#  The histories, as plain text.
#
#  A .cas.h5 needs Fluent and a licence to reopen, and Fluent does not put the
#  residual history in it anyway: reading a case back gives you the converged
#  field, not the road to it.  The record does carry both histories, so the
#  plot can always be redrawn - in the app, or in whatever else the user
#  plots with.  That second one needs a file format that is not ours.
def history_csv(row):
    """One case's residual and Dp histories as a CSV table.

    Both are indexed by the same counter - iteration for a steady run, time
    step for a transient one - so they go in one table joined on it, with a
    blank where one of them has no sample at that index.
    """
    res = list(row.get("residual_history") or [])
    dp = list(row.get("dp_history") or [])
    eqs = []
    for r in res:
        for k in r:
            if k != "iter" and k not in eqs:
                eqs.append(k)
    dpk = []
    for m in dp:
        for k in m:
            if k != "iter" and k not in dpk:
                dpk.append(k)
    by = {}
    for r in res:
        by.setdefault(int(r.get("iter", 0)), {}).update(
            {"res:" + k: r[k] for k in eqs if r.get(k) is not None})
    for m in dp:
        by.setdefault(int(m.get("iter", 0)), {}).update(
            {k: m[k] for k in dpk if m.get(k) is not None})
    cols = ["res:" + e for e in eqs] + list(dpk)
    head = ["time_step" if row.get("transient") else "iter"] + cols
    out = ["# case: %s" % row.get("id", "?"),
           "# %s" % ("transient - the counter is the time step"
                     if row.get("transient")
                     else "steady - the counter is the iteration"),
           ",".join(head)]
    for n in sorted(by):
        r = by[n]
        out.append(",".join([str(n)] + ["" if r.get(c) is None
                                        else repr(float(r[c]))
                                        for c in cols]))
    return "\n".join(out) + "\n"


def write_histories(study, path=None):
    """Every recorded case's histories, one CSV each, into <study>/history/.

    Returns the paths written.  Cases with no history are skipped rather than
    written empty, so what is on disk is what actually ran.
    """
    out_dir = path or os.path.join(study.dir, "history")
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for row in study.results():
        if not (row.get("residual_history") or row.get("dp_history")):
            continue
        name = re.sub(r"[^A-Za-z0-9._-]", "_", str(row.get("id") or "case"))
        p = os.path.join(out_dir, name + ".csv")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(history_csv(row))
        written.append(p)
    return written


def write_report(study, lang="ko", path=None, overview=False):
    body = overview_body(study, lang) if overview else report_body(study, lang)
    path = path or os.path.join(
        study.dir, "overview.html" if overview else "report.html")
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
        force = "--force" in argv
        for s in make_campaign(level):
            try:
                path = s.save(force=force)
            except ValueError as exc:
                #  named and skipped, not silently overwritten and not fatal:
                #  a campaign half of which has run is the normal state, and
                #  the studies that have NOT run should still be rebuilt
                print("skipped %s: %s" % (s.name, exc))
                continue
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
