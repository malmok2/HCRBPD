#!/usr/bin/env python3
"""Every checkable claim in docs/ , checked against the code.

A guide that has drifted from the program is worse than no guide: it is a
confident wrong answer.  So the numbers in docs/ are not typed by hand and
trusted - each one that the code can produce is produced here and compared.

    python3 tests/docs_claims.py

Needs no Fluent and no licence.  Meshes are built in-process, so the ladder
check takes a minute.  Pass --quick to skip it.
"""
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import app as APP                                   # noqa: E402
import correlations as CR                           # noqa: E402
import fluent_case as FC                            # noqa: E402
import mesh_explorer as ME                          # noqa: E402
import study as ST                                  # noqa: E402

FAILED = []


def ck(name, cond, got=""):
    print(("PASS  " if cond else "FAIL  ") + name
          + ("  [%s]" % (got,) if got != "" and not cond else ""))
    if not cond:
        FAILED.append(name)


def az_of(geometry, n_az):
    """Azimuthal cells actually built around one rod for a requested count."""
    c = APP.case_from(geometry, {"nAz": n_az, "nRows": 2, "nCols": 2, "nZ": 1})
    a = ME.quick_az(c)
    return a[1] if isinstance(a, tuple) else a


def docs_text():
    d = os.path.join(ROOT, "docs")
    return {f: io.open(os.path.join(d, f), encoding="utf-8").read()
            for f in sorted(os.listdir(d)) if f.endswith(".md")}


def main(argv):
    quick = "--quick" in argv
    docs = docs_text()
    blob = "\n".join(docs.values())

    print("-- docs/ itself")
    ck("an index plus six numbered documents", len(docs) == 7, sorted(docs))
    ck("no stray Cyrillic (a real typo once)",
       not [c for c in blob if 0x400 <= ord(c) <= 0x4FF])
    #  every relative link resolves
    bad = []
    for name, s in docs.items():
        for m in re.finditer(r"\]\(([^)#][^)]*)\)", s):
            t = m.group(1)
            if t.startswith(("http", "mailto")):
                continue
            q = os.path.normpath(os.path.join(ROOT, "docs", t.split("#")[0]))
            if not os.path.exists(q):
                bad.append((name, t))
    ck("every relative link in docs/ resolves", not bad, bad)

    print("\n-- 02: what nAz actually delivers")
    #  the two tables in 02-격자-직접-지정.md, value for value
    for n, want in ((12, 16), (14, 16), (20, 24), (29, 32), (36, 40),
                    (60, 64), (64, 64), (96, 96), (128, 128), (200, 200)):
        ck("in-line nAz %d -> %d" % (n, want), az_of("rod-inline", n) == want,
           az_of("rod-inline", n))
    for n, want in ((12, 12), (20, 24), (32, 36), (60, 62), (96, 100),
                    (128, 136)):
        ck("staggered nAz %d -> %d" % (n, want),
           az_of("rod-staggered", n) == want, az_of("rod-staggered", n))
    ck("in-line is 8*ceil(nAz/8), as the doc states the rule",
       all(az_of("rod-inline", n) == 8 * -(-n // 8) for n in range(8, 130, 3)))

    print("\n-- 02: the slider ranges the doc quotes")
    h = io.open(os.path.join(ROOT, "mesh_explorer.html"),
                encoding="utf-8").read()
    for pid, lo, hi in (("nAz", "12", "384"), ("nRad", "3", "120"),
                        ("firstLayer", "0.001", "1.5"), ("nZ", "1", "60"),
                        ("nxIn", "2", "240"), ("nxOut", "2", "320")):
        m = re.search(r'id="%s" min="([^"]+)" max="([^"]+)"' % pid, h)
        ck("%s slider %s-%s" % (pid, lo, hi),
           bool(m) and m.group(1) == lo and m.group(2) == hi,
           m.groups() if m else None)
        ck("%s range is quoted in the doc" % pid,
           ("%s –" % lo) in blob or ("%s – %s" % (lo, hi)) in blob
           or ("%s" % hi) in blob)

    print("\n-- 03: the solver settings the doc quotes")
    B = ST.BASE_SETTINGS
    ck("transient by default", B["general"]["steady"] is False)
    ck("k-omega SST", (B["turbulence"]["viscous"],
                       B["turbulence"]["k_omega_variant"]) == ("k-omega", "sst"))
    ck("water 998.2 kg/m3, 1.003e-3 Pa s",
       (B["material"]["density"], B["material"]["viscosity"])
       == (998.2, 0.001003))
    ck("inlet turbulence 5 % / ratio 10",
       (B["inlet"]["intensity"], B["inlet"]["visc_ratio"]) == (5.0, 10.0))
    ck("outlet prevents reverse flow", B["outlet"]["prevent_reverse_flow"] is True)
    ck("rods are walls", B["zones"]["wall_rods"] == "wall")
    ck("all four other boundaries are symmetry",
       all(B["zones"][k] == "symmetry" for k in
           ("wall_side_y0", "wall_side_ymax", "wall_bottom", "wall_top")))
    ck("coupled pressure-velocity", B["methods"]["flow_scheme"] == "Coupled")
    ck("second order throughout",
       all(B["methods"][k].startswith("second-order")
           for k in ("pressure", "momentum", "turb")))
    ck("criterion 1e-5, monitor every 20",
       (B["run"]["residual_criterion"], B["run"]["monitor_every"]) == (1e-5, 20))
    ck("hybrid initialisation", B["run"]["init_method"] == "hybrid")

    print("\n-- 03: the two clocks")
    S = ST.SHEDDING
    ck("St = 0.2", S["strouhal"] == 0.2)
    ck("25 steps per shedding period", S["steps_per_period"] == 25)
    ck("a floor of 20 periods", S["periods"] == 20)
    ck("15 periods averaged", S["average_periods"] == 15)
    ck("2 bundle flow-throughs flushed", S["flush_flowthroughs"] == 2.0)
    ck("12 inner iterations per step", S["max_iter_per_step"] == 12)
    ck("2nd-order implicit is the first choice",
       FC.time_schemes(FC.merge_settings(
           {"general": {"steady": False}}))[0] == "unsteady-2nd-order")
    ck("steady has exactly one spelling",
       FC.time_schemes(FC.merge_settings({"general": {"steady": True}}))
       == ["steady"])
    #  the flush really is 2 sweeps for every case the campaign defines
    flush = []
    for name in ("rod-inline-mesh", "rod-staggered-mesh",
                 "rod-inline-sweep", "rod-staggered-sweep"):
        st = ST.Study.load(name)
        for c in st.cases:
            tc = c.get("transient") or {}
            if not tc.get("bundle_flowthrough_s"):
                continue
            r = st.settings_for(c["id"])["run"]
            flush.append((r["time_steps"] - r["average_last"])
                         * r["time_step"] / tc["bundle_flowthrough_s"])
    ck("all 76 cases carry a bundle flow-through time", len(flush) == 76,
       len(flush))
    ck("and every one flushes ~2 sweeps",
       flush and 1.99 <= min(flush) and max(flush) < 2.5,
       (min(flush), max(flush)) if flush else None)

    print("\n-- 03/04: what is kept, and how much")
    ck("residual history capped at 400", ST.RES_KEEP == 400)
    ck("Dp history capped at 2000 - a 500-step run is whole",
       ST.DP_KEEP == 2000 and len(ST._thin([{"iter": k} for k in range(500)],
                                           ST.DP_KEEP)) == 500)
    ck("cells-per-core turnover quoted as 50 000",
       ST.CELLS_PER_CORE == 50000)

    print("\n-- 05/06: the correlation band")
    ck("four correlations are evaluable",
       set(CR.ENCODED) == {"jakob", "zukauskas", "gunter-shaw", "shen-2024"},
       CR.ENCODED)
    need = sorted(k for k, v in CR.BY_KEY.items()
                  if v.get("status") == "needs-source")
    ck("four still need a source", len(need) == 4, need)
    ck("Gaddis-Gnielinski is one of them", "gaddis-gnielinski" in need)
    s2 = ST.Study.load("rod-inline-mesh")
    _case, fs = ST.case_state("rod-inline", s2.params_for("main-L5"),
                              s2.settings_for("main-L5"))
    for key, want in (("jakob", 0.318), ("zukauskas", 0.333),
                      ("gunter-shaw", 0.405), ("shen-2024", 0.255)):
        got = CR.evaluate(key, fs)["eu_row"]
        ck("%s gives Eu_row %.3f at X=1.5, Re=1e4" % (key, want),
           abs(got - want) < 0.001, got)

    print("\n-- 05: the ladder defect the docs describe")
    ids = ["main-L%d" % i for i in range(1, 6)]
    azs = [az_of("rod-inline", s2.params_for(c)["nAz"]) for c in ids]
    ck("the old ladder's azimuthal counts are 16/24/32/48/64",
       azs == [16, 24, 32, 48, 64], azs)
    fac = ST.refinement_factors(s2, ids)
    ck("the wall-normal direction never refines",
       all(abs(x - 1.0) < 1e-9 for x in fac["wall-normal"]), fac["wall-normal"])
    mean = dict(ST.anisotropy(s2, ids))
    ck("the ladder is detected as anisotropic", bool(mean))
    ck("the fastest direction is ~1.46 per level",
       1.45 < mean.get("__lead__", 0) < 1.47, mean.get("__lead__"))
    ck("radial ~1.20 is named as lagging",
       1.19 < mean.get("radial", 0) < 1.21, mean.get("radial"))

    if quick:
        print("\n-- 02: the worked ladder  (skipped, --quick)")
    else:
        print("\n-- 02: the worked ladder in the doc, meshed and counted")
        p0 = s2.params_for("main-L5")
        for n_az, first, nx_in, want in ((128, 0.0318, 40, 222208),
                                         (168, 0.0245, 52, 315840),
                                         (216, 0.0188, 68, 459648)):
            p = dict(p0, nAz=n_az, firstLayer=first, nxIn=nx_in,
                     nxOut=int(round(nx_in * 1.7)))
            c = APP.case_from("rod-inline", p)
            p["nRad"] = ST.n_rad_for(c.clearance, first, 1.15)
            m = ME.Mesh(APP.case_from("rod-inline", p))
            m.build()
            ck("nAz %d -> %s cells" % (n_az, "{:,}".format(want)),
               len(m.hexes) == want, len(m.hexes))

    print("\n%d failed" % len(FAILED))
    return 1 if FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
