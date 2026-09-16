# -*- coding: utf-8 -*-
"""The correlation library - stage 1 of the pressure-drop study.

WHY THIS FILE EXISTS
--------------------
The point of the study is to put CFD against the published correlations and,
in the end, to fit a better one.  That only means anything if "the published
correlation" is a fixed, citable object rather than a formula somebody typed
in once.  So every correlation lives here, exactly once, with

    * where it comes from - author, year, publication, and which equation,
    * what it is allowed to be asked - Re, X_T, X_L, row count,
    * what it actually returns, in ONE currency (see below),
    * and its status: is the formula here, or only the reference to it?

Nothing in here reads a case, launches anything or knows what a mesh is.  It
takes numbers and returns numbers, which is what makes it testable.

ONE CURRENCY
------------
Sources disagree about what "the friction factor" means: Jakob's f is Fanning-
like and appears as `dp = 2 f N rho u_max^2`, Zukauskas' f is read off a chart
and appears as `dp = N chi f (rho u_max^2 / 2)`.  Comparing the f's directly is
meaningless.  Everything here is therefore converted to the EULER NUMBER PER
ROW,

    Eu_row = dp_bundle / (N_rows * 0.5 * rho * u_max^2)

so Jakob's Eu_row is 4f and Zukauskas' is chi*f, and the two can be put on one
axis.  dp follows from Eu_row and never the other way round.

u_max IS NOT ASSUMED TO BE THE TRANSVERSE GAP
---------------------------------------------
Every correlation here is written on the maximum (gap) velocity, and in a
staggered bank the narrowest gap is not always the transverse one: a rod
sitting between two rods of the row behind may see a tighter DIAGONAL passage,
and then the flow accelerates there instead.  `flow_state` tests both and says
which governed.  Getting this wrong moves both Re and dp (~u_max^2), so it is
decided in one place for every correlation rather than inside each of them.

WHAT IS NOT HERE, AND WHY
-------------------------
Several standard correlations are listed with status "needs-source": the
reference is recorded, the formula is not.  That is deliberate.  A correlation
coefficient written down from memory is worse than an absent one - it runs, it
produces plausible numbers, and nothing ever flags it.  An entry moves from
"needs-source" to "encoded" when its equations have been transcribed from the
publication with the source in hand, not before.
"""

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import zukauskas_charts as ZC          # noqa: E402


# =============================================================================
#  THE FLOW STATE EVERY CORRELATION IS ASKED ABOUT
# =============================================================================
def flow_state(D, ST, SL, n_rows, u_in, rho, nu, staggered, dT=None, dL=None):
    """Reduce a bundle and an approach velocity to what a correlation wants.

    Lengths in metres, `u_in` the mean approach (superficial) velocity in m/s,
    `rho` in kg/m3, `nu` the kinematic viscosity in m2/s.

    `dT`/`dL` are the width and length of the rod's FOOTPRINT in the section.
    For a straight rod both are D.  An inclined tube cuts the section as an
    ellipse, so they differ, and passing them keeps the helical case on the
    same code path instead of a second copy of it.
    """
    dT = D if dT is None else dT
    dL = D if dL is None else dL
    XT, XL = ST / dT, SL / dL
    gapT = ST - dT
    SD = math.hypot(SL, ST / 2.0)
    #  two diagonal passages per transverse pitch, so the area-equivalent
    #  diagonal gap is twice the geometric one - that is the standard test
    gapD = 2.0 * (SD - dT)
    diagonal = bool(staggered and 0.0 < gapD < gapT)
    gap = gapD if diagonal else gapT
    umax = u_in * (ST / gap if gap > 0 else float("inf"))
    return {
        "D": D, "ST": ST, "SL": SL, "dT": dT, "dL": dL,
        "n_rows": int(n_rows), "u_in": u_in, "rho": rho, "nu": nu,
        "staggered": bool(staggered),
        "XT": XT, "XL": XL, "SD": SD,
        "gapT": gapT, "gapD": gapD, "gap": gap, "diagonal": diagonal,
        "umax": umax,
        "Re": umax * D / nu,
        "Re_in": u_in * D / nu,
        "q": 0.5 * rho * umax * umax,          # dynamic head on u_max
    }


def dp_from_eu(eu_row, st):
    """Eu per row -> bundle pressure drop, the only conversion there is."""
    return eu_row * st["n_rows"] * st["q"]


# =============================================================================
#  THE CORRELATIONS
# =============================================================================
def _jakob(st):
    """Jakob (1938), in the form Holman gives it.

        in-line     f = [0.044 + 0.08 X_L / (X_T - 1)^(0.43 + 1.13/X_L)] Re^-0.15
        staggered   f = [0.25 + 0.118 / (X_T - 1)^1.08] Re^-0.16
        dp = 2 f N rho u_max^2           ->   Eu_row = 4 f

    X_T = S_T/D, X_L = S_L/D, Re on u_max.  This is the correlation the
    Geometry tab has always shown; `--check` holds the two implementations
    against each other so they cannot drift.
    """
    XT, XL, Re = st["XT"], st["XL"], max(st["Re"], 1.0)
    if st["staggered"]:
        f = (0.25 + 0.118 / max(XT - 1.0, 1e-3) ** 1.08) * Re ** -0.16
    else:
        f = (0.044 + 0.08 * XL / max(XT - 1.0, 1e-3) ** (0.43 + 1.13 / XL)) \
            * Re ** -0.15
    return 4.0 * f, {"f": f, "f_convention": "dp = 2 f N rho u_max^2"}


def _zukauskas(st):
    """Zukauskas (1972), off the digitised charts - see zukauskas_charts.py.

        dp = N chi f (rho u_max^2 / 2)   ->   Eu_row = chi f
    """
    f, chi, clamped = ZC.friction(st["staggered"], st["XT"], st["XL"],
                                  max(st["Re"], 1.0))
    extra = {"f": f, "chi": chi,
             "f_convention": "dp = N chi f (rho u_max^2 / 2)"}
    if clamped:
        extra["clamped"] = [
            {"axis": a, "asked": v, "lo": lo, "hi": hi} for a, v, lo, hi in clamped]
    return chi * f, extra


def _duct(st, Dh, L):
    """Plain duct friction, for the inlet and outlet boxes.

    Not a bundle correlation - it is the other half of the bookkeeping.  The
    inlet-minus-outlet pressure drop of the CFD spans the boxes as well as the
    bundle, and comparing THAT against a bundle correlation is comparing two
    different quantities.  Either measure the bundle alone (two planes on its
    faces, which is what the report does) or add this to the correlation.

    Darcy f: 64/Re laminar, Blasius 0.316 Re^-0.25 turbulent.
    """
    Re = max(st["u_in"] * Dh / st["nu"], 1.0)
    fD = 64.0 / Re if Re < 2300.0 else 0.316 * Re ** -0.25
    grad = fD / max(Dh, 1e-12) * 0.5 * st["rho"] * st["u_in"] ** 2
    return {"Re": Re, "fD": fD, "grad": grad, "dp": grad * L}


#  ---------------------------------------------------------------------------
#  The registry.  `status`:
#     encoded      - the equations are here, transcribed from the source
#     needs-source - the reference is here, the equations are not.  It will
#                    not be evaluated and will not silently return a number.
#  ---------------------------------------------------------------------------
CORRELATIONS = [
    {
        "key": "jakob",
        "name": "Jakob (1938)",
        "status": "encoded",
        "arrangement": "both",
        "fn": _jakob,
        "source": "M. Jakob, Trans. ASME 60 (1938) 384. In the form given by "
                  "J. P. Holman, Heat Transfer, McGraw-Hill - the in-line and "
                  "staggered friction factors, with dp = 2 f N rho u_max^2.",
        "valid": {"Re": [2000.0, 40000.0], "XT": [1.25, 3.0], "XL": [1.25, 3.0],
                  "n_rows": [10, None]},
        "valid_note": "Holman quotes the correlation for 2000 < Re_max < 40 000. "
                      "The pitch range is the range over which the source data "
                      "were taken and is the weaker of the two limits.",
        "note": "Closed form, and the one the Geometry tab has always shown. "
                "Held against Zukauskas' charts over a square-pitch grid in "
                "the stage-1 report - not asserted here, computed there. The "
                "short version: inside Jakob's own quoted Re range and for "
                "X >= 1.5 the two agree to about 20 % either way; at X = 1.25, "
                "or outside 2e3 < Re < 4e4, they part company by as much as a "
                "factor of 2.2. The tight-pitch corner is where a CFD point is "
                "worth the most.",
        "valid_note_ko": "Holman이 인용한 적용 범위는 2000 < Re_max < 40 "
                         "000 입니다. 피치 범위는 원 데이터가 취해진 "
                         "범위이며, 두 제한 중 더 느슨한 쪽입니다.",
        "note_ko": "닫힌 형태이고, 형상 탭이 계속 보여 온 바로 그 "
                   "상관식입니다. Zukauskas 차트와의 대조는 여기서 "
                   "주장하지 않고 1단계 보고서에서 계산합니다. 요약하면: "
                   "Jakob 자신이 제시한 Re 범위 안에서 X ≥ 1.5 이면 두 "
                   "상관식은 서로 ±20 % 안에 듭니다. 반면 X = 1.25 이거나 "
                   "Re가 2e3–4e4 밖이면 최대 2.2배까지 갈라집니다. 그 "
                   "촘촘한 피치 구석이 CFD 한 점의 가치가 가장 큰 "
                   "곳입니다.",
    },
    {
        "key": "zukauskas",
        "name": "Zukauskas (1972)",
        "status": "encoded",
        "arrangement": "both",
        "fn": _zukauskas,
        "source": "A. Zukauskas, Advances in Heat Transfer 8 (1972) 93-160; "
                  "charts as reprinted in Bergman, Lavine, Incropera & DeWitt, "
                  "Introduction to Heat Transfer, 6th ed., figs. 7.14 / 7.15. "
                  "Read here through the digitisation published in the MIT-"
                  "licensed `ht` library, copied into zukauskas_charts.py.",
        "valid": {"Re": [100.0, 1.0e6], "XT": [1.25, 2.5], "XL": [1.25, 2.5],
                  "n_rows": [10, None]},
        "valid_note": "The four fits carry their own domains, read off their "
                      "knots at import; friction() clamps to them and says so. "
                      "The pitch ratios are drawn only for 1.25 to 2.5.",
        "note": "A chart, not a formula: the chain is chart -> digitiser -> "
                "spline. Where the chart says chi is exactly 1 (square pitch) "
                "the fit returns 1.00 to 1.08, and that 8 % is the honest "
                "floor on what a reading off it is worth. Reproduces both of "
                "`ht`'s own documented examples to the digit, which is what "
                "says the tables came across intact.",
        "valid_note_ko": "네 개의 피팅이 각자의 정의역을 knot에서 직접 "
                         "읽어 가지고 있습니다. friction()이 그 범위로 "
                         "잘라내고 잘라냈다는 사실을 함께 돌려줍니다. "
                         "피치비는 1.25–2.5 에서만 그려져 있습니다.",
        "note_ko": "수식이 아니라 차트입니다: 차트 → 디지타이즈 → "
                   "스플라인 순서로 들어왔습니다. 차트상 chi가 정확히 "
                   "1이어야 하는 정사각 피치에서 이 피팅은 1.00–1.08을 "
                   "내놓으며, 그 8 %가 차트 읽기의 정직한 하한 "
                   "오차입니다. `ht`가 문서에 실어 둔 두 예제를 "
                   "자릿수까지 재현하므로, 표가 손상 없이 옮겨졌다는 것은 "
                   "확인되었습니다.",
    },
    {
        "key": "gaddis-gnielinski",
        "name": "Gaddis & Gnielinski (1985) / VDI",
        "status": "needs-source",
        "arrangement": "both",
        "fn": None,
        "source": "E. S. Gaddis and V. Gnielinski, 'Pressure drop in cross flow "
                  "across tube bundles', International Chemical Engineering "
                  "25 (1) (1985) 1-15; also the VDI Heat Atlas, section L1.4.",
        "valid": {"Re": [1.0, 300000.0], "XT": [1.02, 3.0], "XL": [0.6, 3.0],
                  "n_rows": [5, None]},
        "valid_note": "Quoted range; not verified here.",
        "note": "The best-regarded closed form for this problem, and the one "
                "worth having next: a laminar and a turbulent term superposed, "
                "each a function of the pitch ratios, valid from creeping flow "
                "to Re 3e5, in-line and staggered. The equations are long and "
                "the coefficients are not reproduced in anything reachable "
                "from here, so they are NOT written down from memory. Supply "
                "the source and this becomes an encoded entry.",
        "valid_note_ko": "인용된 범위이며, 여기서 검증한 것은 아닙니다.",
        "note_ko": "이 문제에 대해 가장 신뢰받는 닫힌 형태이고, 다음으로 "
                   "확보할 가치가 가장 큽니다. 층류항과 난류항을 중첩하며 "
                   "각각 피치비의 함수이고, 기어가는 유동부터 Re 3e5까지, "
                   "정렬·엇갈림 모두에 유효합니다. 식이 길고 계수를 "
                   "여기서 확인할 수 있는 경로가 없으므로 기억에 의존해 "
                   "적지 않았습니다. 원문을 주시면 바로 encoded 항목이 "
                   "됩니다.",
    },
    {
        "key": "grimison",
        "name": "Grimison (1937)",
        "status": "needs-source",
        "arrangement": "both",
        "fn": None,
        "source": "E. D. Grimison, Trans. ASME 59 (1937) 583.",
        "valid": {"Re": [2000.0, 40000.0], "XT": [1.25, 3.0], "XL": [1.25, 3.0],
                  "n_rows": [10, None]},
        "valid_note": "Quoted range; not verified here.",
        "note": "The other classical chart set, and the data Jakob's form was "
                "fitted to. Tabulated, not closed form.",
        "valid_note_ko": "인용된 범위이며, 여기서 검증한 것은 아닙니다.",
        "note_ko": "또 하나의 고전적 차트이며, Jakob의 식이 맞춰진 바로 "
                   "그 데이터입니다. 표 형태이고 닫힌 형태가 아닙니다.",
    },
    {
        "key": "idelchik",
        "name": "Idelchik, Handbook of Hydraulic Resistance",
        "status": "needs-source",
        "arrangement": "both",
        "fn": None,
        "source": "I. E. Idelchik, Handbook of Hydraulic Resistance, 3rd ed., "
                  "section 8, diagrams for banks of tubes in cross flow.",
        "valid": {"Re": [3.0, 1.0e5], "XT": [1.1, 4.0], "XL": [1.1, 4.0],
                  "n_rows": [1, None]},
        "valid_note": "Quoted range; not verified here.",
        "note": "Diagram-based, covers the low-Re end better than the others "
                "and handles short banks explicitly.",
        "valid_note_ko": "인용된 범위이며, 여기서 검증한 것은 아닙니다.",
        "note_ko": "선도(diagram) 기반이고, 저 Re 영역을 다른 "
                   "상관식들보다 잘 다루며 짧은 다발을 명시적으로 "
                   "취급합니다.",
    },
    {
        "key": "esdu-79034",
        "name": "ESDU 79034",
        "status": "needs-source",
        "arrangement": "both",
        "fn": None,
        "source": "ESDU 79034, 'Crossflow pressure loss over banks of plain "
                  "tubes in square and equilateral triangular arrays'.",
        "valid": {"Re": [10.0, 1.0e6], "XT": [1.2, 4.0], "XL": [1.2, 4.0],
                  "n_rows": [4, None]},
        "valid_note": "Quoted range; not verified here.",
        "note": "The engineering-practice standard, behind a licence. Listed so "
                "that it is visibly absent rather than invisibly missing.",
        "valid_note_ko": "인용된 범위이며, 여기서 검증한 것은 아닙니다.",
        "note_ko": "실무 표준이지만 라이선스가 필요합니다. 없다는 사실이 보이도록 목록에 넣어 둡니다.",
    },
]

BY_KEY = {c["key"]: c for c in CORRELATIONS}
ENCODED = [c["key"] for c in CORRELATIONS if c["status"] == "encoded"]


# =============================================================================
#  EVALUATION
# =============================================================================
def in_range(c, st):
    """Which of a correlation's stated limits the state falls outside."""
    out = []
    v = c.get("valid") or {}
    for key, label, value in (("Re", "Re_max", st["Re"]),
                              ("XT", "X_T", st["XT"]),
                              ("XL", "X_L", st["XL"]),
                              ("n_rows", "rows", st["n_rows"])):
        lim = v.get(key)
        if not lim:
            continue
        lo, hi = lim
        if lo is not None and value < lo:
            out.append("%s = %.4g is below the quoted %.4g" % (label, value, lo))
        if hi is not None and value > hi:
            out.append("%s = %.4g is above the quoted %.4g" % (label, value, hi))
    return out


def evaluate(key, st):
    """One correlation against one flow state.

    Always returns a record.  A correlation that cannot answer says so in
    `status` and leaves `eu_row` and `dp` as None - it never guesses.
    """
    c = BY_KEY.get(key)
    if c is None:
        raise KeyError("no correlation named %r" % key)
    rec = {"key": key, "name": c["name"], "status": c["status"],
           "eu_row": None, "dp": None, "outside": [], "extra": {}}
    if c["fn"] is None:
        rec["why"] = c["note"]
        return rec
    eu, extra = c["fn"](st)
    rec["eu_row"] = eu
    rec["dp"] = dp_from_eu(eu, st)
    rec["dp_per_row"] = rec["dp"] / max(1, st["n_rows"])
    rec["extra"] = extra
    rec["outside"] = in_range(c, st)
    return rec


def compare(st, keys=None):
    """Every correlation against one flow state, in registry order."""
    return [evaluate(k, st) for k in (keys or [c["key"] for c in CORRELATIONS])]


# =============================================================================
#  THE FORM STAGE 4 WILL FIT
# =============================================================================
#  Stage 4 is "a correlation of our own", and a correlation of our own needs a
#  shape before it needs data.  Every entry above reduces to Eu per row as a
#  power of Re times a function of the pitch ratios, so that is the shape:
#
#      Eu_row = C * Re^(-m) * (X_T - 1)^(-p) * X_L^(q)
#
#  four coefficients, fitted per arrangement, reducing to Jakob's staggered
#  branch at q = 0.  The helical case adds the lean: X_T and X_L are already
#  footprint ratios in flow_state, so the extra freedom a coil needs is one
#  more factor in the same product rather than a different function.
FIT_FORM = {
    "name": "power law in Re and the pitch ratios",
    "expr": "Eu_row = C * Re^(-m) * (X_T-1)^(-p) * X_L^q",
    "coefficients": ["C", "m", "p", "q"],
    "helical_extension": "Eu_row *= (cos alpha)^(-r); alpha the helix angle, "
                         "r the fifth coefficient. At alpha = 0 it is 1, so "
                         "the straight-rod fit is the coil fit's own limit.",
}


def fit_eu(coeff, st):
    """Evaluate the fitted form.  `coeff` is {C, m, p, q[, r]}."""
    eu = (coeff["C"] * max(st["Re"], 1.0) ** (-coeff["m"])
          * max(st["XT"] - 1.0, 1e-3) ** (-coeff["p"])
          * st["XL"] ** coeff["q"])
    r = coeff.get("r")
    if r:
        cos_a = st.get("cos_alpha", 1.0)
        eu *= max(cos_a, 1e-6) ** (-r)
    return eu


def fit_power_law(rows, with_helix=False):
    """Least squares for C, m, p, q on ln(Eu_row), in plain Python.

    ln Eu = ln C - m ln Re - p ln(X_T-1) + q ln X_L   is linear in the
    coefficients, so this is an ordinary normal-equation solve and needs no
    numpy.  `rows` are dicts with Re, XT, XL, eu_row (and cos_alpha if
    with_helix).  Returns the coefficients and the fit quality.
    """
    A, b = [], []
    for r in rows:
        if not r.get("eu_row"):
            continue
        row = [1.0, -math.log(max(r["Re"], 1.0)),
               -math.log(max(r["XT"] - 1.0, 1e-3)), math.log(r["XL"])]
        if with_helix:
            row.append(-math.log(max(r.get("cos_alpha", 1.0), 1e-6)))
        A.append(row)
        b.append(math.log(r["eu_row"]))
    n = len(A[0]) if A else 0
    if len(A) < n:
        raise ValueError("%d usable points cannot fit %d coefficients"
                         % (len(A), n))
    #  normal equations, Gauss-Jordan with partial pivoting
    M = [[sum(A[k][i] * A[k][j] for k in range(len(A))) for j in range(n)]
         + [sum(A[k][i] * b[k] for k in range(len(A)))] for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[p][i]) < 1e-14:
            raise ValueError("the design matrix is singular - the cases do not "
                             "vary every coefficient independently")
        M[i], M[p] = M[p], M[i]
        piv = M[i][i]
        M[i] = [v / piv for v in M[i]]
        for r in range(n):
            if r != i and M[r][i]:
                fac = M[r][i]
                M[r] = [v - fac * w for v, w in zip(M[r], M[i])]
    x = [M[i][n] for i in range(n)]
    names = ["C", "m", "p", "q"] + (["r"] if with_helix else [])
    coeff = dict(zip(names, x))
    coeff["C"] = math.exp(coeff["C"])
    #  how well it fits, in the currency a reader cares about
    errs = []
    for r in rows:
        if not r.get("eu_row"):
            continue
        pred = fit_eu(coeff, {"Re": r["Re"], "XT": r["XT"], "XL": r["XL"],
                              "cos_alpha": r.get("cos_alpha", 1.0)})
        errs.append(abs(pred - r["eu_row"]) / r["eu_row"])
    coeff["n_points"] = len(errs)
    coeff["mean_abs_error"] = sum(errs) / len(errs) if errs else None
    coeff["max_abs_error"] = max(errs) if errs else None
    return coeff


# =============================================================================
#  THE STAGE-1 REPORT
# =============================================================================
def _agreement_grid(keys=("jakob", "zukauskas")):
    """Every encoded correlation against every other, over the pitch/Re grid
    they share.  This is what makes the claims in `note` checkable instead of
    asserted: it is computed, here, from the same code the study will use."""
    out = []
    D, rho, nu = 0.010, 998.2, 1.004e-6
    for staggered in (False, True):
        for X in (1.25, 1.5, 2.0, 2.5):
            for Re in (1e3, 3e3, 1e4, 3e4, 1e5):
                #  the grid is stated in Re, so solve the approach velocity
                #  that lands on it: u_max = u_in * S_T/gap, Re = u_max D/nu
                probe = flow_state(D=D, ST=D * X, SL=D * X, n_rows=10, u_in=1.0,
                                   rho=rho, nu=nu, staggered=staggered)
                st = flow_state(D=D, ST=D * X, SL=D * X, n_rows=10,
                                u_in=Re * nu / D / probe["umax"],
                                rho=rho, nu=nu, staggered=staggered)
                out.append({"staggered": staggered, "X": X, "Re": st["Re"],
                            "eu": {k: evaluate(k, st)["eu_row"] for k in keys}})
    return out


def _fmt(v, n=4):
    return "—" if v is None else ("%.*g" % (n, v))


def report_body(lang="ko"):
    """The stage-1 report as an HTML fragment.

    A fragment, not a document: the tab drops it into the page's own `.paper`
    and the save button wraps it in the page's own stylesheet, so the report on
    screen and the report on disk are the same thing rendered twice.
    """
    ko = (lang == "ko")
    h = []

    def p(s):
        h.append(s)

    t = (lambda a, b: a if ko else b)
    #  the citations stay as published; the prose around them does not
    local = (lambda c, field: c.get(field + "_ko", c[field]) if ko else c[field])

    p('<h1>%s</h1>' % t("1단계 · 상관식 수집 및 정리",
                        "Stage 1 - the correlations, collected"))
    p('<div class="sub">%s</div>' % t(
        "직선 rod 다발 횡유동 압력강하 상관식 · 출처, 적용 범위, 사용 가능 여부",
        "Cross-flow pressure drop over a bank of straight rods - source, "
        "range, and whether it can actually be evaluated"))

    n_enc = len(ENCODED)
    p('<div class="tiles">')
    p('<div class="tile"><div class="k">%s</div><div class="v">%d</div></div>'
      % (t("수록된 상관식", "correlations listed"), len(CORRELATIONS)))
    p('<div class="tile"><div class="k">%s</div><div class="v">%d</div></div>'
      % (t("계산 가능", "evaluable here"), n_enc))
    p('<div class="tile"><div class="k">%s</div><div class="v">%d</div></div>'
      % (t("출처 필요", "awaiting the source"), len(CORRELATIONS) - n_enc))
    p("</div>")

    p("<h2>%s</h2>" % t("공통 정의", "The common definitions"))
    p("<p>%s</p>" % t(
        "상관식마다 마찰계수의 정의가 다르므로 (Jakob은 <code>dp = 2 f N rho u_max^2</code>, "
        "Zukauskas는 <code>dp = N chi f (rho u_max^2/2)</code>) f를 직접 비교하는 것은 "
        "의미가 없습니다. 이 라이브러리는 모두 <b>행당 Euler 수</b>로 환산합니다.",
        "Sources disagree about what f means - Jakob's appears as "
        "<code>dp = 2 f N rho u_max^2</code>, Zukauskas' as "
        "<code>dp = N chi f (rho u_max^2/2)</code> - so the f's are not "
        "comparable. Everything here is converted to the <b>Euler number per "
        "row</b>."))
    p("<table>")
    p("<tr><th>Eu<sub>row</sub></th><td class=\"n\">dp<sub>bundle</sub> / "
      "(N<sub>rows</sub> · &frac12; rho u<sub>max</sub><sup>2</sup>)</td></tr>")
    p("<tr><th>Re<sub>max</sub></th><td class=\"n\">u<sub>max</sub> D / nu</td></tr>")
    p("<tr><th>u<sub>max</sub></th><td class=\"n\">%s</td></tr>" % t(
        "가장 좁은 통로 기준. 엇갈림 배열에서는 횡방향 간극 (S_T − D)과 "
        "대각 간극 2(S_D − D)를 비교하여 좁은 쪽을 사용합니다.",
        "on the narrowest passage. In a staggered bank the transverse gap "
        "(S_T − D) is tested against the diagonal 2(S_D − D) and the tighter "
        "one governs."))
    p("<tr><th>X<sub>T</sub>, X<sub>L</sub></th><td class=\"n\">S_T/D, S_L/D</td></tr>")
    p("</table>")

    p("<h2>%s</h2>" % t("상관식 목록", "The list"))
    for c in CORRELATIONS:
        ok = c["status"] == "encoded"
        p('<h3>%s %s</h3>' % (
            c["name"],
            '' if ok else '<span class="flag">%s</span>'
            % t("출처 필요", "source needed")))
        p("<table>")
        p('<tr><th>%s</th><td>%s</td></tr>' % (t("출처", "source"), c["source"]))
        v = c["valid"]
        rng = "Re %s–%s · X_T %s–%s · X_L %s–%s · %s %s" % (
            _fmt(v["Re"][0]), _fmt(v["Re"][1]), _fmt(v["XT"][0]), _fmt(v["XT"][1]),
            _fmt(v["XL"][0]), _fmt(v["XL"][1]),
            t("행수 ≥", "rows >="), _fmt(v["n_rows"][0]))
        p('<tr><th>%s</th><td class="n">%s</td></tr>' % (t("적용 범위", "range"), rng))
        p('<tr><th>%s</th><td>%s</td></tr>' % (
            t("범위 주석", "on that range"), local(c, "valid_note")))
        p('<tr><th>%s</th><td>%s</td></tr>' % (t("비고", "note"),
                                               local(c, "note")))
        p("</table>")

    p("<h2>%s</h2>" % t("계산 가능한 두 상관식의 상호 대조",
                        "The two evaluable correlations, against each other"))
    p("<p>%s</p>" % t(
        "아래 표는 이 파일의 코드로 <i>지금</i> 계산한 값입니다. 정사각 피치 "
        "(X_T = X_L = X), 10행, 물 기준. 두 상관식은 서로 독립적인 출처이므로 "
        "이 차이가 '기존 상관식'이라는 기준선 자체의 폭입니다.",
        "Computed by this file's own code, now. Square pitch (X_T = X_L = X), "
        "10 rows, water. The two come from independent sources, so the spread "
        "between them is the width of the baseline itself."))
    grid = _agreement_grid()
    for staggered in (False, True):
        p("<h3>%s</h3>" % (t("엇갈림 (staggered)", "staggered") if staggered
                           else t("정렬 (in-line)", "in-line")))
        p("<table><tr><th>X</th><th>Re<sub>max</sub></th>"
          "<th>Eu<sub>row</sub> Jakob</th><th>Eu<sub>row</sub> Zukauskas</th>"
          "<th>Jakob / Zukauskas</th></tr>")
        for g in grid:
            if g["staggered"] != staggered:
                continue
            a, b = g["eu"]["jakob"], g["eu"]["zukauskas"]
            p('<tr><td class="n">%.2f</td><td class="n">%.3g</td>'
              '<td class="n">%.4f</td><td class="n">%.4f</td>'
              '<td class="n">%.2f</td></tr>' % (g["X"], g["Re"], a, b, a / b))
        p("</table>")
    #  The summary, restricted as well as unrestricted.  The unrestricted
    #  spread is dominated by the corners of the grid that Jakob was never
    #  quoted for, and reporting only that would make the baseline look
    #  useless; reporting only the restricted one would make it look better
    #  than it is.  Both, and the restriction named.
    p("<h3>%s</h3>" % t("요약", "In summary"))
    p("<table><tr><th>%s</th><th>%s</th><th>%s</th></tr>"
      % (t("범위", "over"), t("정렬", "in-line"), t("엇갈림", "staggered")))

    def band(pred):
        cells = []
        for staggered in (False, True):
            rs = [g["eu"]["jakob"] / g["eu"]["zukauskas"] for g in grid
                  if g["staggered"] == staggered and pred(g)]
            cells.append('<td class="n">%.2f – %.2f (%s %.2f)</td>'
                         % (min(rs), max(rs), t("평균", "mean"),
                            sum(rs) / len(rs)))
        return "".join(cells)

    p("<tr><th>%s</th>%s</tr>" % (
        t("격자 전체", "the whole grid"), band(lambda g: True)))
    p("<tr><th>%s</th>%s</tr>" % (
        t("Jakob 적용 Re (2e3–4e4)", "Jakob's quoted Re, 2e3–4e4"),
        band(lambda g: 2000.0 <= g["Re"] <= 40000.0)))
    p("<tr><th>%s</th>%s</tr>" % (
        t("그리고 X ≥ 1.5", "... and X >= 1.5"),
        band(lambda g: 2000.0 <= g["Re"] <= 40000.0 and g["X"] >= 1.5)))
    p("</table>")
    p("<p>%s</p>" % t(
        "읽는 법: 두 상관식은 Jakob이 제시한 Re 범위 안에서, 피치비가 1.5 이상이면 "
        "서로 ±20 % 안에 듭니다. 즉 그 영역에서는 '기존 상관식'이 하나의 값에 가깝고, "
        "CFD가 그 밖으로 크게 벗어나면 CFD를 의심해야 합니다. 반대로 <b>X = 1.25의 "
        "촘촘한 피치</b>와 <b>Re가 범위 밖</b>인 곳에서는 두 상관식이 최대 2.2배까지 "
        "갈라집니다 — 기준선 자체가 없는 영역이고, CFD 한 점의 가치가 가장 큰 곳입니다. "
        "3단계의 해석 조건은 이 두 영역을 모두 지나도록 잡습니다.",
        "How to read it: inside Jakob's quoted Re range and for pitch ratios "
        "of 1.5 and up, the two agree to within about 20 % either way. There, "
        "'the existing correlation' is close to a single number and CFD that "
        "misses it by much is CFD to be suspicious of. At <b>X = 1.25</b> and "
        "<b>outside that Re range</b> they part by up to a factor of 2.2 - "
        "there is no baseline there, and that is where one CFD point is worth "
        "the most. The stage-3 matrix is laid out to cross both regions."))

    p("<h2>%s</h2>" % t("4단계에서 맞출 형태", "The form stage 4 will fit"))
    p("<p><code>%s</code></p>" % FIT_FORM["expr"])
    p("<p>%s</p>" % t(
        "위 상관식들은 모두 '(Re의 거듭제곱) × (피치비의 함수)' 형태로 환원되므로, "
        "우리 상관식도 같은 형태에서 출발합니다. 나선 코일은 여기에 경사각 항 "
        "하나가 곱해지며, 경사각 0에서 직선 rod 결과로 정확히 돌아갑니다.",
        "Every correlation above reduces to a power of Re times a function of "
        "the pitch ratios, so ours starts from the same shape. The coil adds "
        "one factor in the helix angle, which is 1 at zero lean - so the "
        "straight-rod fit is the coil fit's own limit rather than a different "
        "correlation."))
    p("<p>%s</p>" % FIT_FORM["helical_extension"])

    p("<h2>%s</h2>" % t("이 다발에 대한 주의", "What to watch for in THIS bundle"))
    for s in (
        t("상관식들은 <b>완전발달</b> 상태, 즉 행수가 10 이상인 다발에 대해 "
          "제시되었습니다. 앞쪽 1–2 행은 압력강하가 다르므로, 행수가 적은 해석을 "
          "행당 값으로 비교할 때 이 차이가 그대로 나타납니다.",
          "Every correlation here is quoted for a FULLY DEVELOPED bank, "
          "meaning about ten rows or more. The first row or two do not drop "
          "the same pressure, so a short bundle compared per row carries that "
          "difference straight into the ratio."),
        t("상관식은 <b>측벽이 없는</b> 무한 폭 다발을 가정합니다. 이 해석 영역은 "
          "양쪽에 벽이 있고 (반쪽 rod 옵션 포함), 벽 마찰과 벽 근처 유로가 "
          "추가됩니다. 열 수가 적을수록 이 효과가 큽니다.",
          "The correlations are for a bank with no side walls. This domain has "
          "them - with or without the half rods - and they add wall friction "
          "and a passage of their own. The fewer the columns, the more it "
          "matters."),
        t("CFD의 입구−출구 압력차는 입·출구 박스를 <b>포함</b>합니다. 상관식은 "
          "다발만을 예측하므로, 다발 앞뒤 면에 평면을 세워 측정한 값과 비교해야 "
          "합니다 (보고서 탭이 이미 이렇게 합니다).",
          "The CFD's inlet-minus-outlet spans the inlet and outlet boxes as "
          "well. The correlation predicts the bundle alone, so the comparison "
          "has to be against two planes on the bundle faces - which is what "
          "the report tab already measures."),
        t("y+ 목표가 난류 모델과 맞아야 합니다. k-omega SST에 y+ 30을 쓰면 첫 셀이 "
          "완충층에 놓여 벽면 전단이 틀리고, 압력강하가 그만큼 틀립니다.",
          "The y+ target has to suit the turbulence model. k-omega SST with a "
          "target of 30 puts the first cell in the buffer layer, the wall "
          "shear comes out wrong, and the pressure drop with it."),
    ):
        p("<p class=\"muted\">· %s</p>" % s)

    p('<div class="foot">%s</div>' % t(
        "correlations.py 에서 생성 · 표의 모든 수치는 같은 코드로 계산됨",
        "generated by correlations.py - every number in the tables computed by "
        "the same code the study uses"))
    return "\n".join(h)


def as_json():
    """The registry over the wire, without the function objects."""
    out = []
    for c in CORRELATIONS:
        d = {k: v for k, v in c.items() if k != "fn"}
        d["evaluable"] = c["fn"] is not None
        out.append(d)
    return {"correlations": out, "encoded": ENCODED, "fit_form": FIT_FORM}


# =============================================================================
#  COMMAND LINE
# =============================================================================
def _check_against_js():
    """Jakob here against Jakob in the browser.

    The Geometry tab computes the same correlation in JavaScript and must go
    on doing so - it works with no server at all, which is half the point of
    the page.  Two implementations of one formula is exactly the arrangement
    that drifts, so they are held against each other here the way the mesh and
    its Python twin are.
    """
    import subprocess
    import tempfile
    cases = []
    for staggered in (False, True):
        for XT in (1.25, 1.5, 2.0, 2.5):
            for XL in (1.25, 1.5, 2.0, 2.5):
                for u in (0.1, 0.5, 2.0):
                    cases.append({"staggered": staggered, "XT": XT, "XL": XL,
                                  "u": u})
    js = """
'use strict';
const cases = %s;
const D = 0.010, rho = 998.2, nu = 1.004e-6, N = 10;
const out = cases.map(c => {
  const ST = D*c.XT, SL = D*c.XL, dT = D;
  const XT = ST/dT, XL = SL/D;
  const gapT = ST - dT;
  const SD = Math.hypot(SL, ST/2);
  const gapD = 2*(SD - dT);
  const diagonal = c.staggered && gapD > 0 && gapD < gapT;
  const gap = diagonal ? gapD : gapT;
  const umax = c.u*(gap > 0 ? ST/gap : 1);
  const ReMax = umax*D/nu;
  const f = c.staggered
    ? (0.25 + 0.118/Math.pow(Math.max(XT-1,1e-3), 1.08))
      * Math.pow(Math.max(ReMax,1), -0.16)
    : (0.044 + 0.08*XL/Math.pow(Math.max(XT-1,1e-3), 0.43+1.13/XL))
      * Math.pow(Math.max(ReMax,1), -0.15);
  return {umax, ReMax, eu: 4*f, dp: 2*f*N*rho*umax*umax};
});
console.log(JSON.stringify(out));
""" % json.dumps(cases)
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False) as fh:
        fh.write(js)
        path = fh.name
    try:
        raw = subprocess.check_output(["node", path]).decode("utf-8")
    finally:
        os.unlink(path)
    theirs = json.loads(raw)
    worst = 0.0
    worst_at = None
    for c, jsv in zip(cases, theirs):
        st = flow_state(D=0.010, ST=0.010 * c["XT"], SL=0.010 * c["XL"],
                        n_rows=10, u_in=c["u"], rho=998.2, nu=1.004e-6,
                        staggered=c["staggered"])
        rec = evaluate("jakob", st)
        for what, a, b in (("u_max", st["umax"], jsv["umax"]),
                           ("Re", st["Re"], jsv["ReMax"]),
                           ("Eu_row", rec["eu_row"], jsv["eu"]),
                           ("dp", rec["dp"], jsv["dp"])):
            e = abs(a - b) / max(abs(b), 1e-30)
            if e > worst:
                worst, worst_at = e, (what, c, a, b)
    print("Jakob, python vs the browser: %d cases x 4 quantities" % len(cases))
    print("  worst relative difference %.3e%s"
          % (worst, "" if worst < 1e-12 else "   at %r" % (worst_at,)))
    ok = worst < 1e-12
    print("  " + ("the two implementations agree" if ok else "THEY DISAGREE"))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--check" in argv:
        return _check_against_js()
    if "--json" in argv:
        print(json.dumps(as_json(), ensure_ascii=False, indent=1))
        return 0
    if "--report" in argv:
        i = argv.index("--report")
        out = argv[i + 1] if len(argv) > i + 1 and not argv[i + 1].startswith("-") \
            else os.path.join(HERE, "studies", "stage1-correlations", "report.html")
        lang = "en" if "--en" in argv else "ko"
        os.makedirs(os.path.dirname(out), exist_ok=True)
        body = report_body(lang)
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(_standalone(body, lang))
        print("wrote %s (%.1f kB)" % (out, os.path.getsize(out) / 1e3))
        return 0
    #  default: the registry, and one worked state
    for c in CORRELATIONS:
        print("%-20s %-13s %s" % (c["key"], c["status"], c["name"]))
    st = flow_state(D=0.010, ST=0.020, SL=0.020, n_rows=4, u_in=0.5,
                    rho=998.2, nu=1.004e-6, staggered=False)
    print("\nexample: D 10 mm, S_T = S_L = 20 mm, 4 rows, water at 0.5 m/s")
    print("  u_max %.4f m/s   Re_max %.4g   X_T %.2f   X_L %.2f"
          % (st["umax"], st["Re"], st["XT"], st["XL"]))
    for rec in compare(st):
        if rec["eu_row"] is None:
            print("  %-20s -  %s" % (rec["key"], rec["status"]))
        else:
            print("  %-20s Eu_row %.4f   dp %.1f Pa%s"
                  % (rec["key"], rec["eu_row"], rec["dp"],
                     "   OUTSIDE: " + "; ".join(rec["outside"])
                     if rec["outside"] else ""))
    return 0


def _standalone(body, lang):
    """The report as a file that opens anywhere.

    The page's own stylesheet is not reachable from the command line, so this
    carries a small one of its own.  When the tab saves the same report it uses
    the page's, exactly as the Report tab does.
    """
    css = """
body{margin:0;background:#eef1f5;padding:26px 16px;
  font-family:system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans KR",sans-serif}
.paper{max-width:900px;margin:0 auto;background:#fff;padding:38px 44px;
  border-radius:10px;box-shadow:0 2px 18px rgba(20,28,40,.10);color:#181f2a;
  font-size:13.5px;line-height:1.62}
.paper h1{font-size:23px;margin:0 0 4px}
.paper h2{font-size:16px;margin:30px 0 10px;padding-bottom:5px;
  border-bottom:1px solid #dde3ea}
.paper h3{font-size:13.5px;margin:20px 0 7px;color:#33415a}
.paper .sub{color:#5a6678;font-size:12px;margin-bottom:18px}
.paper table{border-collapse:collapse;width:100%;margin:9px 0;font-size:12.5px}
.paper th,.paper td{border:1px solid #e2e7ee;padding:5px 9px;text-align:left;
  vertical-align:top}
.paper th{background:#f6f8fb;font-weight:640;width:24%;color:#33415a}
.paper td.n{font-variant-numeric:tabular-nums}
.paper tr>th:first-child:empty{width:auto}
.paper .tiles{display:flex;gap:10px;margin:16px 0 4px;flex-wrap:wrap}
.paper .tile{flex:1 1 130px;background:#f6f8fb;border:1px solid #e2e7ee;
  border-radius:7px;padding:10px 13px}
.paper .tile .k{font-size:11px;color:#5a6678}
.paper .tile .v{font-size:19px;font-weight:660;font-variant-numeric:tabular-nums}
.paper .muted{color:#5a6678;font-size:12px;margin:5px 0}
.paper .flag{background:#fdeaea;color:#9b1c1c;border-radius:4px;padding:1px 7px;
  font-size:11px;font-weight:640;margin-left:7px}
.paper code{background:#f2f5f9;border-radius:3px;padding:1px 5px;font-size:12px}
.paper .foot{margin-top:30px;padding-top:12px;border-top:1px solid #dde3ea;
  color:#7a8595;font-size:11.5px}
@media print{body{background:#fff;padding:0}
  .paper{box-shadow:none;max-width:none;padding:0}}
"""
    return ('<!doctype html>\n<html lang="%s"><head><meta charset="utf-8">\n'
            '<title>Stage 1 - correlations</title>\n<style>%s</style>\n'
            '</head>\n<body>\n<div class="paper">\n%s\n</div>\n</body></html>\n'
            % (lang, css, body))


if __name__ == "__main__":
    raise SystemExit(main())
