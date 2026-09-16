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
def flow_state(D, ST, SL, n_rows, u_in, rho, nu, staggered, dT=None, dL=None,
               helix=0.0, coil_alternating=False):
    """Reduce a bundle and an approach velocity to what a correlation wants.

    Lengths in metres, `u_in` the mean approach (superficial) velocity in m/s,
    `rho` in kg/m3, `nu` the kinematic viscosity in m2/s.

    `dT`/`dL` are the width and length of the rod's FOOTPRINT in the section.
    For a straight rod both are D.  An inclined tube cuts the section as an
    ellipse, so they differ, and passing them keeps the helical case on the
    same code path instead of a second copy of it.

    `helix` is the lean in DEGREES and `coil_alternating` says whether adjacent
    radial layers are wound the opposite way.  Both are zero/False for a
    straight bundle and are only read by the correlations that were fitted on
    coils; they are here so that stage 4 does not need a second flow state.

    TWO PITCH RATIOS, ON PURPOSE.  `XT`/`XL` divide by the FOOTPRINT, which is
    what sets the gap and therefore u_max.  `XT_d`/`XL_d` divide by the TUBE
    DIAMETER, which is what a published S/d means.  For a straight rod they
    are the same number; for a leaning tube they are not, and a correlation
    has to be asked with the one its author used.
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
    #  Gunter & Shaw's volumetric hydraulic diameter: four times the free
    #  volume over the friction surface.  Per lattice cell that is
    #  4(S_T S_L - pi D^2/4) / (pi D), and it reproduces the 0.1334 m the
    #  KAERI CHX paper states for D = 27.2 mm, a = 2.65, b = 1.75 to 6e-6 m.
    Dv = 4.0 * (ST * SL - math.pi * D * D / 4.0) / (math.pi * D)
    return {
        "D": D, "ST": ST, "SL": SL, "dT": dT, "dL": dL,
        "n_rows": int(n_rows), "u_in": u_in, "rho": rho, "nu": nu,
        "staggered": bool(staggered),
        "helix": float(helix), "coil_K": 1 if coil_alternating else 0,
        "cos_eps": math.cos(math.radians(helix * (1.0 - helix / 90.0))),
        "XT": XT, "XL": XL, "XT_d": ST / D, "XL_d": SL / D, "SD": SD, "Dv": Dv,
        "gapT": gapT, "gapD": gapD, "gap": gap, "diagonal": diagonal,
        "umax": umax,
        "Re": umax * D / nu,
        "Re_v": umax * Dv / nu,                # Gunter & Shaw's Reynolds number
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


def _gunter_shaw(st):
    """Gunter & Shaw (1945), on the volumetric hydraulic diameter.

        Re_v = u_max D_v / nu,      D_v = 4 (S_T S_L - pi D^2/4) / (pi D)
        f/2  = 90 / Re_v                      Re_v <= 200
        f/2  = 0.96 Re_v^-0.145               Re_v >  200
        dp   = (f/2) (G^2 L)/(rho D_v) (mu/mu_w)^0.14 (D_v/S_T)^0.4 (S_L/S_T)^0.6

    with G = rho u_max and L = N S_L the depth of the bank.  Reducing that to
    the library's currency, the mass flux and the row count both cancel:

        Eu_row = 2 (f/2) (S_L/D_v) (D_v/S_T)^0.4 (S_L/S_T)^0.6

    THREE THINGS WERE CHECKED, because the leading factor is the one place
    two sources could be read differently:

      * D_v reproduces the 0.1334 m the KAERI CHX paper states for
        D = 27.2 mm, a = 2.65, b = 1.75, to 6e-6 m.
      * the two branches meet at the stated transition of Re_v = 200 to 1.1 %,
        which is what a correctly transcribed pair does and what a misread
        factor of two would destroy.
      * the magnitude lands between the in-line and staggered branches of
        Zukauskas at the same pitch, which is what a correlation fitted across
        arrangements should do.  Read as `2f` instead of `f/2` it would come
        out four times larger than either.

    The viscosity ratio is exactly 1 here: these cases are isothermal.  Which
    way up that ratio goes was NOT settled from the two sources to hand, so it
    is not applied - and saying so is cheaper than a term that is silently
    upside down the first time somebody runs a heated case.
    """
    Re = max(st["Re_v"], 1.0)
    f_half = 90.0 / Re if Re <= 200.0 else 0.96 * Re ** -0.145
    eu = (2.0 * f_half * (st["SL"] / st["Dv"])
          * (st["Dv"] / st["ST"]) ** 0.4 * (st["SL"] / st["ST"]) ** 0.6)
    return eu, {"f_half": f_half, "Re_v": Re, "Dv": st["Dv"],
                "f_convention": "dp = (f/2) G^2 L / (rho D_v) x pitch factors",
                "isothermal": "the (mu/mu_w)^0.14 term is 1 and is not applied"}


def _shen2024(st):
    """Shen et al. (2024), eq. 12 - liquid metal across a HELICAL bundle.

        f = (209.8/Re + 0.598/Re^0.037) (a b)^-0.69 (cos eps)^-(4.2K+3)
        eps = beta (1 - beta/90),  a = S_T/d,  b = S_L/d
        K = 0 same coiling direction, 1 alternating

    Their f is defined as 2 dp / (rho u_max^2 z) with z the row count - which
    is this library's Eu_row exactly, so no conversion is needed.  That is
    worth noticing: the paper and this project already speak one currency.

    Checked against the paper's own quoted numbers: the pitch term gives
    S/d 1.4 as 10.0 % above S/d 1.5, where the paper says 9.9 %, and 20.2 %
    above S/d 1.6 where the paper says 18.6 %.  That is the exponent
    confirmed, sign included.

    Two limits worth being honest about.  At beta = 0 the helix factor is
    exactly 1 by construction and the expression becomes a straight-bundle
    correlation - but the fit never saw a straight bundle, its range starts at
    2 degrees, so a straight-rod comparison against it is an extrapolation and
    is flagged as one.  And it was fitted on lead-bismuth at 1.4 <= S/d <= 1.6;
    the friction factor of a forced flow should not care about Pr, but the
    pitch range is narrow and the study will cross both edges of it.
    """
    Re = max(st["Re"], 1.0)
    a, b = st["XT_d"], st["XL_d"]
    K = st.get("coil_K", 0)
    cos_eps = max(st.get("cos_eps", 1.0), 1e-6)
    eu = ((209.8 / Re + 0.598 / Re ** 0.037) * (a * b) ** -0.69
          * cos_eps ** -(4.2 * K + 3.0))
    return eu, {"a": a, "b": b, "K": K, "cos_eps": cos_eps,
                "helix": st.get("helix", 0.0),
                "f_convention": "f = 2 dp / (rho u_max^2 z) - this is Eu_row"}


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
                "Held against the others over a square-pitch grid in the "
                "stage-1 report - not asserted here, computed there. Against "
                "ZUKAUSKAS ALONE it sits within about 20 % inside its own Re "
                "range at X >= 1.5, and parts by up to 2.2x at X = 1.25 or "
                "outside 2e3 < Re < 4e4. Adding Gunter & Shaw WIDENS that, "
                "which is the more important result: the three applicable "
                "correlations differ from one another by 1.0x to 2.2x over the "
                "grid, typically about 1.45x, and there is no condition at "
                "which they agree closely.",
        "valid_note_ko": "Holman이 인용한 적용 범위는 2000 < Re_max < 40 "
                         "000 입니다. 피치 범위는 원 데이터가 취해진 "
                         "범위이며, 두 제한 중 더 느슨한 쪽입니다.",
        "note_ko": "닫힌 형태이고, 형상 탭이 계속 보여 온 바로 그 "
                   "상관식입니다. 다른 상관식들과의 대조는 여기서 주장하지 "
                   "않고 1단계 보고서에서 계산합니다. Zukauskas 하나만 놓고 "
                   "보면 Jakob 자신의 Re 범위 안, X ≥ 1.5 에서 ±20 % 안에 "
                   "들고, X = 1.25 이거나 Re 가 범위 밖이면 최대 2.2배까지 "
                   "갈라집니다. 그런데 Gunter & Shaw 를 더하면 그 폭이 오히려 "
                   "더 넓어집니다. 그것이 더 중요한 결과입니다 — 적용 가능한 "
                   "세 상관식은 격자 전체에서 서로 1.0–2.2배, 보통 1.45배쯤 "
                   "차이가 나며, 서로 바짝 일치하는 조건은 없습니다.",
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
        "key": "gunter-shaw",
        "name": "Gunter & Shaw (1945)",
        "status": "encoded",
        "arrangement": "both",
        "fn": _gunter_shaw,
        "source": "A. Y. Gunter and W. A. Shaw, 'A general correlation of "
                  "friction factors for various types of surfaces in cross "
                  "flow', Trans. ASME 67 (1945) 643-660. Transcribed from two "
                  "independent secondary sources that agree: H. G. Noh, "
                  "J. Eoh, D. E. Kim and M. H. Kim, Trans. KNS Spring Meeting "
                  "2018, 18S-013 (the KAERI CHX study), and E. K. Kim, "
                  "Y. S. Kim and Y. S. Sim, Trans. KNS Spring Meeting 2000 "
                  "(COMMIX-HSG), which give the same f/2 branches and the same "
                  "pressure-drop form.",
        "valid": {"Re": [1.0, 1.0e6], "XT": [1.25, 5.0], "XL": [1.25, 5.0],
                  "n_rows": [4, None]},
        "valid_note": "Gunter & Shaw fitted tubes of 0.02 to 2 inches at "
                      "transverse and longitudinal pitches of 1.25 to 5 "
                      "diameters, bare and extended surfaces together. The Re "
                      "limits are on Re_v, built on the volumetric hydraulic "
                      "diameter, which for a typical bank is an order of "
                      "magnitude larger than Re on the tube diameter.",
        "valid_note_ko": "Gunter & Shaw 는 관경 0.02–2 inch, 횡·종 피치비 "
                         "1.25–5 의 나관과 확장면 데이터를 함께 맞췄습니다. "
                         "Re 한계는 체적 수력직경 기준 Re_v 이며, 보통의 다발에서 "
                         "관경 기준 Re 보다 한 자릿수 큽니다.",
        "note": "The only correlation here that does NOT separate in-line from "
                "staggered: it collapses both onto one curve through the "
                "volumetric hydraulic diameter and two pitch-ratio factors, "
                "and lands between Zukauskas' two branches at the same pitch. "
                "That is its appeal and its limitation. It is also the "
                "correlation the KAERI CHX study found came CLOSEST to CFD on "
                "a real helical bundle - 45 % out where Zukauskas was 62 % "
                "out, which is the measurement this whole project exists to "
                "improve on.",
        "note_ko": "여기서 유일하게 정렬/엇갈림을 구분하지 않는 상관식입니다. "
                   "체적 수력직경과 두 개의 피치비 인자로 두 배열을 하나의 곡선에 "
                   "모으며, 같은 피치에서 Zukauskas 의 두 분기 사이에 놓입니다. "
                   "그것이 장점이자 한계입니다. 또한 KAERI CHX 연구에서 실제 나선 "
                   "다발 CFD 에 가장 가까웠던 상관식이기도 합니다 — Zukauskas 가 "
                   "62 % 빗나갈 때 45 % 였습니다. 바로 이 수치가 이 과제의 존재 "
                   "이유입니다.",
    },
    {
        "key": "shen-2024",
        "name": "Shen et al. (2024), helical bundle",
        "status": "encoded",
        "arrangement": "both",
        "fn": _shen2024,
        "source": "C. Shen, M. Liu, L. Liu, Z. Xu, C. Zeng, L. Liu and H. Gu, "
                  "'Development of friction factor and heat transfer "
                  "correlation of liquid metal flow in helical tube bundles', "
                  "Annals of Nuclear Energy 201 (2024) 110442, eq. (12). The "
                  "deviation angle eps = beta(1 - beta/90) is Gilli's (1965).",
        "valid": {"Re": [2500.0, 120000.0], "XT": [1.4, 1.6], "XL": [1.4, 1.6],
                  "n_rows": [4, None]},
        "valid_note": "Fitted on lead-bismuth CFD at helix angles of 2 to 15 "
                      "degrees, S/d 1.4 to 1.6, Re 2500 to 120 000, with "
                      "prediction error inside 10 %. A straight bundle is "
                      "OUTSIDE it: the helix factor is exactly 1 at beta = 0 "
                      "by construction, but the fit never saw beta = 0.",
        "valid_note_ko": "납-비스무트 CFD 로, 나선각 2–15°, S/d 1.4–1.6, "
                         "Re 2500–120 000 범위에서 맞췄고 예측 오차는 10 % 이내 "
                         "입니다. 직선 다발은 이 범위 밖입니다 — 나선 인자는 "
                         "beta = 0 에서 구조상 정확히 1이지만, 적합 과정에서 "
                         "beta = 0 을 본 적은 없습니다.",
        "note": "The only entry here written for a COIL, and the closest thing "
                "in the literature to what stage 4 is meant to produce. Its f "
                "is defined as 2 dp/(rho u_max^2 z), which is this library's "
                "Eu_row exactly - the paper and this project already speak one "
                "currency. Its shape is also the shape stage 4 should start "
                "from: a laminar plus turbulent Re function, one factor in the "
                "pitch product, one in the helix angle that is unity when the "
                "bundle is straight. Two of its own quoted numbers were "
                "reproduced from the formula as a transcription check.",
        "note_ko": "여기서 유일하게 나선 코일을 위해 쓰인 항목이고, 4단계가 "
                   "만들려는 것에 문헌상 가장 가까운 식입니다. f 의 정의가 "
                   "2 dp/(rho u_max^2 z) 로 이 라이브러리의 Eu_row 와 정확히 "
                   "같습니다 — 논문과 이 과제가 이미 같은 단위를 쓰고 있습니다. "
                   "형태 또한 4단계가 출발해야 할 형태입니다: 층류항 + 난류항의 "
                   "Re 함수, 피치 곱 인자 하나, 그리고 직선일 때 1이 되는 나선각 "
                   "인자 하나. 논문이 인용한 수치 두 개를 식에서 재현해 전사 "
                   "오류가 없음을 확인했습니다.",
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
#  shape before it needs data.  The shape is NOT invented here: Shen et al.
#  (2024) fitted a helical liquid-metal bundle with
#
#      Eu_row = (209.8/Re + 0.598 Re^-0.037) (X_T X_L)^-0.69 (cos eps)^-(4.2K+3)
#
#  and every other entry in the registry is a special case of that skeleton -
#  a laminar term plus a turbulent term in Re, one factor in the pitch, and a
#  factor in the lean that is exactly 1 when the bundle is straight.  Starting
#  from a published shape and re-fitting its coefficients on our own data is a
#  smaller and far more defensible claim than proposing a new functional form,
#  and it makes the straight-rod fit and the coil fit the SAME correlation at
#  two values of one angle rather than two correlations that have to be
#  reconciled afterwards.
#
#  ONE GENERALISATION, AND THE REASON FOR IT.  Shen's pitch factor is
#  (X_T X_L)^-0.69 - the two ratios only ever appear as a product.  That is
#  not a physical claim, it is a consequence of their dataset: every case they
#  ran had S_T = S_L, so no data of theirs could tell the two exponents apart.
#  It is the same collinearity this project's stage-3 matrix was laid out to
#  avoid.  The form here gives them separate exponents,
#
#      X_T^-p X_L^-q,      which is Shen's when p = q,
#
#  so our data can answer a question theirs could not, and `tie_pitch` refits
#  with p forced equal to q so the two can be compared on the same data.
FIT_FORM = {
    "name": "Shen-form generalised: laminar + turbulent in Re, "
            "separate pitch exponents, helix lean",
    "expr": "Eu_row = (A/Re + B Re^-m) X_T^-p X_L^-q (cos eps)^-(cK K + c0)",
    "coefficients": ["A", "B", "m", "p", "q", "c0", "cK"],
    "reference": "the shape of Shen et al. (2024) eq. 12, with the pitch "
                 "product split, refitted",
    "helical_extension": "eps = beta(1 - beta/90) is Gilli's deviation angle "
                         "and K is 0 for one coiling direction, 1 for "
                         "alternating. At beta = 0 the last factor is exactly "
                         "1, so the straight-rod fit is the coil fit's own "
                         "limit and c0/cK are not identifiable from "
                         "straight-rod data alone - which is why they are held "
                         "at Shen's values until coil cases exist.",
    "shen_values": {"A": 209.8, "B": 0.598, "m": 0.037, "p": 0.69,
                    "q": 0.69, "c0": 3.0, "cK": 4.2},
}


def fit_eu(coeff, st):
    """Evaluate the fitted form.  `st` needs Re, XT_d, XL_d and, for a coil,
    cos_eps and coil_K."""
    Re = max(st.get("Re", 1.0), 1.0)
    a = st.get("XT_d", st.get("XT", 1.0))
    b = st.get("XL_d", st.get("XL", 1.0))
    eu = ((coeff["A"] / Re + coeff["B"] * Re ** -coeff["m"])
          * a ** -coeff["p"] * b ** -coeff.get("q", coeff["p"]))
    cos_eps = max(st.get("cos_eps", 1.0), 1e-6)
    if cos_eps < 1.0:
        K = st.get("coil_K", 0)
        eu *= cos_eps ** -(coeff.get("cK", 4.2) * K + coeff.get("c0", 3.0))
    return eu


def _solve(M):
    """Gauss-Jordan with partial pivoting, on an augmented matrix."""
    n = len(M)
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[piv][i]) < 1e-300:
            raise ValueError("singular: the cases do not vary every "
                             "coefficient independently")
        M[i], M[piv] = M[piv], M[i]
        d = M[i][i]
        M[i] = [v / d for v in M[i]]
        for r in range(n):
            if r != i and M[r][i]:
                f = M[r][i]
                M[r] = [v - f * w for v, w in zip(M[r], M[i])]
    return [M[i][n] for i in range(n)]


def _nelder_mead(f, x0, step, iters=4000, tol=1e-10):
    """Nelder-Mead, written out.

    The Shen form is linear in A and B once the exponents are fixed, so the
    fit is nested: a simplex search over the three or five EXPONENTS, with an
    ordinary least-squares solve for A and B inside it.  That keeps the
    non-linear part to a handful of well-scaled parameters, which is what
    makes a plain simplex adequate and a dependency unnecessary.
    """
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += step[i]
        pts.append(p)
        vals = None
    vals = [f(p) for p in pts]
    for _ in range(iters):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) <= tol * (abs(vals[0]) + tol):
            break
        cen = [sum(p[i] for p in pts[:-1]) / n for i in range(n)]
        ref = [cen[i] + (cen[i] - pts[-1][i]) for i in range(n)]
        fr = f(ref)
        if fr < vals[0]:
            exp = [cen[i] + 2.0 * (cen[i] - pts[-1][i]) for i in range(n)]
            fe = f(exp)
            pts[-1], vals[-1] = (exp, fe) if fe < fr else (ref, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = ref, fr
        else:
            con = [cen[i] + 0.5 * (pts[-1][i] - cen[i]) for i in range(n)]
            fc = f(con)
            if fc < vals[-1]:
                pts[-1], vals[-1] = con, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = [pts[0][j] + 0.5 * (pts[i][j] - pts[0][j])
                              for j in range(n)]
                    vals[i] = f(pts[i])
    i = min(range(n + 1), key=lambda i: vals[i])
    return pts[i], vals[i]


def fit_form(rows, with_helix=False, tie_pitch=False):
    """Fit the Shen-shaped form to measured Eu_row.

    `rows` are dicts with Re, XT_d (or XT), XL_d (or XL), eu_row and, for
    coils, cos_eps and coil_K.  Returns the coefficients and how well it fits,
    in relative error, which is the currency a reader cares about.

    The residual is on ln(Eu), not on Eu: the data span a decade and a half
    and a least-squares fit on the raw value would be decided entirely by the
    tightest, slowest cases.
    """
    pts = [r for r in rows if r.get("eu_row")]
    if len(pts) < 6:
        raise ValueError("%d usable points is not enough to fit this form"
                         % len(pts))
    helical = with_helix and any(r.get("cos_eps", 1.0) < 1.0 for r in pts)

    def linear_part(exps):
        """Given the exponents, the best A and B, and the residual."""
        m, p = exps[0], exps[1]
        q = p if tie_pitch else exps[2]
        rest = exps[2 if tie_pitch else 3:]
        c0, cK = (rest[0], rest[1]) if helical else (3.0, 4.2)
        basis, y = [], []
        for r in pts:
            Re = max(r["Re"], 1.0)
            a = r.get("XT_d", r.get("XT"))
            b = r.get("XL_d", r.get("XL"))
            g = a ** -p * b ** -q
            if helical:
                g *= max(r.get("cos_eps", 1.0), 1e-6) ** -(cK * r.get("coil_K", 0) + c0)
            basis.append((g / Re, g * Re ** -m))
            y.append(r["eu_row"])
        #  least squares on ln is non-linear in A,B, so solve the linear
        #  problem on the value and then report the error on ln - the two
        #  agree to second order and this keeps the inner solve closed form
        S = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        for (u, v), yy in zip(basis, y):
            w = 1.0 / max(yy * yy, 1e-300)          # relative, not absolute
            S[0][0] += w * u * u
            S[0][1] += w * u * v
            S[1][0] += w * u * v
            S[1][1] += w * v * v
            S[0][2] += w * u * yy
            S[1][2] += w * v * yy
        #  A and B are held NON-NEGATIVE.  A is the laminar term and B the
        #  turbulent one; a negative laminar coefficient has no meaning, and
        #  left free the solve will happily use one to fake a steeper Re
        #  dependence than the form can otherwise produce - which is exactly
        #  what it did the first time, returning A = -472 against Shen's
        #  +209.8 and a 27 % fit.  With two unknowns the constrained solution
        #  is one of three candidates, so it is found exactly rather than
        #  iterated: the free solve if it is feasible, else A = 0, else B = 0.
        def resid(A, B):
            e = 0.0
            for (u, v), yy in zip(basis, y):
                pred = A * u + B * v
                if pred <= 0.0:
                    return 1e9
                e += (math.log(pred / yy)) ** 2
            return e / len(pts)

        cands = []
        try:
            A, B = _solve([row[:] for row in S])
            if A >= 0.0 and B >= 0.0:
                cands.append((A, B))
        except ValueError:
            pass
        if S[1][1] > 0:
            cands.append((0.0, max(S[1][2] / S[1][1], 0.0)))      # A pinned
        if S[0][0] > 0:
            cands.append((max(S[0][2] / S[0][0], 0.0), 0.0))      # B pinned
        best = min(cands, key=lambda ab: resid(*ab)) if cands else (0.0, 0.0)
        return resid(*best), best

    x0 = [0.037, 0.69] + ([] if tie_pitch else [0.69]) \
        + ([3.0, 4.2] if helical else [])
    step = [0.02, 0.15] + ([] if tie_pitch else [0.15]) \
        + ([0.6, 0.8] if helical else [])
    best, _ = _nelder_mead(lambda x: linear_part(x)[0], x0, step)
    val, (A, B) = linear_part(best)
    rest = best[2 if tie_pitch else 3:]
    coeff = {"A": A, "B": B, "m": best[0], "p": best[1],
             "q": best[1] if tie_pitch else best[2],
             "c0": rest[0] if helical else 3.0,
             "cK": rest[1] if helical else 4.2,
             "helical": helical, "tie_pitch": tie_pitch}
    errs = []
    for r in pts:
        pred = fit_eu(coeff, r)
        errs.append(abs(pred - r["eu_row"]) / r["eu_row"])
    coeff["n_points"] = len(errs)
    coeff["mean_abs_error"] = sum(errs) / len(errs)
    coeff["max_abs_error"] = max(errs)
    coeff["rms_log"] = math.sqrt(val)
    if coeff["A"] <= 0.0:
        coeff["laminar_dropped"] = (
            "the laminar term went to zero: over the Re range that was run, "
            "a single power of Re describes the data and the 1/Re term is not "
            "needed. Add cases below Re ~ 1000 before claiming it is absent")
        coeff["laminar_dropped_ko"] = (
            "층류항이 0으로 떨어졌습니다. 해석한 Re 범위 안에서는 Re 의 단일 "
            "거듭제곱만으로 데이터가 설명되고 1/Re 항이 필요하지 않다는 뜻입니다. "
            "이 항이 정말 없다고 말하려면 Re 1000 이하 케이스를 더해야 합니다.")
    coeff["fixed"] = (["q"] if tie_pitch else [])
    if not helical:
        coeff["fixed"] = coeff["fixed"] + ["c0", "cK"]
        coeff["fixed_note"] = ("cos eps is 1 on every case, so the helix "
                               "exponents cannot be seen by this data; they "
                               "are held at Shen et al.'s values")
        coeff["fixed_note_ko"] = ("모든 케이스에서 cos eps 가 1 이므로 이 "
                                  "데이터로는 나선각 지수를 볼 수 없습니다. "
                                  "Shen 등의 값으로 고정해 두었습니다.")
    if tie_pitch:
        coeff["fixed_note"] = ((coeff.get("fixed_note") or "") +
                               "; q is tied to p, which is Shen's own "
                               "(X_T X_L)^-p").lstrip("; ")
        coeff["fixed_note_ko"] = ((coeff.get("fixed_note_ko") or "") +
                                  " q 를 p 에 묶었습니다 — Shen 식의 "
                                  "(X_T X_L)^-p 형태입니다.").strip()
    return coeff


#  kept as the simpler, always-identifiable baseline: a single power of Re.
#  It cannot represent the low-Re end, which is why it is not the main form,
#  but it is linear in logs and so has no optimiser and no starting guess.
def fit_power_law(rows, with_helix=False):
    """Least squares for C, m, p, q on ln(Eu_row), in plain Python.

        ln Eu = ln C - m ln Re - p ln(X_T-1) + q ln X_L
    """
    A, b = [], []
    for r in rows:
        if not r.get("eu_row"):
            continue
        XT = r.get("XT_d", r.get("XT"))
        XL = r.get("XL_d", r.get("XL"))
        row = [1.0, -math.log(max(r["Re"], 1.0)),
               -math.log(max(XT - 1.0, 1e-3)), math.log(XL)]
        if with_helix:
            row.append(-math.log(max(r.get("cos_eps", 1.0), 1e-6)))
        A.append(row)
        b.append(math.log(r["eu_row"]))
    n = len(A[0]) if A else 0
    if len(A) < n:
        raise ValueError("%d usable points cannot fit %d coefficients"
                         % (len(A), n))
    M = [[sum(A[k][i] * A[k][j] for k in range(len(A))) for j in range(n)]
         + [sum(A[k][i] * b[k] for k in range(len(A)))] for i in range(n)]
    x = _solve(M)
    names = ["C", "m", "p", "q"] + (["r"] if with_helix else [])
    coeff = dict(zip(names, x))
    coeff["C"] = math.exp(coeff["C"])
    errs = []
    for r in rows:
        if not r.get("eu_row"):
            continue
        XT = r.get("XT_d", r.get("XT"))
        XL = r.get("XL_d", r.get("XL"))
        pred = (coeff["C"] * max(r["Re"], 1.0) ** -coeff["m"]
                * max(XT - 1.0, 1e-3) ** -coeff["p"] * XL ** coeff["q"])
        if with_helix and coeff.get("r"):
            pred *= max(r.get("cos_eps", 1.0), 1e-6) ** -coeff["r"]
        errs.append(abs(pred - r["eu_row"]) / r["eu_row"])
    coeff["n_points"] = len(errs)
    coeff["mean_abs_error"] = sum(errs) / len(errs) if errs else None
    coeff["max_abs_error"] = max(errs) if errs else None
    return coeff


# =============================================================================
#  THE STAGE-1 REPORT
# =============================================================================
def _agreement_grid(keys=None):
    keys = keys or ENCODED
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
                recs = {k: evaluate(k, st) for k in keys}
                out.append({"staggered": staggered, "X": X, "Re": st["Re"],
                            "eu": {k: r["eu_row"] for k, r in recs.items()},
                            "outside": {k: r["outside"]
                                        for k, r in recs.items()}})
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

    p("<h2>%s</h2>" % t("이 과제가 왜 필요한가", "Why this project exists"))
    p("<p>%s</p>" % t(
        "실제 나선 코일 열교환기 형상에서 기존 상관식이 얼마나 빗나가는지는 이미 "
        "측정되어 있습니다. KAERI 의 CHX (나선 코일 배열의 sodium-to-sodium "
        "열교환기, 관경 27.2 mm, P_T/D 2.65, P_L/D 1.75, 4행 × 11열) 에 대한 "
        "CFX 해석과 상관식 비교에서, 100 % 출력 조건의 쉘측 압력강하는 CFD "
        "288.2 Pa 였고 <b>Zukauskas 는 61.9 %, Gunter–Shaw 는 45.2 % 벗어났습니다</b>. "
        "저자들의 결론은 '기존 상관식으로는 예측할 수 없으므로 새 상관식을 "
        "실험적으로 개발할 필요가 있다' 였습니다.",
        "How far the existing correlations miss on a real helical-coil "
        "exchanger has already been measured. For KAERI's CHX - a "
        "sodium-to-sodium exchanger with a helically-coiled arrangement, "
        "27.2 mm tubes, P_T/D 2.65, P_L/D 1.75, 4 rows by 11 columns - the "
        "shell-side pressure drop at full power was 288.2 Pa by CFD, and "
        "<b>Zukauskas was 61.9 % out, Gunter-Shaw 45.2 % out</b>. The authors "
        "concluded that a new correlation had to be developed."))
    p('<p class="muted">%s</p>' % t(
        "출처: H. G. Noh, J. Eoh, D. E. Kim, M. H. Kim, "
        "Trans. Korean Nuclear Society Spring Meeting 2018, 18S-013. "
        "그 논문이 쓴 두 상관식이 모두 이 라이브러리에 들어 있습니다.",
        "Source: H. G. Noh, J. Eoh, D. E. Kim and M. H. Kim, Trans. Korean "
        "Nuclear Society Spring Meeting 2018, 18S-013. Both correlations that "
        "paper used are in this library."))
    p("<p>%s</p>" % t(
        "그래서 순서는 이렇습니다: 먼저 <b>직선 rod</b> 에서 우리 CFD 절차가 "
        "기존 상관식을 재현하는지 확인하고 (그것이 절차의 검증입니다), 그 다음 "
        "나선 코일로 옮겨 가서 상관식이 벗어나는 지점을 측정하고, 마지막으로 "
        "형상 정보를 반영한 상관식을 만듭니다.",
        "Hence the order: establish on <b>straight rods</b> that our CFD "
        "procedure reproduces the published correlations - that is what "
        "validates the procedure - then move to the coil and measure where the "
        "correlations leave off, then fit one that carries the geometry."))

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
        #  a sentence, not a number: `n` right-aligns it in the page's sheet
        p('<tr><th>%s</th><td>%s</td></tr>' % (t("적용 범위", "range"), rng))
        p('<tr><th>%s</th><td>%s</td></tr>' % (
            t("범위 주석", "on that range"), local(c, "valid_note")))
        p('<tr><th>%s</th><td>%s</td></tr>' % (t("비고", "note"),
                                               local(c, "note")))
        p("</table>")

    p("<h2>%s</h2>" % t("계산 가능한 상관식들의 상호 대조",
                        "The evaluable correlations, against each other"))
    p("<p>%s</p>" % t(
        "아래 표는 이 파일의 코드로 <i>지금</i> 계산한 값입니다. 정사각 피치 "
        "(X_T = X_L = X), 10행, 물 기준, 직선 rod. 서로 독립적인 출처이므로 "
        "이들 사이의 폭이 '기존 상관식'이라는 기준선 자체의 폭입니다.",
        "Computed by this file's own code, now. Square pitch (X_T = X_L = X), "
        "10 rows, water, straight rods. They come from independent sources, so "
        "the spread between them is the width of the baseline itself."))
    grid = _agreement_grid()
    keys = ENCODED
    names = {k: BY_KEY[k]["name"] for k in keys}
    for staggered in (False, True):
        p("<h3>%s</h3>" % (t("엇갈림 (staggered)", "staggered") if staggered
                           else t("정렬 (in-line)", "in-line")))
        p('<table class="grid"><tr><th>X</th><th>Re<sub>max</sub></th>'
          + "".join("<th>%s</th>" % esc_name for esc_name in
                    (names[k] for k in keys)) + "</tr>")
        for g in grid:
            if g["staggered"] != staggered:
                continue
            cells = ""
            for k in keys:
                v = g["eu"].get(k)
                out = g["outside"].get(k)
                cells += ('<td class="n">%s%s</td>'
                          % (_fmt(v, 4) if v else "—",
                             ' <i title="%s">·</i>' % _fmt(len(out or []))
                             if out else ""))
            p('<tr><td class="n">%.2f</td><td class="n">%.3g</td>%s</tr>'
              % (g["X"], g["Re"], cells))
        p("</table>")
    p('<p class="muted">%s</p>' % t(
        "· 표시는 그 점이 해당 상관식의 인용 범위 밖이라는 뜻입니다. "
        "Shen 식은 직선 다발(나선각 0)에 대해서는 원래 범위 밖이며, "
        "여기서는 그 외삽이 얼마나 벌어지는지 보려고 함께 계산했습니다.",
        "A dot means that point is outside that correlation's quoted range. "
        "Shen's is outside it for a STRAIGHT bundle by definition; it is shown "
        "so the size of that extrapolation is visible rather than assumed."))

    #  The summary, restricted as well as unrestricted.  The unrestricted
    #  spread is dominated by the corners of the grid that Jakob was never
    #  quoted for, and reporting only that would make the baseline look
    #  useless; reporting only the restricted one would make it look better
    #  than it is.  Both, and the restriction named.
    p("<h3>%s</h3>" % t("요약 — 가장 벌어진 두 상관식의 비",
                        "In summary - the widest pair, as a ratio"))
    p('<table class="grid"><tr><th>%s</th><th>%s</th><th>%s</th></tr>'
      % (t("범위", "over"), t("정렬", "in-line"), t("엇갈림", "staggered")))

    def band(pred, in_range_only=True):
        """The highest correlation over the lowest, at each point.

        Counting only the ones that are IN RANGE there, because a correlation
        being asked outside its own stated limits is not evidence about the
        baseline - it is evidence about the extrapolation.  Shen's is outside
        its range for a straight bundle at every single point, so including it
        would widen every cell for a reason that has nothing to do with how
        well the literature agrees.
        """
        cells = []
        for staggered in (False, True):
            worst = []
            for g in grid:
                if g["staggered"] != staggered or not pred(g):
                    continue
                vals = [v for k, v in g["eu"].items()
                        if v and (not in_range_only or not g["outside"].get(k))]
                if len(vals) >= 2:
                    worst.append(max(vals) / min(vals))
            cells.append('<td class="n">%s</td>' % (
                "—" if not worst else "%.2f – %.2f (%s %.2f)"
                % (min(worst), max(worst), t("평균", "mean"),
                   sum(worst) / len(worst))))
        return "".join(cells)

    p("<tr><th>%s</th>%s</tr>" % (
        t("격자 전체", "the whole grid"), band(lambda g: True)))
    p("<tr><th>%s</th>%s</tr>" % (
        t("Jakob 적용 Re (2e3–4e4)", "Jakob's quoted Re, 2e3–4e4"),
        band(lambda g: 2000.0 <= g["Re"] <= 40000.0)))
    p("<tr><th>%s</th>%s</tr>" % (
        t("그리고 X ≥ 1.5", "... and X >= 1.5"),
        band(lambda g: 2000.0 <= g["Re"] <= 40000.0 and g["X"] >= 1.5)))
    p("<tr><th>%s</th>%s</tr>" % (
        t("(범위 밖 포함)", "(including out-of-range)"),
        band(lambda g: True, in_range_only=False)))
    p("</table>")
    #  and where the coil correlation lands when it is dragged to zero lean,
    #  which is the number stage 4 has to beat
    sh = [g["eu"]["shen-2024"] / (sum(v for k, v in g["eu"].items()
                                      if k != "shen-2024" and v)
                                  / max(1, len([1 for k, v in g["eu"].items()
                                                if k != "shen-2024" and v])))
          for g in grid if g["eu"].get("shen-2024")]
    p("<p>%s</p>" % t(
        "나선 코일용인 Shen 식을 나선각 0으로 외삽하면, 나머지 세 상관식의 "
        "평균 대비 %.2f–%.2f 배 (평균 %.2f) 로 <b>일관되게 낮게</b> 나옵니다. "
        "직선 다발은 그 식의 적합 범위 밖이므로 이것은 결함이 아니라 "
        "외삽의 크기이고, 4단계가 메워야 할 간격입니다."
        % (min(sh), max(sh), sum(sh) / len(sh)),
        "Dragged to zero lean, Shen's coil correlation comes out "
        "<b>consistently low</b> against the mean of the other three - "
        "%.2f to %.2f times it, mean %.2f. A straight bundle is outside that "
        "fit's range, so this is not a fault: it is the size of the "
        "extrapolation, and the gap stage 4 has to close."
        % (min(sh), max(sh), sum(sh) / len(sh))))
    p("<p>%s</p>" % t(
        "읽는 법: 각 칸은 그 조건에서 <b>가장 높은 상관식 ÷ 가장 낮은 상관식</b>"
        "입니다. 1.0이면 전부 일치, 2.0이면 두 배 차이입니다.",
        "How to read it: each cell is the <b>highest correlation divided by "
        "the lowest</b> at that condition. 1.0 is unanimity, 2.0 is a factor "
        "of two."))
    p("<p>%s</p>" % t(
        "그리고 여기서 가장 중요한 결과가 나옵니다: <b>기존 상관식들 사이에 "
        "'정답'이라 부를 만큼 좁은 영역이 없습니다.</b> 적용 범위 안의 세 상관식 "
        "(Jakob, Zukauskas, Gunter–Shaw) 조차 격자 전체에서 서로 평균 1.4–1.5배, "
        "최대 2.2배 차이가 납니다. 두 개만 놓고 보면 좁아 보이지만 세 번째를 "
        "더하면 넓어집니다 — 앞서의 '±20 %' 는 Jakob–Zukauskas 쌍의 성질이지 "
        "문헌의 합의가 아니었습니다.",
        "And here is the result that matters most: <b>there is no region where "
        "the existing correlations are tight enough to call one of them the "
        "answer.</b> Even the three that are in range - Jakob, Zukauskas, "
        "Gunter & Shaw - differ from one another by about 1.4 to 1.5 times on "
        "average and up to 2.2 times over the grid. Two of them look close; "
        "adding a third widens it. The earlier '20 %' was a property of the "
        "Jakob-Zukauskas pair, not a consensus in the literature."))
    p("<p>%s</p>" % t(
        "3단계에 대해 이것이 말해 주는 것은 두 가지입니다. 첫째, CFD를 상관식으로 "
        "'검증'할 수 있는 정밀도의 상한이 이 폭입니다 — 어느 한 상관식의 20 % "
        "안에 든다고 맞는 것도 아니고, 40 % 벗어난다고 틀린 것도 아닙니다. "
        "판정은 <b>세 상관식이 만드는 띠 안에 드는가</b>와 <b>기울기(Re 의존성, "
        "피치 의존성)가 같은가</b>로 해야 합니다. 둘째, 그만큼 새 데이터의 "
        "가치가 큽니다.",
        "Two things follow for stage 3. First, this width is the ceiling on "
        "how precisely CFD can be 'validated' against a correlation at all: "
        "landing inside 20 % of one of them does not make the CFD right, and "
        "missing one by 40 % does not make it wrong. The test has to be "
        "<b>does it fall inside the band the three of them span</b> and "
        "<b>does it have the same slopes</b>, in Re and in pitch. Second, that "
        "is exactly how much a new measurement is worth."))

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
def _extract_js(name, text):
    """One top-level `function name(...)` out of the page, by brace matching."""
    key = "function %s(" % name
    i = text.index(key)
    j = text.index("{", i)
    depth, k = 0, j
    while k < len(text):
        c = text[k]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return text[i:k + 1]
        k += 1
    raise ValueError("unterminated %s" % name)


def _check_against_js():
    """Jakob here against Jakob in the browser - the browser's ACTUAL code.

    The Geometry tab computes the same correlation in JavaScript and must go
    on doing so: it works with no server at all, which is half the point of
    the page.  Two implementations of one formula is exactly the arrangement
    that drifts, so they are held against each other the way the mesh and its
    Python twin are.

    The two functions are LIFTED OUT OF mesh_explorer.html and run, rather
    than transcribed into this file.  A transcription would be a third copy
    and would pass happily while the page said something else - which is the
    failure this check exists to catch.
    """
    import subprocess
    import tempfile
    page = os.path.join(HERE, "mesh_explorer.html")
    with open(page, encoding="utf-8") as fh:
        text = fh.read()
    lifted = "\n".join(_extract_js(n, text) for n in ("gapVelocity", "lossModel"))

    cases = []
    for staggered in (False, True):
        for XT in (1.25, 1.5, 2.0, 2.5):
            for XL in (1.25, 1.5, 2.0, 2.5):
                for u in (0.1, 0.5, 2.0):
                    cases.append({"staggered": staggered, "XT": XT, "XL": XL,
                                  "u": u})
    harness = """
'use strict';
const cases = %s;
//  the page's own globals, reduced to what these two functions read
let P = {}, G = {};
const toM = () => 1e-3;                       // the study works in mm
const isStg = () => !!P._stg;
const smooth01 = t => t;                      // lossAt only; unused here
%s
const out = cases.map(c => {
  P = {D:10, ST:10*c.XT, SL:10*c.XL, nRows:10, vel:c.u, rho:998.2,
       nu:1.004, lIn:60, lOut:120, H:2, nCols:4, _stg:c.staggered};
  G = {bE:P.D/2, aE:P.D/2, W:P.nCols*P.ST, X0:P.lIn, X1:P.lIn+P.nRows*P.SL};
  const L = lossModel();
  return {umax:L.umax, ReMax:L.ReMax, eu:L.Eu, dp:L.dpBundle,
          diagonal:!!L.diagonal};
});
console.log(JSON.stringify(out));
""" % (json.dumps(cases), lifted)
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False,
                                     encoding="utf-8") as fh:
        fh.write(harness)
        path = fh.name
    try:
        raw = subprocess.check_output(["node", path]).decode("utf-8")
    finally:
        os.unlink(path)
    theirs = json.loads(raw)
    worst, worst_at = 0.0, None
    for c, jsv in zip(cases, theirs):
        st = flow_state(D=0.010, ST=0.010 * c["XT"], SL=0.010 * c["XL"],
                        n_rows=10, u_in=c["u"], rho=998.2, nu=1.004e-6,
                        staggered=c["staggered"])
        rec = evaluate("jakob", st)
        if bool(st["diagonal"]) != bool(jsv["diagonal"]):
            print("  the two disagree about which gap governs at %r" % (c,))
            return 1
        for what, a, b in (("u_max", st["umax"], jsv["umax"]),
                           ("Re", st["Re"], jsv["ReMax"]),
                           ("Eu_row", rec["eu_row"], jsv["eu"]),
                           ("dp", rec["dp"], jsv["dp"])):
            e = abs(a - b) / max(abs(b), 1e-30)
            if e > worst:
                worst, worst_at = e, (what, c, a, b)
    print("Jakob, python vs the page's own gapVelocity + lossModel:")
    print("  %d cases x 4 quantities, and which gap governs on every one"
          % len(cases))
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
            fh.write(standalone(body, "Stage 1 - correlations", lang))
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


def standalone(body, title="report", lang="ko"):
    """A report body as a file that opens anywhere.

    The page's own stylesheet is not reachable from the command line, so this
    carries a small one of its own.  When the tab saves the same report it uses
    the page's, exactly as the Report tab does - which is why every class used
    here also exists in mesh_explorer.html's `.paper` rules.
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
.paper table.grid th{width:auto;white-space:nowrap}
.paper table.grid td{white-space:nowrap}
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
.paper .tile .u{font-size:12px;color:#5a6678;margin-left:3px}
.paper .warn{background:#fff8e6;border-left:3px solid #d97706;padding:7px 11px;
  margin:9px 0;font-size:12.5px;color:#7a4a06}
.paper .bad{background:#fdeaea;border-left:3px solid #9b1c1c;padding:7px 11px;
  margin:9px 0;font-size:12.5px;color:#9b1c1c;font-weight:600}
.paper ul{margin:8px 0 8px 18px;padding:0;font-size:12.5px;color:#33415a}
.paper li{margin:3px 0}
.paper svg.plot{display:block;max-width:660px;margin:14px 0 6px;
  background:#fff;border:1px solid #eef1f6;border-radius:6px}
.paper .ovgrid{display:grid;
  grid-template-columns:repeat(auto-fill,minmax(236px,1fr));gap:11px;
  margin:14px 0 6px}
.paper .ovcard{border:1px solid #e4e9f0;border-left:3px solid #94a3b8;
  border-radius:8px;padding:9px 11px 7px;background:#fbfcfe}
.paper .ovhead{display:flex;justify-content:space-between;align-items:baseline;
  gap:8px;font-size:12px}
.paper .ovhead b{font-weight:660;color:#1b2431}
.paper .ovst{font-size:10.5px;font-weight:650;white-space:nowrap}
.paper .ovnums{font-size:10.5px;color:#5a6678;margin:3px 0 5px;
  font-variant-numeric:tabular-nums}
.paper .ovplot{margin-bottom:4px}
.paper .ovplot>span{display:block;font-size:9.5px;color:#7a8595;
  margin-bottom:1px}
.paper svg.spark{display:block;width:100%;height:auto;background:#fff;
  border-radius:4px}
.paper .ovwide svg.plot{max-width:none;width:100%}
@media print{body{background:#fff;padding:0}
  .paper{box-shadow:none;max-width:none;padding:0}
  .paper svg.plot{break-inside:avoid}}
"""
    return ('<!doctype html>\n<html lang="%s"><head><meta charset="utf-8">\n'
            '<title>%s</title>\n<style>%s</style>\n'
            '</head>\n<body>\n<div class="paper">\n%s\n</div>\n</body></html>\n'
            % (lang, title, css, body))


if __name__ == "__main__":
    raise SystemExit(main())
