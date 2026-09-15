#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-phase Fluent case : settings contract, PyFluent driver, offline mock
===========================================================================
This is the layer between the browser front end and ANSYS Fluent.  It holds
three things.

1. ``SETTINGS`` - the schema.  One declaration of every knob the Settings tab
   offers: its type, default, choices, and the settings-API path it drives.
   The browser BUILDS ITS PANEL FROM THIS, served over the API, so the control
   the user sees and the path the driver writes cannot drift apart.  Adding a
   setting means adding one entry here, nothing else.

2. ``FluentDriver`` - applies that schema to a real solver session through
   ``ansys-fluent-core`` and reads results back.

3. ``MockDriver`` - the same interface with no Fluent at all.  It builds the
   real mesh and returns a SYNTHETIC field, so the server, the wire format,
   the contour rendering and the report panel can all be exercised end to end
   on a machine with no licence.  Everything it returns is tagged
   ``mock: true`` and the UI labels it loudly.  It is a plumbing test, not a
   flow solution, and must never be read as one.

On paths: PyFluent generates its settings API from Fluent itself, but it also
ships a static copy per release (v242 ... v271).  ``audit_paths()`` walks every
path named below against those static trees, so a path that a Fluent upgrade
moves or renames is caught here rather than half way through a solver run.
Run it with ``python3 fluent_case.py --audit``.

What CANNOT be checked without Fluent: the string values of the enumerated
settings (``viscous.model = "k-omega"`` and friends).  Fluent supplies those at
runtime.  They are written once in ``_ENUM`` below, and the live driver calls
``allowed_values()`` on connect and reports any that the running Fluent does
not recognise, rather than failing deep inside the setup.
"""

from __future__ import annotations

import importlib
import json
import math
import os
import sys
import time

# the mesh generator is the twin of the browser's; the server builds the mesh
# with it so the browser never has to upload a hundred-megabyte file
import mesh_explorer as ME

PATCHES = ME.PATCH_ORDER          # inlet, outlet, wall_rods, wall_side_*, wall_bottom/top
WALL_PATCHES = [p for p in PATCHES if p.startswith("wall")]


# =============================================================================
#  ENUMERATED VALUES
# =============================================================================
#  These are the strings Fluent expects.  They are NOT in the static settings
#  tree - Fluent serves them at runtime - so they cannot be audited offline.
#  FluentDriver.check_enums() compares them against the live allowed_values()
#  and reports the difference instead of guessing.
_ENUM = {
    "viscous": ["laminar", "k-epsilon", "k-omega", "spalart-allmaras",
                "reynolds-stress", "les", "inviscid"],
    "k_omega_variant": ["standard", "sst", "bsl", "geko"],
    "k_epsilon_variant": ["standard", "realizable", "rng"],
    "wall_treatment": ["standard-wall-fn", "scalable-wall-fn",
                       "non-equilibrium-wall-fn", "enhanced-wall-treatment",
                       "menter-lechner", "user-defined-wall-functions"],
    "flow_scheme": ["SIMPLE", "SIMPLEC", "PISO", "Coupled"],
    "gradient": ["green-gauss-cell-based", "green-gauss-node-based",
                 "least-squares-cell-based"],
    "pressure_disc": ["standard", "second-order", "presto!", "linear",
                      "body-force-weighted"],
    "upwind": ["first-order-upwind", "second-order-upwind", "quick",
               "third-order-muscl", "power-law"],
    "zone_type": ["wall", "symmetry", "periodic"],
    "turb_spec": ["Intensity and Viscosity Ratio", "Intensity and Hydraulic Diameter",
                  "K and Omega", "K and Epsilon"],
}


# =============================================================================
#  SETTINGS PATHS, WITH PER-RELEASE ALTERNATES
# =============================================================================
#  The settings API moves between Fluent releases.  Measured against the trees
#  PyFluent ships: 2024 R2 keeps the discretisation schemes and the surface
#  integrals at a shallower depth, and 2027 R1 renames setup.models.viscous to
#  setup.models.turbulence.  Each entry below is therefore a "|"-separated list
#  of candidates, most recent spelling first; the driver writes the first one
#  the live session actually has, and audit_paths() requires that at least one
#  resolves in every release we claim to support.  One table, used by the
#  schema, the driver, the audit and the generated journal alike.
PATHS = {
    "solver_time":     "setup.general.solver.time",
    "op_pressure":     "setup.general.operating_conditions.operating_pressure",
    "energy":          "setup.models.energy.enabled",

    "viscous_model":   "setup.models.viscous.model|setup.models.turbulence.model",
    "k_omega_model":   "setup.models.viscous.k_omega_model|"
                       "setup.models.turbulence.k_omega_model",
    "k_epsilon_model": "setup.models.viscous.k_epsilon_model|"
                       "setup.models.turbulence.k_epsilon_model",
    "wall_treatment":  "setup.models.viscous.near_wall_treatment.wall_treatment|"
                       "setup.models.turbulence.near_wall_treatment.wall_treatment",

    "materials":       "setup.materials.fluid",
    "density_option":  "setup.materials.fluid.density.option",
    "density_value":   "setup.materials.fluid.density.value",
    "visc_option":     "setup.materials.fluid.viscosity.option",
    "visc_value":      "setup.materials.fluid.viscosity.value",

    "inlet_spec":      "setup.boundary_conditions.velocity_inlet.momentum."
                       "velocity_specification_method",
    "inlet_magnitude": "setup.boundary_conditions.velocity_inlet.momentum."
                       "velocity_magnitude|"
                       "setup.boundary_conditions.velocity_inlet.momentum.velocity",
    "inlet_components": "setup.boundary_conditions.velocity_inlet.momentum."
                        "velocity_components",
    "inlet_turb_spec": "setup.boundary_conditions.velocity_inlet.turbulence."
                       "turbulence_specification",
    "inlet_intensity": "setup.boundary_conditions.velocity_inlet.turbulence."
                       "turbulent_intensity",
    "inlet_visc_ratio": "setup.boundary_conditions.velocity_inlet.turbulence."
                        "turbulent_viscosity_ratio",
    "inlet_hyd_diam":  "setup.boundary_conditions.velocity_inlet.turbulence."
                       "hydraulic_diameter",

    "outlet_pressure": "setup.boundary_conditions.pressure_outlet.momentum."
                       "gauge_pressure",
    "outlet_no_reverse": "setup.boundary_conditions.pressure_outlet.momentum."
                         "prevent_reverse_flow",
    "outlet_bf_intensity": "setup.boundary_conditions.pressure_outlet.turbulence."
                           "backflow_turbulent_intensity",
    "outlet_bf_visc_ratio": "setup.boundary_conditions.pressure_outlet.turbulence."
                            "backflow_turbulent_viscosity_ratio",

    "set_zone_type":   "setup.boundary_conditions.set_zone_type",

    "flow_scheme":     "solution.methods.p_v_coupling.flow_scheme",
    "gradient":        "solution.methods.spatial_discretization.gradient_scheme|"
                       "solution.methods.gradient_scheme",
    "discretization":  "solution.methods.spatial_discretization.discretization_scheme|"
                       "solution.methods.discretization_scheme",
    "pseudo_time":     "solution.methods.pseudo_time_method.formulation.coupled_solver",

    "under_relaxation": "solution.controls.under_relaxation",
    "residual_eqs":    "solution.monitor.residual.equations",

    "init_type":       "solution.initialization.initialization_type",
    "hybrid_init":     "solution.initialization.hybrid_initialize",
    "standard_init":   "solution.initialization.standard_initialize",
    "iterate":         "solution.run_calculation.iterate",
    "interrupt":       "solution.run_calculation.interrupt",

    "read_mesh":       "file.read_mesh",
    "write_case_data": "file.write_case_data",

    "si_area_avg":     "results.report.surface_integrals.get_area_weighted_avg|"
                       "results.report.surface_integrals.area_weighted_avg",
    "si_mass_flow":    "results.report.surface_integrals.get_mass_flow_rate|"
                       "results.report.surface_integrals.mass_flow_rate",
    "si_facet_min":    "results.report.surface_integrals.get_facet_min|"
                       "results.report.surface_integrals.facet_min",
    "si_facet_max":    "results.report.surface_integrals.get_facet_max|"
                       "results.report.surface_integrals.facet_max",
    "si_area":         "results.report.surface_integrals.get_area|"
                       "results.report.surface_integrals.area",
    "plane_surface":   "results.surfaces.plane_surface",
}


def resolve_obj(root, spec):
    """First alternate of `spec` that this session actually has -> the object."""
    for path in spec.split("|"):
        obj, ok = root, True
        for part in path.split("."):
            if not hasattr(obj, part):
                ok = False
                break
            obj = getattr(obj, part)
        if ok:
            return obj
    raise DriverError("none of these settings paths exist here: " + spec)


def set_path(root, spec, value):
    """Write `value` at the first alternate of `spec` that exists."""
    for path in spec.split("|"):
        parts = path.split(".")
        obj, ok = root, True
        for part in parts[:-1]:
            if not hasattr(obj, part):
                ok = False
                break
            obj = getattr(obj, part)
        if ok and hasattr(obj, parts[-1]):
            setattr(obj, parts[-1], value)
            return path
    raise DriverError("none of these settings paths exist here: " + spec)


def _c(*pairs):
    """choices: ('value', 'korean', 'english') triples -> schema form"""
    return [{"v": v, "ko": ko, "en": en} for (v, ko, en) in pairs]


# =============================================================================
#  THE SETTINGS SCHEMA
# =============================================================================
#  field:  id, kind (number|choice|bool|text), label, default, and `paths`:
#          the settings-API paths that field drives, for audit_paths().
#          `show_if` hides a field when another field makes it meaningless.
SETTINGS = [
    {
        "id": "launch", "ko": "실행", "en": "Launch",
        "fields": [
            {"id": "processors", "kind": "number", "default": 4, "min": 1, "max": 256,
             "step": 1, "int": True, "ko": "코어 수", "en": "Processor count"},
            {"id": "precision", "kind": "choice", "default": "double",
             "ko": "정밀도", "en": "Precision",
             "choices": _c(("double", "배정밀도", "Double"),
                           ("single", "단정밀도", "Single"))},
            {"id": "version", "kind": "text", "default": "",
             "ko": "Fluent 버전 (빈칸 = 자동)", "en": "Fluent version (blank = auto)",
             "hint_ko": "예: 26.1.0", "hint_en": "e.g. 26.1.0"},
            {"id": "ui_mode", "kind": "choice", "default": "no_gui",
             "ko": "UI 모드", "en": "UI mode",
             "choices": _c(("no_gui", "GUI 없음", "No GUI"),
                           ("gui", "GUI 표시", "Show GUI"))},
        ],
    },
    {
        "id": "general", "ko": "일반", "en": "General",
        "fields": [
            {"id": "steady", "kind": "bool", "default": True,
             "ko": "정상 상태 (steady)", "en": "Steady state",
             "paths": ["solver_time"]},
            {"id": "operating_pressure", "kind": "number", "default": 101325.0,
             "min": 0, "max": 5e7, "step": 100, "unit": "Pa",
             "ko": "작동 압력", "en": "Operating pressure",
             "paths": ["op_pressure"]},
            {"id": "energy", "kind": "bool", "default": False,
             "ko": "에너지 방정식", "en": "Energy equation",
             "paths": ["energy"]},
        ],
    },
    {
        "id": "turbulence", "ko": "난류 모델", "en": "Turbulence",
        "fields": [
            {"id": "viscous", "kind": "choice", "default": "k-omega",
             "ko": "점성 모델", "en": "Viscous model",
             "choices": _c(("laminar", "층류", "Laminar"),
                           ("k-omega", "k-ω", "k-omega"),
                           ("k-epsilon", "k-ε", "k-epsilon"),
                           ("spalart-allmaras", "Spalart-Allmaras", "Spalart-Allmaras")),
             "paths": ["viscous_model"]},
            {"id": "k_omega_variant", "kind": "choice", "default": "sst",
             "ko": "k-ω 변형", "en": "k-omega variant",
             "show_if": {"viscous": ["k-omega"]},
             "choices": _c(("sst", "SST", "SST"), ("standard", "표준", "Standard"),
                           ("bsl", "BSL", "BSL"), ("geko", "GEKO", "GEKO")),
             "paths": ["k_omega_model"]},
            {"id": "k_epsilon_variant", "kind": "choice", "default": "realizable",
             "ko": "k-ε 변형", "en": "k-epsilon variant",
             "show_if": {"viscous": ["k-epsilon"]},
             "choices": _c(("realizable", "Realizable", "Realizable"),
                           ("standard", "표준", "Standard"), ("rng", "RNG", "RNG")),
             "paths": ["k_epsilon_model"]},
            {"id": "wall_treatment", "kind": "choice", "default": "",
             "ko": "벽면 처리 (빈칸 = Fluent 기본)", "en": "Near-wall treatment (blank = Fluent default)",
             "show_if": {"viscous": ["k-epsilon"]},
             "choices": _c(("", "기본값", "Default"),
                           ("standard-wall-fn", "표준 벽함수", "Standard wall functions"),
                           ("scalable-wall-fn", "Scalable 벽함수", "Scalable wall functions"),
                           ("enhanced-wall-treatment", "Enhanced wall treatment",
                            "Enhanced wall treatment")),
             "paths": ["wall_treatment"]},
        ],
    },
    {
        "id": "material", "ko": "물성치", "en": "Material",
        "fields": [
            {"id": "name", "kind": "text", "default": "water-liquid",
             "ko": "유체 이름", "en": "Fluid name",
             "hint_ko": "Fluent 데이터베이스 이름", "hint_en": "Fluent database name",
             "paths": ["materials"]},
            {"id": "density", "kind": "number", "default": 998.2, "min": 1e-3,
             "max": 2e4, "step": 0.1, "unit": "kg/m³",
             "ko": "밀도 ρ", "en": "Density",
             "paths": ["density_option", "density_value"]},
            {"id": "viscosity", "kind": "number", "default": 1.003e-3, "min": 1e-8,
             "max": 1e3, "step": 1e-6, "unit": "kg/m·s",
             "ko": "점성계수 μ", "en": "Dynamic viscosity",
             "paths": ["visc_option", "visc_value"]},
        ],
    },
    {
        "id": "inlet", "ko": "입구 조건", "en": "Inlet",
        "fields": [
            {"id": "spec", "kind": "choice", "default": "components",
             "ko": "속도 지정 방식", "en": "Velocity specification",
             "choices": _c(("components", "성분 (U,0,0)", "Components (U,0,0)"),
                           ("magnitude-normal", "크기 · 면 법선", "Magnitude, normal to boundary")),
             "paths": ["inlet_spec"],
             "note_ko": "이 메시는 입구면이 유동에 정확히 수직이므로 두 방식이 같은 결과를 줍니다. "
                        "성분 지정이 기울기에 영향을 받지 않아 더 안전합니다.",
             "note_en": "This mesh keeps the inlet exactly perpendicular to the flow, so both "
                        "give the same thing. Components is the safer of the two."},
            {"id": "velocity", "kind": "number", "default": 2.0, "min": 1e-4, "max": 500,
             "step": 0.05, "unit": "m/s", "ko": "유입 속도 U", "en": "Inlet velocity",
             "paths": ["inlet_magnitude", "inlet_components"]},
            {"id": "turb_spec", "kind": "choice", "default": "Intensity and Viscosity Ratio",
             "ko": "난류 지정", "en": "Turbulence specification",
             "show_if": {"~viscous": ["laminar"]},
             "choices": _c(("Intensity and Viscosity Ratio", "강도 · 점성비",
                            "Intensity and viscosity ratio"),
                           ("Intensity and Hydraulic Diameter", "강도 · 수력직경",
                            "Intensity and hydraulic diameter")),
             "paths": ["inlet_turb_spec"]},
            {"id": "intensity", "kind": "number", "default": 5.0, "min": 0.01, "max": 100,
             "step": 0.1, "unit": "%", "ko": "난류 강도", "en": "Turbulent intensity",
             "show_if": {"~viscous": ["laminar"]},
             "paths": ["inlet_intensity"]},
            {"id": "visc_ratio", "kind": "number", "default": 10.0, "min": 1e-3, "max": 1e5,
             "step": 1, "ko": "난류 점성비", "en": "Turbulent viscosity ratio",
             "show_if": {"~viscous": ["laminar"], "turb_spec": ["Intensity and Viscosity Ratio"]},
             "paths": ["inlet_visc_ratio"]},
            {"id": "hydraulic_diameter", "kind": "number", "default": 0.0, "min": 0, "max": 100,
             "step": 1e-4, "unit": "m", "ko": "수력직경 (0 = 자동 계산)",
             "en": "Hydraulic diameter (0 = computed)",
             "show_if": {"~viscous": ["laminar"],
                         "turb_spec": ["Intensity and Hydraulic Diameter"]},
             "paths": ["inlet_hyd_diam"]},
        ],
    },
    {
        "id": "outlet", "ko": "출구 조건", "en": "Outlet",
        "fields": [
            {"id": "gauge_pressure", "kind": "number", "default": 0.0, "min": -1e7,
             "max": 1e7, "step": 10, "unit": "Pa", "ko": "게이지 압력", "en": "Gauge pressure",
             "paths": ["outlet_pressure"]},
            {"id": "prevent_reverse_flow", "kind": "bool", "default": True,
             "ko": "역류 방지", "en": "Prevent reverse flow",
             "paths": ["outlet_no_reverse"]},
            {"id": "backflow_intensity", "kind": "number", "default": 5.0, "min": 0.01,
             "max": 100, "step": 0.1, "unit": "%",
             "ko": "역류 난류 강도", "en": "Backflow turbulent intensity",
             "show_if": {"~viscous": ["laminar"]},
             "paths": ["outlet_bf_intensity"]},
            {"id": "backflow_visc_ratio", "kind": "number", "default": 10.0, "min": 1e-3,
             "max": 1e5, "step": 1, "ko": "역류 난류 점성비",
             "en": "Backflow turbulent viscosity ratio",
             "show_if": {"~viscous": ["laminar"]},
             "paths": ["outlet_bf_visc_ratio"]},
        ],
    },
    {
        "id": "zones", "ko": "벽면 zone 유형", "en": "Wall zone types",
        "note_ko": "측벽과 상·하면을 symmetry로 두면 무한 배열의 한 구획을 푸는 것이 됩니다. "
                   "반쪽 rod 옵션을 켠 메시에서 특히 의미가 있습니다. periodic은 Fluent에서 "
                   "짝을 직접 지정해야 하므로 여기서는 제공하지 않습니다.",
        "note_en": "Setting the sides and the caps to symmetry solves one bay of an infinite "
                   "array, which is what the half-rod option is for. Periodic is left out: it "
                   "needs its pairing set up inside Fluent.",
        "fields": [
            {"id": "wall_rods", "kind": "choice", "default": "wall",
             "ko": "rod 벽면", "en": "Rod walls",
             "choices": _c(("wall", "wall (점착)", "wall (no slip)")),
             "paths": ["set_zone_type"]},
            {"id": "wall_side_y0", "kind": "choice", "default": "wall",
             "ko": "측벽 y = 0", "en": "Side wall y = 0",
             "choices": _c(("wall", "wall", "wall"), ("symmetry", "symmetry", "symmetry")),
             "paths": ["set_zone_type"]},
            {"id": "wall_side_ymax", "kind": "choice", "default": "wall",
             "ko": "측벽 y = W", "en": "Side wall y = W",
             "choices": _c(("wall", "wall", "wall"), ("symmetry", "symmetry", "symmetry")),
             "paths": ["set_zone_type"]},
            {"id": "wall_bottom", "kind": "choice", "default": "wall",
             "ko": "하면 z = 0", "en": "Bottom z = 0",
             "choices": _c(("wall", "wall", "wall"), ("symmetry", "symmetry", "symmetry")),
             "paths": ["set_zone_type"]},
            {"id": "wall_top", "kind": "choice", "default": "wall",
             "ko": "상면 z = H", "en": "Top z = H",
             "choices": _c(("wall", "wall", "wall"), ("symmetry", "symmetry", "symmetry")),
             "paths": ["set_zone_type"]},
        ],
    },
    {
        "id": "methods", "ko": "해석 기법", "en": "Solution methods",
        "fields": [
            {"id": "flow_scheme", "kind": "choice", "default": "Coupled",
             "ko": "압력-속도 결합", "en": "Pressure-velocity coupling",
             "choices": _c(("Coupled", "Coupled", "Coupled"), ("SIMPLE", "SIMPLE", "SIMPLE"),
                           ("SIMPLEC", "SIMPLEC", "SIMPLEC"), ("PISO", "PISO", "PISO")),
             "paths": ["flow_scheme"]},
            {"id": "gradient", "kind": "choice", "default": "least-squares-cell-based",
             "ko": "구배 계산", "en": "Gradient",
             "choices": _c(("least-squares-cell-based", "Least squares cell based",
                            "Least squares cell based"),
                           ("green-gauss-node-based", "Green-Gauss node based",
                            "Green-Gauss node based"),
                           ("green-gauss-cell-based", "Green-Gauss cell based",
                            "Green-Gauss cell based")),
             "paths": ["gradient"]},
            {"id": "pressure", "kind": "choice", "default": "second-order",
             "ko": "압력 이산화", "en": "Pressure discretisation",
             "choices": _c(("second-order", "2차", "Second order"),
                           ("standard", "표준", "Standard"),
                           ("presto!", "PRESTO!", "PRESTO!"),
                           ("body-force-weighted", "Body force weighted",
                            "Body force weighted")),
             "paths": ["discretization"]},
            {"id": "momentum", "kind": "choice", "default": "second-order-upwind",
             "ko": "운동량 이산화", "en": "Momentum discretisation",
             "choices": _c(("second-order-upwind", "2차 풍상", "Second order upwind"),
                           ("first-order-upwind", "1차 풍상", "First order upwind"),
                           ("quick", "QUICK", "QUICK"),
                           ("third-order-muscl", "3차 MUSCL", "Third order MUSCL")),
             "paths": ["discretization"]},
            {"id": "turb", "kind": "choice", "default": "second-order-upwind",
             "ko": "난류량 이산화", "en": "Turbulence discretisation",
             "show_if": {"~viscous": ["laminar"]},
             "choices": _c(("second-order-upwind", "2차 풍상", "Second order upwind"),
                           ("first-order-upwind", "1차 풍상", "First order upwind")),
             "paths": ["discretization"]},
            {"id": "pseudo_time", "kind": "bool", "default": True,
             "ko": "Pseudo transient", "en": "Pseudo transient",
             "paths": ["pseudo_time"]},
        ],
    },
    {
        "id": "controls", "ko": "완화 계수", "en": "Under-relaxation",
        "note_ko": "Coupled 를 쓰면 압력·운동량 완화는 Fluent가 다르게 다룹니다. "
                   "빈칸(0)으로 두면 Fluent 기본값을 그대로 씁니다.",
        "note_en": "Under Coupled, Fluent treats the pressure and momentum factors "
                   "differently. Leave a factor at 0 to keep Fluent's own default.",
        "fields": [
            {"id": "urf_pressure", "kind": "number", "default": 0.0, "min": 0, "max": 1,
             "step": 0.05, "ko": "압력", "en": "Pressure",
             "paths": ["under_relaxation"]},
            {"id": "urf_momentum", "kind": "number", "default": 0.0, "min": 0, "max": 1,
             "step": 0.05, "ko": "운동량", "en": "Momentum",
             "paths": ["under_relaxation"]},
            {"id": "urf_k", "kind": "number", "default": 0.0, "min": 0, "max": 1,
             "step": 0.05, "ko": "난류 운동에너지 k", "en": "Turbulent kinetic energy",
             "show_if": {"~viscous": ["laminar"]},
             "paths": ["under_relaxation"]},
            {"id": "urf_omega", "kind": "number", "default": 0.0, "min": 0, "max": 1,
             "step": 0.05, "ko": "비산일률 ω / ε", "en": "Specific dissipation / epsilon",
             "show_if": {"~viscous": ["laminar"]},
             "paths": ["under_relaxation"]},
        ],
    },
    {
        "id": "run", "ko": "초기화 · 반복", "en": "Initialisation & iteration",
        "fields": [
            {"id": "init_method", "kind": "choice", "default": "hybrid",
             "ko": "초기화", "en": "Initialisation",
             "choices": _c(("hybrid", "Hybrid", "Hybrid"),
                           ("standard", "Standard (입구 값)", "Standard (from inlet)")),
             "paths": ["init_type", "hybrid_init", "standard_init"]},
            {"id": "iterations", "kind": "number", "default": 300, "min": 1, "max": 100000,
             "step": 10, "int": True, "ko": "반복 횟수", "en": "Iterations",
             "paths": ["iterate"]},
            {"id": "residual_criterion", "kind": "number", "default": 1e-4, "min": 1e-12,
             "max": 1e-1, "step": 1e-5, "ko": "수렴 판정 잔차", "en": "Residual criterion",
             "paths": ["residual_eqs"]},
            {"id": "write_case", "kind": "bool", "default": True,
             "ko": "끝나면 case+data 저장", "en": "Write case and data when finished",
             "paths": ["write_case_data"]},
        ],
    },
]


def default_settings():
    """The settings dict the UI starts from."""
    out = {}
    for g in SETTINGS:
        out[g["id"]] = {f["id"]: f["default"] for f in g["fields"]}
    return out


def merge_settings(user):
    """Defaults overlaid with whatever the browser sent, ignoring unknown keys."""
    s = default_settings()
    for gid, grp in (user or {}).items():
        if gid in s and isinstance(grp, dict):
            for fid, v in grp.items():
                if fid in s[gid]:
                    s[gid][fid] = v
    return s


def validate(s):
    """Problems that would waste a solver run.  Returns a list of messages."""
    msgs = []
    if s["material"]["density"] <= 0:
        msgs.append("density must be positive")
    if s["material"]["viscosity"] <= 0:
        msgs.append("viscosity must be positive")
    if s["inlet"]["velocity"] <= 0:
        msgs.append("inlet velocity must be positive")
    if int(s["run"]["iterations"]) < 1:
        msgs.append("iterations must be at least 1")
    if s["zones"]["wall_side_y0"] != s["zones"]["wall_side_ymax"]:
        msgs.append("the two side walls have different zone types; that is rarely "
                    "intended and is not symmetric")
    if (s["turbulence"]["viscous"] != "laminar"
            and s["inlet"]["turb_spec"] == "Intensity and Hydraulic Diameter"
            and s["inlet"]["hydraulic_diameter"] <= 0):
        msgs.append("hydraulic diameter is zero but the inlet turbulence "
                    "specification asks for it")
    return msgs


# =============================================================================
#  OFFLINE AUDIT OF EVERY SETTINGS PATH
# =============================================================================
def _resolve(root, path):
    cls = root
    for part in path.split("."):
        nxt = (getattr(cls, "_child_classes", {}) or {}).get(part)
        if nxt is None:
            return None, part
        cls = nxt
        if getattr(cls, "child_object_type", None) is not None:
            cls = cls.child_object_type
    return cls, None


def _resolve_spec(root, spec):
    """True if ANY alternate of the spec resolves; else the failure detail."""
    fails = []
    for path in spec.split("|"):
        cls, missing = _resolve(root, path)
        if cls is not None:
            return True, path
        fails.append("%s (stops at %r)" % (path, missing))
    return False, "; ".join(fails)


def audit_paths(versions=("242", "251", "252", "261", "271")):
    """Check every path the driver can write against the shipped settings trees.

    This is the part of the Fluent integration that CAN be verified without a
    licence, so it runs in the test suite.  A spec passes a release if at least
    one of its alternates resolves there - that is what the alternates are for.
    Returns (n_specs, {version: [(key, detail), ...]}) and also, per version,
    which alternate won, so a silent move between releases is visible.
    """
    used = {}
    for g in SETTINGS:
        for f in g["fields"]:
            for key in f.get("paths", ()):
                used.setdefault(key, set()).add(g["id"] + "." + f["id"])
    for key in ("read_mesh", "write_case_data", "set_zone_type", "iterate",
                "interrupt", "si_area_avg", "si_mass_flow", "si_facet_min",
                "si_facet_max", "si_area", "plane_surface", "inlet_components",
                "hybrid_init", "standard_init", "init_type", "residual_eqs",
                "under_relaxation", "materials"):
        used.setdefault(key, set()).add("driver")

    unknown = sorted(k for k in used if k not in PATHS)
    report, chosen = {}, {}
    for v in versions:
        try:
            mod = importlib.import_module(
                "ansys.fluent.core.generated.solver.settings_" + v)
        except Exception as exc:                       # noqa: BLE001
            report[v] = [("<module>", str(exc))]
            continue
        bad, picks = [], {}
        for key in sorted(used):
            if key not in PATHS:
                bad.append((key, "not declared in PATHS (used by %s)"
                            % ", ".join(sorted(used[key]))))
                continue
            ok, detail = _resolve_spec(mod.root, PATHS[key])
            if ok:
                picks[key] = detail
            else:
                bad.append((key, detail))
        report[v] = bad
        chosen[v] = picks
    return len(used), report, chosen, unknown


# =============================================================================
#  DRIVERS
# =============================================================================
class DriverError(RuntimeError):
    pass


class BaseDriver(object):
    """What the server needs from a solver, real or mock."""

    mock = False

    def __init__(self, case, settings, mesh_path, log):
        self.case = case              # mesh_explorer.Case
        self.s = settings
        self.mesh_path = mesh_path
        self.log = log                # callable(str)
        self.residuals = []           # list of dicts: {"iter": n, <eq>: value}
        self.stopping = False

    # lifecycle ----------------------------------------------------------
    def launch(self):        raise NotImplementedError
    def setup(self):         raise NotImplementedError
    def initialize(self):    raise NotImplementedError
    def iterate(self, n):    raise NotImplementedError
    def close(self):         pass
    # results ------------------------------------------------------------
    def surfaces(self):                       raise NotImplementedError
    def variables(self):                      raise NotImplementedError
    def field(self, surface, variable):       raise NotImplementedError
    def report(self, kind, surfaces, variable): raise NotImplementedError


# --------------------------------------------------------------------------
#  The real thing
# --------------------------------------------------------------------------
class FluentDriver(BaseDriver):
    """Drives ANSYS Fluent through ansys-fluent-core.

    Every settings write below is a path listed in SETTINGS and checked by
    audit_paths(), so a release that moves one is caught offline.
    """

    def __init__(self, *a, **kw):
        BaseDriver.__init__(self, *a, **kw)
        self.solver = None
        self.enum_warnings = []

    # -- lifecycle --------------------------------------------------------
    def launch(self):
        import ansys.fluent.core as pf
        L = self.s["launch"]
        kw = dict(precision=L.get("precision", "double"),
                  processor_count=int(L.get("processors", 4)),
                  mode="solver", dimension=3,
                  ui_mode=("gui" if L.get("ui_mode") == "gui" else "no_gui"))
        if str(L.get("version", "")).strip():
            kw["product_version"] = str(L["version"]).strip()
        self.log("launch_fluent(%s)" % ", ".join("%s=%r" % kv for kv in sorted(kw.items())))
        try:
            self.solver = pf.launch_fluent(**kw)
        except Exception as exc:                       # noqa: BLE001
            raise DriverError(
                "Fluent would not start: %s\n"
                "  If Fluent is installed, check that AWP_ROOT<version> points at it "
                "(PyFluent finds it that way) and that a licence is reachable.\n"
                "  To work on the app itself without Fluent, restart with "
                "--backend mock." % exc)
        self.log("connected: %s" % getattr(self.solver, "get_fluent_version", lambda: "?")())
        self.check_enums()

    def check_enums(self):
        """Ask the live Fluent which strings it accepts, and say so if ours differ.

        The static settings tree carries no allowed values, so this is the only
        place the enumerated strings can be confirmed.  A mismatch is reported,
        not raised: Fluent may legitimately offer a different set per release
        or per enabled model.
        """
        checks = [("viscous_model", _ENUM["viscous"]),
                  ("flow_scheme", _ENUM["flow_scheme"])]
        for path, ours in checks:
            try:
                obj = self._obj(path)
                live = obj.allowed_values()
            except Exception as exc:                   # noqa: BLE001
                self.log("  (could not read allowed values for %s: %s)" % (path, exc))
                continue
            if not live:
                continue
            unknown = [v for v in ours if v not in live]
            if unknown:
                msg = ("%s: this Fluent does not list %s; it offers %s"
                       % (path, unknown, list(live)))
                self.enum_warnings.append(msg)
                self.log("  WARNING " + msg)

    # every settings write goes through these two, so the alternates in PATHS
    # are honoured everywhere and nothing hardcodes one release's spelling
    def _obj(self, key):
        return resolve_obj(self.solver.settings, PATHS[key])

    def _set(self, key, value):
        used = set_path(self.solver.settings, PATHS[key], value)
        if used != PATHS[key].split("|")[0]:
            self.log("  (%s via the fallback path %s)" % (key, used))
        return used

    def _try(self, key, value):
        """Optional setting: log and carry on if this release will not take it."""
        try:
            self._set(key, value)
            return True
        except Exception as exc:                       # noqa: BLE001
            self.log("  %s not settable here (%s)" % (key, exc))
            return False

    # -- setup ------------------------------------------------------------
    def setup(self):
        s, S = self.s, self.solver.settings
        self.log("reading mesh: %s" % self.mesh_path)
        self._obj("read_mesh")(file_name=self.mesh_path)

        g = s["general"]
        self._set("solver_time", "steady" if g["steady"] else "unsteady")
        self._set("op_pressure", float(g["operating_pressure"]))
        self._set("energy", bool(g["energy"]))

        # --- turbulence ---
        t = s["turbulence"]
        self._set("viscous_model", t["viscous"])
        if t["viscous"] == "k-omega":
            self._set("k_omega_model", t["k_omega_variant"])
        elif t["viscous"] == "k-epsilon":
            self._set("k_epsilon_model", t["k_epsilon_variant"])
            if t["wall_treatment"]:
                self._try("wall_treatment", t["wall_treatment"])

        # --- material ---
        m = s["material"]
        name = m["name"]
        fl = self._obj("materials")
        if name not in fl:                       # not in the case yet: copy from the db
            try:
                S.setup.materials.database.copy_by_name(type="fluid", name=name)
            except Exception as exc:             # noqa: BLE001
                self.log("  material %r not in the database (%s); editing the "
                         "existing fluid instead" % (name, exc))
                name = list(fl.keys())[0]
        mat = fl[name]
        mat.density.option = "constant"
        mat.density.value = float(m["density"])
        mat.viscosity.option = "constant"
        mat.viscosity.value = float(m["viscosity"])
        self.log("material %s: rho=%g, mu=%g" % (name, m["density"], m["viscosity"]))
        for z in S.setup.cell_zone_conditions.fluid:
            S.setup.cell_zone_conditions.fluid[z].material = name

        # --- zone types, before the boundary conditions are written ---
        self.apply_zone_types()
        self.apply_boundary_conditions()
        self.apply_methods()
        self.apply_controls()
        self.apply_residuals()

    def apply_zone_types(self):
        z, set_type = self.s["zones"], self._obj("set_zone_type")
        for patch in WALL_PATCHES:
            want = z.get(patch, "wall")
            if want != "wall":
                self.log("zone %s -> %s" % (patch, want))
                set_type(zone_list=[patch], new_type=want)

    def apply_boundary_conditions(self):
        S, s = self.solver.settings, self.s
        lam = s["turbulence"]["viscous"] == "laminar"

        inl = S.setup.boundary_conditions.velocity_inlet["inlet"]
        i = s["inlet"]
        if i["spec"] == "components":
            inl.momentum.velocity_specification_method = "Components"
            # the flat streamwise direction maps exactly to +X for every
            # geometry this tool makes, rolled or not
            comp = inl.momentum.velocity_components
            comp[0] = float(i["velocity"])
            comp[1] = 0.0
            comp[2] = 0.0
        else:
            inl.momentum.velocity_specification_method = "Magnitude, Normal to Boundary"
            inl.momentum.velocity_magnitude = float(i["velocity"])
        if not lam:
            inl.turbulence.turbulence_specification = i["turb_spec"]
            inl.turbulence.turbulent_intensity = float(i["intensity"]) / 100.0
            if i["turb_spec"] == "Intensity and Viscosity Ratio":
                inl.turbulence.turbulent_viscosity_ratio = float(i["visc_ratio"])
            else:
                inl.turbulence.hydraulic_diameter = float(i["hydraulic_diameter"])

        out = S.setup.boundary_conditions.pressure_outlet["outlet"]
        o = s["outlet"]
        out.momentum.gauge_pressure = float(o["gauge_pressure"])
        out.momentum.prevent_reverse_flow = bool(o["prevent_reverse_flow"])
        if not lam:
            out.turbulence.backflow_turbulent_intensity = \
                float(o["backflow_intensity"]) / 100.0
            out.turbulence.backflow_turbulent_viscosity_ratio = \
                float(o["backflow_visc_ratio"])

    def apply_methods(self):
        m = self.s["methods"]
        self._set("flow_scheme", m["flow_scheme"])
        self._set("gradient", m["gradient"])
        # the discretisation schemes are a map keyed by equation name; which
        # keys exist depends on the models that are on, so each is optional
        disc = self._obj("discretization")
        wanted = {"pressure": m["pressure"], "mom": m["momentum"]}
        if self.s["turbulence"]["viscous"] != "laminar":
            wanted["k"] = m["turb"]
            wanted["omega"] = m["turb"]
            wanted["epsilon"] = m["turb"]
        for key, val in wanted.items():
            try:
                disc[key] = val
            except Exception as exc:                   # noqa: BLE001
                self.log("  discretisation %r not settable here (%s)" % (key, exc))
        self._try("pseudo_time", bool(m["pseudo_time"]))

    def apply_controls(self):
        c = self.s["controls"]
        urf = self._obj("under_relaxation")
        for key, sid in (("pressure", "urf_pressure"), ("mom", "urf_momentum"),
                         ("k", "urf_k"), ("omega", "urf_omega"), ("epsilon", "urf_omega")):
            v = float(c.get(sid, 0) or 0)
            if v <= 0:                                  # 0 means "leave Fluent's default"
                continue
            try:
                urf[key] = v
            except Exception as exc:                   # noqa: BLE001
                self.log("  under-relaxation %r not settable (%s)" % (key, exc))

    def apply_residuals(self):
        crit = float(self.s["run"]["residual_criterion"])
        try:
            eqs = self._obj("residual_eqs")
            for name in list(eqs.keys()):
                eqs[name].absolute_criteria = crit
        except Exception as exc:                       # noqa: BLE001
            self.log("  residual criteria not settable (%s)" % exc)

    # -- run --------------------------------------------------------------
    def initialize(self):
        r = self.s["run"]
        if r["init_method"] == "hybrid":
            self.log("hybrid initialisation")
            self._try("init_type", "hybrid")
            self._obj("hybrid_init")()
        else:
            self.log("standard initialisation from the inlet")
            self._try("init_type", "standard")
            try:
                self.solver.settings.solution.initialization.\
                    compute_defaults.velocity_inlet["inlet"]()
            except Exception:                          # noqa: BLE001
                pass
            self._obj("standard_init")()

    def iterate(self, n):
        """Run n iterations, refreshing the residual history as they arrive.

        Fluent owns the residual history, so each iteration event just re-reads
        the monitor rather than trying to accumulate values itself.  If the
        event stream is unavailable the history is still read once at the end,
        so the plot is never empty on account of the events.
        """
        import ansys.fluent.core as pf

        def on_iter(session=None, event_info=None):    # noqa: ARG001
            self.collect_residuals(quiet=True)

        handle = None
        try:
            handle = self.solver.events.register_callback(
                pf.SolverEvent.ITERATION_ENDED, on_iter)
        except Exception as exc:                       # noqa: BLE001
            self.log("  live residuals unavailable (%s); they arrive at the end" % exc)

        self.log("iterating %d" % n)
        try:
            self._obj("iterate")(iter_count=int(n))
        finally:
            if handle is not None:
                try:
                    self.solver.events.unregister_callback(handle)
                except Exception:                      # noqa: BLE001
                    pass
        self.collect_residuals()
        if self.s["run"]["write_case"]:
            base = os.path.splitext(self.mesh_path)[0]
            self.log("writing %s.cas.h5" % base)
            self._obj("write_case_data")(file_name=base + ".cas.h5")

    def collect_residuals(self, quiet=False):
        """Read the residual history out of Fluent."""
        try:
            data = self.solver.monitors.get_monitor_set_data("residual")
        except Exception as exc:                        # noqa: BLE001
            if not quiet:
                self.log("  residual history unavailable (%s)" % exc)
            return
        if not data:
            return
        try:
            xs, series = data[0], data[1]
            rows = []
            for k in range(len(xs)):
                row = {"iter": int(xs[k])}
                for name, vals in series.items():
                    try:
                        row[str(name)] = float(vals[k])
                    except Exception:                   # noqa: BLE001
                        pass
                rows.append(row)
            self.residuals = rows
        except Exception as exc:                        # noqa: BLE001
            if not quiet:
                self.log("  could not read the residual monitor (%s)" % exc)

    def interrupt(self):
        try:
            self._obj("interrupt")()
        except Exception as exc:                        # noqa: BLE001
            self.log("  interrupt failed (%s)" % exc)

    def close(self):
        if self.solver is not None:
            try:
                self.solver.exit()
            except Exception:                           # noqa: BLE001
                pass
            self.solver = None

    # -- results ----------------------------------------------------------
    def surfaces(self):
        """Our own patches first, in mesh order, then anything else Fluent has
        (planes the user made, for instance)."""
        names = list(PATCHES)
        try:
            live = list(self.solver.fields.field_data.get_surface_ids().keys())
        except Exception:                               # noqa: BLE001
            try:
                live = list(self.solver.field_info.get_surfaces_info().keys())
            except Exception:                           # noqa: BLE001
                live = []
        for n in live:
            if n not in names:
                names.append(n)
        return names

    def variables(self):
        return list(VARIABLES)

    def field(self, surface, variable):
        import ansys.fluent.core as pf
        fd = self.solver.fields.field_data
        geo = fd.get_field_data(pf.SurfaceFieldDataRequest(
            surfaces=[surface],
            data_types=[pf.SurfaceDataType.Vertices, pf.SurfaceDataType.FacesConnectivity]))
        sc = fd.get_field_data(pf.ScalarFieldDataRequest(
            surfaces=[surface], field_name=VARIABLES[variable]["fluent"],
            node_value=True, boundary_value=True))
        verts = _as_list(getattr(geo[surface], "vertices", geo[surface]))
        conn = _as_list(getattr(geo[surface], "connectivity", None))
        vals = _as_list(sc[surface] if not hasattr(sc[surface], "scalar_field")
                        else sc[surface].scalar_field)
        return {"surface": surface, "variable": variable,
                "vertices": verts, "faces": conn, "values": vals, "mock": False}

    def report(self, kind, surfaces, variable):
        fld = VARIABLES[variable]["fluent"]
        key = {"area-weighted-avg": "si_area_avg", "mass-flow-rate": "si_mass_flow",
               "facet-min": "si_facet_min", "facet-max": "si_facet_max",
               "area": "si_area"}.get(kind)
        if key is None:
            raise DriverError("unknown report %r" % kind)
        fn = self._obj(key)
        if kind in ("area",):
            return float(fn(surface_names=list(surfaces)))
        return float(fn(surface_names=list(surfaces), report_of=fld))


def _as_list(x):
    if x is None:
        return []
    try:
        return x.tolist()
    except AttributeError:
        return list(x)


#  variable id -> {fluent field name, label}.  The ids are what the wire
#  protocol and the UI use; the Fluent names are what the solver knows.
VARIABLES = {
    "pressure":      {"fluent": "pressure",           "ko": "정압",        "en": "Static pressure",   "unit": "Pa"},
    "total-pressure": {"fluent": "total-pressure",    "ko": "전압",        "en": "Total pressure",    "unit": "Pa"},
    "velocity-magnitude": {"fluent": "velocity-magnitude", "ko": "속도 크기", "en": "Velocity magnitude", "unit": "m/s"},
    "x-velocity":    {"fluent": "x-velocity",         "ko": "x 속도",      "en": "X velocity",        "unit": "m/s"},
    "y-velocity":    {"fluent": "y-velocity",         "ko": "y 속도",      "en": "Y velocity",        "unit": "m/s"},
    "z-velocity":    {"fluent": "z-velocity",         "ko": "z 속도",      "en": "Z velocity",        "unit": "m/s"},
    "turb-kinetic-energy": {"fluent": "turb-kinetic-energy", "ko": "난류 운동에너지 k", "en": "Turbulent kinetic energy", "unit": "m²/s²"},
    "wall-shear":    {"fluent": "wall-shear",         "ko": "벽면 전단응력", "en": "Wall shear stress", "unit": "Pa"},
}


# --------------------------------------------------------------------------
#  The mock
# --------------------------------------------------------------------------
class MockDriver(BaseDriver):
    """No Fluent.  Real mesh, invented field.

    Its whole job is to let the server, the wire format, the contour renderer
    and the report panel be exercised where there is no licence.  Everything it
    produces carries mock:true and the UI shows a badge; do not read a number
    out of it.
    """

    mock = True

    def __init__(self, *a, **kw):
        BaseDriver.__init__(self, *a, **kw)
        self.mesh = None
        self._surf = {}

    def launch(self):
        self.log("MOCK backend - no Fluent is being launched")
        time.sleep(0.05)

    def setup(self):
        self.log("MOCK: building the mesh in-process instead of reading %s"
                 % os.path.basename(self.mesh_path))
        self.mesh = ME.Mesh(self.case)
        self.mesh.build()
        self.log("MOCK: %d cells, %d boundary faces"
                 % (len(self.mesh.hexes), len(self.mesh.boundary)))

    def initialize(self):
        self.log("MOCK: initialised")

    def iterate(self, n):
        """A residual history that falls off like a real one, so the plot,
        the axes and the stop button all get exercised."""
        n = int(n)
        eqs = ["continuity", "x-velocity", "y-velocity", "z-velocity"]
        if self.s["turbulence"]["viscous"] != "laminar":
            eqs += ["k", "omega"]
        crit = float(self.s["run"]["residual_criterion"])
        for k in range(1, n + 1):
            if self.stopping:
                self.log("MOCK: interrupted at iteration %d" % k)
                break
            row = {"iter": k}
            for j, e in enumerate(eqs):
                start = 1.0 * (0.4 ** j)
                row[e] = start * math.exp(-3.5 * k / max(n, 1)) * (
                    1.0 + 0.25 * math.sin(k * (0.7 + 0.11 * j)))
            self.residuals.append(row)
            if k % max(1, n // 40) == 0:
                time.sleep(0.01)                        # let the UI see it stream
            worst = max(row[e] for e in eqs)
            if worst < crit:
                self.log("MOCK: residuals below %g at iteration %d" % (crit, k))
                break
        self.log("MOCK: %d iterations recorded" % len(self.residuals))

    def interrupt(self):
        self.stopping = True

    # -- results ----------------------------------------------------------
    def surfaces(self):
        return list(PATCHES)

    def variables(self):
        return list(VARIABLES)

    def _patch_geometry(self, surface):
        """Vertices and quad connectivity of one patch, in exported metres."""
        if surface in self._surf:
            return self._surf[surface]
        m, c = self.mesh, self.case
        idx, verts, faces = {}, [], []
        for e in m.patches[surface]:
            f = []
            for i in e[0]:
                j = idx.get(i)
                if j is None:
                    j = len(verts)
                    idx[i] = j
                    verts.append(list(c.XP(m.points[i])))
                f.append(j)
            faces.append(f)
        self._surf[surface] = (verts, faces)
        return verts, faces

    def _value(self, variable, p):
        """An invented but smoothly varying field, so contours look like
        something and the colour scale has a real range."""
        s = self.s
        U = float(s["inlet"]["velocity"])
        rho = float(s["material"]["density"])
        x, y, z = p
        c = self.case
        k = c.export_scale
        L = max(c.l_tot * k, 1e-9)
        xn = min(max(x / L, 0.0), 1.0)
        # a gap-velocity bump across the bundle, and a monotonic pressure drop
        bump = 1.0 + 0.9 * math.exp(-((xn - 0.5) ** 2) / 0.05) * (
            0.5 + 0.5 * math.cos(6.0 * math.pi * y / max(c.W * k, 1e-9)))
        umag = U * bump
        dp = -2.2 * rho * U * U * xn
        table = {
            "velocity-magnitude": umag,
            "x-velocity": umag * 0.97,
            "y-velocity": 0.15 * U * math.sin(9.0 * math.pi * y / max(c.W * k, 1e-9)),
            "z-velocity": 0.05 * U * math.cos(5.0 * math.pi * z / max(c.H * k, 1e-9)),
            "pressure": dp,
            "total-pressure": dp + 0.5 * rho * umag * umag,
            "turb-kinetic-energy": 1.5 * (0.05 * umag) ** 2,
            "wall-shear": 0.5 * rho * umag * umag * 0.0135,
        }
        return table.get(variable, 0.0)

    def field(self, surface, variable):
        verts, faces = self._patch_geometry(surface)
        vals = [self._value(variable, v) for v in verts]
        return {"surface": surface, "variable": variable, "vertices": verts,
                "faces": faces, "values": vals, "mock": True}

    def report(self, kind, surfaces, variable):
        tot_a, acc, lo, hi = 0.0, 0.0, float("inf"), float("-inf")
        mdot = 0.0
        rho = float(self.s["material"]["density"])
        for surface in surfaces:
            verts, faces = self._patch_geometry(surface)
            for f in faces:
                p = [verts[i] for i in f]
                a = _quad_area(p)
                cen = [sum(q[i] for q in p) / len(p) for i in range(3)]
                v = self._value(variable, cen)
                tot_a += a
                acc += a * v
                lo = min(lo, v)
                hi = max(hi, v)
                nx = _quad_normal(p)
                mdot += rho * a * self._value("x-velocity", cen) * nx[0]
        if kind == "area":
            return tot_a
        if kind == "facet-min":
            return lo if lo < float("inf") else 0.0
        if kind == "facet-max":
            return hi if hi > float("-inf") else 0.0
        if kind == "mass-flow-rate":
            return mdot
        return acc / tot_a if tot_a > 0 else 0.0


def _quad_normal(p):
    ax = [p[2][i] - p[0][i] for i in range(3)]
    bx = [p[3][i] - p[1][i] for i in range(3)] if len(p) > 3 else \
         [p[2][i] - p[1][i] for i in range(3)]
    n = [0.5 * (ax[1] * bx[2] - ax[2] * bx[1]),
         0.5 * (ax[2] * bx[0] - ax[0] * bx[2]),
         0.5 * (ax[0] * bx[1] - ax[1] * bx[0])]
    m = math.sqrt(sum(v * v for v in n)) or 1.0
    return [v / m for v in n]


def _quad_area(p):
    ax = [p[2][i] - p[0][i] for i in range(3)]
    bx = [p[3][i] - p[1][i] for i in range(3)] if len(p) > 3 else \
         [p[2][i] - p[1][i] for i in range(3)]
    n = [0.5 * (ax[1] * bx[2] - ax[2] * bx[1]),
         0.5 * (ax[2] * bx[0] - ax[0] * bx[2]),
         0.5 * (ax[0] * bx[1] - ax[1] * bx[0])]
    return math.sqrt(sum(v * v for v in n))


def fluent_available():
    """Is the PyFluent module importable?  That is NOT the same as having
    Fluent - see fluent_installed()."""
    try:
        importlib.import_module("ansys.fluent.core")
        return True
    except Exception:                                   # noqa: BLE001
        return False


def fluent_installed():
    """Is there an actual Fluent for PyFluent to launch?

    ansys-fluent-core installs happily from PyPI on a machine with no Ansys on
    it at all, so importing it proves nothing.  PyFluent finds the real thing
    through AWP_ROOT<ver> environment variables; ask it, and report what it
    found, so the app can say up front that it will have to use the mock rather
    than discovering it half a minute into a run.

    Returns (ok, detail).
    """
    if not fluent_available():
        return False, "ansys-fluent-core is not installed (pip install ansys-fluent-core)"
    try:
        from ansys.fluent.core.utils.fluent_version import FluentVersion
        v = FluentVersion.get_latest_installed()
        return True, str(v)
    except Exception as exc:                            # noqa: BLE001
        return False, str(exc).split("\n")[0]


def make_driver(backend, case, settings, mesh_path, log):
    """backend: 'fluent' | 'mock' | 'auto'"""
    if backend == "mock":
        return MockDriver(case, settings, mesh_path, log)
    if backend == "fluent":
        ok, detail = fluent_installed()
        if not ok:
            raise DriverError("cannot use Fluent: " + detail)
        return FluentDriver(case, settings, mesh_path, log)
    ok, detail = fluent_installed()
    if not ok:
        log("no Fluent found (%s) - falling back to the MOCK backend" % detail)
    return (FluentDriver if ok else MockDriver)(case, settings, mesh_path, log)


# =============================================================================
#  THE SAME CASE AS A STANDALONE SCRIPT
# =============================================================================
_VERSION_ALIAS = {"24.2": "242", "25.1": "251", "25.2": "252",
                  "26.1": "261", "27.1": "271",
                  "2024r2": "242", "2025r1": "251", "2025r2": "252",
                  "2026r1": "261", "2027r1": "271"}
JOURNAL_DEFAULT_VERSION = "261"


def normalise_version(v):
    """'2026 R1', '26.1.0', 'v261' -> '261'.  Unknown -> the default."""
    if not v:
        return JOURNAL_DEFAULT_VERSION
    s = str(v).strip().lower().replace(" ", "").replace("-", "").lstrip("v")
    if s in _VERSION_ALIAS:
        return _VERSION_ALIAS[s]
    digits = "".join(ch for ch in s if ch.isdigit())
    for cand in ("242", "251", "252", "261", "271"):
        if digits.startswith(cand):
            return cand
    if len(digits) >= 3 and digits[:3] in ("242", "251", "252", "261", "271"):
        return digits[:3]
    return JOURNAL_DEFAULT_VERSION


def path_for(version, key):
    """The spelling of `key` that this release actually has.

    The journal is meant to be read and kept, so it gets one plain path rather
    than a runtime alternates helper.  Which one is decided here, against the
    settings tree PyFluent ships for that release.
    """
    spec = PATHS[key]
    try:
        mod = importlib.import_module(
            "ansys.fluent.core.generated.solver.settings_" + version)
    except Exception:                                   # noqa: BLE001
        return spec.split("|")[0]
    ok, detail = _resolve_spec(mod.root, spec)
    return detail if ok else spec.split("|")[0]


def journal(geometry, params, settings, mesh_path, version=None):
    """Emit the equivalent standalone PyFluent script.

    The app is convenient; a script is reproducible, diffable and belongs in a
    paper's supplementary material.  It does what the driver does, in the same
    order, written for one named Fluent release so every line is a plain path.
    """
    s = settings
    ver = normalise_version(version or s["launch"].get("version"))
    P_ = lambda key: path_for(ver, key)                 # noqa: E731
    lam = s["turbulence"]["viscous"] == "laminar"
    L, g, t, m, i, o, me, c, r = (s["launch"], s["general"], s["turbulence"],
                                  s["material"], s["inlet"], s["outlet"],
                                  s["methods"], s["controls"], s["run"])
    out = []
    w = out.append
    w('"""Generated by the HCRBPD bundle mesh explorer.')
    w("")
    w("Geometry : %s" % geometry)
    w("Mesh     : %s" % os.path.basename(mesh_path))
    w("Written for Fluent settings API v%s.  On a different release a few of"
      % ver)
    w("these paths move; regenerate this file from the app with that version")
    w("selected and it will use that release's spelling.")
    w('"""')
    w("import ansys.fluent.core as pyfluent")
    w("")
    w("solver = pyfluent.launch_fluent(")
    w("    precision=%r, processor_count=%d," % (L["precision"], int(L["processors"])))
    w("    mode='solver', dimension=3, ui_mode=%r,"
      % ("gui" if L["ui_mode"] == "gui" else "no_gui"))
    if str(L.get("version", "")).strip():
        w("    product_version=%r," % str(L["version"]).strip())
    w(")")
    w("S = solver.settings")
    w("")
    w("S.%s(file_name=%r)" % (P_("read_mesh"), mesh_path))
    w("S.%s = %r" % (P_("solver_time"), "steady" if g["steady"] else "unsteady"))
    w("S.%s = %g" % (P_("op_pressure"), float(g["operating_pressure"])))
    w("S.%s = %r" % (P_("energy"), bool(g["energy"])))
    w("S.%s = %r" % (P_("viscous_model"), t["viscous"]))
    if t["viscous"] == "k-omega":
        w("S.%s = %r" % (P_("k_omega_model"), t["k_omega_variant"]))
    elif t["viscous"] == "k-epsilon":
        w("S.%s = %r" % (P_("k_epsilon_model"), t["k_epsilon_variant"]))
        if t["wall_treatment"]:
            w("S.%s = %r" % (P_("wall_treatment"), t["wall_treatment"]))
    w("")
    w("mat = S.%s[%r]" % (P_("materials"), m["name"]))
    w("mat.density.option = 'constant'")
    w("mat.density.value = %g" % float(m["density"]))
    w("mat.viscosity.option = 'constant'")
    w("mat.viscosity.value = %g" % float(m["viscosity"]))
    w("")
    zoned = [(pt, s["zones"].get(pt, "wall")) for pt in WALL_PATCHES
             if s["zones"].get(pt, "wall") != "wall"]
    if zoned:
        w("# one bay of an infinite array: the sides and/or caps are not walls")
        for patch, want in zoned:
            w("S.%s(zone_list=[%r], new_type=%r)" % (P_("set_zone_type"), patch, want))
        w("")
    w("inlet = S.setup.boundary_conditions.velocity_inlet['inlet']")
    if i["spec"] == "components":
        w("inlet.momentum.velocity_specification_method = 'Components'")
        w("#  the mesh keeps the inlet plane perpendicular to the flow, and the")
        w("#  streamwise direction maps exactly to +X, rolled or not")
        w("inlet.momentum.velocity_components[0] = %g" % float(i["velocity"]))
        w("inlet.momentum.velocity_components[1] = 0.0")
        w("inlet.momentum.velocity_components[2] = 0.0")
    else:
        w("inlet.momentum.velocity_specification_method = 'Magnitude, Normal to Boundary'")
        w("S.%s = %g" % (P_("inlet_magnitude"), float(i["velocity"])))
    if not lam:
        w("inlet.turbulence.turbulence_specification = %r" % i["turb_spec"])
        w("inlet.turbulence.turbulent_intensity = %g" % (float(i["intensity"]) / 100.0))
        if i["turb_spec"] == "Intensity and Viscosity Ratio":
            w("inlet.turbulence.turbulent_viscosity_ratio = %g" % float(i["visc_ratio"]))
        else:
            w("inlet.turbulence.hydraulic_diameter = %g" % float(i["hydraulic_diameter"]))
    w("")
    w("outlet = S.setup.boundary_conditions.pressure_outlet['outlet']")
    w("outlet.momentum.gauge_pressure = %g" % float(o["gauge_pressure"]))
    w("outlet.momentum.prevent_reverse_flow = %r" % bool(o["prevent_reverse_flow"]))
    if not lam:
        w("outlet.turbulence.backflow_turbulent_intensity = %g"
          % (float(o["backflow_intensity"]) / 100.0))
        w("outlet.turbulence.backflow_turbulent_viscosity_ratio = %g"
          % float(o["backflow_visc_ratio"]))
    w("")
    w("S.%s = %r" % (P_("flow_scheme"), me["flow_scheme"]))
    w("S.%s = %r" % (P_("gradient"), me["gradient"]))
    w("disc = S.%s" % P_("discretization"))
    w("disc['pressure'] = %r" % me["pressure"])
    w("disc['mom'] = %r" % me["momentum"])
    if not lam:
        for key in ("k", "omega", "epsilon"):
            w("disc[%r] = %r" % (key, me["turb"]))
    w("S.%s = %r" % (P_("pseudo_time"), bool(me["pseudo_time"])))
    urfs = [(key, float(c.get(sid, 0) or 0))
            for key, sid in (("pressure", "urf_pressure"), ("mom", "urf_momentum"),
                             ("k", "urf_k"), ("omega", "urf_omega"))]
    urfs = [(k, v) for k, v in urfs if v > 0]
    if urfs:
        w("urf = S.%s" % P_("under_relaxation"))
        for key, v in urfs:
            w("urf[%r] = %g" % (key, v))
    w("")
    w("for eq in S.%s:" % P_("residual_eqs"))
    w("    S.%s[eq].absolute_criteria = %g"
      % (P_("residual_eqs"), float(r["residual_criterion"])))
    w("")
    if r["init_method"] == "hybrid":
        w("S.%s = 'hybrid'" % P_("init_type"))
        w("S.%s()" % P_("hybrid_init"))
    else:
        w("S.%s = 'standard'" % P_("init_type"))
        w("S.%s()" % P_("standard_init"))
    w("S.%s(iter_count=%d)" % (P_("iterate"), int(r["iterations"])))
    if r["write_case"]:
        w("S.%s(file_name=%r)"
          % (P_("write_case_data"), os.path.splitext(mesh_path)[0] + ".cas.h5"))
    w("")
    w("avg = S.%s" % P_("si_area_avg"))
    w("dp = (avg(surface_names=['inlet'], report_of='pressure')")
    w("      - avg(surface_names=['outlet'], report_of='pressure'))")
    w("print('bundle pressure drop: %.4g Pa' % dp)")
    w("solver.exit()")
    return "\n".join(out) + "\n"


# =============================================================================
if __name__ == "__main__":
    if "--audit" in sys.argv:
        n, rep, chosen, unknown = audit_paths()
        bad = len(unknown)
        print("auditing %d settings paths against the shipped settings trees\n" % n)
        for u in unknown:
            print("  PATHS has no entry for %r" % u)
        for v in ("242", "251", "252", "261", "271"):
            problems = rep.get(v, [])
            if not problems:
                alt = sum(1 for k, p in chosen.get(v, {}).items()
                          if p != PATHS[k].split("|")[0])
                print("  v%s   all %d resolve"
                      % (v, n) + ("   (%d on a fallback spelling)" % alt if alt else ""))
            else:
                bad += len(problems)
                print("  v%s   %d PROBLEM(S)" % (v, len(problems)))
                for owner, msg in problems:
                    print("        %-22s %s" % (owner, msg))
        if "-v" in sys.argv:
            print("\nwhich alternate each release uses where they differ:")
            for k in sorted(PATHS):
                picks = {v: chosen.get(v, {}).get(k) for v in chosen}
                if len(set(picks.values())) > 1:
                    print("  %s" % k)
                    for v in sorted(picks):
                        print("     v%s  %s" % (v, picks[v]))
        sys.exit(1 if bad else 0)
    if "--schema" in sys.argv:
        print(json.dumps({"settings": SETTINGS, "defaults": default_settings(),
                          "variables": VARIABLES}, ensure_ascii=False, indent=1))
        sys.exit(0)
    print(__doc__)
