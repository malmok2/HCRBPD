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

The enumerated VALUES are audited the same way, by ``audit_choices()``: a path
that resolves can still be handed a string the setting will not take, and that
is a separate failure.  About half of them the shipped trees do publish, and
those are checked offline; the rest the release keeps to itself, so the live
driver asks the running Fluent - once the mesh is read and the objects are
active - about every value THIS run will send, before it sends any of them.
One run then names every bad string instead of dying on the first.

The schema in ``SETTINGS`` is the only place these strings are written down.
A second copy that nothing compares against is how ``least-squares-cell-based``
survived for a scheme Fluent calls ``least-square-cell-based``.
"""

from __future__ import annotations

import difflib
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
    "outlet_turb_spec": "setup.boundary_conditions.pressure_outlet.turbulence."
                        "turbulence_specification",
    "outlet_bf_hyd_diam": "setup.boundary_conditions.pressure_outlet.turbulence."
                          "backflow_hydraulic_diameter",
    #  the bare .material on a cell zone is a deprecated alias PyFluent still
    #  accepts with a warning; .general.material is the real one, and it is the
    #  real one in every release from 2024 R2 on
    "zone_material":   "setup.cell_zone_conditions.fluid.general.material",

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
    "read_case_data":  "file.read_case_data",

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
                           ("scalable-wall-functions", "Scalable 벽함수",
                            "Scalable wall functions"),
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
             #  these two keys are OURS - each drives a different pair of
             #  writes, so the driver cannot just pass the choice through.
             #  api_values names what Fluent is actually given, so the choice
             #  audit checks those strings rather than skipping the field.
             "api_values": {"components": "Components",
                            "magnitude-normal": "Magnitude, Normal to Boundary"},
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
             "note_ko": "퍼센트로 입력하십시오. Fluent settings API는 분율을 받으므로 100으로 "
                        "나누어 전달하고, 적용 후 실제 값을 해석 로그에 되읽어 남깁니다.",
             "note_en": "Enter a percentage. The settings API takes a fraction, so this is "
                        "divided by 100 on the way in, and the value Fluent ends up holding "
                        "is read back into the run log.",
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
             "note_ko": "켜 두면 역류가 없으므로 Fluent가 아래 역류 난류 입력을 "
                        "비활성화합니다. 그래서 이 항목들도 함께 숨겨집니다.",
             "note_en": "With this on there is no backflow, so Fluent deactivates the "
                        "backflow turbulence inputs below - which is why they are hidden.",
             "paths": ["outlet_no_reverse"]},
            #  hidden, and not sent, while reverse flow is prevented: Fluent
            #  refuses every write to the backflow group in that state
            {"id": "backflow_intensity", "kind": "number", "default": 5.0, "min": 0.01,
             "max": 100, "step": 0.1, "unit": "%",
             "ko": "역류 난류 강도", "en": "Backflow turbulent intensity",
             "show_if": {"~viscous": ["laminar"], "~prevent_reverse_flow": [True]},
             "paths": ["outlet_bf_intensity"]},
            {"id": "backflow_visc_ratio", "kind": "number", "default": 10.0, "min": 1e-3,
             "max": 1e5, "step": 1, "ko": "역류 난류 점성비",
             "en": "Backflow turbulent viscosity ratio",
             "show_if": {"~viscous": ["laminar"], "~prevent_reverse_flow": [True]},
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
            #  Fluent spells it "least-SQUARE-cell-based", singular - the
            #  plural is rejected outright.  audit_choices checks it now.
            {"id": "gradient", "kind": "choice", "default": "least-square-cell-based",
             "ko": "구배 계산", "en": "Gradient",
             "choices": _c(("least-square-cell-based", "Least squares cell based",
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


def merge_settings(user, notes=None):
    """Defaults overlaid with whatever the browser sent, ignoring unknown keys.

    A choice the schema does not offer is dropped back to the default rather
    than passed through.  The browser keeps the last case in local storage, so
    a value that WAS valid survives a correction to the schema and gets sent
    again long after the code stopped offering it - which is how
    'least-squares-cell-based' reached Fluent once more after it was fixed.
    Coercions are appended to `notes` if one is given: silently changing a
    setting behind the user's back would be worse than the stale value.
    """
    s = default_settings()
    allowed = {}
    for g in SETTINGS:
        for f in g["fields"]:
            if f.get("kind") == "choice":
                allowed[(g["id"], f["id"])] = ([c["v"] for c in f["choices"]],
                                               f["default"])
    for gid, grp in (user or {}).items():
        if gid not in s or not isinstance(grp, dict):
            continue
        for fid, v in grp.items():
            if fid not in s[gid]:
                continue
            ok, dflt = allowed.get((gid, fid), (None, None))
            if ok is not None and v not in ok:
                if notes is not None:
                    notes.append("%s.%s: %r is not one of this version's "
                                 "choices; using %r" % (gid, fid, v, dflt))
                continue                       # leave the default in place
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


def field(group_id, field_id):
    """One field out of the schema, by group and id."""
    for g in SETTINGS:
        if g["id"] != group_id:
            continue
        for f in g["fields"]:
            if f["id"] == field_id:
                return f
    raise KeyError("%s.%s" % (group_id, field_id))


def api_values(group_id, field_id):
    """What Fluent is given for each of a field's choices.

    Most choices ARE the Fluent string.  A few are ours - the inlet
    specification picks between two different sets of writes, not between two
    values of one setting - and those declare an api_values map.  Going
    through here means the string the audit checks is the string the driver
    sends, rather than the two agreeing by inspection.
    """
    f = field(group_id, field_id)
    m = f.get("api_values")
    return dict(m) if m else {c["v"]: c["v"] for c in f["choices"]}


#  Enum values the driver writes as literals, outside any choice field, with
#  the PATHS key they are written to.  audit_choices checks these too, so a
#  spelling like "unsteady" - which is not an allowed value, "transient" is -
#  is caught here rather than 70 seconds into a run.
LITERAL_ENUMS = [
    ("solver_time", "steady"),          # solver_time() returns one of these
    ("solver_time", "transient"),
    ("init_type", "standard"),
    ("init_type", "hybrid"),
]


def solver_time(settings):
    """The value setup.general.solver.time takes for this case.

    "unsteady" is NOT one of them - the list is steady / transient /
    unsteady-1st-order / unsteady-2nd-order / unsteady-2nd-order-bounded.
    Driver, journal and audit all read it from here so they cannot drift.
    """
    return "steady" if settings["general"]["steady"] else "transient"


def setting_value(settings, field_id):
    """A field's current value, found by id across the groups (ids are unique)."""
    for g in SETTINGS:
        for f in g["fields"]:
            if f["id"] == field_id:
                return settings.get(g["id"], {}).get(field_id, f.get("default"))
    return None


def visible(field_def, settings):
    """The Python twin of the page's fieldVisible().

    show_if is {field: [values]} to show only for those values, or
    {"~field": [values]} to hide for them.  Both sides evaluate it, so a field
    the panel hides is also a field this does not send or check.
    """
    cond = field_def.get("show_if")
    if not cond:
        return True
    for key, allowed in cond.items():
        neg = key.startswith("~")
        v = setting_value(settings, key[1:] if neg else key)
        hit = v in allowed
        if (hit if neg else not hit):
            return False
    return True


def planned_enums(settings):
    """Every enumerated string one run will send, with the path it goes to.

    Offline, audit_choices checks the whole schema against the shipped trees.
    This is the same question asked of a run: only the values it will actually
    use, so a live session can be asked about exactly those.
    """
    out = []
    for g in SETTINGS:
        got = settings.get(g["id"], {})
        for f in g["fields"]:
            if f.get("kind") != "choice" or not f.get("paths"):
                continue
            v = got.get(f["id"], f.get("default"))
            v = api_values(g["id"], f["id"]).get(v, v)
            if v == "" or v is None:
                continue
            if not visible(f, settings):
                continue
            out.append(("%s.%s" % (g["id"], f["id"]), f["paths"][0], v))
    #  the only enum the driver derives rather than passing through
    out.append(("general.steady", "solver_time", solver_time(settings)))
    return out


def audit_choices(versions=("242", "251", "252", "261", "271")):
    """Check every enum STRING the app can send against the shipped trees.

    audit_paths proves the setting exists; this proves the value is one the
    setting will take.  They fail differently and they failed separately in
    practice - the paths audit was green while the gradient scheme was spelled
    'least-squares-cell-based' for a scheme Fluent calls 'least-square-...'.

    Not every release publishes its allowed values statically: where a setting
    does not carry them, it is reported as unchecked rather than passed.
    Returns (n_checked, unchecked_set, {version: [(what, detail)]}).
    """
    wanted = []                       # (label, paths-key, value)
    for g in SETTINGS:
        for f in g["fields"]:
            if f.get("kind") != "choice" or not f.get("paths"):
                continue
            amap = api_values(g["id"], f["id"])
            for c in f["choices"]:
                v = amap.get(c["v"], c["v"])
                if v == "":           # blank means "leave Fluent's default"
                    continue
                wanted.append(("%s.%s" % (g["id"], f["id"]), f["paths"][0], v))
    for key, v in LITERAL_ENUMS:
        wanted.append(("driver." + key, key, v))

    report, checked, unchecked = {}, set(), set()
    for ver in versions:
        try:
            mod = importlib.import_module(
                "ansys.fluent.core.generated.solver.settings_" + ver)
        except Exception as exc:                       # noqa: BLE001
            report[ver] = [("<module>", str(exc))]
            continue
        bad = []
        for label, key, val in wanted:
            cls = None
            for alt in PATHS[key].split("|"):
                cls, _ = _resolve(mod.root, alt)
                if cls is not None:
                    break
            allowed = getattr(cls, "_allowed_values", None) if cls is not None else None
            if not allowed:
                unchecked.add((label, val))
                continue
            checked.add((label, val))
            if val not in allowed:
                near = difflib.get_close_matches(val, allowed, n=1, cutoff=0.6)
                bad.append((label, "%r is not allowed%s; allowed: %s"
                            % (val, (" (did you mean %r?)" % near[0]) if near else "",
                               ", ".join(allowed))))
        report[ver] = bad
    return len(checked), unchecked - checked, report


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
    for key in ("read_mesh", "write_case_data", "read_case_data",
                "set_zone_type", "iterate",
                "interrupt", "si_area_avg", "si_mass_flow", "si_facet_min",
                "si_facet_max", "si_area", "plane_surface", "inlet_components",
                "hybrid_init", "standard_init", "init_type", "residual_eqs",
                "under_relaxation", "materials", "zone_material",
                "outlet_turb_spec", "outlet_bf_hyd_diam"):
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
    def load_case(self, path):                raise NotImplementedError
    def bbox(self):                           raise NotImplementedError
    def make_plane(self, name, axis, value):  raise NotImplementedError
    def drop_plane(self, name):               raise NotImplementedError
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
        self._surface_listing_errors = []
        self._said_no_listing = False

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

    def check_enums(self):
        """Ask the live Fluent about every enum string this run will send.

        audit_choices settles offline what the shipped settings trees publish,
        but they publish only about half of it, and a value can move between
        releases.  This asks the session in front of us, about exactly the
        values it is about to be given.

        Two things matter about WHEN.  The objects have to be active, which
        means after the mesh is read - asked any earlier, Fluent answers
        "setup.models is currently inactive" and nothing is checked at all.
        And every value is checked before any is written, so one run reports
        every bad string instead of dying on the first and hiding the rest.

        Reported, not raised: Fluent legitimately offers different sets per
        release and per enabled model, and the write that follows will fail
        on its own if the value really is wrong - by which point the log
        already says so, with the spelling Fluent wants.
        """
        rows, seen = planned_enums(self.s), {}
        for label, key, val in rows:
            if key not in PATHS:
                continue
            if key not in seen:
                try:
                    seen[key] = list(self._obj(key).allowed_values() or [])
                except Exception as exc:               # noqa: BLE001
                    seen[key] = None
                    self.log("  (no allowed values for %s: %s)" % (key, exc))
            live = seen[key]
            if not live or val in live:
                continue
            near = difflib.get_close_matches(val, live, n=1, cutoff=0.6)
            msg = ("%s: this Fluent does not accept %r%s (it offers %s)"
                   % (label, val, (" - did you mean %r?" % near[0]) if near else "",
                      ", ".join(live)))
            self.enum_warnings.append(msg)
            self.log("  WARNING " + msg)
        n = sum(1 for k in seen if seen[k])
        self.log("  checked %d enum value(s) against this Fluent across %d setting(s)%s"
                 % (sum(1 for l, k, v in rows if seen.get(k)), n,
                    "" if not self.enum_warnings
                    else " - %d PROBLEM(S) above" % len(self.enum_warnings)))

    # every settings write goes through these two, so the alternates in PATHS
    # are honoured everywhere and nothing hardcodes one release's spelling
    def _obj(self, key):
        return resolve_obj(self.solver.settings, PATHS[key])

    def _set(self, key, value):
        used = set_path(self.solver.settings, PATHS[key], value)
        if used != PATHS[key].split("|")[0]:
            self.log("  (%s via the fallback path %s)" % (key, used))
        return used

    def _soft(self, what, fn):
        """Run a settings write that may legitimately be refused.

        Fluent deactivates an input when the model state makes it meaningless -
        backflow turbulence before a specification method is chosen, say.  That
        is worth reporting but is not a reason to abandon a set-up that is
        otherwise complete.
        """
        try:
            fn()
            return True
        except Exception as exc:                       # noqa: BLE001
            self.log("  %s was refused (%s)" % (what, exc))
            return False

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
        #  only now are the model and method objects active enough to answer
        self.check_enums()

        g = s["general"]
        self._set("solver_time", solver_time(s))
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
        zones = S.setup.cell_zone_conditions.fluid
        for z in zones:
            zones[z].general.material = name

        # --- zone types, before the boundary conditions are written ---
        self.apply_zone_types()
        self.apply_boundary_conditions()
        self.apply_methods()
        self.apply_controls()
        self.apply_residuals()
        self.verify_setup()

    def verify_setup(self):
        """Read the important settings back out of Fluent and log them.

        Writing a setting is not the same as Fluent holding what you meant, and
        a wrong turbulent intensity is invisible in the residuals and obvious
        only in the answer.

        The intensity unit is the case in point.  The settings API takes it as
        a FRACTION while the panel asks for a percentage, so this divides by
        100 - and that is documented in neither the API nor the shipped
        examples.  It was settled by measurement instead: writing 0.05 through
        the API put 5 in the Turbulent Intensity box of Fluent 2025 R1.  The
        read-back below is checked against what was sent, so if a release ever
        changes the convention this says so rather than quietly running a case
        at a hundredth of the intended turbulence.
        """
        S, s = self.solver.settings, self.s
        lam = s["turbulence"]["viscous"] == "laminar"
        rows = []

        def read(label, fn):
            try:
                rows.append((label, fn()))
            except Exception as exc:                   # noqa: BLE001
                rows.append((label, "<unreadable: %s>" % exc))

        inl = S.setup.boundary_conditions.velocity_inlet["inlet"]
        out = S.setup.boundary_conditions.pressure_outlet["outlet"]
        read("inlet spec", lambda: inl.momentum.velocity_specification_method())
        if s["inlet"]["spec"] == "components":
            read("inlet u,v,w", lambda: [inl.momentum.velocity_components[k]()
                                         for k in range(3)])
        else:
            read("inlet |u|", lambda: resolve_obj(S, PATHS["inlet_magnitude"])())
        read("outlet gauge p", lambda: out.momentum.gauge_pressure())
        if not lam:
            read("inlet turb spec", lambda: inl.turbulence.turbulence_specification())
            read("inlet intensity  (sent %g%% as %g)"
                 % (float(s["inlet"]["intensity"]),
                    float(s["inlet"]["intensity"]) / 100.0),
                 lambda: inl.turbulence.turbulent_intensity())
            if not s["outlet"]["prevent_reverse_flow"]:
                read("outlet turb spec",
                     lambda: out.turbulence.turbulence_specification())
                read("outlet backflow intensity  (sent %g%% as %g)"
                     % (float(s["outlet"]["backflow_intensity"]),
                        float(s["outlet"]["backflow_intensity"]) / 100.0),
                     lambda: out.turbulence.backflow_turbulent_intensity())
        read("viscous model", lambda: resolve_obj(S, PATHS["viscous_model"])())
        zones = S.setup.cell_zone_conditions.fluid
        read("cell zone material",
             lambda: [zones[z].general.material() for z in zones])

        self.log("--- what Fluent holds after set-up ---")
        for label, val in rows:
            self.log("    %-46s %s" % (label, val))
        if not lam:
            self.check_intensity_unit(rows)

    def check_intensity_unit(self, rows):
        """Did Fluent keep the intensity we sent, in the unit we sent it in?

        Verified against Fluent 2025 R1: 0.05 written through the settings API
        shows as 5 in the panel's Turbulent Intensity box, so the API unit is
        the fraction and the /100 here is right.  Nothing in the API states
        that, though, so it is worth one comparison per run: if a release ever
        switched to percent the read-back would come back a hundred times the
        value sent, and a case would silently run at 0.05 % turbulence.
        """
        bad, seen = [], 0
        for label, val in rows:
            if "(sent " not in label or not isinstance(val, (int, float)):
                continue
            seen += 1
            sent = float(label.rsplit(" as ", 1)[1].rstrip(")"))
            if sent <= 0:
                continue
            if abs(val - sent) <= 1e-9 + 1e-6 * sent:
                continue                                # holding what we sent
            bad.append((label.split("  (sent")[0].strip(), sent, val,
                        "looks like PERCENT, not the fraction"
                        if abs(val - 100.0 * sent) <= 1e-6 * 100.0 * sent
                        else "neither the fraction nor the percentage"))
        if not seen:
            #  nothing came back to compare against - say that, rather than
            #  reporting a check that never happened
            self.log("    intensity unit NOT checked: Fluent returned no value")
            return
        if not bad:
            self.log("    intensity unit checked: Fluent is holding the "
                     "fraction that was sent (0.05 = 5 %)")
            return
        for what, sent, got, why in bad:
            self.log("  [WARNING] %s: sent %g, Fluent holds %g - %s."
                     % (what, sent, got, why))
        self.log("  [WARNING] the turbulence level of this run is NOT what the "
                 "Settings tab asked for; treat the result as suspect.")

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
        spec = api_values("inlet", "spec")          # the audited strings
        if i["spec"] == "components":
            inl.momentum.velocity_specification_method = spec["components"]
            # the flat streamwise direction maps exactly to +X for every
            # geometry this tool makes, rolled or not
            comp = inl.momentum.velocity_components
            comp[0] = float(i["velocity"])
            comp[1] = 0.0
            comp[2] = 0.0
        else:
            inl.momentum.velocity_specification_method = spec["magnitude-normal"]
            inl.momentum.velocity_magnitude = float(i["velocity"])
        if not lam:
            #  specification first, then the inputs it activates
            self._soft("inlet turbulence specification", lambda: setattr(
                inl.turbulence, "turbulence_specification", i["turb_spec"]))
            #  a FRACTION, not a percentage: 0.05 written here shows as 5 in
            #  the panel's Turbulent Intensity box (measured on 2025 R1).  The
            #  API says nothing about the unit, so verify_setup re-reads it.
            self._soft("inlet turbulent intensity", lambda: setattr(
                inl.turbulence, "turbulent_intensity", float(i["intensity"]) / 100.0))
            if i["turb_spec"] == "Intensity and Viscosity Ratio":
                self._soft("inlet turbulent viscosity ratio", lambda: setattr(
                    inl.turbulence, "turbulent_viscosity_ratio",
                    float(i["visc_ratio"])))
            else:
                self._soft("inlet hydraulic diameter", lambda: setattr(
                    inl.turbulence, "hydraulic_diameter",
                    float(i["hydraulic_diameter"])))

        out = S.setup.boundary_conditions.pressure_outlet["outlet"]
        o = s["outlet"]
        out.momentum.gauge_pressure = float(o["gauge_pressure"])
        no_backflow = bool(o["prevent_reverse_flow"])
        out.momentum.prevent_reverse_flow = no_backflow
        if not lam and no_backflow:
            #  With reverse flow prevented there IS no backflow, so Fluent
            #  deactivates the whole backflow turbulence group - specification
            #  method included - and refuses every write to it.  Writing them
            #  anyway produced three "the object is not active" errors per run
            #  for values that could not have had any effect.
            self.log("  outlet: reverse flow is prevented, so Fluent has no "
                     "backflow turbulence to set; skipping those three inputs")
        if not lam and not no_backflow:
            #  The specification method has to be chosen FIRST: until it is,
            #  Fluent keeps the backflow inputs inactive and rejects a write to
            #  them with "the object is not active".  The outlet follows the
            #  inlet's choice - a case that specifies intensity and hydraulic
            #  diameter going in and viscosity ratio coming back would be odd.
            self._soft("outlet turbulence specification", lambda: setattr(
                out.turbulence, "turbulence_specification", i["turb_spec"]))
            self._soft("backflow turbulent intensity", lambda: setattr(
                out.turbulence, "backflow_turbulent_intensity",
                float(o["backflow_intensity"]) / 100.0))   # fraction, as above
            if i["turb_spec"] == "Intensity and Viscosity Ratio":
                self._soft("backflow turbulent viscosity ratio", lambda: setattr(
                    out.turbulence, "backflow_turbulent_viscosity_ratio",
                    float(o["backflow_visc_ratio"])))
            else:
                self._soft("backflow hydraulic diameter", lambda: setattr(
                    out.turbulence, "backflow_hydraulic_diameter",
                    float(i["hydraulic_diameter"])))

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
    #  Accessors that list the surfaces a live session has, newest spelling
    #  first.  get_surface_ids() is NOT one of them: it takes a list of names
    #  and returns ids, so calling it bare raised every time and the result
    #  was silently no surfaces at all - which is why a plane that had been
    #  created never appeared.  Each entry returns an iterable of names.
    _SURFACE_LISTERS = (
        ("field_data._allowed_surface_names",
         lambda S: S.fields.field_data._allowed_surface_names()),
        ("settings.results.surfaces.plane_surface",
         lambda S: list(S.settings.results.surfaces.plane_surface)),
        ("field_info.get_surfaces_info",
         lambda S: list(S.fields.field_info.get_surfaces_info())),
    )

    def surfaces(self):
        """Our patches in mesh order, then everything else this session has.

        The boundaries are known from the mesh we wrote, so they never depend
        on Fluent being able to enumerate anything.  Whatever else is there -
        a plane made here, or one the user made in Fluent - comes from the
        listing below, and if none of the accessors work that is said once
        rather than quietly returning a short list.
        """
        names = list(PATCHES)
        live, how = [], None
        for label, fn in self._SURFACE_LISTERS:
            try:
                live = [str(n) for n in fn(self.solver)]
                how = label
                break
            except Exception as exc:                    # noqa: BLE001
                self._surface_listing_errors.append("%s: %s" % (label, exc))
        if how is None and not self._said_no_listing:
            self._said_no_listing = True
            self.log("  (this Fluent lists no surfaces through PyFluent; the "
                     "boundaries are known from the mesh and planes made here "
                     "are tracked by the app, so nothing is lost unless you "
                     "made a surface inside Fluent itself)")
            for e in self._surface_listing_errors[:3]:
                self.log("      tried %s" % e)
        for n in live:
            if n not in names:
                names.append(n)
        return names

    #  method name per axis: the plane is the one the OTHER two axes span.
    #  The same spelling in every release from 2024 R2 to 2027 R1.
    PLANE_METHOD = {"x": "yz-plane", "y": "zx-plane", "z": "xy-plane"}

    def make_plane(self, name, axis, value):
        """A constant-x, -y or -z cut through the domain, as a named surface.

        Fluent then treats it exactly like a boundary patch, so the contour,
        the surface integrals and the field request need no special case for
        it - a plane is just another name in the surface list.
        """
        if axis not in self.PLANE_METHOD:
            raise DriverError("axis must be x, y or z, not %r" % axis)
        ps = self._obj("plane_surface")
        ps[name] = {}                          # create, then configure
        pl = ps[name]
        self._soft("plane method",
                   lambda: setattr(pl, "method", self.PLANE_METHOD[axis]))
        setattr(pl, axis, float(value))
        self.log("plane %s: %s = %g" % (name, axis, value))
        return name

    def drop_plane(self, name):
        ps = self._obj("plane_surface")
        try:
            del ps[name]
            return True
        except Exception as exc:                        # noqa: BLE001
            self.log("  could not delete the plane %s (%s)" % (name, exc))
            return False

    def bbox(self):
        """The domain's bounding box in Fluent's own coordinates.

        Every boundary except the rods together wrap the domain, so their
        vertices bound it exactly - and they are small patches, which is why
        this asks for those rather than for the cells.
        """
        import ansys.fluent.core as pf
        want = [n for n in PATCHES if n != "wall_rods"]
        lo = [1e30] * 3
        hi = [-1e30] * 3
        fd = self.solver.fields.field_data
        for name in want:
            try:
                geo = fd.get_field_data(pf.SurfaceFieldDataRequest(
                    surfaces=[name], data_types=[pf.SurfaceDataType.Vertices]))
                verts = _as_list(getattr(geo[name], "vertices", None))
            except Exception as exc:                    # noqa: BLE001
                self.log("  (no vertices for %s: %s)" % (name, exc))
                continue
            for q in verts:
                for k in range(3):
                    if q[k] is None:
                        continue
                    if q[k] < lo[k]:
                        lo[k] = q[k]
                    if q[k] > hi[k]:
                        hi[k] = q[k]
        if lo[0] > hi[0]:
            raise DriverError("could not measure the domain: no boundary "
                              "vertices came back from Fluent")
        return [lo, hi]

    def load_case(self, path):
        """Reopen a written case+data instead of meshing and iterating again.

        Everything downstream - surfaces, planes, fields, reports - works off
        the loaded solution exactly as it does off a fresh one, because none
        of it knows or cares how the data got into Fluent.
        """
        self.log("reading case and data: %s" % path)
        self._obj("read_case_data")(file_name=path)
        self.log("loaded; surfaces: %s" % ", ".join(self.surfaces()[:8]))
        self.collect_residuals(quiet=True)
        return True

    def variables(self):
        return list(VARIABLES)

    def field(self, surface, variable):
        import ansys.fluent.core as pf
        fd = self.solver.fields.field_data
        types = [pf.SurfaceDataType.Vertices, pf.SurfaceDataType.FacesConnectivity]
        #  The structured connectivity is deprecated in favour of the flat one
        #  and warns on every call; ask for flat where the request takes it.
        #  Either shape is decoded below, so the fallback is not a second
        #  code path so much as a second spelling of the same request.
        flat = True
        try:
            req = pf.SurfaceFieldDataRequest(surfaces=[surface], data_types=types,
                                             flatten_connectivity=True)
        except TypeError:
            flat = False
            req = pf.SurfaceFieldDataRequest(surfaces=[surface], data_types=types)
        geo = fd.get_field_data(req)
        sc = fd.get_field_data(pf.ScalarFieldDataRequest(
            surfaces=[surface], field_name=VARIABLES[variable]["fluent"],
            node_value=True, boundary_value=True))
        sd = geo[surface]
        verts = _as_list(getattr(sd, "vertices", None))
        faces = _faces(getattr(sd, "connectivity", None), flat)
        vals = _as_list(getattr(sc[surface], "scalar_field", sc[surface]))
        return {"surface": surface, "variable": variable,
                "vertices": verts, "faces": faces, "values": vals, "mock": False}

    def report(self, kind, surfaces, variable):
        fld = VARIABLES[variable]["fluent"]
        key = {"area-weighted-avg": "si_area_avg", "mass-flow-rate": "si_mass_flow",
               "facet-min": "si_facet_min", "facet-max": "si_facet_max",
               "area": "si_area"}.get(kind)
        if key is None:
            raise DriverError("unknown report %r" % kind)
        fn = self._obj(key)
        if kind in ("area",):
            return _scalar(fn(surface_names=list(surfaces)), kind)
        return _scalar(fn(surface_names=list(surfaces), report_of=fld), kind)


def _scalar(x, what):
    """One number out of whatever a surface integral hands back.

    Some releases return the bare float, some a one-element array, some a
    dict keyed by surface.  float() on the array forms raises a numpy message
    that says nothing about which report failed, so unwrap first and name the
    report if there is still no single number in there.
    """
    v = _plain(x)
    while isinstance(v, (list, tuple)) and len(v) == 1:
        v = v[0]
    if isinstance(v, dict) and len(v) == 1:
        v = _plain(list(v.values())[0])
        while isinstance(v, (list, tuple)) and len(v) == 1:
            v = v[0]
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise DriverError("the %s report came back as %s, not a number: %r"
                          % (what, type(x).__name__, x))
    return float(v)


def _faces(conn, flat):
    """Face connectivity as a list of vertex-index lists, from either shape.

    Flat is one array per surface: [n, v0..vn-1, n, v0..., ...].  Structured
    is a list of arrays, one per face.  Telling them apart by looking is
    unreliable - a structured surface of triangles and a flat array both come
    out as a sequence of ints - so the caller says which was asked for, and
    the flat decode still checks its own counts rather than trusting them.
    """
    a = _as_list(conn)
    if not a:
        return []
    if not flat:
        return [list(f) for f in a]
    out, i, n = [], 0, len(a)
    while i < n:
        k = int(a[i])
        if k <= 0 or i + k >= n:
            raise DriverError(
                "face connectivity is malformed at index %d: a face of %r "
                "vertices does not fit in the %d values left" % (i, a[i], n - i - 1))
        out.append([int(v) for v in a[i + 1:i + 1 + k]])
        i += k + 1
    return out


def _as_list(x):
    """Whatever PyFluent hands back, as plain JSON-able Python.

    It is not one shape.  Vertices come back as an (N,3) ndarray, a scalar
    field as a 1-D ndarray, but face connectivity as a LIST of 1-D ndarrays,
    one per face - and a list has no .tolist(), so a shallow conversion left
    the inner arrays as numpy and the error surfaced two layers away, in
    json.dumps, as "only 0-dimensional arrays can be converted to Python
    scalars".  Recursing costs nothing and the shape stops mattering.

    Non-finite values become None: a diverged run should draw a hole and say
    so, not abort the whole payload on a single NaN.
    """
    if x is None:
        return []
    return _plain(x)


def _plain(x):
    if x is None:
        return None
    tolist = getattr(x, "tolist", None)         # ndarray, and numpy scalars
    if tolist is not None:
        return _plain(tolist())
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, (int, str, bool)) or x is None:
        return x
    item = getattr(x, "item", None)             # any other numpy-ish scalar
    if item is not None:
        return _plain(item())
    return x


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
        self._planes = {}

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
        #  Only the boundaries, deliberately: a real session may not be able
        #  to enumerate the planes it holds, and the app has to survive that
        #  by tracking its own.  Listing them here would hide the case that
        #  actually broke - a plane created in Fluent, never shown.
        return list(PATCHES)

    def variables(self):
        return list(VARIABLES)

    #  planes the mock has been asked for: name -> (axis, exported value)
    def make_plane(self, name, axis, value):
        if axis not in ("x", "y", "z"):
            raise DriverError("axis must be x, y or z, not %r" % axis)
        self._planes[name] = (axis, float(value))
        self._surf.pop(name, None)
        self.log("MOCK plane %s: %s = %g" % (name, axis, value))
        return name

    def drop_plane(self, name):
        self._planes.pop(name, None)
        self._surf.pop(name, None)
        return True

    def load_case(self, path):
        self.log("MOCK: pretending to read %s" % os.path.basename(path))
        if self.mesh is None:
            self.setup()
        return True

    def bbox(self):
        """The domain's exported bounding box, which the plane sliders span."""
        if self.mesh is None:
            self.setup()
        c = self.case
        lo = [1e30] * 3
        hi = [-1e30] * 3
        for q in self.mesh.points:
            e = c.XP(q)
            for k in range(3):
                if e[k] < lo[k]:
                    lo[k] = e[k]
                if e[k] > hi[k]:
                    hi[k] = e[k]
        return [lo, hi]

    def _flat_axis(self, axis, want):
        """The FLAT coordinate whose exported image is `want` on this axis.

        The mock cuts the mesh in its own flat space, but the UI works in the
        exported coordinates Fluent would use, so the two have to be tied
        together.  For the rod family the map is a scale and this is exact;
        for a rolled coil the export mixes the axes and this is the value at
        the middle of the other two, which is why a mock plane on a coil is
        indicative rather than exact.  It is a mock.
        """
        c = self.case
        k = "xyz".index(axis)
        span = {"x": c.l_tot, "y": c.W, "z": c.H}[axis]
        mid = [c.l_tot * 0.5, c.W * 0.5, c.H * 0.5]

        def exported(t):
            q = list(mid)
            q[k] = t
            return c.XP(q)[k]

        lo, hi = 0.0, span
        f_lo, f_hi = exported(lo), exported(hi)
        if f_hi < f_lo:
            lo, hi, f_lo, f_hi = hi, lo, f_hi, f_lo
        if want <= f_lo:
            return lo
        if want >= f_hi:
            return hi
        for _ in range(60):                       # monotonic: plain bisection
            mid_t = 0.5 * (lo + hi)
            if exported(mid_t) < want:
                lo = mid_t
            else:
                hi = mid_t
        return 0.5 * (lo + hi)

    def _plane_geometry(self, name):
        """A structured cut across the domain, with the rods punched out.

        Built in flat space and then exported, so the rod footprint is a plain
        ellipse test and the cut of a rolled coil comes out curved, the way
        the real one does.
        """
        axis, want = self._planes[name]
        c, m = self.case, self.mesh
        k = "xyz".index(axis)
        t = self._flat_axis(axis, want)
        free = [i for i in range(3) if i != k]
        span = [c.l_tot, c.W, c.H]
        n = [132, 88]
        idx, verts, faces = {}, [], []

        def node(a, b):
            key = (a, b)
            j = idx.get(key)
            if j is None:
                q = [0.0, 0.0, 0.0]
                q[k] = t
                q[free[0]] = span[free[0]] * a / n[0]
                q[free[1]] = span[free[1]] * b / n[1]
                j = len(verts)
                idx[key] = j
                verts.append(list(c.XP(q)))
            return j

        def in_rod(a, b):
            """Is the centre of this cell inside a rod footprint?"""
            q = [0.0, 0.0, 0.0]
            q[k] = t
            q[free[0]] = span[free[0]] * (a + 0.5) / n[0]
            q[free[1]] = span[free[1]] * (b + 0.5) / n[1]
            for (cx, cy) in m.centres:
                if ((q[0] - cx) / c.aE) ** 2 + ((q[1] - cy) / c.bE) ** 2 < 1.0:
                    return True
            return False

        for a in range(n[0]):
            for b in range(n[1]):
                if in_rod(a, b):
                    continue                    # a hole, not a cell
                faces.append([node(a, b), node(a + 1, b),
                              node(a + 1, b + 1), node(a, b + 1)])
        return verts, faces

    def _patch_geometry(self, surface):
        """Vertices and quad connectivity of one surface, in exported metres.

        A plane the user made is a surface like any other from here on, which
        is what lets the field, the contour and the reports treat it the same
        way Fluent does.
        """
        if surface in self._planes:
            if surface not in self._surf:
                self._surf[surface] = self._plane_geometry(surface)
            return self._surf[surface]
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


# =============================================================================
#  FIELD SNAPSHOTS
# =============================================================================
#  A .cas.h5 can only be reopened by Fluent, which means a licence and a
#  minute of waiting every time you want to look at a picture you already
#  made.  A snapshot is the other half: the sampled field itself - the
#  vertices, the faces and one value per NODE per variable, for the surfaces
#  and planes you were looking at.  It reloads with no solver at all.
#
#  It is a record of what was sampled, not a substitute for the case: it
#  carries the surfaces that were in it and nothing else, so it says which
#  those were and when, and the app labels a reloaded view accordingly.
SNAPSHOT_VERSION = 1
SNAPSHOT_EXT = ".fields.json"


def snapshot(driver, surfaces, variables, meta=None):
    """Sample `variables` on `surfaces` and return the snapshot dict.

    The geometry of a surface is fetched once and shared by every variable on
    it - it is the same mesh - so adding a variable costs one array, not a
    whole copy of the surface.
    """
    out = {"snapshot_version": SNAPSHOT_VERSION,
           "created": time.strftime("%Y-%m-%d %H:%M:%S"),
           "variables": list(variables), "surfaces": {},
           "mock": bool(getattr(driver, "mock", False))}
    out.update(meta or {})
    for name in surfaces:
        first = driver.field(name, variables[0])
        entry = {"vertices": first["vertices"], "faces": first["faces"],
                 "values": {variables[0]: first["values"]}}
        for v in variables[1:]:
            entry["values"][v] = driver.field(name, v)["values"]
        out["surfaces"][name] = entry
    return out


def write_snapshot(path, data):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(data, fh, separators=(",", ":"), allow_nan=False)
    return path


def read_snapshot(path):
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    v = data.get("snapshot_version")
    if v != SNAPSHOT_VERSION:
        raise DriverError("%s was written by a different version of this tool "
                          "(snapshot_version %r, this one reads %d)"
                          % (os.path.basename(path), v, SNAPSHOT_VERSION))
    if not isinstance(data.get("surfaces"), dict) or not data["surfaces"]:
        raise DriverError("%s carries no surfaces" % os.path.basename(path))
    return data


class SnapshotDriver(BaseDriver):
    """Serves a saved snapshot.  No Fluent, no mesh, no solver.

    It answers the same three questions the Results tab asks - which surfaces,
    which variables, and the field on one of them - and refuses the rest,
    because a snapshot is a record of what was sampled and cannot be asked for
    anything that was not.
    """

    mock = False

    def __init__(self, data, log):
        self.data = data
        self.log = log
        self.residuals = []
        self.snapshot = True

    def launch(self):
        n = len(self.data["surfaces"])
        self.log("snapshot from %s: %d surface(s), %d variable(s)"
                 % (self.data.get("created", "?"), n, len(self.data["variables"])))
        if self.data.get("mock"):
            self.log("  this snapshot was taken from the MOCK backend - "
                     "nothing in it is a result")

    def setup(self):
        pass

    def initialize(self):
        pass

    def iterate(self, n):
        pass

    def interrupt(self):
        pass

    def close(self):
        pass

    def surfaces(self):
        return list(self.data["surfaces"])

    def variables(self):
        return list(self.data["variables"])

    def bbox(self):
        lo = [1e30] * 3
        hi = [-1e30] * 3
        for e in self.data["surfaces"].values():
            for q in e["vertices"]:
                for k in range(3):
                    if q[k] is None:
                        continue
                    lo[k] = min(lo[k], q[k])
                    hi[k] = max(hi[k], q[k])
        if lo[0] > hi[0]:
            raise DriverError("the snapshot has no vertices to measure")
        return [lo, hi]

    def field(self, surface, variable):
        e = self.data["surfaces"].get(surface)
        if e is None:
            raise DriverError("the snapshot has no surface %r; it has %s"
                              % (surface, ", ".join(self.surfaces())))
        vals = e["values"].get(variable)
        if vals is None:
            raise DriverError(
                "the snapshot does not carry %r on %s - it was saved with %s. "
                "Reopen the case to sample a variable it does not have."
                % (variable, surface, ", ".join(self.data["variables"])))
        return {"surface": surface, "variable": variable,
                "vertices": e["vertices"], "faces": e["faces"],
                "values": vals, "mock": bool(self.data.get("mock")),
                "snapshot": True}

    def make_plane(self, name, axis, value):
        raise DriverError(
            "a snapshot cannot be cut: it holds the surfaces that were saved "
            "with it. Reopen the case to make a new plane.")

    def drop_plane(self, name):
        return False

    def load_case(self, path):
        raise DriverError("a snapshot is not a case")

    def report(self, kind, surfaces, variable):
        """The same integrals, computed here from the saved facets."""
        tot_a, acc, lo, hi = 0.0, 0.0, float("inf"), float("-inf")
        for name in surfaces:
            f = self.field(name, variable)
            verts, faces, vals = f["vertices"], f["faces"], f["values"]
            for face in faces:
                pts = [verts[i] for i in face]
                vs = [vals[i] for i in face]
                if any(v is None for v in vs):
                    continue
                a = _quad_area(pts)
                m = sum(vs) / len(vs)
                tot_a += a
                acc += a * m
                lo = min(lo, m)
                hi = max(hi, m)
        if kind == "area":
            return tot_a
        if kind == "facet-min":
            return lo if lo < float("inf") else 0.0
        if kind == "facet-max":
            return hi if hi > float("-inf") else 0.0
        if kind == "mass-flow-rate":
            raise DriverError(
                "a snapshot cannot give a mass flow: it holds one variable at "
                "a time on a surface, not the velocity vector and the normal "
                "together. Reopen the case for that.")
        return acc / tot_a if tot_a > 0 else 0.0


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
    w("S.%s = %r" % (P_("solver_time"), solver_time(s)))
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
    w("for z in S.setup.cell_zone_conditions.fluid:")
    w("    S.%s[z].general.material = %r"
      % (P_("zone_material").rsplit(".general.material", 1)[0]
         .replace("setup.cell_zone_conditions.fluid",
                  "setup.cell_zone_conditions.fluid"), m["name"]))
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
        w("inlet.momentum.velocity_specification_method = %r"
          % api_values("inlet", "spec")["components"])
        w("#  the mesh keeps the inlet plane perpendicular to the flow, and the")
        w("#  streamwise direction maps exactly to +X, rolled or not")
        w("inlet.momentum.velocity_components[0] = %g" % float(i["velocity"]))
        w("inlet.momentum.velocity_components[1] = 0.0")
        w("inlet.momentum.velocity_components[2] = 0.0")
    else:
        w("inlet.momentum.velocity_specification_method = %r"
          % api_values("inlet", "spec")["magnitude-normal"])
        w("S.%s = %g" % (P_("inlet_magnitude"), float(i["velocity"])))
    if not lam:
        w("#  the specification method activates the inputs below it, so it")
        w("#  has to be set first or Fluent refuses them as 'not active'")
        w("inlet.turbulence.turbulence_specification = %r" % i["turb_spec"])
        w("#  intensity is a FRACTION here: 0.05 is the 5 % the panel shows")
        w("inlet.turbulence.turbulent_intensity = %g   # %g %%"
          % (float(i["intensity"]) / 100.0, float(i["intensity"])))
        if i["turb_spec"] == "Intensity and Viscosity Ratio":
            w("inlet.turbulence.turbulent_viscosity_ratio = %g" % float(i["visc_ratio"]))
        else:
            w("inlet.turbulence.hydraulic_diameter = %g" % float(i["hydraulic_diameter"]))
    w("")
    w("outlet = S.setup.boundary_conditions.pressure_outlet['outlet']")
    w("outlet.momentum.gauge_pressure = %g" % float(o["gauge_pressure"]))
    w("outlet.momentum.prevent_reverse_flow = %r" % bool(o["prevent_reverse_flow"]))
    if not lam and o["prevent_reverse_flow"]:
        w("#  reverse flow is prevented, so Fluent deactivates the backflow")
        w("#  turbulence inputs and refuses every write to them")
    if not lam and not o["prevent_reverse_flow"]:
        w("outlet.turbulence.turbulence_specification = %r" % i["turb_spec"])
        w("outlet.turbulence.backflow_turbulent_intensity = %g   # %g %%"
          % (float(o["backflow_intensity"]) / 100.0, float(o["backflow_intensity"])))
        if i["turb_spec"] == "Intensity and Viscosity Ratio":
            w("outlet.turbulence.backflow_turbulent_viscosity_ratio = %g"
              % float(o["backflow_visc_ratio"]))
        else:
            w("outlet.turbulence.backflow_hydraulic_diameter = %g"
              % float(i["hydraulic_diameter"]))
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

        #  a path that resolves can still be handed a value it will not take,
        #  so the two audits run together and both have to pass
        nc, unchecked, crep = audit_choices()
        print("\nauditing the enum values those settings are given")
        print("  %d value(s) checked, %d not published by any release"
              % (nc, len(unchecked)))
        if "-v" in sys.argv and unchecked:
            print("  no release publishes an allowed list for these, so they "
                  "stand unverified:")
            for label in sorted({l for l, _ in unchecked}):
                vals = sorted(v for l, v in unchecked if l == label)
                print("     %-24s %s" % (label, ", ".join(vals)))
        for v in ("242", "251", "252", "261", "271"):
            problems = crep.get(v, [])
            if not problems:
                continue
            bad += len(problems)
            print("  v%s   %d PROBLEM(S)" % (v, len(problems)))
            for owner, msg in problems:
                print("        %-22s %s" % (owner, msg))
        if not any(crep.get(v) for v in crep):
            print("  every release that publishes its allowed values agrees")
        sys.exit(1 if bad else 0)
    if "--schema" in sys.argv:
        print(json.dumps({"settings": SETTINGS, "defaults": default_settings(),
                          "variables": VARIABLES}, ensure_ascii=False, indent=1))
        sys.exit(0)
    print(__doc__)
