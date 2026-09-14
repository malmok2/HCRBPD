# HCRBPD — helical coil / rod bundle pressure drop

Parametric CAD + block-structured hex mesh explorer for cross-flow pressure-drop
studies. One self-contained HTML file: it draws the 3-D geometry, previews the
mesh with quality colouring, estimates the pre-CFD velocity and pressure field
on that mesh, and exports solver-ready files for ANSYS Fluent / CFX, OpenFOAM,
ParaView and STL — all in the browser, no mesher.

| File | What it is |
|---|---|
| `mesh_explorer.html` | the tool — open it directly in a browser, nothing to install |
| `mesh_explorer.py` | the standalone twin — same mesh, separate code, for scripting and cross-checking |

## Four geometries, one tool

Pick the geometry from the panel at the top left. They differ in two
independent traits:

| | in-line | staggered |
|---|---|---|
| **rod** — straight circular rods, channel height `H` set directly | `rod-inline` | `rod-staggered` |
| **helical** — tubes inclined by the helix angle, meshed flat and rolled back onto the coil | `helical-inline` | `helical-staggered` |

Switching geometry keeps the parameters, so the same case can be compared
in-line against staggered, or flat against rolled, without retyping it.
`Reset` restores that geometry's own defaults.

For the helical family the section footprint is an ellipse `2·R/cos α` long by
`2·R` wide, and the channel height is the arc length `R_mid·θ` rather than a
free parameter. Setting `α = 0` and turning the roll off collapses every
helical term back to the rod case — which is why a single code path serves
both.

```
python3 mesh_explorer.py --geometry helical-staggered --helix 20 --sector 15
python3 mesh_explorer.py --geometry rod-inline --json --no-write   # numbers only
```

## The standard

The interface, controls, checks and export formats are a house standard shared
by every geometry; **geometry knowledge lives in exactly five functions**
(`derived`, `lattice`, `cellPolygon`, `sideDiv`, `warp3`), each labelled in
both files. That standard — and the hard-won failure modes behind it — is
written down as a Claude Code skill at
[`.claude/skills/parametric-mesh-explorer/`](.claude/skills/parametric-mesh-explorer/SKILL.md).
Start there before adding a geometry or touching anything the user sees.

This used to be four copied sibling files, one per geometry, which meant every
fix to shared machinery had to be applied four times. A new geometry is now an
entry in the `GEOM` registry, not a new copy.

The helix lean is ramped back to zero across the inlet and outlet boxes, so
those two planes stay exactly perpendicular to the flow: the inlet patch normal
is `(-1,0,0)` at every face, rolled or not, and Fluent's default
"Magnitude, Normal to Boundary" velocity inlet is simply correct. Leaning the
whole patch instead would tilt the inlet by the full helix angle. Give the
inlet and outlet boxes a non-zero length — with `L_in = 0` there is nothing to
ramp through and the tool warns that that end is left oblique.

## Verification

Every export is gated: the mesh must be watertight (0 cracks), every shared
node bit-identical across blocks, every cell volume positive when rebuilt from
its faces, and the total volume must match the analytic volume of the faceted
domain. The Fluent face-orientation convention was settled empirically against
ANSYS Fluent 2026 R1 and is documented in the skill's `mesh-core.md`.

For a rolled coil the analytic target is closed form too: by Pappus the volume
is the section's first moment about the coil axis times the swept angle, and
the chorded cells come out short by exactly `sin(Δθ)` per layer — so that case
is gated at 1e-6 like any other, rather than being skipped as "indicative" the
way the old helical explorers had to.

The Python twin builds the same mesh independently. As of the merge, browser
and twin agree on the **order-sensitive checksums of both the node and the quad
tables**, on every node/face/patch count, and on total volume to 1e-15, across
all four geometries at default and coarse settings. `mesh_explorer.py --json`
exists for exactly that comparison; repeat it after any change to the mesh core.

## Generated files

`*.msh`, `*.vtu`, `*.stl`, `*_foam/` are ignored by git. Each is reproducible
in seconds from the explorer or its Python twin.
