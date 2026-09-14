# HCRBPD — helical coil / rod bundle pressure drop

Parametric CAD + block-structured hex mesh explorers for cross-flow pressure-drop
studies. Each geometry is one self-contained HTML file: it draws the 3-D
geometry, previews the mesh with quality colouring, estimates the pre-CFD
velocity and pressure field on that mesh, and exports solver-ready files for
ANSYS Fluent / CFX, OpenFOAM, ParaView and STL — all in the browser, no mesher.

| Folder | Geometry | Explorer | Python twin |
|---|---|---|---|
| `rod_bundle_lattice/` | in-line rod lattice (**reference implementation**) | `rod_lattice_explorer.html` | `rod_lattice_mesh.py` |
| `rod_bundle_staggered/` | staggered rod bundle | `rod_bundle_explorer.html` | `rod_bundle_mesh.py` |
| `helical_coil_bundle/` | helical coil, in-line | `helical_coil_explorer.html` | `helical_coil_mesh.py` |
| `helical_coil_bundle_staggered/` | helical coil, staggered | `helical_coil_staggered_explorer.html` | `helical_coil_staggered_mesh.py` |

Open any explorer directly in a browser. Nothing to install.

## The standard

The interface, controls, checks and export formats are a house standard shared
by every geometry; only the geometry code differs. That standard — and the
hard-won failure modes behind it — is written down as a Claude Code skill at
[`.claude/skills/parametric-mesh-explorer/`](.claude/skills/parametric-mesh-explorer/SKILL.md).
Start there before adding a geometry or touching anything the user sees.

The explorers are siblings by copy, not by shared library. A fix to the common
machinery must be applied to every folder; `grep` across the repository to be
sure it was.

## Verification

Every export is gated: the mesh must be watertight (0 cracks), every shared node
bit-identical across blocks, every cell volume positive when rebuilt from its
faces, and the total volume must match the analytic volume of the faceted
domain. The Fluent face-orientation convention was settled empirically against
ANSYS Fluent 2026 R1 and is documented in the skill's `mesh-core.md`.

The Python twin builds the same mesh independently; comparing node and face
counts, patch sizes and total volume between the two is the strongest check
available and is worth repeating after any change to the mesh core.

## Generated files

`*.msh`, `*.vtu`, `*.stl`, `*_foam/` are ignored by git. Each is reproducible in
seconds from the explorer or its Python twin.
