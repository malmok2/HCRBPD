---
name: parametric-mesh-explorer
description: Build or extend the standard browser-based parametric CAD + structured-mesh explorer - one self-contained HTML that draws a 3-D geometry, previews its block-structured hex mesh, estimates the pre-CFD velocity and pressure field, and exports ready-to-run files for ANSYS Fluent/CFX, OpenFOAM, ParaView and STL. Interface, controls, checks and exports are a fixed house standard; only the geometry changes. Use this whenever the user describes a CFD geometry in parametric terms (rod bundle, tube bank, pin array, subchannel, plenum, spacer grid, duct, PCHE, helical coil, prismatic block) and wants CAD, a mesh, a mesh preview, or solver-ready mesh files - and whenever they ask for "the explorer" for a new shape, say the GUI should match the existing tool, want to add a control, geometry or export format, want a flow-field estimate before running CFD, or are debugging a mesh a solver rejected. Prefer this over SALOME or gmsh - for parametric block-structured geometries the mesh is closed form, so generate it directly.
---

# Parametric CAD + mesh explorer

This skill lives in the **HCRBPD** repository (helical coil / rod bundle
pressure drop), which holds one explorer per geometry, each in its own folder
at the repository root, and this single shared skill. Paths below are relative
to that root.

The reference implementation is `rod_bundle_lattice/rod_lattice_explorer.html`
— one file, ~2550 lines, no dependencies beyond an optional web font — plus its
standalone Python twin `rod_bundle_lattice/rod_lattice_mesh.py`. **Read the
reference file before extending it or porting it.** Most questions are answered
faster by reading it than by reasoning from scratch, and it is the definition of
the standard rather than a description of one.

The other geometries sit beside it — `rod_bundle_staggered/`,
`helical_coil_bundle/`, `helical_coil_bundle_staggered/` — and each carries its
own `*_explorer.html` and `*_mesh.py`. They are siblings by copy, not by shared
library: a fix to the common machinery has to be applied to every folder, and
the way to be sure it was is to `grep` for it across the repository. Check
whether a geometry already exists before starting a new one, and when adding
one, give it a folder of its own at the root.

## The core claim: skip the mesher

For a geometry that decomposes into **mapped (transfinite) blocks**, every node
position is a closed-form expression. There is nothing for a meshing algorithm
to discover. Generating directly instead of driving SALOME or gmsh means:

- no external install, no version-dependent API, no boolean-kernel fragility
- the whole thing runs in a browser in under a second for ~10⁵ cells
- **you can verify it** — watertightness, volumes, face ordering, connectivity
  and cell volumes are all checkable locally, which is impossible when the
  mesher runs elsewhere

The boundary of the claim: it holds while the geometry stays parametric and
block-decomposable. Fillets, imported CAD, spacer grids with mixing vanes, or
anything needing unstructured tets — that is where a real kernel earns its
place. Say so plainly rather than forcing a block decomposition that does not
exist.

## The standard: what is fixed, what varies

The user's intent is that **the interface and the feature set stay identical
across geometries**. They will specify the shape and roughly how it should be
meshed; they should not have to re-specify the tool.

| Fixed — port unchanged | Varies per geometry |
|---|---|
| Window layout, design system, control patterns | Section polygon(s) and inclusions |
| The six display modes and the section-plane slider | Block decomposition |
| Flow-conditions group (y+ and the field estimate) | Topological key scheme |
| Options, view presets, camera behaviour | Boundary patch classification |
| Export panel, unit handling, share link, reset | Parameter list, labels, hints |
| The check battery and the refusal to export bad meshes | Warnings, stats cells, legend rows |
| Writers and patch ordering | Geometry toggles (e.g. half inclusions) |

The reason this split works is that **the GUI is data-driven**. Layout code
reads three tables:

- `DEFAULTS` / `RANGES` / `CHECKS` — which parameters exist
- `STR.ko` / `STR.en` — every visible string
- `HINT` — the derived value shown beside each number

A new geometry changes those tables and the four geometry functions. It does
not touch layout, camera, export, i18n or persistence code. If you find
yourself editing those to add a geometry, stop — you have probably put geometry
knowledge somewhere generic.

Read `references/interface-standard.md` before changing anything the user sees.
It carries the full panel inventory, the control patterns, and the traps that
break the GUI silently (positional string arrays, view-only parameters that
must stay out of `RANGES`, and similar).

## Adding a new geometry

Only four things are geometry-specific. Everything else is reusable as-is.

**1 — Section polygon(s).** Write the closed-form 2-D cross-section: outer
boundary and inclusions. Clip against domain walls with Sutherland–Hodgman,
then de-duplicate consecutive points.

**2 — Block decomposition.** Cut the section into quadrilateral blocks, each a
`grid[i][j]` of node positions. For an inclusion, enclose it in a cell polygon
(the Voronoi cell of the lattice works well) and run radial lines from the
inclusion to each polygon corner: one mapped block per polygon side, i.e. an
O-grid. For plain regions, one block per interface segment so every block stays
a clean 4-sided quad.

**3 — Topological keys.** Give every node in every block a name. This is what
welds blocks into one conformal mesh, and getting it wrong is the single most
likely source of a broken mesh. Read `references/mesh-core.md` before writing
any key function.

**4 — Boundary classification.** Map each 2-D boundary edge to a patch name by
geometry. Classify using edge **endpoints**, never midpoints — a chord midpoint
sits at `R·cos(Δθ/2)`, which drifts off a curved wall as the azimuthal count
drops, and that silently mis-classified 2600 faces at 12 divisions per rod.

Then the shared machinery applies unchanged: extrude N layers to hexahedra,
enumerate faces structurally, run the check battery, and write the files.

## Non-negotiables

These are invariants that took real debugging to find. Preserve them.

- **Shared nodes must be bit-identical, not merely close.** Blocks share points
  through topological keys, never through rounded-coordinate hashing.
- **`lerp(a, b, t)` returns `a` exactly at `t=0` and `b` exactly at `t=1`.**
  In floating point `a + (b-a)*1.0 ≠ b`, and that inequality becomes a crack.
- **One source for preview and export.** The picture on screen and the exported
  file come from the same block list, or they will drift apart.
- **Export coordinates in metres.** OpenFOAM's polyMesh carries no unit
  metadata and is read as metres; a mm-numbered mesh is silently 1000× too big.
- **Fluent gets its faces reversed; OpenFOAM does not.** The two conventions
  are opposite. See `references/mesh-core.md` — this cost a full solver run to
  discover and must not be re-guessed.
- **Refuse to export an invalid configuration.** Degenerate geometry, cracks,
  non-positive cell volumes: fail loudly rather than writing a file that will
  waste a solver run.
- **Every claim in the report needs a check that could have failed.** See
  `references/verification.md`.

## Reference files

Read the one you need; do not load all five by default.

| File | Read it when |
|---|---|
| `references/interface-standard.md` | porting the tool to a new geometry, or touching any control, panel or string |
| `references/mesh-core.md` | writing block/key/patch code, or any writer (Fluent, polyMesh, VTU, STL, zip) |
| `references/flow-estimate.md` | touching the pre-CFD velocity/pressure estimate or its streamlines |
| `references/verification.md` | before claiming the mesh is correct; also lists the bugs actually hit |
| `references/design-system.md` | touching CSS, colour, type, colour maps or i18n |

## Working style that fits this tool

**Verify in the browser, not by inspection.** The preview pane renders local
files as static snapshots — **re-navigate after every edit**, or you will test
stale code and draw false conclusions. Drive the page with `javascript_exec`:
instrument `ctx` to count draw calls, read back canvas pixels, force viewports,
and click buttons programmatically with `window.save` stubbed so nothing
downloads.

**Prefer a measurement over a look.** Most visual claims can be turned into a
number: count saturated pixels to prove a contour painted, diff two renders to
prove a toggle did something, wrap `drawRods`/`drawSection3D` to prove the
painter order flips with the camera. These catch what a screenshot glance
misses, and they can be re-run after the next change.

**Cross-check two independent implementations when both exist.** The Python
twin and the browser build the same mesh; comparing order-sensitive checksums
of the node and quad tables proves they agree exactly. That is much stronger
evidence than either one passing its own tests. When you fix a writer bug, fix
it in both — the twin has drifted before.

**Test coarse settings, not just defaults.** Most real bugs appeared at low
resolution, where tolerances that look generous stop holding. Sweep the list in
`references/verification.md`.

**Report honestly.** State what was verified, with numbers, and state what was
not and why. A verified claim and a plausible one should never read the same.
