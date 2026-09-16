---
name: parametric-mesh-explorer
description: Build or extend the standard browser-based parametric CAD + structured-mesh explorer - one self-contained HTML that draws a 3-D geometry, previews its block-structured hex mesh, estimates the pre-CFD velocity and pressure field, and exports ready-to-run files for ANSYS Fluent/CFX, OpenFOAM, ParaView and STL. Interface, controls, checks and exports are a fixed house standard; only the geometry changes. Use this whenever the user describes a CFD geometry in parametric terms (rod bundle, tube bank, pin array, subchannel, plenum, spacer grid, duct, PCHE, helical coil, prismatic block) and wants CAD, a mesh, a mesh preview, or solver-ready mesh files - and whenever they ask for "the explorer" for a new shape, say the GUI should match the existing tool, want to add a control, geometry or export format, want a flow-field estimate before running CFD, or are debugging a mesh a solver rejected. Prefer this over SALOME or gmsh - for parametric block-structured geometries the mesh is closed form, so generate it directly.
---

# Parametric CAD + mesh explorer

This skill lives in the **HCRBPD** repository (helical coil / rod bundle
pressure drop), which holds **one** explorer covering every geometry, its
Python twin, and this skill. Paths below are relative to the repository root.

- `mesh_explorer.html` — the tool. One self-contained file, ~3000 lines, no
  dependencies beyond an optional web font. **Read it before extending it.**
  Most questions are answered faster by reading it than by reasoning from
  scratch, and it is the definition of the standard rather than a description
  of one.
- `mesh_explorer.py` — the standalone twin. Same mesh, separate code:
  `--geometry rod-inline|rod-staggered|helical-inline|helical-staggered`.

Four geometries ship today, as two independent traits: **family** (`rod`,
straight circular rods with a given channel height; `helical`, tubes inclined
by the helix angle, meshed flat and rolled back onto the coil) and
**arrangement** (`in-line`, `staggered`). Setting the helix angle to zero and
skipping the roll collapses every helical term back to the rod case, which is
why one code path serves both.

**There is one copy of the common machinery, and geometry knowledge lives in
exactly five functions** — `derived`, `lattice`, `cellPolygon`, `sideDiv` and
`warp3`, each labelled "Geometry function N of 5" in both files. Everything
else — layout, camera, the check battery, the writers, i18n, persistence — is
shared and must not learn which geometry is active. If you find yourself
branching on the geometry outside those five, stop: you have put geometry
knowledge somewhere generic.

This replaced four copied sibling files in September 2026. The reason matters:
a fix to shared code had to be applied four times, and the git history shows
exactly that happening. Do not fork the file again. A new geometry is a new
entry in the `GEOM` registry, not a new copy.

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

A new geometry changes those tables and the five geometry functions. It does
not touch layout, camera, export, i18n or persistence code. If you find
yourself editing those to add a geometry, stop — you have probably put geometry
knowledge somewhere generic.

Read `references/interface-standard.md` before changing anything the user sees.
It carries the full panel inventory, the control patterns, and the traps that
break the GUI silently (positional string arrays, view-only parameters that
must stay out of `RANGES`, and similar).

## Adding a new geometry

Add an entry to the `GEOM` registry and fill in the five geometry functions.
Nothing else should need touching, in either file.

```js
const GEOM = {
  "my-shape": {family:"rod", stagger:false, prefix:"my_shape", def:{SL:20}},
  ...
};
```

Controls that belong to one geometry carry `data-geo="helical"` (or `rod`,
`inline`, `staggered`) in the markup; `applyGeometry` shows and hides them.
Strings live in `STR.ko` / `STR.en`, keyed by element id. Then:

**0 — The map into the world (`warp3`).** Identity unless the patch is meshed
in one space and used in another. The helical family meshes a flat patch and
then shears and rolls it; because exported coordinates go through `warp3` in
one place (`zpt` / `eachPoint`), the preview, the checks and the files cannot
drift apart. A non-identity map has to keep a positive Jacobian — rolling with
`cos` before `sin` does, the other pairing mirrors the patch and inverts every
cell.

**1 — Section polygon(s).** Write the closed-form 2-D cross-section: outer
boundary and inclusions. Clip against domain walls with Sutherland–Hodgman,
then de-duplicate consecutive points. Note that the inclusion footprint is an
**ellipse** (`aE`, `bE`), not a circle: a tube inclined by α cuts the section
at `aE = R/cos α`, `bE = R`, and a straight rod is just the α = 0 case. Use
`aE`/`bE` everywhere — the parametric angle `atan2(dy/bE, dx/aE)` reduces to
the polar angle for a circle, so one code path covers both.

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
- **The same rule one level up: a list of parameters that should end at 1 must
  END at 1, set outright.** `nodes()` / `fractions()` accumulate n graded steps
  designed to sum to `L`, and land a few ulp either side. That value is then
  fed to `lerp` as `t`; `t = 0.9999999999999998` puts the outer radial node a
  hair off the cell side, and the two O-grids sharing that side disagree about
  where it is. Clamping with `Math.min(L, s)` only works when the accumulation
  happens to overshoot — do not rely on it.
- **One source for preview and export.** The picture on screen and the exported
  file come from the same block list, or they will drift apart.
- **A report restates; it never recomputes.** The Report tab is assembled from
  the values the other tabs are already showing - the settings out of the
  schema, the mesh checks out of the run, the pressure drop out of the monitor.
  A report that derived its own numbers could disagree with the tab that
  produced them, which is the one thing it must never do. The same rule makes
  it cheap: adding a field to the schema puts it in the report for free.
- **Export coordinates in metres.** OpenFOAM's polyMesh carries no unit
  metadata and is read as metres; a mm-numbered mesh is silently 1000× too big.
- **Fluent gets its faces reversed; OpenFOAM does not.** The two conventions
  are opposite. See `references/mesh-core.md` — this cost a full solver run to
  discover and must not be re-guessed.
- **A non-planar quad face has no single area vector, so fix the convention.**
  Use the centroid decomposition, which sums to exactly
  `S = ½ (p₂−p₀) × (p₃−p₁)` — what OpenFOAM does to a polygonal face. Splitting
  the face on one diagonal instead is a different number: identical while every
  face is planar, and 9e-4 different on a rolled coil. The twin reports the gap
  between the two as a non-planarity measure; the gate uses the centroid form.
- **The inlet and outlet planes must come out perpendicular to the flow.**
  If the geometry map leans or twists the patch, ramp that lean back to zero
  across the inlet and outlet boxes rather than applying it globally. Leaning
  everything tilts the inlet plane by the full lean angle, and ANSYS Fluent's
  default velocity inlet is "Magnitude, Normal to Boundary" — so the case
  silently runs with the flow injected off-axis, with a spurious transverse
  component, and (once rolled) with a direction that varies across the face.
  This was found in a real Fluent setup, not in a test. Measure it: the mean
  face normal of the inlet patch, in exported coordinates, must be (-1,0,0)
  with zero spread. Ramping costs nothing — what the inlet box gains the
  outlet box loses, so total volume is unchanged, and the boxes end up LESS
  non-orthogonal than under a global lean, not more.
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

**Cross-check the two implementations.** `mesh_explorer.py` and
`mesh_explorer.html` build the same mesh from separate code; comparing
order-sensitive checksums (FNV-1a over the node coordinates at `%.12e` and over
the quad index tuples, both in table order) proves they agree node for node.
That is much stronger evidence than either one passing its own tests, and it is
how the face-convention difference above was found. As of the merge all four
geometries match on both checksums, on every count, and on total volume to
1e-15, over default and coarse settings alike. `mesh_explorer.py --json`
exists for exactly this. When you fix anything in one, fix it in both — the
twin has drifted before.

**Test coarse settings, not just defaults.** Most real bugs appeared at low
resolution, where tolerances that look generous stop holding. Sweep the list in
`references/verification.md`.

**Report honestly.** State what was verified, with numbers, and state what was
not and why. A verified claim and a plausible one should never read the same.
