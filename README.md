# HCRBPD — helical coil / rod bundle pressure drop

Parametric CAD + block-structured hex mesh explorer for cross-flow pressure-drop
studies. One self-contained HTML file: it draws the 3-D geometry, previews the
mesh with quality colouring, estimates the pre-CFD velocity and pressure field
on that mesh, and exports solver-ready files for ANSYS Fluent / CFX, OpenFOAM,
ParaView and STL — all in the browser, no mesher.

| File | What it is |
|---|---|
| `mesh_explorer.html` | the front end — five tabs: Geometry, Settings, Run, Results, Report |
| `mesh_explorer.py` | the mesh generator, standalone — same mesh as the browser, separate code |
| `fluent_case.py` | the Fluent layer — settings schema, PyFluent driver, offline mock, journal writer |
| `app.py` | the local server that ties them together and launches Fluent |

## Two ways to use it

**As a mesh tool** — open `mesh_explorer.html` in a browser. Nothing to install,
nothing to run. The Geometry tab works exactly as before and every export format
is there. The other three tabs explain why they need a server.

**As a CFD app** — get the repository, install PyFluent, run the server.

Linux / macOS:

```
git clone -b claude/funny-mendel-qzcep2 https://github.com/malmok2/HCRBPD.git
cd HCRBPD
pip install ansys-fluent-core
python3 app.py
```

Windows (PowerShell) — there is no `python3` on Windows, use `python` or `py`:

```powershell
git clone -b claude/funny-mendel-qzcep2 https://github.com/malmok2/HCRBPD.git
cd HCRBPD
pip install ansys-fluent-core
python app.py
```

Note the `-b`: the work is on that branch, not on `main`. The files the app
writes are LF on every platform, including Windows, so a polyMesh written on a
Windows workstation drops straight onto a Linux cluster.

**Two Fluent releases on one machine is the normal case** — a lab licence on
one, a student install on another — so the release is a setting, not a guess.
The Settings tab lists every release the app knows and marks which are actually
here, found the way PyFluent finds them (`AWP_ROOT<nnn>`, then the executable),
with the variable and its value on the tooltip. Picking one launches Fluent at
that release *and* writes the generated script against the same release's
settings-API spelling, which is not cosmetic: 2024 R2 keeps the discretisation
schemes shallower and 2027 R1 calls `models.viscous` `models.turbulence`.
Choosing a release that is not installed is said in the run log before Fluent
is launched, not a minute into it, and the log says so again if a different
release answers than the one asked for.

`app.py` prints what it found before serving — `backend: auto -> Fluent 26.1.0`
if it can drive the real thing, or the reason it cannot. PyFluent locates Fluent
through the `AWP_ROOT<version>` environment variable that the Ansys installer
sets; if that is missing the app says so up front and falls back to the mock
rather than failing part way into a run.

It also prints the revision it is running — `code fb7741f (2026-09-15) ...` —
and that same string heads every run log and sits in the page header. Python
imports its modules once, so a browser refresh picks up a change to the HTML
but **not** to `fluent_case.py` or `mesh_explorer.py`: after a `git pull`,
restart the server. The stamp is there so a stale process is a glance rather
than a guess, and turns orange when the checkout has uncommitted edits.

It serves the same page on `127.0.0.1` (loopback only) and opens it. The
Settings tab configures the Fluent case, the Run tab meshes, launches Fluent and
plots residuals live — and beside them the **pressure drop itself**, sampled
every N iterations, because residuals settling is not the same as the answer
settling and a pressure-drop study cares about the second. The plot calls out
the current Δp and how much it moved over the final stretch; the interval is a
setting, and 0 turns it off. It is sampled by the app rather than kept by
Fluent, so a case file has no history — reopening one reports the converged
value and says that is all it can. A snapshot carries both histories, so
reopening one puts the plots back as they were.

**The Report tab** is one page that says what was run: the headline numbers,
the geometry, the mesh and its checks, the solver set-up straight out of the
settings schema so it cannot drift from the panel, the convergence, and the
fields. Nothing in it recomputes anything — a report that derived its own
numbers could disagree with the tab that produced them. It puts the pressure drop
against the Jakob correlation — and does it honestly: inlet-minus-outlet spans
the whole domain, boxes included, which is not what the correlation predicts,
so the report also puts two planes on the bundle faces and measures the drop
across the bundle alone. That is the number the ratio uses, with the per-row
figure beside it. It says which branch of the correlation was spot-checked and
which was not, warns when a short bundle is being compared against a
fully-developed correlation, and flags a MOCK run in red at the top. `Save as HTML` writes a standalone file with the
figures embedded as data URIs and the stylesheet lifted out of the page's own,
so it opens anywhere with nothing else beside it; `Print · PDF` uses a print
stylesheet that drops the app around it.

The Results tab lists what can be reopened too, so an old result is one click
away from where you would look at it, not a walk back to the Run tab. The
Results tab itself pulls fields back and draws contours, with area-weighted
averages, mass flows and the bundle pressure drop.

**Contours on a plane, not just on the boundary.** The Results tab cuts a
constant-x, -y or -z plane through the domain — the slider moves it in real
coordinates over the measured extent of the case — and Fluent then treats that
plane as a surface like any other, so the contour, the surface values and the
reports all work on it with no special case. Levels are continuous by default;
set a band count and the fill quantises to those bands, optionally with the
contour lines drawn on the band edges, and the colour bar labels them.

**The colour map is a choice, and an honest one.** Ten maps, their control
points sampled out of matplotlib rather than written from memory and reduced to
the fewest stops that reproduce each to under ~3/255 per channel. They are
grouped by what they are: *sequential* (viridis, magma, inferno, plasma,
cividis) rises in lightness one way, so the order survives greyscale and colour
blindness; *diverging* (cool-warm, RdBu, Spectral) is lightest in the middle,
for a quantity with a meaningful zero such as pressure about the outlet; and
*rainbow* (turbo, jet) is neither. Every claim there was measured — CIE L\*
along each map, and again through a deuteranopia simulation — and turbo and
jet are the only two that fail, so they are the only two marked. They are
offered because they are asked for. Levels run to 200, and `Auto` picks a
diverging map for pressure and a sequential one for everything else.

**The y+ target is checked against the turbulence model.** They live in
different tabs — the target is a meshing parameter, the model a solver setting
— so nothing compared them, and the shipped default paired k-ω SST with a
target of 30. That puts the first cell in the buffer layer, where neither the
viscous sublayer nor the log law is resolved: the wall shear comes out wrong
and the pressure drop with it. The default is 1 now, and the pair is checked in
the Geometry tab and again in the report: k-ω and Spalart–Allmaras want y+ of
order 1, k-ε with standard wall functions wants 30–300, and k-ε with enhanced
or scalable treatment takes either end. 5 < y+ < 30 is called out wherever the
model has no blending for it.

**Wall y+** is among the variables, marked as a wall quantity: ask for it on a
plane or an inlet and the tab says so rather than drawing an empty surface.

**Fixed views and a figure you can publish.** Seven view buttons — ISO and the
six axis views — where the axis a button names is the one pointing at the
viewer; each was solved from the projection rather than guessed, and the triad
confirms it live. `Fit` resets zoom and pan without touching the angle, a free
drag leaves the preset behind, and auto-rotate spins on the fast draw. The residual plot and
the Δp plot save the same way, from the same code, with the case name and the
final number in a clean strip above them rather than laid over an axis label —
a plot for a report, not a screenshot of one. `Save as
PNG` renders at 1× to 6× of the on-screen size — real resolution, because the
backing store and the transform go up together, so the line widths and type
scale with it — on white or transparent, with the case name, variable, range,
surfaces and the MOCK flag optionally stamped in the corner. The filename
carries the case, the surfaces and the variable, so a figure stays identifiable
long after the session that made it.

**The CAD silhouette** draws the domain and the rods as a wireframe behind the
field — a cut through a bundle is hard to place without the bundle around it.
The polylines come from the same `Case` and through the same `XP()` as the mesh,
so they register exactly and curve with a rolled coil; a half rod on a side wall
is cut at the wall rather than drawn sticking out of the box. **A triad** in the
corner shows where x, y and z currently point, built from the projection itself.

**Drawn at the resolution it was measured at.** Fluent returns a value at
every *node*, not one per facet, and colouring a facet with the mean of its
corners throws that away — the contour then has exactly as many tiles as the
mesh has cells. Each face is now subdivided and its corner values interpolated
across it, which is what Fluent's own contours do; the note under the controls
says how far it subdivided and how many cells it drew, and the factor adapts so
a coarse patch is smoothed hard and an already-fine plane is left nearly alone.
It is interpolation, not new data. Rotating and zooming drop to the flat draw
and the quality comes back when the pointer stops.

**Where the files went, and getting them back.** The Run tab opens with the
output folder, its path ready to copy, a button that opens it in the file
manager, and everything written so far with its kind, size and time. A case file
carries an Open button: it reopens that case and its data in Fluent and lands in
the Results tab with nothing re-solved, because a solution that took an hour
should not have to be produced twice to be looked at twice.

A **field snapshot** goes further. `Save the chosen surfaces` writes the sampled
field itself — vertices, faces and one value per node per variable, for the
patches and planes you were looking at — as a `.fields.json`. That reopens with
no Fluent at all, on any machine, with no licence: every variable re-colours it,
the surface integrals are recomputed from the saved facets, and the contour
draws the same way. It is a record of what was sampled, so it says so, and it
refuses what it cannot answer — a plane it was not saved with, a variable it
does not carry, a mass flow that needs the velocity vector and the face normal
together — naming what to do instead rather than guessing.
`python3 app.py --backend mock` runs the whole thing with an invented field and
no Fluent, which is how the plumbing is tested; anything it produces is labelled
MOCK on screen and must not be quoted as a result.

### The Fluent settings are declared once

`fluent_case.py` holds one schema: for each setting, its type, default, choices
and the settings-API path it drives. The Settings panel is generated from that
schema over the API and the driver writes those same paths, so a setting cannot
appear in the panel without a path, or be applied to a path the panel never
showed. `python3 fluent_case.py --audit` walks every path against the settings
trees PyFluent ships for Fluent 2024 R2 through 2027 R1, and then the enum
**values** those settings are given — a path that resolves can still be handed
a string the setting will not take, which is how `least-squares-cell-based`
survived for a scheme Fluent calls `least-square-cell-based`. Roughly half the
values the shipped trees publish; the rest the live driver asks the running
Fluent about, after the mesh is read and before it writes any of them, so one
run names every bad string instead of dying on the first. `--audit -v` lists
what no release publishes, so what stands unverified is visible rather than
implied. The API does move
between releases — 2024 R2 keeps the discretisation schemes and the surface
integrals shallower, 2027 R1 renames `models.viscous` to `models.turbulence` —
so each path carries its alternates and the audit requires one to resolve in
every release. `--schema` prints the whole contract as JSON.

A setting Fluent has deactivated is not written at all. Preventing reverse
flow at the outlet means there is no backflow, so Fluent greys out the whole
backflow turbulence group; the panel hides those fields to match, and neither
the driver nor the generated script touches them. And because the browser keeps
the last case in local storage, a saved choice can outlive a correction to the
schema — those are dropped back to the default on load and on arrival at the
server, both of which say which ones, rather than being sent to Fluent as a
string it rejects.

The Settings tab also emits the equivalent standalone PyFluent script, written
for the release you pick, so a case can be reproduced and archived without the
app.

One unit is worth knowing: the settings API takes turbulent intensity as a
**fraction** while the panel asks for a percentage, and this is documented in
neither the API nor the shipped examples. It was settled by measurement —
0.05 written through the API shows as 5 in the Turbulent Intensity box of
Fluent 2025 R1 — so the app divides the panel value by 100. Every run re-reads
it back out of Fluent and compares, and says so if a release ever changes the
convention, rather than quietly running the case at a hundredth of the
intended turbulence.

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
