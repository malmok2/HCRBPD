# HCRBPD — helical coil / rod bundle pressure drop

Parametric CAD + block-structured hex mesh explorer for cross-flow pressure-drop
studies. One self-contained HTML file: it draws the 3-D geometry, previews the
mesh with quality colouring, estimates the pre-CFD velocity and pressure field
on that mesh, and exports solver-ready files for ANSYS Fluent / CFX, OpenFOAM,
ParaView and STL — all in the browser, no mesher.

| File | What it is |
|---|---|
| `mesh_explorer.html` | the front end — six tabs: Geometry, Settings, Run, Results, Report, Campaign |
| `mesh_explorer.py` | the mesh generator, standalone — same mesh as the browser, separate code |
| `fluent_case.py` | the Fluent layer — settings schema, PyFluent driver, offline mock, journal writer |
| `app.py` | the local server that ties them together and launches Fluent |
| `correlations.py` | the correlation library — source, validity range and status for each |
| `zukauskas_charts.py` | Žukauskas' four digitised charts, vendored with their provenance |
| `study.py` | the parametric campaign — case matrices, the runner, GCI, the stage reports |
| `studies/` | one folder per study: what is to be run, what came back, and the report |

## Two ways to use it

**As a mesh tool** — open `mesh_explorer.html` in a browser. Nothing to install,
nothing to run. The Geometry tab works exactly as before and every export format
is there. The other four tabs explain why they need a server.

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

## The campaign

The point of the tool is not one pressure drop. It is to put CFD against the
published correlations over a matrix of arrangements and conditions, and in
the end to fit a better one — first for straight rods, then for the helical
coil. That is four stages, and the **파라메트릭** tab is where they run.

| stage | what it settles | needs Fluent |
|---|---|---|
| 1 | which correlations exist, what each may be asked, and how far apart they already are | no |
| 2 | how fine the mesh has to be before Δp stops moving | yes |
| 3 | Δp over the arrangement/condition matrix, against the correlations | yes |
| 4 | the coil, and a correlation of our own | yes |

A study is a **definition** and its **results**, kept apart:
`studies/<name>/study.json` says what is to be run and is written before
anything is launched, so the matrix is a decision somebody made rather than
whatever happened to get run; `results.json` beside it gets one record per case
as it finishes, so a campaign that stops after nine of fifty picks up at the
tenth. `python3 study.py --make` builds the four straight-rod studies,
`--list` says how far each got, `--report NAME` writes its report.

**Stage 1 is done and needs no solver.** `correlations.py` holds every
correlation with its source, the range it is allowed to be asked, and a
status: the equations are here, or only the reference to them is. Four of the
eight are `needs-source` **on purpose** — a coefficient written down from
memory is worse than an absent one, because it runs, it looks plausible, and
nothing ever flags it. Four are encoded:

* **Jakob (1938)**, Holman's form, the one the Geometry tab has always shown.
  `--check` lifts `gapVelocity` and `lossModel` out of `mesh_explorer.html`
  and runs them against the Python over 96 cases and four quantities; they
  agree to 3e-16. It runs the page's own code rather than a copy of the
  formula, because a third copy would pass while the page said something else.
* **Žukauskas (1972)**, off the charts. There is no closed form to transcribe,
  so the four bicubic fits are copied verbatim from the MIT-licensed
  [`ht`](https://github.com/CalebBell/ht) library — a digitisation of the
  figures as reprinted in Incropera — with the licence in `third_party/` and
  the provenance in the module header. `bisplev` is written out in plain
  Python rather than adding scipy for four fixed tables, and matches scipy to
  1e-14; both of `ht`'s own documented examples reproduce to the digit.

* **Gunter & Shaw (1945)**, on the volumetric hydraulic diameter — the only
  one here that does not separate in-line from staggered. Transcribed from two
  independent secondary sources that agree, and checked three ways: its `D_v`
  reproduces the 0.1334 m the KAERI CHX paper states for its bundle to 6e-6 m,
  its two branches meet at the stated transition of `Re_v = 200` to 1.1 %, and
  its magnitude lands between Žukauskas' two branches.
* **Shen et al. (2024)**, the helical-bundle correlation from *Annals of
  Nuclear Energy* 201 110442. Its friction factor is defined as
  `2Δp/(ρu_max²z)`, which is this library's Euler number per row **exactly** —
  the paper and this project already speak one currency. Its pitch exponent
  reproduces the paper's own quoted numbers, which is the transcription
  confirmed.

Everything is converted to **one currency**, the Euler number per row, because
the friction factors are not comparable: Jakob's appears as
`dp = 2 f N rho u_max²` and Žukauskas' as `dp = N chi f (rho u_max²/2)`.
`u_max` is decided once — including the staggered diagonal-gap test — instead
of inside each correlation, and the browser and the Python now share that rule.

The stage-1 report **computes** the spread between them rather than asserting
it, and the result is the most important thing stage 1 has to say: **there is
no condition at which the published correlations are tight enough to call one
of them the answer.** The three that apply to a straight bundle differ from
one another by 1.0× to 2.2× over the grid, typically about 1.45×. Two of them
look close; adding a third widens it.

That sets a ceiling on what stage 3 can conclude. Landing inside 20 % of one
correlation does not make a CFD result right and missing one by 40 % does not
make it wrong, so the sweep reports whether the CFD falls **inside the band
the in-range correlations span** — with any correlation asked outside its own
limits daggered and left out of the band — and whether it has the same
**slopes**, in Re and in pitch.

It also says why the project exists, in numbers: on the real CHX helical
bundle the CFD gave 288.2 Pa, and Žukauskas was 61.9 % out, Gunter & Shaw
45.2 % out.

**Stage 2 was run once, steady, and every case of it failed.** All eight
stalled at a residual of 7.5e-2 in continuity after 800 iterations, and the
Δp kept falling with every refinement without settling. That is not a mesh
being too coarse. A strictly two-dimensional bank at `Re_max` of 1e4 and above
has **no steady solution to converge to** — a 2-D bluff-body wake is
time-periodic above `Re` of order 200, and these are fifty to five hundred
times that. The domain's own symmetry planes, which are there so the geometry
matches what a correlation describes, remove the spanwise decorrelation that
would otherwise break the vortices up, so this is the most shedding-prone
version of the problem rather than the least. Shen et al. solved the same
problem with URANS and reported time averages.

So the campaign runs **transient**, and the answer of a case is the **time
average** of its bundle pressure drop, not whatever the last instant happened
to be. The time step is **derived, not typed** — the same rule as the first
cell height, which follows from a y+ target rather than from somebody's
judgement:

```
f = St·u_max/D,  T = 1/f,  Δt = T/25,  20 periods run, the last 15 averaged
```

`St = 0.2`. Getting that wrong by 30 % costs 30 % of the run time; getting it
wrong by a factor of ten loses the oscillation entirely, which is what typing
a number would eventually do.

**Stages 2 and 3 need a licence to produce a number.** Three decisions in them
are worth arguing with before you press Run:

*Every wall that is not a rod is a symmetry plane.* The correlations are for a
bank that is infinitely wide and made of infinitely long tubes; a box with four
real walls is not that, and both biases grow as the domain shrinks. With half
rods on the sides the mirror lands on a rod centreline and gives back the bank
the correlation describes. It also makes the case effectively two-dimensional,
which is why sixty of them are affordable at 18k–74k cells each.

*The mesh rules are applied per case rather than typed into a panel.* The first
cell is pinned at y+ = 1, so its height follows `u_max` — which moves with both
the velocity and the pitch — and the radial layer count is then whatever it
takes to reach the gap at a bounded growth ratio. Holding the layer count fixed
instead would let the growth ratio run from 1.2 at the slowest case to 3 at the
fastest, and the fast cases would come back wrong for a reason having nothing
to do with the physics.

*X_T and X_L vary independently.* On the square diagonal the two pitch terms of
any `(X_T-1)^-p X_L^q` form are collinear and neither coefficient can be
identified — a sweep down the diagonal would produce a fit that cannot be
fitted, and the mistake is invisible until the solve is singular.

The mesh study answers with a **grid-convergence index** by the procedure in
Celik et al. (2008) / ASME V&V 20, implicit observed order and all, because a
block mesh refines by integer counts and never lands on the ratio asked for.
Every number of that paper's own worked example reproduces — p 1.534, φ_ext
6.1685, GCI 2.17 %. It names the **coarsest** mesh inside tolerance: the point
of a mesh study is the cheapest adequate mesh, not the finest one that fits in
the night.

**Stage 4 starts from a published shape, not an invented one.** `FIT_FORM` is
the skeleton of Shen et al.'s eq. (12) — a laminar plus a turbulent term in
Re, a pitch factor, and a helix factor that is exactly 1 when the bundle is
straight — refitted on our data rather than proposed from scratch, which is a
far smaller claim and makes the straight-rod fit and the coil fit the *same*
correlation at two values of one angle.

With one generalisation. Shen's pitch factor is `(X_T X_L)^-0.69`: the two
ratios only ever appear as a product. That is not a physical claim but a
consequence of their dataset — every case in that paper had `S_T = S_L`, so no
data of theirs could split the exponents. It is the same collinearity the
stage-3 matrix was laid out to avoid. Given separate exponents and
Jakob-shaped data, the fit error drops from 38 % to 9 % in-line and 18 % to
7 % staggered, with `p` much larger than `q` — the transverse pitch governs
and the longitudinal one barely enters, which is what `u_max` being set by the
transverse gap predicts. The sweep report fits both and prints the comparison,
so the question is settled on the data rather than argued.

The fitter is a nested Nelder–Mead over the exponents with a closed-form
non-negative least squares for `A` and `B` inside it, written out rather than
depended on. `A` and `B` are held non-negative because a negative laminar
coefficient has no meaning and, left free, the solve used one to fake a
steeper Re dependence than the form can otherwise produce. It recovers Shen's
own coefficients from Shen's own formula to 1e-15, which is the test that it
fits rather than merely converges.

**Several cases at once, and it finds out how many.** One case cannot fill a
workstation: the solve is a fraction of each case and the largest case in the
campaign is 44k cells per core, so the way to use the machine is more *cases*,
not more cores per case. What limits that is solver **tasks** in the licence,
which the app cannot know — so it finds out. Workers start one at a time and
each has to get a session up before the next is added; a launch that fails on
anything that reads like a licence caps the ramp there and the queue is
finished with the workers it has. Measured on a 30-case queue: 3.3× at four
workers, 5× at eight, and a licence that refuses past worker 2 still completes
all 30. The first case still runs **alone**, as the smoke test — if the set-up
does not produce a usable solution, running seven more of it in parallel only
wastes the machine faster.

**A case that already ran opens again without a solver.** A run leaves two
different things behind. The *solution* — the field you draw contours of —
lives in a Fluent case file and needs Fluent and a licence to reopen. The
*record* — how the residuals came down, how the pressure drop moved, how big
the mesh was, what the bundle drop measured — is in `results.json` and needs
nothing at all. Clicking a case in the campaign tab puts its histories back on
the Run tab and says in the banner that this is a record, not a run. Mesh
studies write their case files too, so the fields are recoverable there; the
thirty-case sweeps do not, and say so.

**A campaign is read as a wall, not as a stack of cases.** The stage report
argues end to end and draws the residual curve of the one case its argument
turns on. The question asked far more often is "how did the whole batch go?",
and answering that by opening eight cases one at a time is how a bad set-up
survives a night. **한눈에 보기** puts every case on one screen: one panel each,
with its status colour, its numbers, its residual envelope and its Δp trace —
and every panel on the **same axes**, because the comparison is the whole
point. A flat residual curve beside a descending one is obvious; a flat curve
on its own is not. The abscissa is each run's progress from 0 to 1 rather than
its iteration count, so runs of different length still compare by shape, and
the Δp panel is drawn as a percentage of that case's own settled mean, so cases
whose pressure drops differ by orders of magnitude are still on one scale.
Under the wall the same traces are overlaid, coloured by status rather than by
case — with thirty cases and six palette colours a per-case key would put the
same blue against six names, and the question a stack of curves is read for is
which ones went wrong. Clicking a panel replays that case. Nothing here is
computed that the report does not already compute; it is the same record,
arranged for the eye.

**The shipped settings tree is a claim about a release, not about an
installation.** `audit_choices()` was green on all five trees for
`setup.general.solver.time = 'transient'` — the 2026 R1 tree lists `transient`
among its allowed values — and a real 2026 R1 refused it 61 seconds into a run,
naming four values, none of them `transient`. So the value is negotiated with
the session in front of us: the driver asks `allowed_values()`, falls back to
trying candidates in order, logs which one won, and records it on the result
row. The order is **second-order implicit first**, and not because it is newer:
the time step is sized from the shedding period as T/25 to resolve an
oscillation whose amplitude is the measurement, and first-order implicit damps
exactly that. Preferring it everywhere also means the 2025 R1 in the lab and
the 2026 R1 on the student machine advance time the same way, so their answers
can be compared. A study whose cases did not all use one scheme says so in red;
one that fell back to first-order says so in amber. The journal carries the
same negotiation, since it runs without the app.

**Monitoring is not free, and on a transient run it is the whole argument.**
Every surface integral is a round trip to the solver, and on this link a round
trip costs about 0.2 s whatever it asks for — measured, not guessed, by
decomposing two runs of different mesh size at the same iteration count. A
transient case takes 500 time steps, so what happens inside the per-step
callback is multiplied by 500. The Δp trace is sampled **every step and
deliberately not on a cadence**: it is a shedding oscillation resolved at 25
steps per period, and sampling every other step aliases the frequency the run
exists to capture. But it reads the two bundle planes only — the inlet/outlet
pair spans the inlet and outlet boxes and is not what the study reports — and
the residual history, which is read for its shape and polled by the browser
about once a second, goes on a cadence. Five round trips per step became two:
about two minutes of monitoring per case instead of ten.

Do not go looking for the residual history in the `.cas.h5`. Reopening one
needs Fluent and a licence, and what Fluent puts in it is the converged
**field**, not the road the solver took to get there — the app tries anyway
when it loads a case and says in the log which of the two it got. The traces
are in `results.json`, which is why they are kept there: the Δp trace is
capped at 2000 samples so a 500-step transient run is stored whole (thinning
a shedding oscillation at a non-integer stride aliases it, and the frequency
is a thing worth reading off later), the residual trace at 400, which draws
the same shape a six-thousand-iteration history would. **CSV로 내보내기** in
the campaign tab writes them out one file per case, under
`studies/<name>/history/`, for plotting in anything.

**A mock run is not a result, and the tab is built so it cannot become one.**
Every mock row is badged MOCK, the progress counter does not count it, and the
report says at the top that it excluded them and then declines to draw any
conclusion rather than producing a grid-convergence table from invented
numbers. Stages 2 and 3 must be run on a machine with a Fluent licence.

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
