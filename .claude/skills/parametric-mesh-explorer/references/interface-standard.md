# Interface standard

The GUI is a house standard: a new geometry changes what the controls *mean*,
never how the tool looks or behaves. This file is the contract. Read it before
touching anything the user sees.

- [Window layout](#window-layout)
- [The sidebar, group by group](#the-sidebar-group-by-group)
- [Control patterns](#control-patterns)
- [Display modes and the section plane](#display-modes-and-the-section-plane)
- [Viewport furniture](#viewport-furniture)
- [Camera](#camera)
- [Localisation](#localisation)
- [Persistence and sharing](#persistence-and-sharing)
- [Adding a control without breaking anything](#adding-a-control-without-breaking-anything)

## Window layout

```
header  [ title · subtitle · EN/한국어 ]
├── aside .main            fixed-width scrolling sidebar of .grp groups
└── .viewport              flex, position:relative
    ├── #scene             the 3-D canvas, fills the pane
    ├── .inset             section view (x–y), top-right, with a zoom button
    ├── .warn              warning banner, top-left, hidden unless triggered
    ├── .legend            colour key, bottom-left
    └── .stats             8 derived numbers, bottom-right
.zoomwrap                  full-screen section view, opened from the inset
```

Below 900 px the sidebar moves above the viewport and the inset shrinks. There
is one theme (white) and no dark mode — see `design-system.md`.

## The sidebar, group by group

Every group is a `.grp` with an `<h2>`. The order is fixed; geometry groups come
first, then meshing, then physics, then display, then export.

| # | Group | Contents | Geometry-specific? |
|---|---|---|---|
| 1 | Inclusion | its size, channel height `H` | yes |
| 2 | Array | pitches, counts, a "snap to square" button | yes |
| 3 | Inlet & outlet boxes | `lIn`, `lOut` | mostly reusable |
| 4 | Mesh | azimuthal, radial, first-layer, `nZ`, upstream/downstream counts | mostly reusable |
| 5 | Flow conditions | `vel`, fluid preset, `nu`, `rho`, `ypl`, apply-to-first-layer, two read-outs | **reusable as-is** |
| 6 | Mesh display · quality | the six mode buttons, mode note, section-plane slider | **reusable as-is** |
| 7 | Options | display toggles + geometry toggles | partly |
| 8 | View | ISO / Top / Side / Front | **reusable as-is** |
| 9 | (unnamed) Export | unit segment, four format buttons, status note, share / fit / reset | **reusable as-is** |

Groups 5, 6, 8 and 9 are the ones the user means by "same GUI". Port them
verbatim; only the strings change.

### Group 5 — flow conditions

Serves two purposes and both are geometry-independent:

- **y+ sizing.** `flow()` computes gap velocity, `Re`, `u_τ` from Blasius, and
  the first-cell height that lands the wall node at the requested y+. `#physOut`
  shows `Re`, `u_τ`, the required first layer, and the y+ the current mesh
  actually gives. "Apply to first layer" writes it into the mesh parameter.
- **Field estimate.** `#flowOut` shows bundle Δp, Eu per row, solved and 1-D
  `u_max`, the CFL=1 time step, the mass-balance residual, and a standing
  caveat. See `flow-estimate.md`.

Only the gap-velocity expression and the pressure-loss correlation are geometry
knowledge. Everything else stands.

### Group 7 — options

Display toggles (fixed): `showMesh`, `mesh3d`, `showBox`, `showFlow`, `stream`,
`autoRotate`. Geometry toggles are added alongside — `halfRods` in the
reference — and these do change topology, so they call `rebuild()`.

## Control patterns

Three patterns, and nothing else. Consistency here is most of what makes the
tool feel like one instrument.

**Slider + editable number + derived hint.** The canonical parameter control.

```html
<div class="ctrl">
  <div class="row"><label for="ST">Transverse pitch S<sub>T</sub></label><output id="STV"></output></div>
  <input type="range" id="ST" min="..." max="..." step="..." value="...">
</div>
```

At startup a loop over `RANGES` replaces each `<output>` with a text input plus
an `<em>` hint. So one `.ctrl` block plus a `RANGES` entry gets you: a slider, a
typed value, arrow-key nudging (shift = ×10), clamping, integer rounding for
count parameters, a derived hint from `HINT`, share-link participation, and
reset. Do not hand-roll a control that already exists.

The hint is where the derived quantity goes — `P/D`, the growth ratio, the
resulting `Δz`, `u_max`. It is the single highest-value piece of the interface:
it turns a slider into a calculation the user can trust.

**Segmented control** (`.seg`) for small exclusive choices: display mode, view
preset, fluid, unit. Buttons carry `data-*` attributes and the handler toggles
`.on`.

**Checkbox** (`.chk`) for booleans, declared in `CHECKS`.

## Display modes and the section plane

Six modes in `#modeSeg`: `grid`, `skew`, `ar`, `size`, `vel`, `pres`. The first
four are mesh properties, the last two are the flow estimate.

**The semantic that must hold:** the *mode buttons* decide what is drawn; the
*mesh-preview toggles* decide only whether mesh **lines** are drawn on top. A
coloured map — quality or field — is the picture the user asked for, not an
overlay on a mesh preview, so it paints whether or not lines are on. Coupling
the two was a real bug: turning off "grid preview" blanked the velocity field.

Concretely:

- `showFill = canPaint && MODE !== "grid"` — never gated on `showMesh`
- `showWire = canPaint && P.showMesh && MODE === "grid"`
- per-cell outlines and the 3-D wall mesh follow `showMesh` / `mesh3d`

**The section plane.** In 3-D the map is drawn on a horizontal plane at
`z = P.zCut·H` (mid height by default), not on the top cap, so a contour reads
as a cut through the geometry. The rods are drawn in two slabs — below the plane
and above it — with the plane painted between them, so geometry standing above
the cut correctly hides what runs behind it. Which slab comes first is decided
by the depth component of the world `+z` axis, `rotate(0,0,1)[2]`, not by the
sign of pitch: `cos(pitch)` is even, so a camera tilted down is still looking
from above.

`zCut` is a **view** parameter. It is deliberately *not* in `RANGES`, because
every `RANGES` handler calls `rebuild()`, which rebuilds the geometry and throws
away the solved flow field. It has its own handler that only sets `needsDraw`.
Any future view-only control must follow the same rule.

## Viewport furniture

**Section inset.** The x–y section, always live, with a zoom button that opens
the same `paintSection()` full screen. One painter serves the inset, the zoom
view and the 3-D plane, so they can never disagree.

**Stats bar.** Eight cells, one line each: derived geometry, mesh size, worst
skew, worst aspect ratio. Keep it at eight — the bar is a flex row and a ninth
cell starts overlapping the legend. New numbers belong in a `.phys` read-out in
the sidebar instead.

**Warning banner.** One message at a time, worst first, with `.err` for
conditions that make the geometry invalid (overlapping inclusions) and plain
caution for advisory ones (short outlet box, first layer thicker than the gap).
An `.err` state must also block export.

**Legend.** One row per colour actually used in the current drawing.

## Camera

Orthographic, yaw/pitch, drag to rotate, shift/ctrl/middle-drag to pan, wheel to
zoom (0.3–8×), double-click resets to ISO and fits. Four presets: ISO, Top,
Side, Front. `busy()` marks a 180 ms interacting window during which the mesh
and the field are not redrawn — that is what keeps dragging smooth on a 10⁵-cell
mesh, and it is also why anything expensive must check `interacting` before
running.

## Localisation

Korean and English, switched by one header button, persisted in `localStorage`.
Every visible string lives in `STR.ko` / `STR.en`; `T(key)` reads the current
one. Canvas text reads the body font so it follows automatically.

`setLang()` walks the DOM and rewrites text. Several of its lookups are
**positional**, which is the main trap:

| Array | Bound to |
|---|---|
| `h2` | `aside .grp h2` in DOM order |
| `modes` | `#modeSeg button` in DOM order |
| `stats` | `.stats .k` in DOM order |
| `legend` | `.legend .li` in DOM order |
| `dl` | `.dl .btn span` in DOM order |
| `fluids` | `#fluidSeg button` in DOM order |

Add a group, a mode, a stat or a legend row and you must extend the matching
array **in both languages**, or the labels silently shift by one. `setLang` also
takes `document.querySelectorAll("aside .note")[0]` as the mesh note — do not
insert a `.note` element above it.

Labels and checkbox captions are keyed, not positional: `STR.lab[id]` matches
`label[for=id]`, `STR.chk[id]` matches the text node after the checkbox. Those
are safe to extend.

## Persistence and sharing

`SHARE_KEYS = RANGES + ["unit","fluid"] + geometry toggles`. `encodeState()`
builds a `k=v&…` string used both for `localStorage` (restored on load) and for
the copy-link button (written to `location.hash`, which wins over the stored
case). Display toggles and camera state are deliberately excluded — a shared
link should reproduce a *case*, not someone's view.

## Adding a control without breaking anything

A numeric parameter, end to end:

1. `DEFAULTS.newParam = …`
2. add `"newParam"` to `RANGES` (this also puts it in the share link)
3. add the `.ctrl` block with `<output id="newParamV">` in the right group
4. add `STR.ko.lab.newParam` and `STR.en.lab.newParam`
5. optionally `HINT.newParam = d => …` for the derived value
6. if it belongs to a preset (fluid, material), extend that preset map

A checkbox: `DEFAULTS`, `CHECKS`, the `.chk` markup, `STR.*.chk[id]`, and decide
in the `CHECKS` handler whether it needs `rebuild()` (topology) or only
`needsDraw` (appearance).

A view-only control: keep it out of `RANGES`, give it its own handler that sets
`needsDraw`, and add explicit lines to the reset handler and any sync function —
`RANGES.forEach` will not cover it.

Before declaring it done, check the four combinations that break most often:
each display toggle on and off in each mode, both languages, reset, and a
restored share link.
