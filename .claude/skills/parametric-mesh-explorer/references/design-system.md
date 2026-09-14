# Design system

The visual and interaction rules this tool settled on, and the constraints
behind them. The look is deliberately plain: an engineering instrument, not a
consumer app.

## Ground rules

**White, single theme.** Background is pure `#ffffff`, forced. No dark theme, no
`prefers-color-scheme`, no theme switching — the user works next to plots and
papers that are white. Removing the dark path also removes a whole class of
contrast bugs.

**No borders. Depth by shadow.** Nothing in the UI uses `border`. Panels sit on
two-stage shadows (a tight one and a wide soft one) plus a 1 px inset highlight
along the top edge, which reads as a lit surface. Verify with a DOM sweep that
`borderTopWidth` is `0` everywhere.

**No faded text.** Secondary copy is near-black `#22262b` (~15:1) and primary is
`#15181c` (~17.8:1). Hierarchy comes from size and weight only. This was an
explicit preference and it holds up: grey small text on white is the first thing
to fail on a projector or a bad monitor.

**No fluorescent colour.** Accent is deep navy `#0f4c81`. Audit with an HSV
scan: nothing in the cyan band (hue 165–200°) and nothing with saturation > 0.72
at value > 0.82. A teal/cyan accent was rejected for looking neon on white.

**Contrast is measured, not eyeballed.** Compute WCAG ratios for every token
against the background and keep them ≥ 4.5. Current: text 17.8, secondary 15.2,
accent 8.9, ok 6.1, caution 6.4, bad 6.2.

## Type

```
Pretendard Variable, Pretendard, Noto Sans KR, Malgun Gothic, Arial, sans-serif
```

Pretendard is SIL OFL — free for commercial use and redistribution, which makes
it the safe choice for a product that will be registered. Loaded from a CDN as
progressive enhancement; offline it falls back silently to the system Korean
face. Numbers use the same family with `font-variant-numeric: tabular-nums`
rather than a separate monospace, so columns align without a second typeface.

Canvas text reads `getComputedStyle(document.body).fontFamily` so it follows
automatically.

## Avoiding trade dress

This tool is headed for software registration, so the design deliberately avoids
another product's recognisable signature. The distinction that matters: a
*technique* is common property, a *specific combination* is not.

Avoid: vendor font names in the stack (`SF Pro`, `-apple-system`, `Helvetica
Neue`); another OS's exact system colour values; the sliding-pill toggle switch;
a segmented control with a raised light pill on a filled track; full-pill radii
(`980px`); platform emoji as icons.

Fine to use: translucent blurred surfaces (Windows Acrylic, Android, and web
glassmorphism all do this — it is an industry-wide pattern), rounded rectangles,
soft shadows, a circular slider thumb, an accent colour.

What this tool does instead: square-ish 19 px checkbox with a check mark; a
segmented control with a *sunken* track and an accent-tinted active chip; 12–16
px radii; hand-drawn inline SVG line icons; its own navy/neutral palette.

Say plainly that this is a design judgement, not legal advice, and that counsel
should review before registration.

## Tokens

Single `:root` block. Layered translucency plus shadows:

```
--glass  panel fill          --glass2 raised (buttons)     --glass3 hover
--sunk   pressed track       --hl     top inset highlight
--sh-1/2/3   tight / panel / modal shadow stacks
--blur   blur(20px) saturate(120%)
```

Blur only does visible work where something varies behind the surface — the
overlays above the 3-D canvas. On the flat sidebar the same tokens read as a
light panel with a soft shadow, which is fine.

## Components

- **Panels** — 16 px radius, `--glass`, `--sh-2`, inset top highlight
- **Buttons** — 12 px radius, `--glass2`, translate 1 px down on press and drop
  the shadow; primary is the accent fill
- **Checkbox** — 19 px rounded square, accent fill + check when on
- **Segmented** — sunken track, accent-tinted active chip
- **Slider** — 4 px sunken track, 14 px accent thumb with a soft ring; the
  filled portion is painted by JS as a gradient
- **Number + slider pairs** — every parameter has both. The number field accepts
  values outside the slider range (the slider clamps its handle, the value does
  not), supports ↑↓ arrows with Shift for ×10, and never overwrites itself while
  focused. Sliders alone cannot hit "D = 9.5 exactly", which is what a user with
  a spec sheet needs.
- **Status line** — kept as data (`{ok, args}`), not as a formatted string, so it
  re-renders when the language changes

## Localisation

One `STR` object, two language keys, and `setLang(l)` that walks the DOM by
selector — no `data-*` attributes in the markup, so the HTML stays clean and
adding a string never means touching two places.

Dynamic text (warnings, quality summaries, status) lives in the dictionary as
**functions** taking the formatted numbers. Anything transient that is displayed
must be stored structurally so it survives a language switch — a formatted
string left on screen will stay in the old language, which is exactly the bug
that showed up first.

The toggle button shows the language you would switch *to*. Choice persists in
`localStorage` (wrapped in try/catch).

## Interaction

- drag to orbit; **Shift or middle-drag to pan**; wheel to zoom; double-click
  resets the view; a "fit" button clears zoom and pan
- heavy drawing is deferred: while dragging or moving a slider, render the light
  version and restore the full mesh ~180 ms after input stops
- a coordinate triad is pinned in the corner, depth-sorted so the axis pointing
  away is drawn first
- overall dimensions are always shown as offset dimension lines with extension
  lines, arrowheads and a value chip; a dimension that projects to less than
  ~14 px is edge-on and is hidden rather than collapsed to a stray label
- boundary-condition arrows sit **outside** the domain, upstream of the inlet
  plane, and are drawn before or after the solids depending on whether that face
  points at the camera

## Physics in the UI

A geometry tool becomes a meshing tool when it answers the question the user
actually has. First-layer thickness is meaningless on its own; y⁺ is what they
need. Take inlet velocity and fluid, use the **gap velocity**
`U_max = U·S_T/(S_T − D)` for a bundle, and show `Re`, `u_τ`, the required first
layer and the y⁺ the current setting gives — with a button to apply it.

```
Cf = 0.079·Re^-0.25      u_τ = U_max·√(Cf/2)      Δ₁ = 2·y⁺·ν/u_τ
```

The factor 2 puts the first *cell centre* at the target y⁺, which is the finite
volume convention. State the convention in the UI; different codes assume
different things and a silent factor of 2 is a real trap.

## Case sharing

Encode the parameter set into the URL hash (`k=v&…`, ~180 chars) and restore it
on load; fall back to `localStorage` when there is no hash. A colleague gets the
exact geometry from a link, and a parameter study survives a reload.

## Colour maps for fields

Quality modes (skew, aspect ratio, size) use the green-amber-red `ramp()`: it
reads as a verdict, which is what a quality metric is.

Field modes must not reuse it — a velocity or pressure map is data, not a
judgement, and a red patch should not imply "bad".

- **velocity** — sequential, dark blue through teal and green to yellow
  (viridis-like). Zero is the dark end, so stagnation reads as absence.
- **pressure** — diverging, blue through a neutral grey to warm red, spanning
  the actual min and max. The neutral middle keeps the eye on the extremes.

Both are anchor lists interpolated by `cmap(t, C)`. The colour bar always
carries real units, not `min`/`max` — a scale without numbers invites the reader
to over-read the picture.

Streamlines take their colour from the map underneath: near-white over the dark
end of the sequential map, near-black over the diverging one. Neither is the
accent colour; the accent stays reserved for interface state.

## What the toggles mean

Worth restating here because it is a design rule, not just code structure: the
mode buttons decide **what is drawn**, the mesh-preview checkboxes decide
**whether lines are drawn on top of it**. A user who turns off "grid preview"
while looking at a pressure field wants a cleaner field, not a blank pane.
