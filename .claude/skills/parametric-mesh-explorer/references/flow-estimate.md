# Pre-CFD flow estimate

The `vel` and `pres` display modes solve a cheap flow field **on the very mesh
that will be exported**, so the user can see roughly what they are about to
compute before committing a solver run. It is an estimate, and the interface
says so in every place it appears.

Its real value is not the picture. It is that a field solve exercises the mesh
end to end — connectivity, face areas, boundary patches, units — and reports a
mass-balance residual. A mesh that solves cleanly here is a mesh whose topology
and patches are right.

## What it solves

Potential flow on the 2-D section, by finite volume on the welded quad mesh:

```
∇²φ = 0
inlet   ∂φ/∂n = -1        (unit inflow; the solution is u/U, dimensionless)
outlet  φ = 0
walls   ∂φ/∂n = 0          (rods, side walls)
u = ∇φ                     (Green-Gauss reconstruction)
```

Solving with a unit inlet velocity and scaling afterwards keeps units out of the
linear system entirely — `∇φ` comes out as `u/U` and multiplying by `P.vel`
gives m/s. Do not put physical velocity into the boundary condition.

**Discretisation.** Two-point flux `a = L/|d|` on internal faces; a Dirichlet
face uses the normal distance from the owner centroid. The matrix is symmetric
positive definite, assembled in CSR, solved by Jacobi-preconditioned conjugate
gradient started from `φ = x - lTot` (uniform flow). That initial guess already
satisfies the equation everywhere except near the inclusions, which cuts the
iteration count sharply — ~600 iterations and 0.4 s for 12 k cells.

**Reconstruction.** Green-Gauss, `u_P = (1/A) Σ φ_f n̂_f L_f`, with distance-
weighted interpolation on internal faces, `φ_f = φ_P` on no-flux walls, and
`φ_P - d_n` at the inlet.

## Why the pressure needs a correlation on top

Potential flow has no drag — d'Alembert's paradox — so pressure recovers
completely behind every inclusion and the bundle Δp comes out as **zero**. A
pressure picture built from Bernoulli alone is therefore worse than useless for
the quantity engineers actually care about.

The fix is to superpose an empirical streamwise loss:

```
p(x,y) = ½ρ(U² - |u|²) - loss(x)
```

`loss(x)` accumulates duct friction through the inlet and outlet boxes
(Blasius or laminar `64/Re` on the hydraulic diameter) and one smoothstep per
row through the bundle, each carrying `dpBundle / nRows`. The local detail comes
from the potential solution; the total is an engineering number.

For an **in-line** tube bank the reference uses Jakob's friction factor in
Holman's form:

```
f  = [0.044 + 0.08 X_L / (X_T - 1)^(0.43 + 1.13/X_L)] · Re_max^-0.15
Δp = 2 f N ρ u_max²          Eu per row = 4 f
```

The factor 2 matters and is easy to get wrong — the same `f` appears in the
literature with and without it. Spot-check any correlation you adopt against
Zukauskas' charts before trusting it: `X = 1.5` at `Re = 10⁴` should land near
`Eu ≈ 0.35` per row, `X = 2` at `Re = 4·10⁴` near `0.15`. Without the factor 2
this formula gives half those values, which is how the error was caught.

A staggered array, a coil, or a PCHE channel needs its own correlation. Keep it
in one function (`lossModel`) so swapping it is a single edit, and name the
source in a comment.

## Streamlines

Traced by RK2 with an arc-length step over a coarse raster of `u` (320 columns,
holes dilated from neighbours), seeded across the inlet. Integration stops at
the domain edge or when the point enters an inclusion — tested exactly against
the inclusion geometry, not against the raster, so a streamline never crosses a
solid.

They are drawn in both the section view and on the 3-D section plane, between
the two geometry slabs, so anything standing above the cut hides the stretches
behind it. Colour follows the map underneath: light on the dark end of a
sequential map, dark on a diverging one.

## Verification that actually bites

These are identities, not tolerances — they either hold or something is wrong.

| Check | Expected | Catches |
|---|---|---|
| `∫ u_x dA` over a plain slab / `W·L_slab` | `1.000` | reconstruction, face areas, orientation |
| same over the bundle slab | `1.000` | blockage handled consistently |
| area-weighted mean `u_y` over the whole domain | `0` to machine precision | asymmetric assembly |
| outlet flux vs inlet area | equal | solver convergence |
| same geometry entered in mm and in m | identical `u_max`, `Δp`, `Δt` | unit handling everywhere |
| upstream and outlet `|u|/U` | `1.000` | scaling |

The first one is the strong check: 2-D continuity fixes the plane-integral of
`u_x` at `W·U` for every `x`, so integrating over a slab gives `W · L_slab`
exactly. The reference reproduces it to 0.04 %.

Also worth reporting: peak `|u|` from the solve against the 1-D gap value. They
differ — potential flow peaks *at the wall*, so the solved peak runs about 1.5×
the gap mean — and showing both is more honest than showing either alone.

## Cost control

- skip entirely above `FIELD_MAX` (80 k section cells) and say so
- never solve while `interacting` is true
- cache on `G`; `build()` replaces `G`, so a parameter change invalidates it for
  free. `G.field === undefined` means "not attempted", `null` means "attempted,
  unavailable" — keep that distinction, the UI reports different things.

## Honesty rules

No separation, no wake, no turbulence. Past the first row the leeward picture is
wrong, and the front/rear stagnation points come out near-symmetric, which is
the visible signature of that. The peak velocity is also resolution-dependent
(4.51 m/s at 12×3 versus 6.05 m/s at 36×12 in the reference case).

So: label the modes "estimate", keep the standing caveat in the read-out, and
frame the purpose as checking mesh and boundary conditions and getting the
scale — never as a substitute for the solver run. If a user asks for something
this method cannot give — wake structure, separation point, heat transfer — say
so rather than presenting a potential-flow answer.
