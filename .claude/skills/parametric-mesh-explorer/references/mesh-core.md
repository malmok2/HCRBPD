# Mesh core

How the block-structured mesh is built and written. Read the sections you need.

- [Blocks and transfinite grids](#blocks-and-transfinite-grids)
- [Topological node keys](#topological-node-keys) ← the part that breaks if rushed
- [Canonical vertices and exact interpolation](#canonical-vertices-and-exact-interpolation)
- [Cell polygons: in-line and staggered](#cell-polygons-in-line-and-staggered)
- [Extrusion and structural face enumeration](#extrusion-and-structural-face-enumeration)
- [Boundary patches](#boundary-patches)
- [Units](#units)
- [Writers](#writers)

## Blocks and transfinite grids

A block is `grid[i][j]` — an `(ni+1) × (nj+1)` array of 2-D points, plus a
`key(i,j)` function and a `zone` tag. Quads are formed from adjacent grid
entries and flipped to positive area if needed.

For a block bounded by a curve and a straight side, with straight sides joining
their endpoints, the Coons patch reduces exactly to

```
P(u,v) = (1-v)·B(u) + v·T(u)
```

so a simple blend between the two opposite curves is the correct transfinite
map, not an approximation. This is why the browser preview matches what a
mapped mesher would produce.

Spacing along an edge with `n` cells and a constant ratio `q`:

```
first = L(q-1)/(qⁿ-1),  then multiply by q each step
```

`q` is solved by bisection from a requested first-layer height. `scale = q^(n-1)`
is the "last / first" ratio — the same quantity SALOME's `NumberOfSegments`
takes, which is handy when cross-checking against a mesher.

## Topological node keys

Adjacent blocks must agree on their shared nodes. Two ways to do that:

**Coordinate hashing (do not use).** Round positions to a grid and hash. Two
points that should be identical but differ by ~1e-16 can straddle a bucket
boundary and hash apart. The failure rate is tiny per node, which is exactly
what makes it dangerous: it passes for most parameter sets and produces a crack
in a mesh someone ships.

**Topological keys (use this).** Each node gets a *name* derived from what it
is, not where it is. Blocks that share an entity generate the same name, so the
weld is exact and deterministic.

Resolve keys in strict priority order — corner first, then edge, then interior:

```js
key(i,j){
  if (corner)        return "v"+vertexId;                  // shared by 3+ sides
  if (on cell side)  return sideKey(ia, ib, i, n);          // shared by 2 cells
  if (on radial ray) return "r"+cellId+"_"+vertexId+"_"+j;  // shared by 2 sectors
  if (on the arc)    return "a"+cellId+"_"+sector+"_"+i;
  return "o"+cellId+"_"+sector+"_"+i+"_"+j;                 // block-interior
}
```

The corner case is the one that bites. A polygon corner belongs to **two sides
of its own cell and one or more sides of neighbouring cells** — up to three
distinct side keys for one physical point. Naming it by side gives it three
identities and tears the mesh. In the rod bundle this produced 2600 cracked
faces before corners got their own `"v"+id` key.

Edge keys must be direction-free. Walk every side from its lower vertex id:

```js
sideKey(ia, ib, i, n){ return ia<=ib ? `s${ia}_${ib}_${i}` : `s${ib}_${ia}_${n-i}`; }
sidePt (V, ia, ib, i, n){ return ia<=ib ? lerp(V[ia],V[ib], i/n)
                                        : lerp(V[ib],V[ia], (n-i)/n); }
```

Verify the scheme rather than trusting it: when a key is seen twice, compare the
two positions and track the maximum disagreement. **It should be exactly 0.** A
non-zero value means two blocks used the same name for different points; a crack
in the watertightness check means they used different names for the same point.

## Canonical vertices and exact interpolation

Neighbouring cells compute a shared corner from different formulas
(`cx + XE` vs `cx' - XV`), which agree algebraically but not in floating point.
Merge all polygon corners once into a canonical table with a tolerance, then
have every downstream formula read canonical coordinates. After that, shared
points are bit-identical by construction.

Interpolation must be endpoint-exact:

```js
function lerp(a, b, t){
  if (t === 0) return a;          // a + (b-a)*0 is fine, but be explicit
  if (t === 1) return b;          // a + (b-a)*1 is NOT reliably b
  return [a[0]+(b[0]-a[0])*t, a[1]+(b[1]-a[1])*t];
}
```

Points on a curved wall should come from the direction to the canonical vertex
(`c + R·û`), not from `atan2` + `cos/sin`, wherever two sectors must agree:
angle unwrapping (`θ` vs `θ+2π`) gives different trig results for the same
direction.

## Cell polygons: in-line and staggered

The cell polygon is the tile an inclusion gets to itself, and its shape decides
how clean the O-grid around it can be.

**In-line lattice — a rectangle.** The Voronoi cell is simply `SL x ST`. Its
left and right sides already lie on the zone-interface planes, so unlike the
staggered cell it needs no squaring-off, and every block stays a clean quad.
Clip it against the side walls and de-duplicate. This is the easy case; prefer
it when the physics allows.

**Staggered lattice — a hexagon.** For a staggered lattice with transverse pitch `ST` and row spacing `SL`, the
Voronoi cell of each rod is a hexagon with flat top and bottom:

```
half height  HH = ST/2
half length  XV = SL/2 + ST²/(8·SL)      (left/right vertex)
flat edge    XE = SL/2 − ST²/(8·SL)      (where the flat side starts)
vertices: (±XV, 0), (±XE, ±HH)
```

Two structural details that make the decomposition work:

**Square off the first and last rows.** A cell's pointed end sticks past the
zone interface and leaves an uncoverable notch. Replacing the spike with a
vertical edge at the interface plane keeps six sides, tiles the strip exactly,
and lets the neighbouring region stay a clean rectangle. The price is a 90°
corner split by one radial line, so those cells sit at ~45° skew — that is the
floor for this construction, and it is worth stating rather than hiding.

**Rod-to-cell clearance** is `min(ST/2, S_D/2, SL/2) − R` with
`S_D = hypot(SL, ST/2)`. Verified against direct point-to-segment distance
across the parameter space. Refuse to build when it is ≤ 0.

Where an inclusion is absent (walls without half rods), keep the cell and mesh
it as one plain block so the tiling never has holes.

## Extrusion and structural face enumeration

Mesh the section once, then extrude `nZ` layers. Cell index `q*nZ + l` (quad
`q`, layer `l`). Node index `l*N2 + n`.

Do not hash faces — enumerate them structurally:

- **vertical faces** ← 2-D edges: interior 2-D edge → internal face; boundary
  2-D edge → boundary face. `nZ` faces each.
- **horizontal faces** ← stacked quads: `nZ-1` internal per quad, plus the two
  caps.

Totals: `internal = E_int·nZ + Q·(nZ-1)`, `boundary = E_bnd·nZ + 2Q`. This
matched a hash-based Python implementation exactly (334,268 / 34,184), which is
a good sanity check to repeat when porting.

Hexahedron node order: bottom quad CCW seen from +z, then the same four on the
top. Outward faces in that order:

```
(0,3,2,1) −z   (4,5,6,7) +z   (0,1,5,4) −y
(1,2,6,5) +x   (2,3,7,6) +y   (3,0,4,7) −x
```

## Boundary patches

Classify 2-D boundary edges by **endpoint** position. A curved wall's chord
midpoint is at `R·cos(Δθ/2)`, which is 0.4 % off at 36 divisions but 3.4 % off
at 12 — a midpoint test with a 2 % tolerance passes at the default and fails at
coarse settings. Endpoints are exactly on the surface at any resolution.

Patch order is fixed and shared by every writer, so face indices line up:

```
inlet, outlet, wall_rods, wall_side_y0, wall_side_ymax, wall_bottom, wall_top
```

## Units

Model numbers are unitless; the UI declares what they mean (`mm` or `m`) and
**every exported coordinate is scaled to metres**. Route all writers through one
`eachPoint()` / `zpt()` pair so the scale cannot be forgotten in one of them.

- OpenFOAM: no unit metadata, always metres → scaling is mandatory
- Fluent / CFX: the importer asks, but a metre-scaled file needs no thought
- Physics (Re, y⁺) must use the same scale, or it is wrong by 10³

## Writers

**Fluent ASCII `.msh`** — section headers and connectivity are **hexadecimal**;
coordinates are decimal. Node/cell/face indices are 1-based.

```
(2 3)                                   dimension
(10 (0 1 <nnodes> 0 3))                 node declaration
(10 (1 1 <nnodes> 1 3)( x y z ... ))
(12 (0 1 <ncells> 0)) (12 (2 1 <ncells> 1 4))    hexahedra
(13 (0 1 <nfaces> 0))                   face declaration
(13 (<zone> <first> <last> <bc> 4)( n0 n1 n2 n3 c0 c1 ))
(45 (<zone> <type> <name>)())           zone naming
```

BC codes: `2` interior, `3` wall, `5` pressure-outlet, `10` velocity-inlet.
Face zones must be contiguous and sum to the declared face count.

**Orientation — settled, do not re-derive.** Fluent and OpenFOAM use
**opposite** conventions, and this cost a wasted solver run to discover:

| | normal from the node order points |
|---|---|
| OpenFOAM polyMesh | owner to neighbour; outward on a boundary |
| Fluent `.msh` | **into `c0`**; into the domain on a boundary |

The face enumerators (`eachInternalFace`, `eachPatchFace`) hand out the
OpenFOAM orientation, so `foamFiles()` writes them unchanged and
`writeFluent()` emits every face with its **node order reversed** — internal
and boundary alike, keeping `c0` as the owner. Reversing the nodes is the whole
fix; do not also swap `c0`/`c1`. Those two operations cancel on an internal face
and put `0` in the owner slot on a boundary face, which is exactly the bug the
Python twin shipped with.

The evidence: ANSYS Fluent 2026 R1 read a mesh written the OpenFOAM way and
reported **610560 of 610560 cells with non-positive volume**, then read the
reversed file cleanly. One hundred percent of cells failing is the signature of
a global orientation flip — a geometry problem hits some cells, not all of them.

Guard it with the volume check below rather than a flag and a hope.

**Cell volumes, rebuilt from the faces.** Before writing anything, integrate
`x.n dA` over the faces of every cell the way a solver does, and refuse to
export unless all volumes are positive and the total matches the analytic
volume of the faceted domain. On an extruded mesh every face is planar, so
`S = ½(p2-p0)x(p3-p1)` with the four-node centroid is exact and the total lands
within ~1e-13. This is the check that would have caught the orientation bug at
the source; it costs ~0.3 s on 610 k cells.

It cannot validate the `c0`/`c1` convention itself — only a solver knows that —
so it verifies consistency and non-degeneracy under the convention stated above,
and the convention rests on the Fluent run. Say it that way in reports.

**OpenFOAM polyMesh** — four rules, all checkable:

1. `owner < neighbour` for every internal face
2. internal faces sorted by `(owner, neighbour)` — emit per cell in index order,
   taking the vertical-stack neighbour first, then quad neighbours ascending;
   that is already upper-triangular, no sort needed
3. face normal points owner → neighbour; on a boundary, out of the domain
4. boundary faces contiguous per patch, appended after all internal faces

**VTU** — `UnstructuredGrid`, cell type `12`. Carry quality fields (skew,
aspect ratio, volume) as `CellData`; it makes ParaView useful for spotting where
a mesh is bad rather than just that it is.

**STL** — binary, quads split into two triangles. Summing the closed-surface
volume and comparing with the mesh volume is a cheap, strong watertightness
proof of the boundary.

**ZIP** (for polyMesh, which is five files) — store or deflate. Use
`CompressionStream("deflate-raw")` when present (18 MB → 3.7 MB) and fall back
to stored entries. CRC32 is computed over the *uncompressed* bytes either way.

Blob building: push strings into an array and hand the array to `new Blob(...)`.
Never concatenate a 16 MB string.
