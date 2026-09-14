#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helical coil bundle in cross-flow : structured all-hexahedral mesh generator
============================================================================
Standalone.  No SALOME, no CAD kernel, no third-party module - stdlib only.

A local patch of a helical coil steam generator, unrolled.  Wind a tube around
a shell at helix angle ALPHA and unroll that cylindrical surface flat: the tube
becomes a STRAIGHT line inclined by ALPHA.  So an unrolled patch is a STAGGERED
bundle whose tubes lean into the flow - which is exactly what this builds.

  x = shell axis  = cross-flow direction
  z = circumferential (unrolled)     y = radial, i.e. coil layer to coil layer

The tube axis is therefore (sin ALPHA, 0, cos ALPHA); at ALPHA = 0 this reduces
to the plain staggered bundle, which is used as a regression check.

WRAP_TO_COIL then rolls that patch back onto the shell, so the tubes become
real helices:  y -> radius R_IN+y,  z -> angle z/R_MID,  x stays the shell axis.
A tube advances x by z*tan(ALPHA), i.e. R_MID*theta*tan(ALPHA) per radian, so a
layer at radius r sees helix angle atan(tan(ALPHA)*R_MID/r) - the angle falls
off with radius exactly as it does in a real coil wound at constant axial pitch.

Two consequences of the inclination, both handled here:
  * a tube cut by a z = const plane is an ELLIPSE, semi-axes R/cos(ALPHA) along
    the flow and R across it.  Using circles instead would model a squashed
    tube, not a round one.
  * the pattern shifts downstream with height, dx = z*tan(ALPHA).  The WHOLE
    patch is sheared by that, boxes included, so A and C keep a constant
    thickness and their outer faces run parallel to the tube axes.  A uniform
    shear has unit determinant, so the volume check stays exact, and the lean
    no longer has to fit inside the inlet box - large sectors are free.

  flow +X     rod axis +Z     transverse +Y

        y
        ^     |<-- L_in -->|<---- L_bundle ---->|<-- L_out -->|
        |     +------------+--------------------+-------------+
        |     |            |  o     o     o     |             |
   W    |     |  A : inlet |     o     o        |  C : outlet |
        |     |            |  o     o     o     |             |
        |     +------------+--------------------+-------------+  --> x
                                  B : bundle

Every node position is closed form, so no boolean geometry and no meshing
algorithm is involved:

  B : rods sit on a staggered lattice, so the Voronoi cell of each rod is the
      hexagon cut by its six nearest neighbours - half height ST/2, half length
      SL/2 + ST^2/(8 SL) on the axis and SL/2 - ST^2/(8 SL) at the corners.
      Radial lines from the rod to the six corners cut it into one mapped block
      per side (an O-grid).  Radial spacing is geometric -> wall inflation.
  A,C : one mapped block per interface segment, graded towards the bundle.
  The 2-D section is then extruded N_Z layers -> 100 % hexahedra.

Compared with the in-line version the hexagon costs one special case: the first
and last rows are squared off against the A|B and B|C planes, because a hexagon
corner would otherwise poke through an interface that has to stay flat.  What
it buys is skew - a hexagon corner is 117-127 deg and, split by one radial
line, leaves about 60 deg, against the 90 -> 45 deg of a rectangular corner.
Only the squared-off end rows keep 90 deg corners, so they hold the worst cells.

Blocks share nodes through *topological* keys (not coordinate rounding), so
shared points are bit-identical and the mesh is watertight by construction.

Output
  <name>.msh   Fluent/ANSYS ASCII  -> ANSYS Fluent and CFX-Pre
  <name>_foam/ OpenFOAM polyMesh   -> copy into <case>/constant/polyMesh
  <name>.vtu   VTK unstructured    -> ParaView (mesh + quality fields)
  <name>.stl   boundary surface    -> figures / documentation

Run:  python helical_coil_staggered_mesh.py
"""

import math
import os
import struct
import sys

# =============================================================================
#  USER PARAMETERS
# =============================================================================
# --- geometry ---------------------------------------------------------------
D           = 10.0     # rod diameter                        [mm]
ST          = 20.0     # transverse pitch (y)                [mm]
SL          = 20.0     # longitudinal pitch (x)              [mm]  (= ST -> square lattice)

N_ROWS      = 5        # rod rows along the flow
N_COLS      = 4        # full lateral pitches

L_IN        = 60.0     # inlet  box length  (A)              [mm]
L_OUT       = 120.0    # outlet box length  (C)              [mm]

SECTOR_DEG  = 12.0     # HOW FAR THE PATCH RUNS AROUND THE SHELL, degrees.
                       # This is the real handle on "how much coil" is modelled;
                       # the tube length follows from it and the coil radius,
                       # H = radians(SECTOR_DEG) * R_MID, so growing the sector
                       # no longer means inventing a taller channel.
R_MID       = 200.0    # coil mid radius [mm]; the bundle spans R_MID +- W/2

HELIX_ANGLE = 15.0     # tube inclination from the z axis, degrees.
                       # 0 -> plain in-line bundle.  Typical HCSG: 5 .. 25 deg.

WRAP_TO_COIL = True    # True  : roll the patch onto the shell -> real helices
                       # False : keep it flat (unrolled), for comparison

HALF_RODS_AT_WALL = True   # True : the non-offset rows put rods on y = 0 and
                           #        y = W, halved by the walls, so the lattice is
                           #        periodic across them (N_COLS+1 rods per row)
                           # False: those two are left empty, walls stay clear.
                           #        Costs the length-proportional division count
                           #        (see side_div), so True is the better mesh.

FIRST_ROW_OFFSET  = False  # True : row 0 is the offset (mid-pitch) row.
                           # False: row 0 is the full row that touches the walls.

# --- mesh -------------------------------------------------------------------
N_AZ_MIN    = 36       # minimum azimuthal cells around one rod (rounded up to 6*n,
                       # one hexagon side at a time)
N_RAD       = 12       # radial cells, rod wall -> cell boundary
FIRST_LAYER = 0.15     # first cell height on the rod wall   [mm]

N_Z         = 10       # layers along the arc H, i.e. around the sector.
                       # This is the resolution in the tube-axis direction;
                       # SECTOR_DEG / N_Z is the angle each layer spans and
                       # sets how finely the arc is followed (1 = quasi 2-D).

N_X_IN      = 24       # streamwise cells in A
N_X_OUT     = 40       # streamwise cells in C
X_SCALE_IN  = 4.0      # A : last/first cell length ratio (fine at the bundle)
X_SCALE_OUT = 4.0      # C : idem

# --- output -----------------------------------------------------------------
UNIT        = "mm"     # what the numbers above mean: "mm" or "m".
                       # Exported coordinates are always written in METRES:
                       # OpenFOAM's polyMesh carries no unit metadata and is read
                       # as metres, so a mm-numbered mesh would be 1000x too big.
OUT_DIR     = os.path.dirname(os.path.abspath(__file__))
NAME        = "helical_coil_staggered"
WRITE_FLUENT   = True   # .msh   (Fluent + CFX)
WRITE_OPENFOAM = True   # polyMesh
WRITE_VTU      = True   # ParaView
WRITE_STL      = True   # boundary surface for figures

# Fluent face convention: the face normal built from the node order points
# from c0 to c1 (c1 = 0 on a boundary, normal outward).  Set False to swap.
FLUENT_NORMAL_C0_TO_C1 = True

# =============================================================================
#  DERIVED
# =============================================================================
THETA     = math.radians(SECTOR_DEG)                   # azimuthal span
H         = THETA * R_MID                              # arc length at mid radius
R         = 0.5 * D
W         = N_COLS * ST
L_BUND    = N_ROWS * SL
L_TOT     = L_IN + L_BUND + L_OUT
X_B0      = L_IN
X_B1      = L_IN + L_BUND
HH        = 0.5 * ST                                   # cell half height
XV        = 0.5 * SL + ST * ST / (8.0 * SL)            # half length on  y = cy
XE        = 0.5 * SL - ST * ST / (8.0 * SL)            # half length at y = cy +- HH
S_D       = math.hypot(SL, 0.5 * ST)                   # diagonal neighbour distance
ALPHA     = math.radians(HELIX_ANGLE)
TANA      = math.tan(ALPHA)
A_ELL     = R / math.cos(ALPHA)                        # tube footprint: semi-axis
B_ELL     = R                                          # along x / across it


def wall_gap(vx, vy):
    """Distance from the elliptical tube footprint to the Voronoi side that
    faces the neighbour at centre offset (vx, vy).

    That side is the perpendicular bisector of the offset: distance |v|/2 from
    the centre, unit normal v/|v|.  An ellipse reaches sqrt(A^2 nx^2 + B^2 ny^2)
    in direction n, so the gap is exact and closed form.  For a circle it
    collapses to the familiar |v|/2 - R."""
    d = math.hypot(vx, vy)
    return 0.5 * d - math.hypot(A_ELL * vx, B_ELL * vy) / d


#  A staggered cell is bounded by six neighbours: two across the transverse
#  pitch and four on the diagonal.  Which one is tightest depends on the pitch
#  ratio AND on the tube inclination, since the footprint is an ellipse, so both
#  directions are tested rather than assumed.
CLEARANCE = min(wall_gap(0.0, ST), wall_gap(SL, 0.5 * ST))
N_AZ      = max(1, int(math.ceil(N_AZ_MIN / 6.0)))     # divisions on a slanted side
AZ_STEP   = math.hypot(XV - XE, HH) / N_AZ             # target spacing on a cell side
SCALE     = max(L_TOT, W, H)
TOL       = 1.0e-9 * SCALE


def fail(msg):
    print("[ERROR] " + msg)
    sys.exit(1)


if UNIT not in ("mm", "m"):
    fail('UNIT must be "mm" or "m"')
EXPORT_SCALE = {"mm": 1.0e-3, "m": 1.0}[UNIT]     # model units -> metres


def XP(p):
    """A point on its way out of the program: always metres."""
    return (p[0] * EXPORT_SCALE, p[1] * EXPORT_SCALE, p[2] * EXPORT_SCALE)


if CLEARANCE <= 0.0:
    fail("the tube footprint does not fit in its cell : the z-section is an "
         "ellipse %.3f x %.3f, the staggered cell is %.3f x %.3f "
         "(transverse pitch %.3f, diagonal %.3f)"
         % (2 * A_ELL, 2 * B_ELL, 2 * XV, 2 * HH, ST, S_D))
if abs(HELIX_ANGLE) >= 80.0:
    fail("HELIX_ANGLE must stay well below 90 deg")
SHIFT = H * TANA                                  # how far the patch leans over


R_IN   = R_MID - 0.5 * W
R_OUT  = R_MID + 0.5 * W
if WRAP_TO_COIL:
    if R_IN <= 0.0:
        fail("R_MID must exceed W/2 = %.3f so the inner shroud has a radius" % (0.5 * W))
    if THETA >= 2.0 * math.pi:
        fail("H/R_MID = %.3f rad wraps more than a full turn" % THETA)


def shear(x, y, z):
    """Lean the whole patch downstream with height.

    Every z-slice is the same section pushed along the shell axis by
    s = z*tan(ALPHA), so the inlet and outlet boxes keep their thickness and
    their outer faces stay parallel to the tubes.  The map has unit
    determinant, so it does not change any volume."""
    return (x + z * TANA, y, z)


def warp(x, y, z):
    """Flat patch coordinates -> final geometry.

    Shear first (tube inclination), then optionally roll the patch onto the
    shell: y becomes radius, z becomes angle, x stays the shell axis."""
    x, y, z = shear(x, y, z)
    if not WRAP_TO_COIL:
        return (x, y, z)
    #  y -> radius, z -> angle.  Order chosen so the Jacobian is +r/R_MID:
    #  the other pairing mirrors the patch and inverts every cell.
    th = z / R_MID
    r = R_IN + y
    return (x, r * math.cos(th), r * math.sin(th))


# =============================================================================
#  1) LATTICE AND CELLS
# =============================================================================
def lattice():
    """Staggered arrangement: consecutive rows are offset by half a pitch.

    A full row carries N_COLS+1 rods, the outermost two sitting on the walls and
    halved by them; an offset row carries N_COLS whole rods at mid-pitch."""
    cells = []
    for i in range(N_ROWS):
        x = X_B0 + 0.5 * SL + i * SL
        offset = (i % 2 == 0) if FIRST_ROW_OFFSET else (i % 2 == 1)
        if offset:
            for j in range(N_COLS):           # y = ST/2 .. W-ST/2, walls stay clear
                cells.append((x, (j + 0.5) * ST, True, i))
        else:
            for j in range(N_COLS + 1):       # y = 0 .. W, the ends on the walls
                wall = (j == 0 or j == N_COLS)
                cells.append((x, j * ST, HALF_RODS_AT_WALL or not wall, i))
    return cells


def clip_y(poly, c, keep_above):
    out = []
    n = len(poly)
    inside = (lambda p: p[1] >= c - TOL) if keep_above else (lambda p: p[1] <= c + TOL)
    for k in range(n):
        a, b = poly[k], poly[(k + 1) % n]
        ia, ib = inside(a), inside(b)
        if ia:
            out.append(a)
        if ia != ib:
            t = (c - a[1]) / (b[1] - a[1])
            out.append((a[0] + (b[0] - a[0]) * t, c))
    return out


def dedupe(poly):
    out = []
    for p in poly:
        if not out or abs(p[0] - out[-1][0]) > TOL or abs(p[1] - out[-1][1]) > TOL:
            out.append(p)
    if len(out) > 1 and abs(out[0][0] - out[-1][0]) < TOL and abs(out[0][1] - out[-1][1]) < TOL:
        out.pop()
    return out


def cell_polygon(cx, cy, i_row):
    """CCW Voronoi cell of a staggered lattice: the hexagon cut by the six
    nearest neighbours.

    The first and last rows are squared off against the A|B and B|C planes: the
    hexagon's pointed corner would otherwise stick through an interface that has
    to stay a flat plane for the box blocks to butt against."""
    pts = []
    if i_row == N_ROWS - 1:
        pts += [(X_B1, cy - HH), (X_B1, cy), (X_B1, cy + HH)]
    else:
        pts += [(cx + XE, cy - HH), (cx + XV, cy), (cx + XE, cy + HH)]
    if i_row == 0:
        pts += [(X_B0, cy + HH), (X_B0, cy), (X_B0, cy - HH)]
    else:
        pts += [(cx - XE, cy + HH), (cx - XV, cy), (cx - XE, cy - HH)]
    return dedupe(clip_y(clip_y(pts, 0.0, True), W, False))


CELLS = lattice()
#  An empty (rod-free) wall cell is one transfinite block, and a block's
#  opposite sides must carry equal counts - which a four-different-lengths
#  quadrilateral cannot give.  Its presence therefore forces a fixed count per
#  side for the whole mesh; see side_div.
UNIFORM_DIV = any(not rod for (_cx, _cy, rod, _i) in CELLS)
POLY = [(cx, cy, rod, cell_polygon(cx, cy, row)) for (cx, cy, rod, row) in CELLS]
N_RODS = sum(1 for c in CELLS if c[2])

# --- canonical vertex table --------------------------------------------------
#  Neighbouring cells compute a shared corner from different formulas, so the
#  corners are merged once here and every block uses the canonical coordinate.
VERT = []           # canonical (x, y)
_vgrid = {}         # coarse bucket -> [vertex ids]
_CELL = 4.0 * TOL   # bucket size, far below the smallest real spacing


def vertex_id(p):
    bx, by = int(math.floor(p[0] / _CELL)), int(math.floor(p[1] / _CELL))
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for vid in _vgrid.get((bx + dx, by + dy), ()):
                q = VERT[vid]
                if abs(q[0] - p[0]) < TOL and abs(q[1] - p[1]) < TOL:
                    return vid
    vid = len(VERT)
    VERT.append(p)
    _vgrid.setdefault((bx, by), []).append(vid)
    return vid


PCELL = []          # (cx, cy, has_rod, [vertex ids])
for (cx, cy, rod, pg) in POLY:
    PCELL.append((cx, cy, rod, [vertex_id(p) for p in pg]))


# =============================================================================
#  2) 1-D DISTRIBUTIONS
# =============================================================================
def growth_ratio(length, first, n):
    """q with first*(q^n-1)/(q-1) = length."""
    if n < 2 or first <= 0.0:
        return 1.0
    f = lambda q: first * n if abs(q - 1.0) < 1e-12 else first * (q ** n - 1.0) / (q - 1.0)
    lo, hi = 0.2, 5.0
    if f(lo) > length or f(hi) < length:
        return 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(mid) < length:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fractions(q, n):
    """n+1 normalised node positions 0..1 with a constant ratio q."""
    if abs(q - 1.0) < 1e-12:
        return [i / float(n) for i in range(n + 1)]
    out, s, step = [0.0], 0.0, (q - 1.0) / (q ** n - 1.0)
    for i in range(n):
        s += step
        out.append(s)
        step *= q
    out[-1] = 1.0
    return out


Q_RAD = growth_ratio(CLEARANCE, FIRST_LAYER, N_RAD)
FR_RAD = fractions(Q_RAD, N_RAD)
FR_A = fractions(X_SCALE_IN ** (1.0 / max(1, N_X_IN - 1)), N_X_IN)
FR_C = fractions(X_SCALE_OUT ** (1.0 / max(1, N_X_OUT - 1)), N_X_OUT)


def lerp(a, b, t):
    """Endpoint-exact linear interpolation."""
    if t == 0.0:
        return a
    if t == 1.0:
        return b
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)


def side_point(ia, ib, i, n):
    """Point i/n along the segment (ia,ib); identical from either side."""
    if ia <= ib:
        return lerp(VERT[ia], VERT[ib], i / float(n))
    return lerp(VERT[ib], VERT[ia], (n - i) / float(n))


def side_key(ia, ib, i, n):
    return ("side", ia, ib, i) if ia <= ib else ("side", ib, ia, n - i)


def side_div(ia, ib):
    """Divisions on a cell side, proportional to its length.

    A hexagon's flat and slanted sides differ by about a third, and the squared
    -off end rows differ again, so a fixed count per side would bunch nodes on
    part of the tube wall and starve the rest.  Scaling by length keeps the
    spacing even; both cells sharing a side see the same two endpoints, so they
    compute the same count and the mesh stays conformal.

    UNIFORM_DIV disables it when a rod-free wall cell exists - that cell is a
    single transfinite block, whose opposite sides must match, and its four
    sides all differ in length.  A fixed count per side is then the only
    conformal choice."""
    if UNIFORM_DIV:
        return N_AZ
    a, b = VERT[ia], VERT[ib]
    return max(1, int(round(math.hypot(b[0] - a[0], b[1] - a[1]) / AZ_STEP)))


# =============================================================================
#  3) NODES AND QUADS  (2-D section)
# =============================================================================
NODES = []          # (x, y)
NKEY = {}
QUADS = []          # (n0, n1, n2, n3) CCW
QZONE = []          # 0 = A, 1 = B, 2 = C
_maxdev = [0.0]
ROD_AREA = [0.0]    # area actually removed by the faceted tubes
ROD_MY   = [0.0]    # and its first moment about y = 0, for the wrapped volume
AZ_N     = {}       # (cx, cy) -> divisions actually laid around that tube


def node_id(key, xy):
    nid = NKEY.get(key)
    if nid is None:
        nid = len(NODES)
        NKEY[key] = nid
        NODES.append(xy)
        return nid
    q = NODES[nid]                       # same key from another block: must agree
    dev = max(abs(q[0] - xy[0]), abs(q[1] - xy[1]))
    if dev > _maxdev[0]:
        _maxdev[0] = dev
    return nid


def add_block(ni, nj, pos, key, zone):
    ids = [[node_id(key(i, j), pos(i, j)) for j in range(nj + 1)] for i in range(ni + 1)]
    for i in range(ni):
        for j in range(nj):
            a, b = ids[i][j], ids[i + 1][j]
            c, d = ids[i + 1][j + 1], ids[i][j + 1]
            ax, ay = NODES[a]
            area = ((NODES[b][0] - ax) * (NODES[c][1] - ay) -
                    (NODES[b][1] - ay) * (NODES[c][0] - ax))
            QUADS.append((a, b, c, d) if area > 0 else (a, d, c, b))
            QZONE.append(zone)
    return ids


# --- B : one O-grid block per rectangle side --------------------------------
def build_bundle():
    for (cx, cy, rod, vids) in PCELL:
        n = len(vids)
        if n < 3:
            continue
        cid = (round(cx / TOL), round(cy / TOL))

        if not rod:                                    # empty wall cell : one block
            if n != 4:
                fail("empty cell with %d corners" % n)
            v0, v1, v2, v3 = vids

            def pos(i, j, v0=v0, v1=v1, v2=v2, v3=v3):
                return lerp(side_point(v0, v1, i, N_AZ),
                            side_point(v3, v2, i, N_AZ), j / float(N_AZ))

            def key(i, j, v0=v0, v1=v1, v2=v2, v3=v3, cid=cid):
                #  a corner belongs to several sides at once -> its own key
                if i == 0 and j == 0:
                    return ("vert", v0)
                if i == N_AZ and j == 0:
                    return ("vert", v1)
                if i == N_AZ and j == N_AZ:
                    return ("vert", v2)
                if i == 0 and j == N_AZ:
                    return ("vert", v3)
                if j == 0:
                    return side_key(v0, v1, i, N_AZ)
                if j == N_AZ:
                    return side_key(v3, v2, i, N_AZ)
                if i == 0:
                    return side_key(v0, v3, j, N_AZ)
                if i == N_AZ:
                    return side_key(v1, v2, j, N_AZ)
                return ("g", cid, i, j)

            add_block(N_AZ, N_AZ, pos, key, 1)
            continue

        # rod cell : skip the side that lies on the channel wall (inside the rod)
        for k in range(n):
            ia, ib = vids[k], vids[(k + 1) % n]
            A, B = VERT[ia], VERT[ib]
            on_wall = ((abs(cy) < TOL and abs(A[1]) < TOL and abs(B[1]) < TOL) or
                       (abs(cy - W) < TOL and abs(A[1] - W) < TOL and abs(B[1] - W) < TOL))
            if on_wall:
                continue

            def wall_pt(vid, cx=cx, cy=cy):
                """Point on the elliptical tube footprint towards a cell corner."""
                V = VERT[vid]
                th = math.atan2((V[1] - cy) / B_ELL, (V[0] - cx) / A_ELL)
                return (cx + A_ELL * math.cos(th), cy + B_ELL * math.sin(th))

            tA = math.atan2((A[1] - cy) / B_ELL, (A[0] - cx) / A_ELL)
            tB = math.atan2((B[1] - cy) / B_ELL, (B[0] - cx) / A_ELL)
            while tB <= tA:
                tB += 2.0 * math.pi
            ns = side_div(ia, ib)                    # this side's divisions
            AZ_N[(cx, cy)] = AZ_N.get((cx, cy), 0) + ns

            def pos(i, j, ia=ia, ib=ib, tA=tA, tB=tB, cx=cx, cy=cy, ns=ns):
                if i == 0:
                    S = wall_pt(ia)
                elif i == ns:
                    S = wall_pt(ib)
                else:
                    th = tA + (tB - tA) * (i / float(ns))
                    S = (cx + A_ELL * math.cos(th), cy + B_ELL * math.sin(th))
                return lerp(S, side_point(ia, ib, i, ns), FR_RAD[j])

            # the meshed rod is a polygon through the arc nodes, not a circle:
            # accumulate its exact area for the volume check further down
            dth = (tB - tA) / ns                     # fan of ns triangles
            at = 0.5 * A_ELL * B_ELL * math.sin(dth)
            for _q in range(ns):
                t0 = tA + _q * dth
                y0 = cy + B_ELL * math.sin(t0)
                y1 = cy + B_ELL * math.sin(t0 + dth)
                ROD_AREA[0] += at
                ROD_MY[0] += at * (cy + y0 + y1) / 3.0

            def key(i, j, ia=ia, ib=ib, cid=cid, k=k, ns=ns):
                if j == N_RAD:                      # on the cell side
                    if i == 0:                      # corner: shared by 3 sides
                        return ("vert", ia)
                    if i == ns:
                        return ("vert", ib)
                    return side_key(ia, ib, i, ns)
                if i == 0:                          # radial ray, shared with sector k-1
                    return ("ray", cid, ia, j)
                if i == ns:                         # radial ray, shared with sector k+1
                    return ("ray", cid, ib, j)
                if j == 0:
                    return ("arc", cid, k, i)
                return ("o", cid, k, i, j)

            add_block(ns, N_RAD, pos, key, 1)


# --- A and C : one block per interface segment ------------------------------
def interface_sides(x_plane):
    """The vertical cell sides lying on x = x_plane, as (va, vb) with va below."""
    out = []
    for (cx, cy, rod, vids) in PCELL:
        n = len(vids)
        for k in range(n):
            ia, ib = vids[k], vids[(k + 1) % n]
            A, B = VERT[ia], VERT[ib]
            if abs(A[0] - x_plane) < TOL and abs(B[0] - x_plane) < TOL:
                out.append((ia, ib) if A[1] < B[1] else (ib, ia))
    return out


def build_box(x_plane, length, nx, fr, zone, tag):
    if length <= TOL:
        return
    sgn = -1.0 if zone == 0 else 1.0            # A grows towards -x, C towards +x
    for (va, vb) in interface_sides(x_plane):
        ns = side_div(va, vb)

        def pos(i, k, va=va, vb=vb, ns=ns):
            P = side_point(va, vb, i, ns)
            return (x_plane + sgn * length * fr[k], P[1])

        def key(i, k, va=va, vb=vb, ns=ns):
            if k == 0:                              # on the A|B (B|C) interface
                if i == 0:
                    return ("vert", va)
                if i == ns:
                    return ("vert", vb)
                return side_key(va, vb, i, ns)
            if i == 0:                              # streamwise line, shared with the
                return (tag, va, k)                 # neighbouring box block
            if i == ns:
                return (tag, vb, k)
            if k == nx:
                return (tag + "e", va, vb, i)
            return (tag + "i", va, vb, i, k)

        add_block(ns, nx, pos, key, zone)


build_bundle()
build_box(X_B0, L_IN, N_X_IN, FR_A, 0, "a")
build_box(X_B1, L_OUT, N_X_OUT, FR_C, 2, "c")

#  rods on the walls are halved, so only the whole ones say what the azimuthal
#  resolution really is
_azf = [v for (k, v) in AZ_N.items() if abs(k[1]) > TOL and abs(k[1] - W) > TOL]
AZ_FULL = (min(_azf), max(_azf)) if _azf else (0, 0)
AZ_TXT = ("%d" % AZ_FULL[0]) if AZ_FULL[0] == AZ_FULL[1] else ("%d..%d" % AZ_FULL)

print("[2D] %d nodes, %d quads   (shared-node mismatch %.2e)"
      % (len(NODES), len(QUADS), _maxdev[0]))
if _maxdev[0] > TOL:
    fail("blocks disagree on a shared node -> topology keys are wrong")

# --- 2-D watertightness: an edge is shared by two quads unless it is boundary
_CEN = [(cx, cy) for (cx, cy, rod, _v) in PCELL if rod]


def on_boundary_2d(a, b):
    A, B = NODES[a], NODES[b]
    if abs(A[0]) < TOL and abs(B[0]) < TOL:
        return True
    if abs(A[0] - L_TOT) < TOL and abs(B[0] - L_TOT) < TOL:
        return True
    if abs(A[1]) < TOL and abs(B[1]) < TOL:
        return True
    if abs(A[1] - W) < TOL and abs(B[1] - W) < TOL:
        return True
    for (cx, cy) in _CEN:                       # on an elliptical tube footprint
        if (abs(((A[0] - cx) / A_ELL) ** 2 + ((A[1] - cy) / B_ELL) ** 2 - 1.0) < 1e-6 and
                abs(((B[0] - cx) / A_ELL) ** 2 + ((B[1] - cy) / B_ELL) ** 2 - 1.0) < 1e-6):
            return True
    return False


_ec = {}
for q in QUADS:
    for i in range(4):
        a, b = q[i], q[(i + 1) % 4]
        k = (a, b) if a < b else (b, a)
        _ec[k] = _ec.get(k, 0) + 1
_open = [k for k, v in _ec.items() if v == 1]
_crack = [k for k in _open if not on_boundary_2d(*k)]
print("      open edges %d (boundary %d, cracks %d)"
      % (len(_open), len(_open) - len(_crack), len(_crack)))
if _crack:
    for k in _crack[:8]:
        A, B = NODES[k[0]], NODES[k[1]]
        print("        crack (%.4f,%.4f)-(%.4f,%.4f)" % (A[0], A[1], B[0], B[1]))
    fail("the 2-D mesh is not watertight")


# =============================================================================
#  4) EXTRUDE TO HEXAHEDRA
# =============================================================================
NN2 = len(NODES)
ZS = [H * k / float(N_Z) for k in range(N_Z + 1)]
POINTS = [warp(x, y, z) for z in ZS for (x, y) in NODES]


def n3(nid, layer):
    return layer * NN2 + nid


HEXES = []
for (a, b, c, d) in QUADS:
    for l in range(N_Z):
        HEXES.append((n3(a, l), n3(b, l), n3(c, l), n3(d, l),
                      n3(a, l + 1), n3(b, l + 1), n3(c, l + 1), n3(d, l + 1)))

# outward faces of a hexahedron in the node order above
HEX_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))

FACEMAP = {}
for ci, hx in enumerate(HEXES):
    for lf in HEX_FACES:
        f = (hx[lf[0]], hx[lf[1]], hx[lf[2]], hx[lf[3]])
        k = tuple(sorted(f))
        e = FACEMAP.get(k)
        if e is None:
            FACEMAP[k] = [f, ci, -1]
        elif e[2] == -1:
            e[2] = ci
        else:
            fail("a face is shared by more than two cells")

INTERNAL = [e for e in FACEMAP.values() if e[2] >= 0]
BOUNDARY = [e for e in FACEMAP.values() if e[2] < 0]
print("[3D] %d nodes, %d hexahedra, %d faces (%d internal, %d boundary)"
      % (len(POINTS), len(HEXES), len(FACEMAP), len(INTERNAL), len(BOUNDARY)))


# =============================================================================
#  5) BOUNDARY PATCHES
# =============================================================================
def flat(i):
    """Patch coordinates of a node, before shear and wrap."""
    n = NODES[i % NN2]
    return (n[0], n[1], ZS[i // NN2])


def classify(face):
    #  classify in patch coordinates: the walls are planes there whatever the
    #  wrap does to them afterwards
    xs = [flat(i)[0] for i in face]
    ys = [flat(i)[1] for i in face]
    zs = [flat(i)[2] for i in face]
    if max(xs) - min(xs) < TOL:
        if abs(xs[0]) < TOL:
            return "inlet"
        if abs(xs[0] - L_TOT) < TOL:
            return "outlet"
    if max(ys) - min(ys) < TOL:
        if abs(ys[0]) < TOL:
            return "wall_side_y0"
        if abs(ys[0] - W) < TOL:
            return "wall_side_ymax"
    if max(zs) - min(zs) < TOL:
        if abs(zs[0]) < TOL:
            return "wall_bottom"
        if abs(zs[0] - H) < TOL:
            return "wall_top"
    return "wall_rods"


PATCH_ORDER = ["inlet", "outlet", "wall_rods", "wall_side_y0", "wall_side_ymax",
               "wall_bottom", "wall_top"]
PATCHES = dict((p, []) for p in PATCH_ORDER)
for e in BOUNDARY:
    PATCHES[classify(e[0])].append(e)

# every "rod" face must really sit on a rod wall.  Its four corners are exactly
# on the circle (the face centre is not, so testing the centre would fail at a
# coarse azimuthal resolution).
CENTERS = [(cx, cy) for (cx, cy, rod, _v) in PCELL if rod]
bad = 0
for e in PATCHES["wall_rods"]:
    #  test the unwarped 2-D node against the ellipse, not the warped point
    q = [NODES[i % NN2] for i in e[0]]
    if not any(all(abs(((n[0] - cx) / A_ELL) ** 2 + ((n[1] - cy) / B_ELL) ** 2 - 1.0)
                   < 1e-6 for n in q) for (cx, cy) in CENTERS):
        bad += 1
if bad:
    fail("%d faces classified as rod wall are not on any rod" % bad)

for p in PATCH_ORDER:
    print("      %-15s %7d faces" % (p, len(PATCHES[p])))


# =============================================================================
#  6) CHECKS  (volume, watertightness, quality)
# =============================================================================
def hex_volume(h):
    """Volume by the divergence theorem over the six outward faces.

    Each face is fanned about its own centroid rather than split on a diagonal.
    Once the patch is rolled onto the shell some faces are no longer planar, and
    only the centroid fan makes two neighbouring cells agree on the shared face -
    which is also how OpenFOAM and Fluent decompose such a face."""
    v = 0.0
    for lf in HEX_FACES:
        p = [POINTS[h[i]] for i in lf]
        c = (sum(q[0] for q in p) / 4.0,
             sum(q[1] for q in p) / 4.0,
             sum(q[2] for q in p) / 4.0)
        for k in range(4):
            a, b = p[k], p[(k + 1) % 4]
            ux, uy, uz = a[0] - c[0], a[1] - c[1], a[2] - c[2]
            wx, wy, wz = b[0] - c[0], b[1] - c[1], b[2] - c[2]
            nx = uy * wz - uz * wy
            ny = uz * wx - ux * wz
            nz = ux * wy - uy * wx
            v += (c[0] * nx + c[1] * ny + c[2] * nz) / 6.0
    return v


VOL = [hex_volume(h) for h in HEXES]
vmin, vmax, vtot = min(VOL), max(VOL), sum(VOL)
if vmin <= 0.0:
    fail("%d cells have non-positive volume" % sum(1 for v in VOL if v <= 0.0))

#  expected volume of the *meshed* domain: the rods are regular polygons through
#  the arc nodes, so their area is known exactly (not pi*R^2)
if WRAP_TO_COIL:
    #  the roll is not volume preserving: dV = (r/R_MID) dx dy dz
    a_fluid = L_TOT * W - ROD_AREA[0]
    m_fluid = L_TOT * W * W / 2.0 - ROD_MY[0]          # first moment about y = 0
    exact = (H / R_MID) * (R_IN * a_fluid + m_fluid)
    dth = THETA / N_Z                                  # each layer is a chord, not
    expect = exact * (math.sin(dth) / dth)             # an arc: exact deficit
else:
    expect = (L_TOT * W - ROD_AREA[0]) * H
circle = (L_TOT * W - sum(math.pi * A_ELL * B_ELL
                          * (0.5 if (cy < TOL or cy > W - TOL) else 1.0)
                          for (cx, cy) in CENTERS)) * H


def quality(h):
    """In-plane corner skew, plus the lean of the extrusion edge off vertical."""
    dev, lmin, lmax = 0.0, 1e30, 0.0
    e = (POINTS[h[4]][0] - POINTS[h[0]][0], POINTS[h[4]][1] - POINTS[h[0]][1],
         POINTS[h[4]][2] - POINTS[h[0]][2])
    dev = math.degrees(math.atan2(math.hypot(e[0], e[1]), abs(e[2]) or 1e-30))
    for face in ((0, 1, 2, 3), (4, 5, 6, 7)):
        p = [POINTS[h[i]] for i in face]
        for i in range(4):
            a, b, c = p[(i + 3) % 4], p[i], p[(i + 1) % 4]
            v1 = (a[0] - b[0], a[1] - b[1], a[2] - b[2])     # full 3-D: the face
            v2 = (c[0] - b[0], c[1] - b[1], c[2] - b[2])     # is tilted once wrapped
            n1 = math.sqrt(sum(t * t for t in v1)) or 1e-12
            n2 = math.sqrt(sum(t * t for t in v2)) or 1e-12
            cs = max(-1.0, min(1.0, sum(v1[k] * v2[k] for k in range(3)) / (n1 * n2)))
            dev = max(dev, abs(math.degrees(math.acos(cs)) - 90.0))
            lmin, lmax = min(lmin, n2), max(lmax, n2)
    dz = H / float(N_Z)
    return dev, max(lmax, dz) / min(lmin, dz)


QUAL = [quality(h) for h in HEXES]
dev_max = max(q[0] for q in QUAL)
ar_max = max(q[1] for q in QUAL)
bins = [0, 0, 0, 0, 0]
for (dv, _a) in QUAL:
    bins[0 if dv < 15 else 1 if dv <= 30 else 2 if dv < 40 else 3 if dv < 50 else 4] += 1

#  Rolled or flat, the expected volume is known in closed form.  Rolled, each
#  layer chords its arc, and the deficit is exactly sin(dtheta)/dtheta - verified
#  to 6e-15 here and to fall off as dtheta^2 when the layers are refined.
if abs(vtot - expect) > 1e-9 * expect:
    fail("volume %.6f != expected %.6f" % (vtot, expect))
print("[check] volume  %.4f = expected %.4f (err %.1e); cells %.3e .. %.3e"
      % (vtot, expect, abs(vtot - expect) / expect, vmin, vmax))
if WRAP_TO_COIL:
    print("        exact annulus %.4f, chorded by %.3e (= dtheta^2/6)"
          % (exact, (exact - vtot) / exact))
    if math.degrees(dth) > 4.0:
        print("[note] each layer spans %.2f deg; the sector is chorded to %.2f %%."
              " Raise N_Z to sharpen the arc." % (math.degrees(dth),
                                                 100.0 * (exact - vtot) / exact))
if WRAP_TO_COIL:
    print("        rolled onto the shell: sector %.3f rad, r %.1f .. %.1f, "
          "z-facet factor %.9f" % (THETA, R_IN, R_OUT, math.sin(dth) / dth))
print("        faceted tubes (%s-gon) hold %.3f %% more fluid than true ellipses"
      % (AZ_TXT, 100.0 * (vtot - circle) / circle))
print("[check] skew    max %.1f deg   <15 %.1f%% | 15-30 %.1f%% | 30-40 %.1f%% | "
      "40-50 %.1f%% | >50 %.1f%%"
      % ((dev_max,) + tuple(100.0 * b / len(QUAL) for b in bins)))
print("        %d cells sit on a cell corner (~%.0f deg); a hexagon corner split by"
      % (sum(1 for q in QUAL if abs(q[0] - dev_max) < 0.05), dev_max))
print("        one radial line leaves ~60 deg, but the squared-off first and last rows")
print("        keep 90 deg corners, and those four cells per rod are the design floor")
print("[check] aspect  max %.1f  = dz %.3f / first layer %.4f"
      % (ar_max, H / N_Z, CLEARANCE * FR_RAD[1]))


# =============================================================================
#  7) WRITERS
# =============================================================================
def write_vtu(path):
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" '
                'version="0.1" byte_order="LittleEndian">\n <UnstructuredGrid>\n')
        f.write('  <Piece NumberOfPoints="%d" NumberOfCells="%d">\n'
                % (len(POINTS), len(HEXES)))
        f.write('   <Points>\n    <DataArray type="Float64" NumberOfComponents="3" '
                'format="ascii">\n')
        f.write("".join("%.9g %.9g %.9g\n" % XP(p) for p in POINTS))
        f.write('    </DataArray>\n   </Points>\n   <Cells>\n')
        f.write('    <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        f.write("".join("%d %d %d %d %d %d %d %d\n" % h for h in HEXES))
        f.write('    </DataArray>\n    <DataArray type="Int32" Name="offsets" '
                'format="ascii">\n')
        f.write("".join("%d\n" % (8 * (i + 1)) for i in range(len(HEXES))))
        f.write('    </DataArray>\n    <DataArray type="UInt8" Name="types" '
                'format="ascii">\n')
        f.write("12\n" * len(HEXES))
        f.write('    </DataArray>\n   </Cells>\n   <CellData>\n')
        for name, vals, fmt in (("skew_deg", [q[0] for q in QUAL], "%.4g"),
                                ("aspect_ratio", [q[1] for q in QUAL], "%.4g"),
                                ("volume", VOL, "%.6g")):
            f.write('    <DataArray type="Float64" Name="%s" format="ascii">\n' % name)
            f.write("".join((fmt + "\n") % v for v in vals))
            f.write('    </DataArray>\n')
        f.write('   </CellData>\n  </Piece>\n </UnstructuredGrid>\n</VTKFile>\n')


def write_openfoam(folder):
    if not os.path.isdir(folder):
        os.makedirs(folder)

    def head(f, cls, obj):
        f.write("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
                "    class       %s;\n    location    \"constant/polyMesh\";\n"
                "    object      %s;\n}\n\n" % (cls, obj))

    # internal faces in upper-triangular order, then the patches
    ordered, owner, neighbour = [], [], []
    inter = []
    for e in INTERNAL:
        o, n = e[1], e[2]
        inter.append((o, n, e[0]) if o < n else (n, o, e[0][::-1]))
    inter.sort(key=lambda t: (t[0], t[1]))
    for (o, n, fc) in inter:
        ordered.append(fc)
        owner.append(o)
        neighbour.append(n)
    bstart = {}
    for p in PATCH_ORDER:
        bstart[p] = (len(ordered), len(PATCHES[p]))
        for e in PATCHES[p]:
            ordered.append(e[0])
            owner.append(e[1])

    with open(os.path.join(folder, "points"), "w") as f:
        head(f, "vectorField", "points")
        f.write("%d\n(\n" % len(POINTS))
        f.write("".join("(%.10g %.10g %.10g)\n" % XP(p) for p in POINTS))
        f.write(")\n")
    with open(os.path.join(folder, "faces"), "w") as f:
        head(f, "faceList", "faces")
        f.write("%d\n(\n" % len(ordered))
        f.write("".join("4(%d %d %d %d)\n" % tuple(fc) for fc in ordered))
        f.write(")\n")
    with open(os.path.join(folder, "owner"), "w") as f:
        head(f, "labelList", "owner")
        f.write("%d\n(\n" % len(owner))
        f.write("".join("%d\n" % o for o in owner))
        f.write(")\n")
    with open(os.path.join(folder, "neighbour"), "w") as f:
        head(f, "labelList", "neighbour")
        f.write("%d\n(\n" % len(neighbour))
        f.write("".join("%d\n" % n for n in neighbour))
        f.write(")\n")
    with open(os.path.join(folder, "boundary"), "w") as f:
        head(f, "polyBoundaryMesh", "boundary")
        f.write("%d\n(\n" % len(PATCH_ORDER))
        for p in PATCH_ORDER:
            s, n = bstart[p]
            f.write("    %s\n    {\n        type            %s;\n"
                    "        nFaces          %d;\n        startFace       %d;\n    }\n"
                    % (p, "wall" if p.startswith("wall") else "patch", n, s))
        f.write(")\n")


FLUENT_BC = {"inlet": 10, "outlet": 5, "wall_rods": 3, "wall_side_y0": 3,
             "wall_side_ymax": 3, "wall_bottom": 3, "wall_top": 3}


def write_fluent(path):
    with open(path, "w") as f:
        f.write('(0 "helical coil bundle, staggered : '
                'structured hexahedral mesh")\n')
        f.write("(0 \"dimension\")\n(2 3)\n\n")
        f.write('(0 "nodes")\n(10 (0 1 %x 0 3))\n' % len(POINTS))
        f.write("(10 (1 1 %x 1 3)(\n" % len(POINTS))
        f.write("".join("%.10g %.10g %.10g\n" % XP(p) for p in POINTS))
        f.write("))\n\n")
        f.write('(0 "cells")\n(12 (0 1 %x 0))\n' % len(HEXES))
        f.write("(12 (2 1 %x 1 4))\n\n" % len(HEXES))
        f.write('(0 "faces")\n(13 (0 1 %x 0))\n' % len(FACEMAP))

        zid, first = 3, 1

        def emit(items, bc, name, zid, first):
            last = first + len(items) - 1
            f.write("(13 (%x %x %x %x 4)(\n" % (zid, first, last, bc))
            for e in items:
                fc, c0, c1 = e[0], e[1] + 1, (e[2] + 1 if e[2] >= 0 else 0)
                if not FLUENT_NORMAL_C0_TO_C1:
                    fc = fc[::-1]
                    c0, c1 = c1, c0
                f.write("%x %x %x %x %x %x\n" % (fc[0] + 1, fc[1] + 1, fc[2] + 1,
                                                 fc[3] + 1, c0, c1))
            f.write("))\n")
            return last + 1, (zid, name, bc)

        zones = []
        first, z = emit(INTERNAL, 2, "interior", zid, first)
        zones.append(z)
        for p in PATCH_ORDER:
            if not PATCHES[p]:
                continue
            zid += 1
            first, z = emit(PATCHES[p], FLUENT_BC[p], p, zid, first)
            zones.append(z)
        f.write("\n(45 (2 fluid fluid)())\n")
        names = {2: "interior", 3: "wall", 5: "pressure-outlet", 10: "velocity-inlet"}
        for (zi, nm, bc) in zones:
            f.write("(45 (%d %s %s)())\n" % (zi, names[bc], nm))


def write_stl(path):
    tris = []
    for p in PATCH_ORDER:
        for e in PATCHES[p]:
            q = [XP(POINTS[i]) for i in e[0]]
            tris.append((q[0], q[1], q[2]))
            tris.append((q[0], q[2], q[3]))
    with open(path, "wb") as f:
        f.write(b"rod bundle boundary".ljust(80, b"\0"))
        f.write(struct.pack("<I", len(tris)))
        for (a, b, c) in tris:
            ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
            vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
            nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
            m = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
            f.write(struct.pack("<12fH", nx / m, ny / m, nz / m,
                                a[0], a[1], a[2], b[0], b[1], b[2],
                                c[0], c[1], c[2], 0))


def mb(path):
    return os.path.getsize(path) / 1048576.0


print("[write] model unit %s -> exported coordinates in metres (x %g)"
      % (UNIT, EXPORT_SCALE))
if WRITE_FLUENT:
    p = os.path.join(OUT_DIR, NAME + ".msh")
    write_fluent(p)
    print("      %-28s %6.1f MB   (Fluent / CFX-Pre)" % (os.path.basename(p), mb(p)))
if WRITE_OPENFOAM:
    d = os.path.join(OUT_DIR, NAME + "_foam")
    write_openfoam(d)
    print("      %-28s %6.1f MB   (-> case/constant/polyMesh)"
          % (os.path.basename(d) + "/", sum(mb(os.path.join(d, x))
                                            for x in os.listdir(d))))
if WRITE_VTU:
    p = os.path.join(OUT_DIR, NAME + ".vtu")
    write_vtu(p)
    print("      %-28s %6.1f MB   (ParaView, quality fields)"
          % (os.path.basename(p), mb(p)))
if WRITE_STL:
    p = os.path.join(OUT_DIR, NAME + ".stl")
    write_stl(p)
    print("      %-28s %6.1f MB   (boundary surface, figures)"
          % (os.path.basename(p), mb(p)))

print("-" * 68)
print(" helix      : %.1f deg -> tube footprint %.3f x %.3f, patch leans %.3f"
      % (HELIX_ANGLE, 2 * A_ELL, 2 * B_ELL, SHIFT))
if WRAP_TO_COIL:
    print(" coil       : r %.2f .. %.2f (mid %.2f), sector %.2f deg -> arc %.2f, "
          "helix %.2f .. %.2f deg across the layers"
          % (R_IN, R_OUT, R_MID, SECTOR_DEG, H,
             math.degrees(math.atan(TANA * R_MID / R_OUT)),
             math.degrees(math.atan(TANA * R_MID / R_IN))))
print(" lattice    : staggered (%s row offset), pitch %.3f x %.3f (SL x ST), "
      "diagonal %.3f"
      % ("odd" if FIRST_ROW_OFFSET else "even", SL, ST, S_D))
print("              cell %.3f x %.3f (hexagon), tube -> wall %.4f "
      "(transverse %.4f, diagonal %.4f)"
      % (2 * XV, 2 * HH, CLEARANCE, wall_gap(0.0, ST), wall_gap(SL, 0.5 * ST)))
print(" domain %.2f x %.2f x %.2f | rods %d | cells %d"
      % (L_TOT, W, H, N_RODS, len(HEXES)))
print(" around tube %s | radial %d (first %.4f, growth %.3f) | flow A %d / C %d"
      % (AZ_TXT, N_RAD, CLEARANCE * FR_RAD[1], Q_RAD, N_X_IN, N_X_OUT))
if UNIFORM_DIV:
    print("              no rods on the walls -> fixed %d divisions per cell side"
          % N_AZ)
print(" along the arc H : N_Z = %d layers over %.2f deg -> %.3f deg each, step %.3f"
      % (N_Z, SECTOR_DEG, SECTOR_DEG / N_Z, H / N_Z))
print("-" * 68)
