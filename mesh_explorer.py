#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bundle mesh generator : structured all-hexahedral meshes for four geometries
============================================================================
Standalone.  No SALOME, no CAD kernel, no third-party module - stdlib only.

This is the independent twin of ``mesh_explorer.html``.  Both build the same
mesh from the same closed-form rules by separate code, so comparing their node
and quad tables, patch sizes and total volume is a much stronger check than
either one passing its own tests.  Keep them in step: a fix here belongs there
too, and the ``--json`` output exists so a script can compare the two.

  flow +X     rod axis +Z     transverse +Y

        y
        ^     |<-- L_in -->|<---- L_bundle ---->|<-- L_out -->|
        |     +------------+--------------------+-------------+
        |     |            |  o     o     o     |             |
   W    |     |  A : inlet |  o     o     o     |  C : outlet |
        |     |            |  o     o     o     |             |
        |     +------------+--------------------+-------------+  --> x
                                  B : bundle

Four geometries, two independent traits:

  family "rod"      straight circular rods; the channel height H is given.
  family "helical"  tubes inclined by the helix angle, meshed flat and then
                    rolled back onto the coil.  The section footprint is an
                    ellipse 2*aE long by 2*bE wide with aE = R/cos(alpha), and
                    H is the arc length R_mid * theta rather than a free
                    parameter.  Setting alpha = 0 and skipping the roll
                    collapses every helical term back to the rod case, which is
                    why one code path serves both.

  in-line           every row carries rods at the same y; the Voronoi cell is
                    the SL x ST rectangle.
  staggered         consecutive rows are offset by half a transverse pitch; the
                    Voronoi cell is the hexagon cut by the six nearest
                    neighbours, squared off against the A|B and B|C planes.

Everything else is shared.  Geometry knowledge lives in exactly five places:
derived(), lattice(), cell_polygon(), side_div() and warp().

Every node position is closed form, so no boolean geometry and no meshing
algorithm is involved:

  B   : radial lines from each rod to the corners of its Voronoi cell cut the
        cell into one mapped block per side - an O-grid.  Radial spacing is
        geometric, so the wall gets an inflation layer.
  A,C : one mapped block per interface segment, graded towards the bundle.
  The 2-D section is then extruded N_Z layers -> 100 % hexahedra.

Blocks share nodes through *topological* keys (not coordinate rounding), so
shared points are bit-identical and the mesh is watertight by construction.

Output
  <name>.msh   Fluent/ANSYS ASCII  -> ANSYS Fluent and CFX-Pre
  <name>_foam/ OpenFOAM polyMesh   -> copy into <case>/constant/polyMesh
  <name>.vtu   VTK unstructured    -> ParaView (mesh + quality fields)
  <name>.stl   boundary surface    -> figures / documentation

Run:  python3 mesh_explorer.py --geometry rod-inline
      python3 mesh_explorer.py --geometry helical-staggered --helix 20 --sector 15
      python3 mesh_explorer.py --geometry rod-staggered --json --no-write
"""

import argparse
import json
import math
import os
import struct
import sys

TAU = 2.0 * math.pi

# =============================================================================
#  THE GEOMETRY REGISTRY
# =============================================================================
GEOMETRIES = {
    "rod-inline":        dict(family="rod",     stagger=False, prefix="rod_inline",
                              defaults=dict(SL=20.0)),
    "rod-staggered":     dict(family="rod",     stagger=True,  prefix="rod_stg",
                              defaults=dict(SL=17.32, first_offset=True)),
    "helical-inline":    dict(family="helical", stagger=False, prefix="helical_inline",
                              defaults=dict(SL=20.0)),
    "helical-staggered": dict(family="helical", stagger=True,  prefix="helical_stg",
                              defaults=dict(SL=20.0, first_offset=False)),
}

PATCH_ORDER = ["inlet", "outlet", "wall_rods", "wall_side_y0", "wall_side_ymax",
               "wall_bottom", "wall_top"]
FLUENT_BC = {"inlet": 10, "outlet": 5, "wall_rods": 3, "wall_side_y0": 3,
             "wall_side_ymax": 3, "wall_bottom": 3, "wall_top": 3}
FLUENT_TYPE = {2: "interior", 3: "wall", 5: "pressure-outlet", 10: "velocity-inlet"}
UNIT_SCALE = {"mm": 1.0e-3, "m": 1.0}

# outward faces of a hexahedron in the node order used below
HEX_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))


def fail(msg):
    print("[ERROR] " + msg)
    sys.exit(1)


# =============================================================================
#  CASE : parameters and the quantities derived from them
# =============================================================================
class Case(object):
    """One parameter set plus everything derived from it.

    Geometry function 1 of 5 lives here, in ``derived``.
    """

    def __init__(self, geometry, **kw):
        if geometry not in GEOMETRIES:
            fail("unknown geometry %r; pick one of %s"
                 % (geometry, ", ".join(sorted(GEOMETRIES))))
        self.geometry = geometry
        g = GEOMETRIES[geometry]
        self.family = g["family"]
        self.stagger = g["stagger"]
        self.prefix = g["prefix"]
        # defaults, then this geometry's overrides, then whatever was asked for
        p = dict(D=10.0, ST=20.0, SL=20.0, n_rows=5, n_cols=4,
                 l_in=60.0, l_out=120.0, H=40.0,
                 helix=15.0, r_mid=200.0, sector=12.0, wrap=True,
                 n_az=36, n_rad=12, first_layer=0.15, n_z=10,
                 n_x_in=24, n_x_out=40, x_scale=4.0,
                 half_rods=True, first_offset=False, unit="mm")
        p.update(g["defaults"])
        p.update({k: v for k, v in kw.items() if v is not None})
        self.__dict__.update(p)
        if self.unit not in UNIT_SCALE:
            fail('unit must be "mm" or "m"')
        self.export_scale = UNIT_SCALE[self.unit]
        self.derived()

    # -- geometry function 1 of 5 ------------------------------------------
    def derived(self):
        if self.family == "helical":
            self.H = math.radians(self.sector) * self.r_mid   # arc at the mid radius
        self.R = 0.5 * self.D
        self.W = self.n_cols * self.ST
        self.l_bund = self.n_rows * self.SL
        self.l_tot = self.l_in + self.l_bund + self.l_out
        self.x_b0 = self.l_in
        self.x_b1 = self.l_in + self.l_bund
        self.scale = max(self.l_tot, self.W, self.H)
        self.tol = 1.0e-9 * self.scale

        al = math.radians(self.helix) if self.family == "helical" else 0.0
        self.alpha = al
        self.tana = math.tan(al)
        # z-section of an inclined tube: an ellipse stretched along the flow
        self.aE = self.R / math.cos(al)
        self.bE = self.R
        self.shift = self.H * self.tana
        self.r_in = self.r_mid - 0.5 * self.W
        self.r_out = self.r_mid + 0.5 * self.W
        self.theta = math.radians(self.sector)
        self.wrap_ok = (self.family == "helical" and self.r_in > 0.0
                        and self.theta < TAU)
        self.rolled = bool(self.wrap and self.wrap_ok)

        self.HH = 0.5 * self.ST                        # cell half height
        if self.stagger:
            self.XV = 0.5 * self.SL + self.ST * self.ST / (8.0 * self.SL)
            self.XE = 0.5 * self.SL - self.ST * self.ST / (8.0 * self.SL)
            self.S_D = math.hypot(self.SL, 0.5 * self.ST)
            # Distance from the elliptical footprint to the Voronoi side facing
            # the neighbour at centre offset (vx, vy).  That side is the
            # perpendicular bisector - |v|/2 away, normal v/|v| - and an ellipse
            # reaches sqrt(A^2 nx^2 + B^2 ny^2) in direction n, so this is
            # exact.  For a circle it collapses to the familiar |v|/2 - R.  A
            # staggered cell has six sides and which one is tightest depends on
            # the pitch ratio AND on the inclination, so both are measured.
            def gap(vx, vy):
                dd = math.hypot(vx, vy)
                return 0.5 * dd - math.hypot(self.aE * vx, self.bE * vy) / dd
            self.gap_t = gap(0.0, self.ST)
            self.gap_d = gap(self.SL, 0.5 * self.ST)
            self.clearance = min(self.gap_t, self.gap_d)
            self.n_az_s = max(1, int(math.ceil(self.n_az / 6.0)))
            self.az_step = math.hypot(self.XV - self.XE, self.HH) / self.n_az_s
        else:
            self.HL = 0.5 * self.SL
            self.S_D = math.hypot(self.SL, self.ST)    # never binding in-line
            self.clearance = min(0.5 * self.ST - self.bE, 0.5 * self.SL - self.aE)
            # even, so a side halved by a wall splits exactly
            self.n_az_s = 2 * max(1, int(math.ceil(self.n_az / 8.0)))

    # -- geometry function 5 of 5 ------------------------------------------
    def lean_frac(self, x):
        """How much of the lean applies at streamwise position x: 1 across the
        bundle, ramped linearly to 0 at the inlet and outlet planes."""
        if x <= self.x_b0:
            f = x / self.l_in if self.l_in > 1e-12 else 1.0
        elif x >= self.x_b1:
            f = (self.l_tot - x) / self.l_out if self.l_out > 1e-12 else 1.0
        else:
            f = 1.0
        return 0.0 if f < 0.0 else (1.0 if f > 1.0 else f)

    def warp(self, p):
        """The map from the meshed patch into the world.

        Identity for a rod bundle.  For a coil, two steps.

        Step 1, the lean.  The tubes are inclined by alpha, so the bundle has
        to shear by z*tan(alpha).  The lean is ramped linearly back to zero
        across the inlet and outlet boxes, so B shears by the full amount while
        A stretches and C compresses by that same amount, and the inlet and
        outlet planes stay exactly perpendicular to the flow.  Leaning the
        whole patch instead would tilt the inlet plane by exactly alpha, and
        Fluent's default "Magnitude, Normal to Boundary" velocity inlet would
        then inject the flow alpha off the shell axis.  The cost is paid in the
        boxes, whose cells lean instead; total volume is unchanged because what
        A gains C loses.

        Step 2, the roll: y becomes radius, z becomes angle, x stays the shell
        axis.  Cos before sin keeps the Jacobian positive (+r/R); the other
        pairing mirrors the patch and inverts every cell.
        """
        x, y, z = p
        if self.tana:
            x += z * self.tana * self.lean_frac(x)
        if not self.rolled:
            return (x, y, z)
        th = z / self.r_mid
        r = self.r_in + y
        return (x, r * math.cos(th), r * math.sin(th))

    def XP(self, p):
        """A point on its way out of the program: warped, then always metres."""
        w = self.warp(p)
        s = self.export_scale
        return (w[0] * s, w[1] * s, w[2] * s)

    def name(self, az=None):
        """The case name the browser explorer writes too.  `az` is the total
        azimuthal count around a whole tube, which only the built mesh knows,
        so it is passed in; without it the per-side count stands in."""
        n = "%s_D%g_ST%g_SL%g" % (self.prefix, self.D, self.ST, self.SL)
        n += ("_a%g_s%g_R%g" % (self.helix, self.sector, self.r_mid)
              if self.family == "helical" else "_H%g" % self.H)
        n += "_%dx%d_az%dr%dz%d" % (self.n_rows, self.n_cols,
                                    az if az is not None else self.n_az_s,
                                    self.n_rad, self.n_z)
        if not self.half_rods:
            n += "_nohalf"
        if self.stagger and self.first_offset:
            n += "_off"
        if self.family == "helical" and not self.wrap:
            n += "_flat"
        return n.replace(".", "p")

    def validate(self):
        d_t, d_l = 2.0 * self.bE, 2.0 * self.aE
        if self.ST <= d_t:
            fail("rods in a row overlap: ST %.3f <= footprint %.3f" % (self.ST, d_t))
        if not self.stagger and self.SL <= d_l:
            fail("rods overlap along the flow: SL %.3f <= footprint %.3f"
                 % (self.SL, d_l))
        if self.stagger and self.XE <= 0.0:
            fail("the hexagonal cell folds: SL %.3f <= ST/2 %.3f"
                 % (self.SL, 0.5 * self.ST))
        if self.clearance <= 0.0:
            fail("the rod does not fit in its cell: wall gap %.4f" % self.clearance)
        if self.family == "helical" and self.wrap and not self.wrap_ok:
            fail("the coil radius is invalid: R - W/2 = %.3f must be positive and "
                 "the sector must stay under a full turn" % self.r_in)
        if self.family == "helical" and self.tana and (self.l_in <= 0 or self.l_out <= 0):
            print("[warn] %s is zero, so that face keeps the full %.3g deg helix lean "
                  "- there is no box to absorb it. Set the solver's inlet velocity by "
                  "components, not normal-to-boundary."
                  % ("l_in" if self.l_in <= 0 else "l_out", self.helix))


# =============================================================================
#  1) LATTICE AND CELLS   (geometry functions 2 and 3 of 5)
# =============================================================================
def lattice(c):
    """Where the rods sit.

    In-line   - every row carries rods at the same y positions.
    Staggered - consecutive rows are offset by half a transverse pitch.  A full
                row carries n_cols+1 rods, the outermost two sitting on the
                walls and halved by them; an offset row carries n_cols whole
                rods at mid-pitch.
    """
    cells = []
    for i in range(c.n_rows):
        x = c.x_b0 + 0.5 * c.SL + i * c.SL
        if c.stagger:
            off = (i % 2 == 0) if c.first_offset else (i % 2 == 1)
            if off:
                for j in range(c.n_cols):
                    cells.append((x, (j + 0.5) * c.ST, True, i))
            else:
                for j in range(c.n_cols + 1):
                    wall = (j == 0 or j == c.n_cols)
                    cells.append((x, j * c.ST, c.half_rods or not wall, i))
        elif c.half_rods:
            for j in range(c.n_cols + 1):     # y = 0 .. W, the ends halved by walls
                cells.append((x, j * c.ST, True, i))
        else:
            for j in range(c.n_cols):         # y = ST/2 .. W-ST/2, walls stay clear
                cells.append((x, (j + 0.5) * c.ST, True, i))
    return cells


def clip_y(poly, cut, keep_above, tol):
    out = []
    n = len(poly)
    inside = ((lambda p: p[1] >= cut - tol) if keep_above
              else (lambda p: p[1] <= cut + tol))
    for k in range(n):
        a, b = poly[k], poly[(k + 1) % n]
        ia, ib = inside(a), inside(b)
        if ia:
            out.append(a)
        if ia != ib:
            t = (cut - a[1]) / (b[1] - a[1])
            out.append((a[0] + (b[0] - a[0]) * t, cut))
    return out


def dedupe(poly, tol):
    out = []
    for p in poly:
        if not out or abs(p[0] - out[-1][0]) > tol or abs(p[1] - out[-1][1]) > tol:
            out.append(p)
    if (len(out) > 1 and abs(out[0][0] - out[-1][0]) < tol
            and abs(out[0][1] - out[-1][1]) < tol):
        out.pop()
    return out


def cell_polygon(c, cx, cy, row):
    """CCW Voronoi cell enclosing one rod.

    In-line   - the SL x ST rectangle.  No squaring-off is needed at the bundle
                ends: its left and right sides already sit on the A|B and B|C
                planes.
    Staggered - the hexagon cut by the six nearest neighbours.  The first and
                last rows ARE squared off; a hexagon corner would otherwise
                stick through an interface that has to stay flat for the box
                blocks to butt against it.
    """
    if c.stagger:
        pts = []
        if row == c.n_rows - 1:
            pts += [(c.x_b1, cy - c.HH), (c.x_b1, cy), (c.x_b1, cy + c.HH)]
        else:
            pts += [(cx + c.XE, cy - c.HH), (cx + c.XV, cy), (cx + c.XE, cy + c.HH)]
        if row == 0:
            pts += [(c.x_b0, cy + c.HH), (c.x_b0, cy), (c.x_b0, cy - c.HH)]
        else:
            pts += [(cx - c.XE, cy + c.HH), (cx - c.XV, cy), (cx - c.XE, cy - c.HH)]
    else:
        pts = [(cx + c.HL, cy - c.HH), (cx + c.HL, cy + c.HH),
               (cx - c.HL, cy + c.HH), (cx - c.HL, cy - c.HH)]
    return dedupe(clip_y(clip_y(pts, 0.0, True, c.tol), c.W, False, c.tol), c.tol)


# =============================================================================
#  2) 1-D DISTRIBUTIONS
# =============================================================================
def growth_ratio(length, first, n):
    """q with first*(q^n-1)/(q-1) = length."""
    if n < 2 or first <= 0.0:
        return 1.0
    def f(q):
        return first * n if abs(q - 1.0) < 1e-12 else first * (q ** n - 1.0) / (q - 1.0)
    lo, hi = 0.2, 5.0
    if f(lo) > length or f(hi) < length:
        return 1.0
    for _ in range(60):          # same count as the browser twin, to the last bit
        mid = 0.5 * (lo + hi)
        if f(mid) < length:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fractions(q, n):
    """n+1 normalised node positions 0..1 with a constant ratio q.

    The last entry is set to 1.0 outright: accumulating n terms designed to sum
    to 1 lands a few ulp either side, and this list is fed to lerp as t, where
    t = 0.9999999999999998 instead of 1 puts the outer radial node a hair off
    the cell side and the two O-grids sharing it disagree - a crack.  Same
    reasoning as lerp's own endpoint rule, one level up."""
    if abs(q - 1.0) < 1e-9:
        return [i / float(n) for i in range(n)] + [1.0]
    out, s, step = [0.0], 0.0, (q - 1.0) / (q ** n - 1.0)
    for _ in range(n - 1):
        s += step
        out.append(s)
        step *= q
    out.append(1.0)
    return out


def lerp(a, b, t):
    """Endpoint-exact linear interpolation.

    In floating point a + (b-a)*1.0 is not always b, and that inequality is
    exactly what becomes a crack between two blocks.
    """
    if t == 0.0:
        return a
    if t == 1.0:
        return b
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)


# =============================================================================
#  3) THE 2-D SECTION : nodes and quads
# =============================================================================
class Mesh(object):
    def __init__(self, case):
        self.c = case
        self.nodes = []          # (x, y) in the flat patch
        self.nkey = {}
        self.quads = []          # (n0, n1, n2, n3) CCW
        self.qzone = []          # 0 = A, 1 = B, 2 = C
        self.maxdev = 0.0
        self.rod_area = 0.0
        self.vert = []           # canonical cell corners
        self._vgrid = {}
        self.pcell = []          # (cx, cy, has_rod, [vertex ids])

    # -- canonical corner table -------------------------------------------
    def vertex_id(self, p):
        c = self.c
        cell = 4.0 * c.tol
        bx, by = int(math.floor(p[0] / cell)), int(math.floor(p[1] / cell))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for vid in self._vgrid.get((bx + dx, by + dy), ()):
                    q = self.vert[vid]
                    if abs(q[0] - p[0]) < c.tol and abs(q[1] - p[1]) < c.tol:
                        return vid
        vid = len(self.vert)
        self.vert.append(p)
        self._vgrid.setdefault((bx, by), []).append(vid)
        return vid

    def side_point(self, ia, ib, i, n):
        """Point i/n along the segment (ia,ib); identical from either side."""
        if ia <= ib:
            return lerp(self.vert[ia], self.vert[ib], i / float(n))
        return lerp(self.vert[ib], self.vert[ia], (n - i) / float(n))

    @staticmethod
    def side_key(ia, ib, i, n):
        return ("side", ia, ib, i) if ia <= ib else ("side", ib, ia, n - i)

    # -- geometry function 4 of 5 -----------------------------------------
    def side_div(self, ia, ib):
        """Divisions on a cell side, proportional to its length.

        Both cells sharing a side see the same two endpoints, so they compute
        the same count and the mesh stays conformal.

        In-line   - a side cut in half by a channel wall gets half the
                    divisions, so cells there keep the same size.
        Staggered - a hexagon's flat and slanted sides differ by about a third
                    and the squared-off end rows differ again, so a fixed count
                    per side would bunch nodes on part of the tube wall and
                    starve the rest.  A rod-free wall cell is ONE transfinite
                    block whose opposite sides must carry equal counts, which
                    four different lengths cannot give - so its presence forces
                    a fixed count everywhere.
        """
        c = self.c
        a, b = self.vert[ia], self.vert[ib]
        if c.stagger:
            if self.uniform_div:
                return c.n_az_s
            return max(1, int(round(math.hypot(b[0] - a[0], b[1] - a[1]) / c.az_step)))
        dx, dy = abs(b[0] - a[0]), abs(b[1] - a[1])
        full = c.SL if dx > dy else c.ST        # a full side of that orientation
        return max(1, int(round(c.n_az_s * math.hypot(dx, dy) / full)))

    # -- node / block assembly --------------------------------------------
    def node_id(self, key, xy):
        nid = self.nkey.get(key)
        if nid is None:
            nid = len(self.nodes)
            self.nkey[key] = nid
            self.nodes.append(xy)
            return nid
        q = self.nodes[nid]           # same key from another block: must agree
        dev = max(abs(q[0] - xy[0]), abs(q[1] - xy[1]))
        if dev > self.maxdev:
            self.maxdev = dev
        return nid

    def add_block(self, ni, nj, pos, key, zone):
        ids = [[self.node_id(key(i, j), pos(i, j)) for j in range(nj + 1)]
               for i in range(ni + 1)]
        for i in range(ni):
            for j in range(nj):
                a, b = ids[i][j], ids[i + 1][j]
                cc, d = ids[i + 1][j + 1], ids[i][j + 1]
                ax, ay = self.nodes[a]
                area = ((self.nodes[b][0] - ax) * (self.nodes[cc][1] - ay) -
                        (self.nodes[b][1] - ay) * (self.nodes[cc][0] - ax))
                self.quads.append((a, b, cc, d) if area > 0 else (a, d, cc, b))
                self.qzone.append(zone)
        return ids

    # -- B : one O-grid block per cell side --------------------------------
    def build_bundle(self):
        c = self.c
        fr_rad = self.fr_rad
        for (cx, cy, rod, vids) in self.pcell:
            n = len(vids)
            if n < 3:
                continue
            cid = (round(cx / c.tol), round(cy / c.tol))

            if not rod:                              # empty wall cell : one block
                if n != 4:
                    fail("empty cell with %d corners" % n)
                v0, v1, v2, v3 = vids
                na = c.n_az_s

                def pos(i, j, v0=v0, v1=v1, v2=v2, v3=v3, na=na):
                    return lerp(self.side_point(v0, v1, i, na),
                                self.side_point(v3, v2, i, na), j / float(na))

                def key(i, j, v0=v0, v1=v1, v2=v2, v3=v3, cid=cid, na=na):
                    #  a corner belongs to several sides at once -> its own key
                    if i == 0 and j == 0:
                        return ("vert", v0)
                    if i == na and j == 0:
                        return ("vert", v1)
                    if i == na and j == na:
                        return ("vert", v2)
                    if i == 0 and j == na:
                        return ("vert", v3)
                    if j == 0:
                        return self.side_key(v0, v1, i, na)
                    if j == na:
                        return self.side_key(v3, v2, i, na)
                    if i == 0:
                        return self.side_key(v0, v3, j, na)
                    if i == na:
                        return self.side_key(v1, v2, j, na)
                    return ("g", cid, i, j)

                self.add_block(na, na, pos, key, 1)
                continue

            # rod cell : skip the side lying on the channel wall (inside the rod)
            for k in range(n):
                ia, ib = vids[k], vids[(k + 1) % n]
                A, B = self.vert[ia], self.vert[ib]
                on_wall = ((abs(cy) < c.tol and abs(A[1]) < c.tol and abs(B[1]) < c.tol)
                           or (abs(cy - c.W) < c.tol and abs(A[1] - c.W) < c.tol
                               and abs(B[1] - c.W) < c.tol))
                if on_wall:
                    continue

                def wall_pt(vid, cx=cx, cy=cy):
                    """Point on the footprint; for a circle aE = bE and the
                    parametric angle is the polar angle."""
                    V = self.vert[vid]
                    th = math.atan2((V[1] - cy) / c.bE, (V[0] - cx) / c.aE)
                    return (cx + c.aE * math.cos(th), cy + c.bE * math.sin(th))

                tA = math.atan2((A[1] - cy) / c.bE, (A[0] - cx) / c.aE)
                tB = math.atan2((B[1] - cy) / c.bE, (B[0] - cx) / c.aE)
                while tB <= tA:
                    tB += TAU
                ns = self.side_div(ia, ib)

                def pos(i, j, ia=ia, ib=ib, tA=tA, tB=tB, cx=cx, cy=cy, ns=ns):
                    if i == 0:
                        S = wall_pt(ia)              # corner rays are shared with
                    elif i == ns:
                        S = wall_pt(ib)              # the neighbouring sector
                    else:
                        th = tA + (tB - tA) * (i / float(ns))
                        S = (cx + c.aE * math.cos(th), cy + c.bE * math.sin(th))
                    return lerp(S, self.side_point(ia, ib, i, ns), fr_rad[j])

                # the meshed footprint is a polygon through the arc nodes, not
                # an ellipse: accumulate its exact area for the volume check
                self.rod_area += 0.5 * c.aE * c.bE * ns * math.sin((tB - tA) / ns)

                def key(i, j, ia=ia, ib=ib, cid=cid, k=k, ns=ns):
                    if j == c.n_rad:                 # on the cell side
                        if i == 0:                   # corner: shared by 3 sides
                            return ("vert", ia)
                        if i == ns:
                            return ("vert", ib)
                        return self.side_key(ia, ib, i, ns)
                    if i == 0:                       # ray shared with sector k-1
                        return ("ray", cid, ia, j)
                    if i == ns:                      # ray shared with sector k+1
                        return ("ray", cid, ib, j)
                    if j == 0:
                        return ("arc", cid, k, i)
                    return ("o", cid, k, i, j)

                self.add_block(ns, c.n_rad, pos, key, 1)

    # -- A and C : one block per interface segment -------------------------
    def interface_sides(self, x_plane):
        """The vertical cell sides on x = x_plane, as (va, vb) with va below."""
        c = self.c
        out = []
        for (cx, cy, rod, vids) in self.pcell:
            n = len(vids)
            for k in range(n):
                ia, ib = vids[k], vids[(k + 1) % n]
                A, B = self.vert[ia], self.vert[ib]
                if abs(A[0] - x_plane) < c.tol and abs(B[0] - x_plane) < c.tol:
                    out.append((ia, ib) if A[1] < B[1] else (ib, ia))
        return out

    def build_box(self, x_plane, length, nx, fr, zone, tag):
        c = self.c
        if length <= c.tol:
            return
        sgn = -1.0 if zone == 0 else 1.0        # A grows towards -x, C towards +x
        for (va, vb) in self.interface_sides(x_plane):
            ns = self.side_div(va, vb)

            def pos(i, k, va=va, vb=vb, ns=ns):
                Pp = self.side_point(va, vb, i, ns)
                return (x_plane + sgn * length * fr[k], Pp[1])

            def key(i, k, va=va, vb=vb, ns=ns):
                if k == 0:                      # on the A|B (B|C) interface
                    if i == 0:
                        return ("vert", va)
                    if i == ns:
                        return ("vert", vb)
                    return self.side_key(va, vb, i, ns)
                if i == 0:                      # streamwise line, shared with the
                    return (tag, va, k)         # neighbouring box block
                if i == ns:
                    return (tag, vb, k)
                if k == nx:
                    return (tag + "e", va, vb, i)
                return (tag + "i", va, vb, i, k)

            self.add_block(ns, nx, pos, key, zone)

    # -- the whole build ---------------------------------------------------
    def build_section(self):
        """The 2-D part: cells, canonical corners, and the counts that follow
        from them.  Cheap, and enough to NAME the case - which is what anything
        predicting the output filename needs, without meshing the volume."""
        c = self.c
        cells = lattice(c)
        polys = [(cx, cy, rod, cell_polygon(c, cx, cy, row))
                 for (cx, cy, rod, row) in cells]
        self.n_rods = sum(1 for x in cells if x[2])
        for (cx, cy, rod, pg) in polys:
            self.pcell.append((cx, cy, rod, [self.vertex_id(p) for p in pg]))
        #  a rod-free wall cell forces a fixed count per side (see side_div)
        self.uniform_div = c.stagger and any(not x[2] for x in self.pcell)
        #  the azimuthal total depends only on the section, so it is known here
        self.az_full = self.whole_tube_divisions()

    def build(self):
        c = self.c
        self.build_section()
        self.q_rad = growth_ratio(c.clearance, c.first_layer, c.n_rad)
        self.fr_rad = fractions(self.q_rad, c.n_rad)
        fr_a = fractions(c.x_scale ** (1.0 / max(1, c.n_x_in - 1)), c.n_x_in)
        fr_c = fractions(c.x_scale ** (1.0 / max(1, c.n_x_out - 1)), c.n_x_out)

        self.build_bundle()
        self.build_box(c.x_b0, c.l_in, c.n_x_in, fr_a, 0, "a")
        self.build_box(c.x_b1, c.l_out, c.n_x_out, fr_c, 2, "c")

        if self.maxdev > c.tol:
            fail("blocks disagree on a shared node by %.2e -> topology keys are wrong"
                 % self.maxdev)
        self.check_watertight_2d()
        self.extrude()
        self.classify_patches()

    def whole_tube_divisions(self):
        """Divisions around a tube that the walls do NOT halve - what the case
        name and the explorer's hint both report."""
        c = self.c
        lo, hi = None, 0
        for (cx, cy, rod, vids) in self.pcell:
            if not rod or abs(cy) < c.tol or abs(cy - c.W) < c.tol:
                continue
            n = len(vids)
            tot = sum(self.side_div(vids[k], vids[(k + 1) % n]) for k in range(n))
            lo = tot if lo is None else min(lo, tot)
            hi = max(hi, tot)
        return (lo or 0, hi)

    # -- 2-D watertightness ------------------------------------------------
    def on_boundary_2d(self, a, b):
        c = self.c
        A, B = self.nodes[a], self.nodes[b]
        if abs(A[0]) < c.tol and abs(B[0]) < c.tol:
            return True
        if abs(A[0] - c.l_tot) < c.tol and abs(B[0] - c.l_tot) < c.tol:
            return True
        if abs(A[1]) < c.tol and abs(B[1]) < c.tol:
            return True
        if abs(A[1] - c.W) < c.tol and abs(B[1] - c.W) < c.tol:
            return True
        return self.on_rod(A) and self.on_rod(B)

    def on_rod(self, q):
        """Is this point exactly on some rod's footprint?

        Both ends of a wall edge lie on the ellipse; the chord midpoint does
        not, so testing the midpoint would fail at a coarse azimuthal count.
        """
        c = self.c
        for (cx, cy) in self.centres:
            if abs(((q[0] - cx) / c.aE) ** 2 + ((q[1] - cy) / c.bE) ** 2 - 1.0) < 1e-6:
                return True
        return False

    def check_watertight_2d(self):
        self.centres = [(cx, cy) for (cx, cy, rod, _v) in self.pcell if rod]
        ec = {}
        for q in self.quads:
            for i in range(4):
                a, b = q[i], q[(i + 1) % 4]
                k = (a, b) if a < b else (b, a)
                ec[k] = ec.get(k, 0) + 1
        self.open_edges = [k for k, v in ec.items() if v == 1]
        self.cracks = [k for k in self.open_edges if not self.on_boundary_2d(*k)]
        if self.cracks:
            for k in self.cracks[:8]:
                A, B = self.nodes[k[0]], self.nodes[k[1]]
                print("        crack (%.4f,%.4f)-(%.4f,%.4f)" % (A[0], A[1], B[0], B[1]))
            fail("the 2-D mesh is not watertight (%d cracks)" % len(self.cracks))

    # -- extrude -----------------------------------------------------------
    def extrude(self):
        c = self.c
        self.nn2 = len(self.nodes)
        zs = [c.H * k / float(c.n_z) for k in range(c.n_z + 1)]
        #  points stay in the FLAT patch here: the patches are classified on
        #  planes that only exist before the roll.  case.XP applies the map on
        #  the way out, so preview, checks and files share one path.
        self.points = [(x, y, z) for z in zs for (x, y) in self.nodes]
        self.hexes = []
        for (a, b, cc, d) in self.quads:
            for l in range(c.n_z):
                o, o2 = l * self.nn2, (l + 1) * self.nn2
                self.hexes.append((a + o, b + o, cc + o, d + o,
                                   a + o2, b + o2, cc + o2, d + o2))
        facemap = {}
        for ci, hx in enumerate(self.hexes):
            for lf in HEX_FACES:
                f = (hx[lf[0]], hx[lf[1]], hx[lf[2]], hx[lf[3]])
                k = tuple(sorted(f))
                e = facemap.get(k)
                if e is None:
                    facemap[k] = [f, ci, -1]
                elif e[2] == -1:
                    e[2] = ci
                else:
                    fail("a face is shared by more than two cells")
        self.facemap = facemap
        self.internal = [e for e in facemap.values() if e[2] >= 0]
        self.boundary = [e for e in facemap.values() if e[2] < 0]

    # -- boundary patches --------------------------------------------------
    def classify(self, face):
        c = self.c
        xs = [self.points[i][0] for i in face]
        ys = [self.points[i][1] for i in face]
        zs = [self.points[i][2] for i in face]
        if max(xs) - min(xs) < c.tol:
            if abs(xs[0]) < c.tol:
                return "inlet"
            if abs(xs[0] - c.l_tot) < c.tol:
                return "outlet"
        if max(ys) - min(ys) < c.tol:
            if abs(ys[0]) < c.tol:
                return "wall_side_y0"
            if abs(ys[0] - c.W) < c.tol:
                return "wall_side_ymax"
        if max(zs) - min(zs) < c.tol:
            if abs(zs[0]) < c.tol:
                return "wall_bottom"
            if abs(zs[0] - c.H) < c.tol:
                return "wall_top"
        return "wall_rods"

    def classify_patches(self):
        self.patches = dict((p, []) for p in PATCH_ORDER)
        for e in self.boundary:
            self.patches[self.classify(e[0])].append(e)
        # every "rod" face must really sit on a rod wall
        bad = sum(1 for e in self.patches["wall_rods"]
                  if not all(self.on_rod(self.points[i]) for i in e[0]))
        if bad:
            fail("%d faces classified as rod wall are not on any rod" % bad)


def quick_az(case):
    """Azimuthal divisions around a whole tube, without meshing the volume.

    The case name carries that count and the file on disk is called after it,
    so anything that has to predict the filename - the generated journal, for
    one - can get it without paying for the full build.
    """
    m = Mesh(case)
    m.build_section()
    return m.az_full


# =============================================================================
#  4) CHECKS
# =============================================================================
def hex_volume(m, h, split=False):
    """Volume by the divergence theorem over the six outward faces, on the
    EXPORTED points - so the check sees exactly what the solver will read.

    A quad face that is not planar has no single area vector, and the
    convention matters.  Default here is the CENTROID decomposition: fan the
    quad from the average of its four corners, which sums to exactly
    S = 1/2 (p2-p0) x (p3-p1).  That is what OpenFOAM does to a polygonal face,
    and it is what the browser explorer's cellVolumes() computes, so the two
    implementations are comparing the same number.

    split=True instead cuts each face on the 0-2 diagonal.  For a flat mesh
    every face is planar and the two agree to round-off; on a rolled coil they
    differ, and that difference is a useful measure of how non-planar the
    lateral faces have become.  It is reported, not gated.
    """
    v = 0.0
    for lf in HEX_FACES:
        p = [m.c.XP(m.points[h[i]]) for i in lf]
        if split:
            for t in ((0, 1, 2), (0, 2, 3)):
                a, b, cc = p[t[0]], p[t[1]], p[t[2]]
                nx = (b[1] - a[1]) * (cc[2] - a[2]) - (b[2] - a[2]) * (cc[1] - a[1])
                ny = (b[2] - a[2]) * (cc[0] - a[0]) - (b[0] - a[0]) * (cc[2] - a[2])
                nz = (b[0] - a[0]) * (cc[1] - a[1]) - (b[1] - a[1]) * (cc[0] - a[0])
                v += (a[0] * nx + a[1] * ny + a[2] * nz) / 6.0
            continue
        ax, ay, az = (p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2])
        bx, by, bz = (p[3][0] - p[1][0], p[3][1] - p[1][1], p[3][2] - p[1][2])
        sx = 0.5 * (ay * bz - az * by)
        sy = 0.5 * (az * bx - ax * bz)
        sz = 0.5 * (ax * by - ay * bx)
        cx = 0.25 * (p[0][0] + p[1][0] + p[2][0] + p[3][0])
        cy = 0.25 * (p[0][1] + p[1][1] + p[2][1] + p[3][1])
        cz = 0.25 * (p[0][2] + p[1][2] + p[2][2] + p[3][2])
        v += (cx * sx + cy * sy + cz * sz) / 3.0
    return v


def fnv1a(text):
    """32-bit FNV-1a.  Small, stdlib-only, and easy to reproduce in the
    browser, which is the point: an order-sensitive checksum over the node and
    quad tables proves the two implementations agree node for node, not merely
    on totals."""
    h = 0x811c9dc5
    for b in text.encode("utf-8"):
        h = ((h ^ b) * 0x01000193) & 0xFFFFFFFF
    return h


def table_hashes(m):
    nodes = ";".join("%.12e,%.12e" % (x, y) for (x, y) in m.nodes)
    quads = ";".join("%d,%d,%d,%d" % q for q in m.quads)
    return fnv1a(nodes), fnv1a(quads)


def analytic_volume(m):
    """What the cell volumes must add up to.

    Flat: the faceted section area times the height.

    Rolled: the roll has Jacobian r/r_mid, so by Pappus the smooth volume of a
    layer is the section's first moment about the coil axis times the swept
    angle.  The exported cells chord that arc instead of following it, and the
    chorded hexahedron's volume is exactly  int(r dA) * sin(dTheta)  -- not to
    leading order, exactly.  So the target is closed form either way and the
    gate stays at 1e-6 for both.  Area and first moment are taken from the
    exported quads, so the target is built from the same facets as the mesh.
    """
    c = m.c
    a2 = 0.0
    my6 = 0.0
    for q in m.quads:
        for e in range(4):
            a, b = m.nodes[q[e]], m.nodes[q[(e + 1) % 4]]
            cr = a[0] * b[1] - b[0] * a[1]
            a2 += cr
            my6 += (a[1] + b[1]) * cr
    area = abs(a2) / 2.0
    my = math.copysign(1.0, a2) * my6 / 6.0            # int(y dA)
    s3 = c.export_scale ** 3
    if not c.rolled:
        return area * c.H * s3, 0.0
    dth = c.theta / c.n_z
    mom = (c.r_in * area + my) * s3
    return mom * c.n_z * math.sin(dth), 1.0 - math.sin(dth) / dth


def quality(m, h):
    c = m.c
    dev, lmin, lmax = 0.0, 1e30, 0.0
    for face in ((0, 1, 2, 3), (4, 5, 6, 7)):
        p = [m.points[h[i]] for i in face]
        for i in range(4):
            a, b, cc = p[(i + 3) % 4], p[i], p[(i + 1) % 4]
            v1 = (a[0] - b[0], a[1] - b[1])
            v2 = (cc[0] - b[0], cc[1] - b[1])
            n1 = math.hypot(*v1) or 1e-12
            n2 = math.hypot(*v2) or 1e-12
            cs = max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1]) / (n1 * n2)))
            dev = max(dev, abs(math.degrees(math.acos(cs)) - 90.0))
            lmin, lmax = min(lmin, n2), max(lmax, n2)
    dz = c.H / float(c.n_z)
    return dev, max(lmax, dz) / min(lmin, dz)


def run_checks(m, quiet=False):
    c = m.c
    vol = [hex_volume(m, h) for h in m.hexes]
    vmin, vmax, vtot = min(vol), max(vol), sum(vol)
    #  the same volume under the other face convention; equal on a flat mesh,
    #  and a measure of face non-planarity once the patch is rolled
    vsplit = sum(hex_volume(m, h, split=True) for h in m.hexes)
    nonplanar = abs(vsplit - vtot) / max(abs(vtot), 1e-300)
    neg = sum(1 for v in vol if v <= 0.0)
    if neg:
        fail("%d cells have non-positive volume - the face orientation is inverted"
             % neg)
    expect, chord = analytic_volume(m)
    err = abs(vtot - expect) / expect
    if err > 1e-6:
        fail("the cell volumes do not sum to the analytic volume (error %.2e)" % err)
    qual = [quality(m, h) for h in m.hexes]
    dev_max = max(q[0] for q in qual)
    ar_max = max(q[1] for q in qual)
    nh, qh = table_hashes(m)
    #  a sparse sample of the node table itself, so a cross-check script can
    #  compare the numbers and not only their checksum
    step = max(1, len(m.nodes) // 400)
    sample = [[i, m.nodes[i][0], m.nodes[i][1]]
              for i in range(0, len(m.nodes), step)]
    res = dict(node_hash=nh, quad_hash=qh, nodes_sample=sample,
               nodes2d=len(m.nodes), quads=len(m.quads), points=len(m.points),
               cells=len(m.hexes), faces=len(m.facemap),
               internal=len(m.internal), boundary=len(m.boundary),
               cracks=len(m.cracks), maxdev=m.maxdev,
               vmin=vmin, vmax=vmax, vtot=vtot, expect=expect, vol_err=err,
               chord_deficit=chord, nonplanar=nonplanar,
               dev_max=dev_max, ar_max=ar_max,
               rods=m.n_rods, rod_area=m.rod_area, az_full=list(m.az_full),
               patches=dict((p, len(m.patches[p])) for p in PATCH_ORDER))
    m.vol, m.qual, m.report = vol, qual, res
    if quiet:
        return res
    print("[2D] %d nodes, %d quads   (shared-node mismatch %.2e)"
          % (len(m.nodes), len(m.quads), m.maxdev))
    print("      open edges %d (boundary %d, cracks %d)"
          % (len(m.open_edges), len(m.open_edges) - len(m.cracks), len(m.cracks)))
    print("[3D] %d nodes, %d hexahedra, %d faces (%d internal, %d boundary)"
          % (len(m.points), len(m.hexes), len(m.facemap),
             len(m.internal), len(m.boundary)))
    for p in PATCH_ORDER:
        print("      %-15s %7d faces" % (p, len(m.patches[p])))
    print("[check] volume  %.6e = expected %.6e (err %.1e); cells %.3e .. %.3e"
          % (vtot, expect, err, vmin, vmax))
    if chord:
        print("        rolled: each layer chords its arc, so the analytic target "
              "carries sin(dTheta)/dTheta = 1 - %.3e" % chord)
    print("        face convention: centroid vs diagonal split differ by %.1e "
          "(non-planarity of the lateral faces)" % nonplanar)
    print("[check] skew    max %.1f deg      aspect max %.1f" % (dev_max, ar_max))
    return res


# =============================================================================
#  5) WRITERS
# =============================================================================
def write_vtu(m, path):
    c = m.c
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" '
                'version="0.1" byte_order="LittleEndian">\n <UnstructuredGrid>\n')
        f.write('  <Piece NumberOfPoints="%d" NumberOfCells="%d">\n'
                % (len(m.points), len(m.hexes)))
        f.write('   <Points>\n    <DataArray type="Float64" NumberOfComponents="3" '
                'format="ascii">\n')
        f.write("".join("%.9g %.9g %.9g\n" % c.XP(p) for p in m.points))
        f.write('    </DataArray>\n   </Points>\n   <Cells>\n')
        f.write('    <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        f.write("".join("%d %d %d %d %d %d %d %d\n" % h for h in m.hexes))
        f.write('    </DataArray>\n    <DataArray type="Int32" Name="offsets" '
                'format="ascii">\n')
        f.write("".join("%d\n" % (8 * (i + 1)) for i in range(len(m.hexes))))
        f.write('    </DataArray>\n    <DataArray type="UInt8" Name="types" '
                'format="ascii">\n')
        f.write("12\n" * len(m.hexes))
        f.write('    </DataArray>\n   </Cells>\n   <CellData>\n')
        for name, vals, fmt in (("skew_deg", [q[0] for q in m.qual], "%.4g"),
                                ("aspect_ratio", [q[1] for q in m.qual], "%.4g"),
                                ("volume", m.vol, "%.6g")):
            f.write('    <DataArray type="Float64" Name="%s" format="ascii">\n' % name)
            f.write("".join((fmt + "\n") % v for v in vals))
            f.write('    </DataArray>\n')
        f.write('   </CellData>\n  </Piece>\n </UnstructuredGrid>\n</VTKFile>\n')


def write_openfoam(m, folder):
    c = m.c
    if not os.path.isdir(folder):
        os.makedirs(folder)

    def head(f, cls, obj):
        f.write("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
                "    class       %s;\n    location    \"constant/polyMesh\";\n"
                "    object      %s;\n}\n\n" % (cls, obj))

    # internal faces in upper-triangular order, then the patches
    ordered, owner, neighbour = [], [], []
    inter = []
    for e in m.internal:
        o, n = e[1], e[2]
        inter.append((o, n, e[0]) if o < n else (n, o, e[0][::-1]))
    inter.sort(key=lambda t: (t[0], t[1]))
    for (o, n, fc) in inter:
        ordered.append(fc)
        owner.append(o)
        neighbour.append(n)
    bstart = {}
    for p in PATCH_ORDER:
        bstart[p] = (len(ordered), len(m.patches[p]))
        for e in m.patches[p]:
            ordered.append(e[0])
            owner.append(e[1])

    with open(os.path.join(folder, "points"), "w") as f:
        head(f, "vectorField", "points")
        f.write("%d\n(\n" % len(m.points))
        f.write("".join("(%.10g %.10g %.10g)\n" % c.XP(p) for p in m.points))
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


def write_fluent(m, path, reverse=True):
    """Fluent face convention.  The face lists carry the OpenFOAM rule: the
    right-hand-rule normal of the node order runs from the owner into the
    neighbour, and outward on a boundary.  Fluent is the opposite - its normal
    points INTO c0 - so the node order is reversed on the way out.  Settled
    empirically: ANSYS Fluent 2026 R1 read a mesh written the OpenFOAM way and
    reported every cell as non-positive volume.  Pass reverse=False to write it
    the OpenFOAM way instead (which Fluent will reject)."""
    c = m.c
    with open(path, "w") as f:
        f.write('(0 "%s : structured hexahedral mesh")\n' % c.name(m.az_full[1]))
        f.write("(0 \"dimension\")\n(2 3)\n\n")
        f.write('(0 "nodes")\n(10 (0 1 %x 0 3))\n' % len(m.points))
        f.write("(10 (1 1 %x 1 3)(\n" % len(m.points))
        f.write("".join("%.10g %.10g %.10g\n" % c.XP(p) for p in m.points))
        f.write("))\n\n")
        f.write('(0 "cells")\n(12 (0 1 %x 0))\n' % len(m.hexes))
        f.write("(12 (2 1 %x 1 4))\n\n" % len(m.hexes))
        f.write('(0 "faces")\n(13 (0 1 %x 0))\n' % len(m.facemap))

        def emit(items, bc, name, zid, first):
            last = first + len(items) - 1
            f.write("(13 (%x %x %x %x 4)(\n" % (zid, first, last, bc))
            for e in items:
                fc, c0, c1 = e[0], e[1] + 1, (e[2] + 1 if e[2] >= 0 else 0)
                if reverse:
                    fc = fc[::-1]      # c0 stays the owner; only the normal flips
                f.write("%x %x %x %x %x %x\n" % (fc[0] + 1, fc[1] + 1, fc[2] + 1,
                                                 fc[3] + 1, c0, c1))
            f.write("))\n")
            return last + 1, (zid, name, bc)

        zones = []
        zid, first = 3, 1
        first, z = emit(m.internal, 2, "interior", zid, first)
        zones.append(z)
        for p in PATCH_ORDER:
            if not m.patches[p]:
                continue
            zid += 1
            first, z = emit(m.patches[p], FLUENT_BC[p], p, zid, first)
            zones.append(z)
        f.write("\n(45 (2 fluid fluid)())\n")
        for (zi, nm, bc) in zones:
            f.write("(45 (%d %s %s)())\n" % (zi, FLUENT_TYPE[bc], nm))


def write_stl(m, path):
    c = m.c
    tris = []
    for p in PATCH_ORDER:
        for e in m.patches[p]:
            q = [c.XP(m.points[i]) for i in e[0]]
            tris.append((q[0], q[1], q[2]))
            tris.append((q[0], q[2], q[3]))
    with open(path, "wb") as f:
        f.write(c.name(m.az_full[1]).encode("ascii", "replace")[:79].ljust(80, b"\0"))
        f.write(struct.pack("<I", len(tris)))
        for (a, b, cc) in tris:
            ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
            vx, vy, vz = cc[0] - a[0], cc[1] - a[1], cc[2] - a[2]
            nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
            mm = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
            f.write(struct.pack("<12fH", nx / mm, ny / mm, nz / mm,
                                a[0], a[1], a[2], b[0], b[1], b[2],
                                cc[0], cc[1], cc[2], 0))


def mb(path):
    return os.path.getsize(path) / 1048576.0


# =============================================================================
#  6) COMMAND LINE
# =============================================================================
def parse_args(argv):
    ap = argparse.ArgumentParser(
        description="Structured hexahedral mesh for a rod bundle or a helical "
                    "coil bundle, in-line or staggered.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--geometry", "-g", default="rod-inline",
                    choices=sorted(GEOMETRIES),
                    help="which bundle to build")
    gg = ap.add_argument_group("geometry (mm)")
    gg.add_argument("--D", type=float, help="rod / tube diameter")
    gg.add_argument("--ST", type=float, help="transverse pitch (y)")
    gg.add_argument("--SL", type=float, help="longitudinal pitch (x)")
    gg.add_argument("--n-rows", type=int, help="rod rows along the flow")
    gg.add_argument("--n-cols", type=int, help="full lateral pitches")
    gg.add_argument("--l-in", type=float, help="inlet box length (A)")
    gg.add_argument("--l-out", type=float, help="outlet box length (C)")
    gg.add_argument("--H", type=float,
                    help="channel height (rod family only; a coil derives it)")
    gg.add_argument("--half-rods", dest="half_rods", action="store_true", default=None,
                    help="rods on y=0 and y=W, halved by the walls")
    gg.add_argument("--no-half-rods", dest="half_rods", action="store_false",
                    help="rods at mid-pitch, walls stay clear")
    gg.add_argument("--first-offset", dest="first_offset", action="store_true",
                    default=None, help="staggered: make row 0 the offset row")
    gg.add_argument("--no-first-offset", dest="first_offset", action="store_false")
    hg = ap.add_argument_group("helix (helical family only)")
    hg.add_argument("--helix", type=float, help="helix angle alpha [deg]")
    hg.add_argument("--r-mid", type=float, help="coil mid radius [mm]")
    hg.add_argument("--sector", type=float, help="azimuthal sector theta [deg]")
    hg.add_argument("--no-wrap", dest="wrap", action="store_false", default=None,
                    help="keep the patch flat instead of rolling it onto the coil")
    mg = ap.add_argument_group("mesh")
    mg.add_argument("--n-az", type=int, help="minimum azimuthal cells around one rod")
    mg.add_argument("--n-rad", type=int, help="radial cells, rod wall -> cell boundary")
    mg.add_argument("--first-layer", type=float, help="first cell height on the wall")
    mg.add_argument("--n-z", type=int, help="layers in z (1 = quasi 2-D)")
    mg.add_argument("--n-x-in", type=int, help="streamwise cells in A")
    mg.add_argument("--n-x-out", type=int, help="streamwise cells in C")
    og = ap.add_argument_group("output")
    og.add_argument("--unit", choices=sorted(UNIT_SCALE),
                    help="what the numbers above mean; files are always metres")
    og.add_argument("--out-dir", default=None, help="where to write")
    og.add_argument("--name", default=None, help="base name (default: from the case)")
    og.add_argument("--formats", default="fluent,foam,vtu,stl",
                    help="comma-separated subset of fluent,foam,vtu,stl")
    og.add_argument("--no-write", action="store_true", help="check only, write nothing")
    og.add_argument("--json", action="store_true",
                    help="print the check numbers as JSON (for cross-checking "
                         "against the browser explorer)")
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    keys = ("D", "ST", "SL", "n_rows", "n_cols", "l_in", "l_out", "H",
            "helix", "r_mid", "sector", "wrap", "n_az", "n_rad", "first_layer",
            "n_z", "n_x_in", "n_x_out", "half_rods", "first_offset", "unit")
    case = Case(a.geometry, **dict((k, getattr(a, k, None)) for k in keys))
    case.validate()
    quiet = a.json
    if not quiet:
        print("=" * 68)
        print(" %s   (%s, %s)" % (case.name(), case.family,
                                  "staggered" if case.stagger else "in-line"))
        print("=" * 68)
    m = Mesh(case)
    m.build()
    rep = run_checks(m, quiet=quiet)

    if not a.no_write:
        out_dir = a.out_dir or os.path.dirname(os.path.abspath(__file__))
        name = a.name or case.name(m.az_full[1])
        want = set(x.strip() for x in a.formats.split(",") if x.strip())
        if not quiet:
            print("[write] model unit %s -> exported coordinates in metres (x %g)"
                  % (case.unit, case.export_scale))
        if "fluent" in want:
            p = os.path.join(out_dir, name + ".msh")
            write_fluent(m, p)
            if not quiet:
                print("      %-34s %6.1f MB   (Fluent / CFX-Pre)"
                      % (os.path.basename(p), mb(p)))
        if "foam" in want:
            d = os.path.join(out_dir, name + "_foam")
            write_openfoam(m, d)
            if not quiet:
                print("      %-34s %6.1f MB   (-> case/constant/polyMesh)"
                      % (os.path.basename(d) + "/",
                         sum(mb(os.path.join(d, x)) for x in os.listdir(d))))
        if "vtu" in want:
            p = os.path.join(out_dir, name + ".vtu")
            write_vtu(m, p)
            if not quiet:
                print("      %-34s %6.1f MB   (ParaView, quality fields)"
                      % (os.path.basename(p), mb(p)))
        if "stl" in want:
            p = os.path.join(out_dir, name + ".stl")
            write_stl(m, p)
            if not quiet:
                print("      %-34s %6.1f MB   (boundary surface, figures)"
                      % (os.path.basename(p), mb(p)))

    if a.json:
        print(json.dumps(rep, sort_keys=True))
        return 0
    print("-" * 68)
    print(" cell       : %s, clearance %.4f, footprint %.3f x %.3f"
          % ("hexagon" if case.stagger else "%.3f x %.3f rectangle"
             % (case.SL, case.ST), case.clearance, 2 * case.aE, 2 * case.bE))
    print(" domain %.2f x %.2f x %.2f | rods %d | cells %d"
          % (case.l_tot, case.W, case.H, m.n_rods, len(m.hexes)))
    print(" radial %d (first %.4f, growth %.3f) | z %d%s"
          % (case.n_rad, case.clearance * m.fr_rad[1], m.q_rad, case.n_z,
             "  | rolled onto R %.1f, %.2f deg/layer"
             % (case.r_mid, case.sector / case.n_z) if case.rolled else ""))
    print("-" * 68)
    return 0


if __name__ == "__main__":
    sys.exit(main())
