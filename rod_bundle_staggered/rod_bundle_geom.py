# -*- coding: utf-8 -*-
"""
SALOME GEOM + SMESH script : inlet box -> staggered rod bundle -> outlet box
                             fully structured (all-hexahedral) mesh

  flow direction : +X      rod axis : +Z      transverse : +Y

        y
        ^     |<-- L_in -->|<---- L_bundle ---->|<-- L_out -->|
        |     +------------+--------------------+-------------+
        |     |            |  o     o     o     |             |
   W    |     |  A : inlet |     o     o        |  C : outlet |
        |     |            |  o     o     o     |             |
        |     +------------+--------------------+-------------+  --> x
                                  B : bundle

MESH STRATEGY
  A, C : mapped (structured) quadrangles, graded towards the bundle
  B    : each rod is enclosed in its Voronoi hexagon of the staggered lattice
         -> the hexagon is split by radial lines into one quad block per side
         -> mapped quadrangles in every block  = O-grid around each rod
         -> radial distribution is graded, first cell at the rod wall (inflation)
         -> azimuthal cells per rod = 6 * N_AZ_SECT  (>= N_AZ_MIN)
  the whole 2-D section is meshed once (conformal everywhere) and then
  extruded N_Z layers in +Z  ->  100 % hexahedra.

Run:  salome -t rod_bundle_geom.py
 or : SALOME GUI > File > Load Script...
"""

import math
import salome

salome.salome_init()

import GEOM
from salome.geom import geomBuilder
import SMESH
from salome.smesh import smeshBuilder

geompy = geomBuilder.New()
smesh = smeshBuilder.New()

# =============================================================================
#  USER PARAMETERS
# =============================================================================
# --- geometry ----------------------------------------------------------------
D          = 10.0      # rod diameter                       [mm]
ST         = 20.0      # transverse pitch (y-direction)     [mm]
SL         = 17.32     # longitudinal pitch (x, row-to-row) [mm]  (ST*sqrt(3)/2 -> equilateral)

N_ROWS     = 5         # number of rod rows along the flow (x)
N_COLS     = 4         # number of full lateral pitches (y)

L_IN       = 60.0      # inlet  box length  (A, upstream)   [mm]
L_OUT      = 120.0     # outlet box length  (C, downstream) [mm]
H          = 40.0      # channel height (z), = rod length   [mm]

HALF_RODS_AT_WALL = True   # True : aligned rows carry half rods on both side walls
                           #        (classic periodic staggered bank, W = N_COLS*ST)
                           # False: no rod on the walls; the wall cells stay empty

FIRST_ROW_OFFSET  = True   # True : row 0 is the "centered" row (y = (j+0.5)*ST)

# --- mesh --------------------------------------------------------------------
BUILD_MESH  = True     # build the structured mesh (SMESH part)
BUILD_SOLID = True     # also build the 3-D solid + face groups (STEP export, CAD view)

N_AZ_MIN    = 36       # minimum azimuthal cells around one rod (-> 6 per hexagon side)
N_RAD       = 12       # radial cells between the rod wall and the hexagon boundary
FIRST_LAYER = 0.15     # first cell height on the rod wall (inflation)   [mm]

N_Z         = 10       # layers in z  (use 1 for a quasi-2D single-layer mesh)

N_X_IN      = 24       # streamwise cells in the inlet box  A
N_X_OUT     = 40       # streamwise cells in the outlet box C
X_SCALE_IN  = 4.0      # cell growth in A : last/first length ratio (fine at the bundle)
X_SCALE_OUT = 4.0      # cell growth in C : idem

# --- export ------------------------------------------------------------------
EXPORT_STEP = ""       # e.g. r"C:\ClaudeCode\rod_bundle\rod_bundle.step"  ("" = off)
EXPORT_MED  = ""       # e.g. r"C:\ClaudeCode\rod_bundle\rod_bundle.med"   ("" = off)
EXPORT_UNV  = ""       # e.g. r"C:\ClaudeCode\rod_bundle\rod_bundle.unv"   ("" = off)

# =============================================================================
#  DERIVED DIMENSIONS
# =============================================================================
R        = 0.5 * D                     # rod radius
W        = N_COLS * ST                 # channel width (y)
L_BUND   = N_ROWS * SL                 # bundle region length (x)
L_TOT    = L_IN + L_BUND + L_OUT       # total domain length
X_B0     = L_IN                        # bundle region start
X_B1     = L_IN + L_BUND               # bundle region end

HH       = 0.5 * ST                                  # cell half height
XV       = 0.5 * SL + ST * ST / (8.0 * SL)           # cell half length (vertex)
XE       = 0.5 * SL - ST * ST / (8.0 * SL)           # cell flat-edge half length
N_AZ_SECT = max(1, int(math.ceil(N_AZ_MIN / 6.0)))   # azimuthal cells per hexagon side

EXT      = 0.05 * H                    # rod z-overshoot (robust boolean cut)
TOL      = 1.0e-6 * max(L_TOT, W, H)   # geometric tolerance

O   = geompy.MakeVertex(0, 0, 0)
OZ  = geompy.MakeVectorDXDYDZ(0, 0, 1)
OX  = geompy.MakeVectorDXDYDZ(1, 0, 0)


# =============================================================================
#  1) LATTICE : rod centers and their meshing cells
# =============================================================================
def row_x(i):
    return X_B0 + 0.5 * SL + i * SL


def lattice():
    """[(x, y, has_rod), ...] : one entry per meshing cell of the bundle.

    Cells always tile the whole bundle region.  A cell without a rod (only
    possible on the side walls when HALF_RODS_AT_WALL is False) is meshed as a
    single mapped block."""
    cells = []
    for i in range(N_ROWS):
        x = row_x(i)
        offset_row = (i % 2 == 0) if FIRST_ROW_OFFSET else (i % 2 == 1)
        if offset_row:
            for j in range(N_COLS):
                cells.append((x, (j + 0.5) * ST, True))
        else:
            for j in range(N_COLS + 1):
                wall = (j == 0 or j == N_COLS)
                cells.append((x, j * ST, HALF_RODS_AT_WALL or not wall))
    return cells


CELLS   = lattice()
CENTERS = [(x, y) for (x, y, rod) in CELLS if rod]


def cell_polygon(cx, cy, i_row):
    """CCW polygon of the meshing cell around (cx, cy).

    Interior rows  -> Voronoi hexagon of the staggered lattice.
    First/last row -> the pointed end is squared off on the A|B (B|C) plane,
                      keeping 6 sides so that every side carries the same
                      number of azimuthal segments.
    Finally clipped to the side walls (half cells on y = 0 / y = W)."""
    pts = []
    if i_row == N_ROWS - 1:                       # right end -> flat on X_B1
        pts += [(X_B1, cy - HH), (X_B1, cy), (X_B1, cy + HH)]
    else:
        pts += [(cx + XE, cy - HH), (cx + XV, cy), (cx + XE, cy + HH)]
    if i_row == 0:                                # left end -> flat on X_B0
        pts += [(X_B0, cy + HH), (X_B0, cy), (X_B0, cy - HH)]
    else:
        pts += [(cx - XE, cy + HH), (cx - XV, cy), (cx - XE, cy - HH)]
    pts = clip_y(pts, 0.0, True)
    pts = clip_y(pts, W, False)
    return dedupe(pts)


def clip_y(poly, c, keep_above):
    """Sutherland-Hodgman clip of a polygon against the line y = c."""
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


# --- sanity checks -----------------------------------------------------------
#  the rod must stay strictly inside its cell; the three cell walls closest to
#  the rod centre are at ST/2 (flat side), S_D/2 (staggered neighbour) and
#  SL/2 (the squared-off end of the first/last row).
S_D       = math.hypot(SL, 0.5 * ST)                 # diagonal (staggered) pitch
CLEARANCE = min(0.5 * ST, 0.5 * S_D, 0.5 * SL) - R   # rod wall -> cell boundary

_problems = []
if SL <= 0.5 * ST:
    _problems.append("SL <= ST/2 : the staggered Voronoi cell degenerates")
if 0.5 * ST - R <= TOL:
    _problems.append("D >= ST : rods touch across the transverse pitch")
if 0.5 * S_D - R <= TOL:
    _problems.append("D >= S_D = %.3f : rods of neighbouring staggered rows touch" % S_D)
if 0.5 * SL - R <= TOL:
    _problems.append("D >= SL : the first/last row reaches the A|B interface plane")
for _p in _problems:
    print("[ERROR] " + _p)
if not _problems:
    print("[check] rod wall -> cell boundary clearance = %.4f (%.1f %% of D)"
          % (CLEARANCE, 100.0 * CLEARANCE / D))


# =============================================================================
#  2) 3-D SOLID  (CAD model, STEP export, face groups)
# =============================================================================
if BUILD_SOLID:
    channel = geompy.MakeBoxDXDYDZ(L_TOT, W, H)

    rods = [geompy.MakeCylinder(geompy.MakeVertex(xc, yc, -EXT), OZ, R, H + 2.0 * EXT)
            for (xc, yc) in CENTERS]
    rods_cpd = geompy.MakeCompound(rods)

    fluid = geompy.MakeCutList(channel, rods, True)

    trim = 2.0 * max(L_TOT, W, H)
    p_in = geompy.MakePlane(geompy.MakeVertex(X_B0, 0, 0), OX, trim)
    p_out = geompy.MakePlane(geompy.MakeVertex(X_B1, 0, 0), OX, trim)
    fluid = geompy.MakePartition([fluid], [p_in, p_out], [], [],
                                 geompy.ShapeType["SOLID"], 0, [], 0)
    geompy.addToStudy(fluid, "fluid_solid")

    def new_group(parent, stype, ids, name):
        if not ids:
            return None
        g = geompy.CreateGroup(parent, geompy.ShapeType[stype])
        geompy.UnionIDs(g, ids)
        geompy.addToStudyInFather(parent, g, name)
        return g

    def flat_at(bb, axis, value):
        lo, hi = bb[2 * axis], bb[2 * axis + 1]
        return abs(hi - lo) < TOL and abs(lo - value) < TOL

    ids = {"inlet": [], "outlet": [], "side_y0": [], "side_ymax": [],
           "bottom": [], "top": [], "rods": [], "interface": []}
    for f in geompy.SubShapeAll(fluid, geompy.ShapeType["FACE"]):
        fid = geompy.GetSubShapeID(fluid, f)
        if geompy.KindOfShape(f)[0] == geompy.kind.CYLINDER2D:
            ids["rods"].append(fid)
            continue
        bb = geompy.BoundingBox(f)
        if flat_at(bb, 0, 0.0):        ids["inlet"].append(fid)
        elif flat_at(bb, 0, L_TOT):    ids["outlet"].append(fid)
        elif flat_at(bb, 1, 0.0):      ids["side_y0"].append(fid)
        elif flat_at(bb, 1, W):        ids["side_ymax"].append(fid)
        elif flat_at(bb, 2, 0.0):      ids["bottom"].append(fid)
        elif flat_at(bb, 2, H):        ids["top"].append(fid)
        else:                          ids["interface"].append(fid)
    for key in ids:
        new_group(fluid, "FACE", ids[key], key)

    zone = {"zone_inlet": [], "zone_bundle": [], "zone_outlet": []}
    for s in geompy.SubShapeAll(fluid, geompy.ShapeType["SOLID"]):
        sid = geompy.GetSubShapeID(fluid, s)
        xc = geompy.PointCoordinates(geompy.MakeCDG(s))[0]
        key = "zone_inlet" if xc < X_B0 else ("zone_bundle" if xc < X_B1 else "zone_outlet")
        zone[key].append(sid)
    for k, v in zone.items():
        new_group(fluid, "SOLID", v, k)

    if EXPORT_STEP:
        geompy.ExportSTEP(fluid, EXPORT_STEP, GEOM.LU_MILLIMETER)
        print(" STEP exported : %s" % EXPORT_STEP)


# =============================================================================
#  3) 2-D SECTION : hexagonal cells + radial cuts  (the meshing geometry)
# =============================================================================
if BUILD_MESH:
    def V(x, y):
        return geompy.MakeVertex(x, y, 0.0)

    # --- outer rectangle minus the rod disks ---------------------------------
    rect = geompy.MakeFaceWires(
        [geompy.MakePolyline([V(0, 0), V(L_TOT, 0), V(L_TOT, W), V(0, W)], True)], True)
    disks = [geompy.MakeDiskPntVecR(V(cx, cy), OZ, R) for (cx, cy) in CENTERS]
    section = geompy.MakeCutList(rect, disks, True) if disks else rect

    # --- cutting tools -------------------------------------------------------
    tools = []
    seen_edges = set()

    def add_edge(p, q):
        """Add a cutting edge once (adjacent cells share their sides)."""
        a = (round(p[0] / TOL), round(p[1] / TOL))
        b = (round(q[0] / TOL), round(q[1] / TOL))
        key = (a, b) if a <= b else (b, a)
        if key in seen_edges:
            return
        seen_edges.add(key)
        tools.append(geompy.MakeEdge(V(p[0], p[1]), V(q[0], q[1])))

    POLY = []          # (cx, cy, has_rod, polygon) for every cell, in lattice order
    for i in range(N_ROWS):
        x = row_x(i)
        row_cells = [c for c in CELLS if abs(c[0] - x) < TOL]
        for (cx, cy, has_rod) in row_cells:
            pg = cell_polygon(cx, cy, i)
            POLY.append((cx, cy, has_rod, pg))
            for k in range(len(pg)):                       # cell sides
                add_edge(pg[k], pg[(k + 1) % len(pg)])
            if not has_rod:
                continue
            for (vx, vy) in pg:                            # radial cuts
                on_wall = (abs(cy) < TOL and abs(vy) < TOL) or \
                          (abs(cy - W) < TOL and abs(vy - W) < TOL)
                if on_wall:                                # already a boundary edge
                    continue
                # start inside the rod so the tool crosses the rod wall cleanly
                # (never starts on the rod centre, which lies on the side wall
                #  for the half rods)
                ln = math.hypot(vx - cx, vy - cy)
                add_edge((cx + 0.5 * R * (vx - cx) / ln,
                          cy + 0.5 * R * (vy - cy) / ln), (vx, vy))

    # --- A and C : split at every cell breakpoint so each block is a clean quad
    n_break = int(round(W / HH))
    for m in range(1, n_break):
        y = m * HH
        if L_IN > TOL:
            add_edge((0.0, y), (X_B0, y))
        if L_OUT > TOL:
            add_edge((X_B1, y), (L_TOT, y))

    sect_p = geompy.MakePartition([section], tools, [], [],
                                  geompy.ShapeType["FACE"], 0, [], 0)
    sect_entry = geompy.addToStudy(sect_p, "section_2D")

    n_faces_geom = len(geompy.SubShapeAll(sect_p, geompy.ShapeType["FACE"]))
    print("[section] %d cells (%d with a rod) -> %d partitioned faces"
          % (len(POLY), len(CENTERS), n_faces_geom))


# =============================================================================
#  4) EDGE CLASSIFICATION  (radial / axial-A / axial-C ; the rest is azimuthal)
# =============================================================================
if BUILD_MESH:
    EDGE = geompy.ShapeType["EDGE"]

    def pt(v):
        c = geompy.PointCoordinates(v)
        return (c[0], c[1])

    def on_rod(p):
        """True if p lies on a rod wall."""
        for (cx, cy) in CENTERS:
            if abs(math.hypot(p[0] - cx, p[1] - cy) - R) < 1.0e-4 * D:
                return True
        return False

    rad_edges, rad_rev = [], []
    ax_in, ax_in_rev = [], []
    ax_out, ax_out_rev = [], []
    n_arc = n_vert = n_side = 0

    for e in geompy.SubShapeAll(sect_p, EDGE):
        p0 = pt(geompy.GetVertexByIndex(e, 0))
        p1 = pt(geompy.GetVertexByIndex(e, 1))
        w0, w1 = on_rod(p0), on_rod(p1)
        if w0 and w1:                                   # arc on a rod  -> azimuthal
            n_arc += 1
            continue
        eid = geompy.GetSubShapeID(sect_p, e)
        if w0 or w1:                                    # radial (or wall) segment
            rad_edges.append(e)
            if not w0:                                  # must start on the rod wall
                rad_rev.append(eid)
            continue
        if abs(p1[0] - p0[0]) < TOL:                    # vertical -> azimuthal
            n_vert += 1
            continue
        xm = 0.5 * (p0[0] + p1[0])
        if xm < X_B0 - TOL:
            ax_in.append(e)
            if p0[0] < p1[0]:                           # start at the bundle side
                ax_in_rev.append(eid)
        elif xm > X_B1 + TOL:
            ax_out.append(e)
            if p0[0] > p1[0]:                           # start at the bundle side
                ax_out_rev.append(eid)
        else:                                           # cell side inside B
            n_side += 1

    def edge_group(edges, name):
        if not edges:
            return None
        g = geompy.CreateGroup(sect_p, EDGE)
        geompy.UnionList(g, edges)
        geompy.addToStudyInFather(sect_p, g, name)
        return g

    grp_rad = edge_group(rad_edges, "radial_edges")
    grp_ain = edge_group(ax_in, "axial_edges_A")
    grp_aout = edge_group(ax_out, "axial_edges_C")

    print("[edges] radial %d (reversed %d) | axial A %d / C %d | arcs %d | "
          "cell sides %d | vertical %d"
          % (len(rad_edges), len(rad_rev), len(ax_in), len(ax_out),
             n_arc, n_side, n_vert))


# =============================================================================
#  5) RADIAL GRADING  (first cell at the rod wall = FIRST_LAYER)
# =============================================================================
def growth_ratio(length, first, n):
    """q such that first*(q^n - 1)/(q - 1) = length   (q = 1 -> uniform)."""
    if n < 2 or first <= 0.0:
        return 1.0
    if abs(first * n - length) < 1.0e-12:
        return 1.0
    lo, hi = 0.2, 5.0
    f = lambda q: first * (q ** n - 1.0) / (q - 1.0) if abs(q - 1.0) > 1e-9 else first * n
    if f(lo) > length or f(hi) < length:
        return 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) < length:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if BUILD_MESH:
    # the radial edges are not all the same length : the shortest one (the tightest
    # gap, where the wall shear is highest) is used as the reference, so FIRST_LAYER
    # is the first cell height there.  All edges share the same growth ratio, hence
    # the first cell of a longer edge is larger in the same proportion.
    L_RAD = CLEARANCE                                 # shortest wall -> cell distance
    L_RAD_MAX = L_RAD
    for (cx, cy, has_rod, pg) in POLY:
        if not has_rod:
            continue
        for (vx, vy) in pg:
            L_RAD_MAX = max(L_RAD_MAX, math.hypot(vx - cx, vy - cy) - R)

    Q_RAD = growth_ratio(L_RAD, FIRST_LAYER, N_RAD)
    RAD_SCALE = Q_RAD ** (N_RAD - 1)                  # SALOME "scale" = last/first
    y1_min = L_RAD * (Q_RAD - 1.0) / (Q_RAD ** N_RAD - 1.0) if abs(Q_RAD - 1.0) > 1e-9 \
        else L_RAD / N_RAD
    print("[inflation] %d radial cells, growth %.3f (scale %.2f)" % (N_RAD, Q_RAD, RAD_SCALE))
    print("            radial span %.3f .. %.3f  ->  first cell %.4f .. %.4f"
          % (L_RAD, L_RAD_MAX, y1_min, y1_min * L_RAD_MAX / L_RAD))
    if Q_RAD < 1.0:
        print("[note] FIRST_LAYER exceeds the uniform size -> cells shrink outwards;"
              " reduce FIRST_LAYER or N_RAD")


# =============================================================================
#  6) 2-D STRUCTURED MESH
# =============================================================================
if BUILD_MESH:
    mesh = smesh.Mesh(sect_p, "rod_bundle")

    mesh.Segment().NumberOfSegments(N_AZ_SECT)          # azimuthal / cell sides
    mesh.Quadrangle(algo=smeshBuilder.QUADRANGLE)       # mapped -> structured quads

    def sub_segments(grp, n, scale, reversed_ids):
        if grp is None:
            return None
        hyp = mesh.Segment(geom=grp).NumberOfSegments(n, scale) if scale != 1.0 \
            else mesh.Segment(geom=grp).NumberOfSegments(n)
        if reversed_ids:
            try:
                hyp.SetObjectEntry(sect_entry)
                hyp.SetReversedEdges(reversed_ids)
            except Exception as ex:
                print("[warn] reversed edges not applied (%s)" % ex)
        return hyp

    sub_segments(grp_rad, N_RAD, RAD_SCALE, rad_rev)
    sub_segments(grp_ain, N_X_IN, X_SCALE_IN, ax_in_rev)
    sub_segments(grp_aout, N_X_OUT, X_SCALE_OUT, ax_out_rev)

    ok2d = mesh.Compute()
    n_tri = mesh.NbTriangles()
    print("[2D] %s : %d quadrangles, %d triangles, %d nodes"
          % ("OK" if ok2d else "FAILED", mesh.NbQuadrangles(), n_tri, mesh.NbNodes()))
    if n_tri:
        print("[warn] %d triangles -> some block is not mappable; check the cell shapes"
              % n_tri)


# =============================================================================
#  7) EXTRUSION -> HEXAHEDRA  +  BOUNDARY GROUPS
# =============================================================================
if BUILD_MESH and ok2d:
    faces2d = mesh.GetElementsByType(SMESH.FACE)
    mesh.ExtrusionSweep(faces2d, [0.0, 0.0, H / float(N_Z)], N_Z)
    print("[3D] extruded %d layers -> %d hexahedra, %d nodes"
          % (N_Z, mesh.NbHexas(), mesh.NbNodes()))

    # drop the leftover 2-D / 1-D elements of the section, then rebuild the skin
    mesh.RemoveElements(mesh.GetElementsByType(SMESH.FACE))
    mesh.RemoveElements(mesh.GetElementsByType(SMESH.EDGE))

    skin = None
    try:
        skin = mesh.MakeBoundaryElements(SMESH.BND_2DFROM3D, "skin")[1]
    except Exception as ex:
        print("[warn] boundary faces not created (%s) -> groups may be empty" % ex)

    big = 10.0 * max(L_TOT, W, H)

    def plane_group(name, origin, normal):
        pl = geompy.MakePlane(geompy.MakeVertex(*origin),
                              geompy.MakeVectorDXDYDZ(*normal), big)
        flt = smesh.GetFilter(SMESH.FACE, SMESH.FT_BelongToPlane, SMESH.FT_Undefined,
                              pl, SMESH.FT_Undefined, SMESH.FT_Undefined, 1.0e-4)
        return mesh.GroupOnFilter(SMESH.FACE, name, flt)

    g_in = plane_group("inlet", (0, 0, 0), (1, 0, 0))
    g_out = plane_group("outlet", (L_TOT, 0, 0), (1, 0, 0))
    g_y0 = plane_group("side_y0", (0, 0, 0), (0, 1, 0))
    g_yw = plane_group("side_ymax", (0, W, 0), (0, 1, 0))
    g_bot = plane_group("bottom", (0, 0, 0), (0, 0, 1))
    g_top = plane_group("top", (0, 0, H), (0, 0, 1))

    if skin is not None:
        try:
            g_rod = mesh.CutListOfGroups([skin], [g_in, g_out, g_y0, g_yw, g_bot, g_top],
                                         "rod_walls")
            mesh.UnionListOfGroups([g_y0, g_yw, g_bot, g_top, g_rod], "walls_all")
            print("[groups] inlet %d, outlet %d, rod_walls %d, sides %d/%d, bot/top %d/%d"
                  % (g_in.Size(), g_out.Size(), g_rod.Size(), g_y0.Size(), g_yw.Size(),
                     g_bot.Size(), g_top.Size()))
        except Exception as ex:
            print("[warn] rod_walls group failed (%s)" % ex)

    if EXPORT_MED:
        mesh.ExportMED(EXPORT_MED)
        print(" MED exported : %s" % EXPORT_MED)
    if EXPORT_UNV:
        mesh.ExportUNV(EXPORT_UNV)
        print(" UNV exported : %s" % EXPORT_UNV)


# =============================================================================
#  8) REPORT
# =============================================================================
n_half = sum(1 for (x, y) in CENTERS if y < TOL or y > W - TOL)

print("-" * 66)
print(" domain     : %.3f (x) x %.3f (y) x %.3f (z)" % (L_TOT, W, H))
print(" zones      : A 0..%.3f | B %.3f..%.3f | C %.3f..%.3f"
      % (X_B0, X_B0, X_B1, X_B1, L_TOT))
print(" rods       : %d  (%d full, %d half on the walls)"
      % (len(CENTERS), len(CENTERS) - n_half, n_half))
print(" D=%.3f  ST=%.3f  SL=%.3f  P/D=%.3f  S_D=%.3f  gap=%.3f"
      % (D, ST, SL, ST / D, S_D, ST - D))
print(" cell       : hexagon  half-height %.3f, half-length %.3f (flat %.3f)"
      % (HH, XV, XE))
if BUILD_MESH:
    print(" azimuthal  : %d per rod (%d per hexagon side, min asked %d)"
          % (6 * N_AZ_SECT, N_AZ_SECT, N_AZ_MIN))
    print(" radial     : %d cells, growth %.3f, first %.4f" % (N_RAD, Q_RAD, y1_min))
    print(" axial      : A %d / C %d cells, z %d layers" % (N_X_IN, N_X_OUT, N_Z))
print("-" * 66)

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()
