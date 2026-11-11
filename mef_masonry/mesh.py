"""Meshes: 2D horizontal section of the prism (triangles and quadrilaterals) and 3D
prism by extrusion: 6-node wedges (W6), 8-node hexahedra (H8) or 15-node
quadratic wedges (W15).

The section is constant over the height (holes aligned between courses), so
the 3D prism is obtained by extruding the 2D mesh. The 2D mesh is structured and
the nodes on the hole boundaries are projected onto the circle, which gives smooth
holes and a conforming mesh without an external mesh generator.
"""
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .parameters import Geometry

# Materials
AIR, CERAMIC, JOINT, CAPPING, STEEL, MORTAR, GROUT = 0, 1, 2, 3, 4, 5, 6
NAMES = {AIR: "air (holes)", CERAMIC: "unit", JOINT: "dry joint",
         CAPPING: "capping", STEEL: "steel (plates)", MORTAR: "mortar", GROUT: "grout"}

# Markers of the section triangles
SOLID, HOLE = -1, 0  # head-joint bands: 1, 2, ...


@dataclass
class Mesh2D:
    xy: np.ndarray          # (n, 2) coordinates, m
    tri: np.ndarray         # (ne, 3) connectivity, counter-clockwise
    mark: np.ndarray       # (ne,) SOLID, HOLE or band index (1, 2, ...)
    bands: Tuple[float, ...]
    length: float
    model_width: float   # modelled width (half if there is symmetry)
    width: float
    symmetry: bool
    quad: Optional[np.ndarray] = None        # (nq, 4) source quadrilaterals (counter-clockwise)
    mark_quad: Optional[np.ndarray] = None  # (nq,)
    band_width: float = 0.002             # thickness of the head joints (and set-back of the cells)
    rectangular: bool = False                 # rectangular cells (holes may change from course to course)
    vertical: bool = False                   # vertical section of one course (horizontal holes): the first
                                             # coordinate is the height in the course, the second the thickness

    @property
    def n_nodes(self):
        return len(self.xy)

    def areas(self):
        p = self.xy[self.tri]
        return 0.5 * ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                      - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))

    def boundary_edges(self):
        """Edges that belong to a single triangle: (a, b)."""
        e = np.vstack([self.tri[:, [0, 1]], self.tri[:, [1, 2]], self.tri[:, [2, 0]]])
        es = np.sort(e, axis=1)
        _, inv, cnt = np.unique(es, axis=0, return_inverse=True, return_counts=True)
        return e[cnt[inv.ravel()] == 1]


def _coordinates(a, b, h, required=()):
    n = max(1, int(np.ceil((b - a) / h - 1e-9)))
    c = np.linspace(a, b, n + 1)
    req = np.asarray(sorted(required), float)
    if req.size:
        keep = np.ones(c.size, bool)
        for o in req:
            keep &= np.abs(c - o) > 0.35 * h
        keep[0] = keep[-1] = True
        c = np.unique(np.round(np.concatenate([c[keep], req]), 12))
    return c


def _quads_to_triangles(q):
    """q: (nq, 4) quadrilateral nodes; splits them alternating the diagonal."""
    alt = (np.arange(len(q)) % 2 == 0)[:, None]
    t1 = np.where(alt, q[:, [0, 1, 2]], q[:, [0, 1, 3]])
    t2 = np.where(alt, q[:, [0, 2, 3]], q[:, [1, 2, 3]])
    return np.vstack([t1, t2])


def _grid(points, h):
    """Coordinates that pass through all given points, with a spacing <= h between them."""
    p = np.unique(np.round(np.asarray(points, float), 9))
    out = [p[0]]
    for a, b in zip(p[:-1], p[1:]):
        n = max(1, int(np.ceil((b - a) / h - 1e-9)))
        out.extend(a + (b - a) * np.arange(1, n + 1) / n)
    return np.asarray(out)


def rect_section_mesh(geo: Geometry, h: float = 0.005, symmetry: bool = True,
                           bands: Optional[Sequence[float]] = None,
                           t_joint: Optional[float] = None) -> Mesh2D:
    """Section (plan) of units with rectangular cells (concrete, structural
    clay) or with horizontal holes (infill: the grid follows the edges of the
    columns of holes, and the extrusion opens the holes in the rows).

    Structured grid aligned with all cell edges of all courses
    (with running bond the cells change position from one course to the next) and
    with the head-joint bands. The HOLE mark follows course 0 (2D section);
    the extrusion classifies each element course by course (course_states).
    """
    L, W = geo.length, geo.width
    tj = geo.joint_thickness if t_joint is None else t_joint
    bands = tuple(geo.joint_positions if bands is None else bands)
    px = [0.0, L] + [p for x in bands for p in (x - tj / 2, x + tj / 2)]
    py = [0.0, W / 2, W] + [v for f in geo.holes_h_y() for v in f]
    for i in range(geo.n_courses):
        for c in geo.cells(i, tj):
            px += [c[0], c[1]]
            py += [c[2], c[3]]
    xs, ys = _grid(px, h), _grid(py, h)
    if symmetry:
        ys = ys[ys <= W / 2 + 1e-12]
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    nx, ny = X.shape
    xy = np.c_[X.ravel(), Y.ravel()]
    idx = np.arange(nx * ny).reshape(nx, ny)
    i, j = (a.ravel() for a in np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), indexing="ij"))
    n00, n10, n11, n01 = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
    alt = ((i + j) % 2 == 0)[:, None]
    tri = np.vstack([np.where(alt, np.c_[n00, n10, n11], np.c_[n00, n10, n01]),
                     np.where(alt, np.c_[n00, n11, n01], np.c_[n10, n11, n01])])
    quad = np.c_[n00, n10, n11, n01]

    def tag(pol):
        cnt = xy[pol].mean(1)
        mk = np.full(len(pol), SOLID)
        for k, xf in enumerate(bands, start=1):
            mk[np.abs(cnt[:, 0] - xf) < tj / 2] = k
        for c in geo.cells(0, tj):
            mk[(cnt[:, 0] > c[0]) & (cnt[:, 0] < c[1]) & (cnt[:, 1] > c[2]) & (cnt[:, 1] < c[3])] = HOLE
        return mk
    return Mesh2D(xy=xy, tri=tri, mark=tag(tri), bands=bands, length=L,
                   model_width=W / 2 if symmetry else W, width=W, symmetry=symmetry,
                   quad=quad, mark_quad=tag(quad), band_width=tj, rectangular=True)


def vertical_section_mesh(geo: Geometry, h: float = 0.005, t_joint: Optional[float] = None) -> Mesh2D:
    """Vertical section of one course of units with horizontal holes (thickness ×
    height, centre to centre), for the 2D thermal analysis and the Wall mode. The first
    coordinate is the height in the course (0 to course_height) and the second the thickness
    (y = 0 is the exposed face). Half a joint below and half above (marks 1 and 2);
    holes with the HOLE mark. As the holes run through the unit, the section holds
    over the whole length of the wall."""
    W, hf = geo.width, geo.course_height
    tj = geo.joint_thickness if t_joint is None else t_joint
    z_lines = [0.0, tj / 2, hf - tj / 2, hf] + [v for f in geo.holes_h_z(hf / 2, tj) for v in f]
    ys = _grid([0.0, W / 2, W] + [v for f in geo.holes_h_y() for v in f], h)
    zs = _grid(z_lines, h)
    Z, Y = np.meshgrid(zs, ys, indexing="ij")
    nz, ny = Z.shape
    xy = np.c_[Z.ravel(), Y.ravel()]
    idx = np.arange(nz * ny).reshape(nz, ny)
    i, j = (a.ravel() for a in np.meshgrid(np.arange(nz - 1), np.arange(ny - 1), indexing="ij"))
    n00, n10, n11, n01 = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
    alt = ((i + j) % 2 == 0)[:, None]
    tri = np.vstack([np.where(alt, np.c_[n00, n10, n11], np.c_[n00, n10, n01]),
                     np.where(alt, np.c_[n00, n11, n01], np.c_[n10, n11, n01])])
    quad = np.c_[n00, n10, n11, n01]

    def tag(pol):
        cnt = xy[pol].mean(1)
        mk = np.full(len(pol), SOLID)
        mk[cnt[:, 0] < tj / 2] = 1
        mk[cnt[:, 0] > hf - tj / 2] = 2
        for z0, z1 in geo.holes_h_z(hf / 2, tj):
            for y0, y1 in geo.holes_h_y():
                mk[(cnt[:, 0] > z0) & (cnt[:, 0] < z1) & (cnt[:, 1] > y0) & (cnt[:, 1] < y1)] = HOLE
        return mk
    return Mesh2D(xy=xy, tri=tri, mark=tag(tri), bands=(0.0, hf), length=hf, model_width=W,
                   width=W, symmetry=False, quad=quad, mark_quad=tag(quad), band_width=tj, vertical=True)


def course_states(m2d: Mesh2D, geo: Geometry, i: int, quad: bool = False):
    """Hole and grout per 2D element (triangle, or quadrilateral if quad) in course i.

    Round holes: aligned in all courses (HOLE mark); the holes of the
    chosen columns are grouted. Rectangular cells: classified by the
    cells of course i; a main cell is grouted if it covers more than
    half of a grouted column (defined in course 0)."""
    pol = m2d.quad if quad else m2d.tri
    mark = m2d.mark_quad if quad else m2d.mark
    cnt = m2d.xy[pol].mean(1)
    cols = geo.grouted_columns(m2d.band_width)
    if m2d.rectangular and geo.rectangular:
        hole = np.zeros(len(pol), bool)
        grout = np.zeros(len(pol), bool)
        for c in geo.cells(i, m2d.band_width):
            inside = (cnt[:, 0] > c[0]) & (cnt[:, 0] < c[1]) & (cnt[:, 1] > c[2]) & (cnt[:, 1] < c[3])
            hole |= inside
            if c[4] and any(min(c[1], b) - max(c[0], a) > 0.5 * (c[1] - c[0]) for a, b in cols):
                grout |= inside
        return hole, grout
    hole = mark == HOLE
    grout = np.zeros(len(pol), bool)
    for a, b in cols:
        grout |= hole & (cnt[:, 0] > a - 1e-9) & (cnt[:, 0] < b + 1e-9)
    return hole, grout


def section_mesh(geo: Geometry = Geometry(), h: float = 0.005, symmetry: bool = True,
                bands: Optional[Sequence[float]] = (0.14, 0.28), band_width: float = 0.002,
                holes: Optional[Sequence[Tuple[float, float]]] = None,
                radius: Optional[float] = None) -> Mesh2D:
    """Mesh of the horizontal section of the one-and-a-half-unit prism (420 x 140 mm).

    Outside the holes the mesh is rectangular and structured; around each hole there is
    a square box with an "O" mesh (concentric rings), which follows the
    circle exactly. The inside of the hole (air) is also meshed,
    for the thermal analysis.

    h        : element size (m)
    symmetry : models only y <= W/2 (symmetry plane through the hole centres)
    bands    : x positions of the head joints (thin bands where a joint may exist)
    holes    : centres (x, y) of the holes, all with the same y; default = BCEPQ
    radius   : hole radius; default = hole_diameter/2 (0 = solid section)
    """
    if geo.horizontal and holes is None and radius is None:
        return vertical_section_mesh(geo, h)
    if geo.rectangular and holes is None and radius is None:
        return rect_section_mesh(geo, h, symmetry, bands or (),
                                      band_width if bands else geo.joint_thickness)
    L, W = geo.length, geo.width
    R = geo.hole_diameter / 2 if radius is None else radius
    if holes is None:
        holes = [(x, W / 2) for x in geo.hole_centres]
    holes = sorted(holes) if R > 0 else []
    bands = tuple(bands or ())
    borders = [p for x in bands for p in (x - band_width / 2, x + band_width / 2)]

    if holes:
        cy = holes[0][1]
        if any(abs(f[1] - cy) > 1e-12 for f in holes):
            raise ValueError("The holes must be aligned (same y).")
        xs_f = [f[0] for f in holes]
        gap = [min(cy, W - cy), xs_f[0], L - xs_f[-1]]
        if len(xs_f) > 1:
            gap.append(float(np.min(np.diff(xs_f))) / 2)
        a = R + (min(gap) - R) / 2               # half-width of the O box
        n = max(8, 2 * int(round(a / h)))          # segments per side (even)
        side = lambda c: c - a + 2 * a * np.arange(n + 1) / n
    # x and y coordinates of the rectangular grid
    parts, ini = [], 0.0
    for cx, _ in holes:
        parts += [_coordinates(ini, cx - a, h, [b for b in borders if ini < b < cx - a]), side(cx)]
        ini = cx + a
    parts.append(_coordinates(ini, L, h, [b for b in borders if ini < b < L]))
    xs = np.unique(np.round(np.concatenate(parts), 12))
    if holes:
        ys = np.unique(np.round(np.concatenate(
            [_coordinates(0.0, cy - a, h), side(cy), _coordinates(cy + a, W, h)]), 12))
    else:
        ys = _coordinates(0.0, W, h, [W / 2])

    X, Y = np.meshgrid(xs, ys, indexing="ij")
    nx, ny = X.shape
    pts = [np.c_[X.ravel(), Y.ravel()]]
    nn = nx * ny
    idx = np.arange(nn).reshape(nx, ny)
    i, j = (v.ravel() for v in np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), indexing="ij"))
    cxr = 0.5 * (xs[i] + xs[i + 1]); cyr = 0.5 * (ys[j] + ys[j + 1])
    outside = np.ones(len(i), bool)
    for cx, _ in holes:
        outside &= ~((np.abs(cxr - cx) < a) & (np.abs(cyr - cy) < a))
    i, j = i[outside], j[outside]
    q = np.c_[idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]]
    tris = [_quads_to_triangles(q)]
    marks = [np.full(len(tris[0]), SOLID)]
    quads, marks_q = [q], [np.full(len(q), SOLID)]

    for cx, _ in holes:
        ix0 = int(np.argmin(np.abs(xs - (cx - a)))); iy0 = int(np.argmin(np.abs(ys - (cy - a))))
        # boundary of the box, counter-clockwise from the bottom-left corner
        per = np.concatenate([idx[ix0:ix0 + n, iy0], idx[ix0 + n, iy0:iy0 + n],
                              idx[ix0 + n:ix0:-1, iy0 + n], idx[ix0, iy0 + n:iy0:-1]])
        Q = pts[0][per]
        c = np.array([cx, cy])
        d = Q - c
        C = c + R * d / np.linalg.norm(d, axis=1)[:, None]
        b = 0.45 * R
        I = c + d * (b / a)
        m1 = max(2, int(round(((a - R) + (a * np.sqrt(2) - R)) / 2 / h)))
        m2 = max(2, int(round((R - b * 1.2) / h)))
        # rings: inside (I) -> circle (C) -> box (Q)
        rings = [I + (C - I) * s for s in np.arange(m2) / m2] + \
                [C + (Q - C) * s for s in np.arange(m1) / m1]
        ring_ids = []
        for ring in rings:
            ring_ids.append(nn + np.arange(len(ring))); pts.append(ring); nn += len(ring)
        ring_ids.append(per)
        for k in range(len(ring_ids) - 1):
            p0, p1 = ring_ids[k], ring_ids[k + 1]
            qq = np.c_[p0, np.roll(p0, -1), np.roll(p1, -1), p1]
            t = _quads_to_triangles(qq)
            tris.append(t)
            marks.append(np.full(len(t), HOLE if k < m2 else SOLID))
            quads.append(qq); marks_q.append(np.full(len(qq), HOLE if k < m2 else SOLID))
        # structured inner square (air)
        g = -np.ones((n + 1, n + 1), int)
        ring0 = ring_ids[0]
        pos = np.concatenate([np.c_[np.arange(n), np.zeros(n, int)],
                              np.c_[np.full(n, n), np.arange(n)],
                              np.c_[np.arange(n, 0, -1), np.full(n, n)],
                              np.c_[np.zeros(n, int), np.arange(n, 0, -1)]])
        g[pos[:, 0], pos[:, 1]] = ring0
        ii, jj = np.nonzero(g < 0)
        new_ones = c + np.c_[-b + 2 * b * ii / n, -b + 2 * b * jj / n]
        g[ii, jj] = nn + np.arange(len(ii)); pts.append(new_ones); nn += len(ii)
        a_, b_ = (v.ravel() for v in np.meshgrid(np.arange(n), np.arange(n), indexing="ij"))
        qq = np.c_[g[a_, b_], g[a_ + 1, b_], g[a_ + 1, b_ + 1], g[a_, b_ + 1]]
        t = _quads_to_triangles(qq)
        tris.append(t); marks.append(np.full(len(t), HOLE))
        quads.append(qq); marks_q.append(np.full(len(qq), HOLE))

    xy = np.vstack(pts)
    tri = np.vstack(tris)
    mark = np.concatenate(marks)
    quad = np.vstack(quads)
    mark_q = np.concatenate(marks_q)
    cnt = xy[tri].mean(axis=1)
    cenq = xy[quad].mean(axis=1)
    for k, xf in enumerate(bands, start=1):
        mark[(mark == SOLID) & (np.abs(cnt[:, 0] - xf) < band_width / 2)] = k
        mark_q[(mark_q == SOLID) & (np.abs(cenq[:, 0] - xf) < band_width / 2)] = k
    if symmetry:
        keep = cnt[:, 1] < W / 2
        tri, mark = tri[keep], mark[keep]
        mq = cenq[:, 1] < W / 2
        quad, mark_q = quad[mq], mark_q[mq]
    # removes orphan nodes and renumbers
    used = np.unique(tri)
    new = -np.ones(len(xy), int); new[used] = np.arange(len(used))
    xy, tri, quad = xy[used], new[tri], new[quad]
    # counter-clockwise orientation
    p = xy[tri]
    air = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
    tri[air < 0] = tri[air < 0][:, [0, 2, 1]]
    pq = xy[quad]
    fname = 0.5 * np.sum(pq[:, :, 0] * np.roll(pq[:, :, 1], -1, axis=1) - np.roll(pq[:, :, 0], -1, axis=1) * pq[:, :, 1], axis=1)
    quad[fname < 0] = quad[fname < 0][:, ::-1]
    m = Mesh2D(xy=xy, tri=tri, mark=mark, bands=bands, length=L,
                model_width=W / 2 if symmetry else W, width=W, symmetry=symmetry,
                quad=quad, mark_quad=mark_q, band_width=band_width)
    if (np.abs(m.areas()) < 1e-14).any():
        raise RuntimeError("2D mesh with degenerate elements.")
    return m


def rectangle_mesh(Lx: float, Ly: float, h: float, hy: Optional[float] = None) -> Mesh2D:
    """Solid rectangle (used in the verifications); hy allows an anisotropic mesh."""
    xs = _coordinates(0.0, Lx, h)
    ys = _coordinates(0.0, Ly, hy or h)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    nx, ny = X.shape
    xy = np.c_[X.ravel(), Y.ravel()]
    idx = np.arange(nx * ny).reshape(nx, ny)
    i, j = (a.ravel() for a in np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), indexing="ij"))
    n00, n10, n11, n01 = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
    alt = ((i + j) % 2 == 0)[:, None]
    t1 = np.where(alt, np.c_[n00, n10, n11], np.c_[n00, n10, n01])
    t2 = np.where(alt, np.c_[n00, n11, n01], np.c_[n10, n11, n01])
    tri = np.vstack([t1, t2])
    quad = np.c_[n00, n10, n11, n01]
    return Mesh2D(xy=xy, tri=tri, mark=np.full(len(tri), SOLID), bands=(),
                   length=Lx, model_width=Ly, width=Ly, symmetry=False,
                   quad=quad, mark_quad=np.full(len(quad), SOLID))


# ----------------------------------------------------------------------------
# 3D extrusion
# ----------------------------------------------------------------------------
@dataclass
class Layer:
    z0: float
    z1: float
    kind: str            # 'plate', 'capping', 'joint', 'unit'
    course: int = -1      # 'unit': the course; 'joint': the course below; 'capping': the neighbouring course
    hole_h: int = -1     # 'unit' layer of a unit with horizontal holes: row of holes in this sublayer (-1 = shell)
    material: int = CERAMIC  # used when kind == 'uniform'


@dataclass
class Mesh3D:
    xyz: np.ndarray        # (n, 3)
    conn3: np.ndarray        # (ne, 6): 3 bottom nodes + 3 top nodes, same order
    mat: np.ndarray        # (ne,) material
    layer: np.ndarray     # (ne,) sublayer index
    tri_of: np.ndarray     # (ne,) source 2D triangle (quadrilateral, for H8)
    zs: np.ndarray         # heights of the node planes
    m2d: Mesh2D
    layers: List[Layer] = field(default_factory=list)
    kind: str = "W6"       # element: W6, H8 or W15

    @property
    def n_nodes(self):
        return len(self.xyz)

    @property
    def elem(self):
        from .elements3d import kind
        return kind(self.kind)

    @property
    def n2(self):
        return self.m2d.n_nodes


def prism_layers(geo: Geometry = Geometry(), t_joint: float = 0.002,
                   capping: bool = True, plates: bool = True,
                   t_capping: Optional[float] = None) -> List[Layer]:
    """Layers of the prism: bottom plate, capping, courses with joints
    between them, capping and top plate. The joints are centred on the
    interfaces between courses (z = 60, 120, ... mm). The capping has the joint
    thickness, unless t_capping is given."""
    hf, n = geo.course_height, geo.n_courses
    H = hf * n
    tc = t_joint if t_capping is None else t_capping
    lay: List[Layer] = []
    if plates:
        lay.append(Layer(-geo.plate_thickness, 0.0, "plate"))
    z = 0.0
    if capping:
        lay.append(Layer(0.0, tc, "capping", course=0))
        z = tc
    for i in range(n):
        top = (i + 1) * hf - (t_joint / 2 if i < n - 1 else (tc if capping else 0.0))
        # units with horizontal holes: the unit layer is split at the edges of the rows of holes
        cuts = [(z0, z1, k) for k, (z0, z1) in enumerate(geo.holes_h_z((i + 0.5) * hf, t_joint))]
        for z0, z1, k in cuts:
            if z0 > z + 1e-9:
                lay.append(Layer(z, z0, "unit", course=i))
            lay.append(Layer(z0, z1, "unit", course=i, hole_h=k))
            z = z1
        lay.append(Layer(z, top, "unit", course=i))
        z = top
        if i < n - 1:
            lay.append(Layer(z, z + t_joint, "joint", course=i))
            z += t_joint
    if capping:
        lay.append(Layer(z, H, "capping", course=n - 1))
    if plates:
        lay.append(Layer(H, H + geo.plate_thickness, "plate"))
    return lay


def _tri_edges(tri):
    """Unique edges of the triangles and, per triangle, the edge indices (0-1, 1-2, 2-0)."""
    e = np.stack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]], axis=1)      # (nt, 3, 2)
    es = np.sort(e, axis=2).reshape(-1, 2)
    uniq, inv = np.unique(es, axis=0, return_inverse=True)
    return uniq, inv.reshape(-1, 3)


def extrude(m2d: Mesh2D, layers: Sequence[Layer], dz: float = 0.010,
             geo: Geometry = Geometry(), kind: str = "W6") -> Mesh3D:
    """Extrudes the section over the given layers (subdividing each with a step <= dz).
    kind: 'W6' (triangles -> wedges), 'H8' (quadrilaterals -> hexahedra) or
    'W15' (triangles -> quadratic wedges, with mid-edge nodes)."""
    zs = [layers[0].z0]
    info = []  # (layer, index)
    for ic, c in enumerate(layers):
        nsub = max(1, int(np.ceil((c.z1 - c.z0) / dz - 1e-9)))
        for s in range(nsub):
            zs.append(c.z0 + (c.z1 - c.z0) * (s + 1) / nsub)
            info.append(ic)
    zs = np.asarray(zs)
    n2 = m2d.n_nodes
    nz = len(zs) - 1
    if kind == "H8":
        if m2d.quad is None:
            raise ValueError("The section has no quadrilaterals for H8.")
        base2d, mark2d = m2d.quad, m2d.mark_quad
    else:
        base2d, mark2d = m2d.tri, m2d.mark
    nt = len(base2d)
    k = np.repeat(np.arange(nz), nt)
    t = np.tile(base2d, (nz, 1))
    if kind == "W15":
        edges, eid = _tri_edges(m2d.tri)
        ne2 = len(edges)
        n2q = n2 + ne2
        xy2 = np.vstack([m2d.xy, 0.5 * (m2d.xy[edges[:, 0]] + m2d.xy[edges[:, 1]])])
        zm = 0.5 * (zs[:-1] + zs[1:])
        xyz = np.vstack([np.c_[np.tile(xy2, (nz + 1, 1)), np.repeat(zs, n2q)],
                         np.c_[np.tile(m2d.xy, (nz, 1)), np.repeat(zm, n2)]])
        e = np.tile(eid, (nz, 1))
        kb, kt = (k * n2q)[:, None], ((k + 1) * n2q)[:, None]
        km = ((nz + 1) * n2q + k * n2)[:, None]
        conn3 = np.c_[t + kb, t + kt, n2 + e + kb, n2 + e + kt, t + km]
    else:
        xyz = np.c_[np.tile(m2d.xy, (nz + 1, 1)), np.repeat(zs, n2)]
        conn3 = np.c_[t + (k * n2)[:, None], t + ((k + 1) * n2)[:, None]]
    mat = np.empty(nz * nt, dtype=np.int8)
    joints = geo.joints_per_course
    mortar = geo.joint_type == "mortar"
    mat_vertical = (MORTAR if geo.filled_head_joints else AIR) if mortar else JOINT
    # face-shell bedding: mortar only on the bands of the face shells (|y - W/2| >= W/2 - shell)
    yc = m2d.xy[base2d].mean(1)[:, 1]
    bed = np.ones(len(base2d), bool)
    if mortar and geo.bedding == "lateral" and geo.has_holes:
        bed = np.abs(yc - geo.width / 2) >= geo.width / 2 - geo.shell_thickness - 1e-9
    # state of each 2D element in each course: 0 solid, 1 hole (air), 2 grouted hole
    nf = max(1, geo.n_courses)
    states = []
    for i in range(nf):
        f, g = course_states(m2d, geo, i, quad=(kind == "H8"))
        states.append(np.where(g, 2, np.where(f, 1, 0)))

    def sty(i):
        return states[min(max(i, 0), nf - 1)]
    for kk in range(nz):
        c = layers[info[kk]]
        sl = slice(kk * nt, (kk + 1) * nt)
        mk = mark2d
        if c.kind == "plate":
            mat[sl] = STEEL
        elif c.kind == "uniform":
            mat[sl] = np.where(mk == HOLE, AIR, c.material)
        elif c.kind == "capping":
            e = sty(c.course)
            mat[sl] = np.where(e == 1, AIR, CAPPING)
        elif c.kind == "joint":
            # joint between courses i and i+1: grout where there is grout and no empty hole
            # on either side; mortar (or dry joint) where it is solid on both sides
            e0, e1 = sty(c.course), sty(c.course + 1 if c.course >= 0 else 0)
            empty = (e0 == 1) | (e1 == 1)
            grt = ~empty & ((e0 == 2) | (e1 == 2))
            sol = (e0 == 0) & (e1 == 0)
            base = np.where(sol & bed, MORTAR, AIR) if mortar else np.where(sol, JOINT, AIR)
            mat[sl] = np.where(grt, GROUT, base)
        else:  # unit
            e = sty(c.course)
            m = np.where(e == 2, GROUT, np.where(e == 1, AIR, CERAMIC))
            band = np.zeros(len(base2d), bool)            # head joints of this course
            if m2d.bands and 0 <= c.course < len(joints):
                for xj in joints[c.course]:
                    kf = int(np.argmin([abs(xj - f) for f in m2d.bands])) + 1
                    band |= mk == kf
                m = np.where(band & (e == 0), mat_vertical, m)
            if c.hole_h >= 0:                              # row of horizontal holes
                col = np.zeros(len(base2d), bool)
                for y0, y1 in geo.holes_h_y():
                    col |= (yc > y0) & (yc < y1)
                m = np.where(col & ~band, AIR, m)
            mat[sl] = m
    return Mesh3D(xyz=xyz, conn3=conn3, mat=mat, layer=np.repeat(np.asarray(info), nt),
                   tri_of=np.tile(np.arange(nt), nz), zs=zs, m2d=m2d, layers=list(layers), kind=kind)


def prism_mesh(geo: Geometry = Geometry(), h: float = 0.006, dz: float = 0.010,
                 t_joint: float = 0.002, symmetry: bool = True, plates: bool = True,
                 capping: bool = True, element: str = "W6") -> Mesh3D:
    """Shortcut: 3D mesh of the prism (default: 5 courses of one and a half units, dry joints)."""
    problems = geo.validate()
    if problems:
        raise ValueError(" ".join(problems))
    if geo.rectangular or geo.horizontal:
        m2d = rect_section_mesh(geo, h=h, symmetry=symmetry, bands=geo.joint_positions, t_joint=t_joint)
    else:
        m2d = section_mesh(geo, h=h, symmetry=symmetry, band_width=t_joint, bands=geo.joint_positions)
    t_cap = t_joint if geo.joint_type == "dry" else min(t_joint, 0.003)
    lay = prism_layers(geo, t_joint=t_joint, capping=capping, plates=plates, t_capping=t_cap)
    return extrude(m2d, lay, dz=dz, geo=geo, kind=element)


def nearest_node(xyz: np.ndarray, point) -> int:
    return int(np.argmin(np.linalg.norm(xyz - np.asarray(point)[None, :], axis=1)))


def summary(m) -> str:
    if isinstance(m, Mesh3D):
        cont = {NAMES[k]: int((m.mat == k).sum()) for k in NAMES if (m.mat == k).any()}
        name = {"W6": "W6 wedges", "H8": "H8 hexahedra", "W15": "W15 wedges"}.get(m.kind, m.kind)
        return (f"{m.n_nodes} nodes, {len(m.conn3)} {name}, {len(m.zs) - 1} layers; "
                + ", ".join(f"{k}: {v}" for k, v in cont.items()))
    a = m.areas()
    return f"{m.n_nodes} nodes, {len(m.tri)} triangles; min/median area = {a.min() / np.median(a):.2f}"
