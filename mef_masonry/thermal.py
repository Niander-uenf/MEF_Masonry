"""Transient heat conduction (2D or 3D) with convection and radiation on the boundary.

Formulation: Galerkin with lumped capacity and backward Euler.
Properties and radiation coefficients are evaluated at the temperature of the
previous step (semi-implicit); the step dt controls the accuracy (the verification
example shows the convergence).
"""
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .mesh import AIR, CERAMIC, JOINT, CAPPING, STEEL, MORTAR, GROUT, Mesh2D, Mesh3D
from .parameters import Ceramic, Mortar, Grout, DryJoint, Steel, HoleAir, SIGMA
from . import elements as el


@dataclass
class Materials:
    ceramic: Ceramic = field(default_factory=Ceramic)
    mortar: Mortar = field(default_factory=Mortar)
    grout: Grout = field(default_factory=Grout)
    joint: DryJoint = field(default_factory=DryJoint)
    capping: DryJoint = field(default_factory=DryJoint)
    steel: Steel = field(default_factory=Steel)
    air: HoleAir = field(default_factory=HoleAir)
    hole_diameter: float = 0.084

    def k(self, mat, Te):
        k = np.empty(len(mat))
        cer = self.ceramic
        for m, sel in _groups(mat):
            T = Te[sel]
            if m == CERAMIC:
                k[sel] = cer.k_T(T)
            elif m == MORTAR:
                k[sel] = self.mortar.k_T(T)
            elif m == GROUT:
                k[sel] = self.grout.k_T(T)
            elif m == JOINT:
                k[sel] = cer.k_T(T) * self.joint.k_factor
            elif m == CAPPING:
                k[sel] = cer.k_T(T) * self.capping.k_factor
            elif m == STEEL:
                k[sel] = self.steel.k
            elif m == AIR:
                TK = T + 273.15
                eps = self.air.wall_emissivity
                k[sel] = (0.0241 * (TK / 273.15) ** 0.81
                          + 4 * SIGMA * TK ** 3 * self.hole_diameter * eps / (2 - eps))
        return k

    @property
    def variable_capacity(self):
        return (self.ceramic.variable_capacity or self.mortar.variable_capacity
                or self.grout.variable_capacity)

    def rho_c(self, mat, Te=None):
        """Volumetric heat capacity per element (J/m³.K); with Te, uses the
        c(T) and rho(T) tables (apparent heat capacity method)."""
        rc = np.empty(len(mat))
        for m, sel in _groups(mat):
            if m in (CERAMIC, JOINT, CAPPING, MORTAR, GROUT):
                mt = {MORTAR: self.mortar, GROUT: self.grout}.get(m, self.ceramic)
                if Te is None:
                    rc[sel] = mt.rho * mt.c
                else:
                    T = Te[sel]
                    rc[sel] = mt.rho_T(T) * mt.c_T(T)
            elif m == STEEL:
                rc[sel] = self.steel.rho * self.steel.c
            else:
                rc[sel] = self.air.rho_c
        return rc


def _groups(mat):
    for m in np.unique(mat):
        yield int(m), mat == m


@dataclass
class Boundary:
    """Group of faces exposed to a gas at temperature T_gas(t) (°C).
    node_area: area (3D, m²) or length (2D, m) associated with each node.
    faces: optional, (face connectivity (nf, k), mass matrices of the
    faces (nf, k, k)) for the consistent heat exchange (quadratic elements)."""
    node_area: np.ndarray
    temperature: Callable[[float], float]
    h_convection: float
    emissivity: float
    name: str = ""
    faces: Optional[tuple] = None


@dataclass
class ThermalResult:
    t: np.ndarray
    T_gas: np.ndarray
    probes: Dict[str, np.ndarray]
    T_final: np.ndarray
    fields: Dict[float, np.ndarray]
    stopped_at: Optional[float] = None


def csr_structure(conn, n):
    """Sparsity pattern and map element -> position in data (for fast assembly)."""
    nen = conn.shape[1]
    r = np.repeat(conn, nen, axis=1).ravel().astype(np.int64)
    c = np.tile(conn, (1, nen)).ravel().astype(np.int64)
    key = r * n + c
    uniq, emap = np.unique(key, return_inverse=True)
    lines = uniq // n
    indices = (uniq % n).astype(np.int32)
    indptr = np.searchsorted(lines, np.arange(n + 1)).astype(np.int32)
    diag = np.searchsorted(uniq, np.arange(n, dtype=np.int64) * (n + 1))
    return indices, indptr, emap.ravel(), diag, len(uniq)


def _cg(A, b, x0, rtol=1e-10):
    d = A.diagonal()
    M = spla.LinearOperator(A.shape, matvec=lambda v: v / d)
    try:
        x, info = spla.cg(A, b, x0=x0, rtol=rtol, atol=0.0, M=M, maxiter=5000)
    except TypeError:  # old SciPy
        x, info = spla.cg(A, b, x0=x0, tol=rtol, M=M, maxiter=5000)
    if info != 0:
        x = spla.spsolve(A.tocsc(), b)
    return x


def solve_(conn: np.ndarray, Ke1: np.ndarray, volumes: np.ndarray, mat: np.ndarray,
             n_nodes: int, M: Materials, boundaries: Sequence[Boundary], t_final: float,
             dt: float = 30.0, T_initial: float = 30.0, probes: Optional[Dict] = None,
             record_every: float = 60.0, field_instants: Sequence[float] = (),
             stop_if: Optional[Callable[[float, np.ndarray], bool]] = None,
             verbose: bool = False, cap_weights: Optional[np.ndarray] = None,
             cap_elem: Optional[np.ndarray] = None) -> ThermalResult:
    """Integrates the transient conduction from t = 0 to t_final.

    conn, Ke1, volumes : connectivity, conduction matrices with k = 1 and
                         volumes (or areas, in 2D) of the elements
    probes    : {name: node index or list of nodes (mean)}
    field_instants : times at which the full field is stored
    stop_if   : function (t, T) -> bool to stop early (e.g. insulation criterion)
    cap_weights : fractions (ne, nen) of the lumped capacity at each node (default:
                equal shares)
    cap_elem  : consistent capacity matrices per unit of rho c,
                (ne, nen, nen); if given, they replace the lumped capacity
    Boundaries with .faces use the consistent heat-exchange matrix
    (h integrated with N_i N_j on each face, h at the mean temperature of the face).
    """
    nen = conn.shape[1]
    indices, indptr, emap, diag, nnz = csr_structure(conn, n_nodes)
    if cap_weights is None:
        cap_weights = np.full(conn.shape, 1.0 / nen)

    def capacity(rc):
        """Capacity (consistent: CSR data and matrix; lumped: nodal vector)."""
        if cap_elem is not None:
            dC = np.bincount(emap, weights=(rc[:, None, None] * cap_elem).ravel(), minlength=nnz)
            return dC, sp.csr_matrix((dC, indices, indptr), shape=(n_nodes, n_nodes)), None
        return None, None, np.bincount(conn.ravel(), weights=((rc * volumes)[:, None] * cap_weights).ravel(),
                                       minlength=n_nodes)

    variable = M.variable_capacity
    data_C, Cmat, C = capacity(M.rho_c(mat))
    pairs = {}
    for i, g in enumerate(boundaries):
        if g.faces is not None:
            fc = g.faces[0]
            k = fc.shape[1]
            pairs[i] = (np.repeat(fc, k, 1).ravel(), np.tile(fc, (1, k)).ravel())
    T = np.full(n_nodes, float(T_initial))
    probes = probes or {}
    probes = {k: np.atleast_1d(v) for k, v in probes.items()}
    hist_t, hist_g, hist_s = [], [], {k: [] for k in probes}
    fields = {}
    pend = sorted(field_instants)

    def record(t):
        hist_t.append(t)
        hist_g.append(boundaries[0].temperature(t) if boundaries else np.nan)
        for k, v in probes.items():
            hist_s[k].append(T[v].mean())

    record(0.0)
    t, nxt = 0.0, record_every
    n_steps = int(np.ceil(t_final / dt - 1e-9))
    stopped = None
    for step in range(n_steps):
        tn = min(t_final, t + dt)
        h = tn - t
        Te = T[conn].mean(axis=1)
        data = np.bincount(emap, weights=(M.k(mat, Te)[:, None, None] * Ke1).ravel(), minlength=nnz)
        if variable:                     # apparent capacity at the temperature of the previous step
            data_C, Cmat, C = capacity(M.rho_c(mat, Te))
        if C is None:
            data += data_C / h
            diagonal = np.zeros(n_nodes)
            b = (Cmat @ T) / h
        else:
            diagonal = C / h
            b = C / h * T
        TK = T + 273.15
        extra = None
        for i, g in enumerate(boundaries):
            Tg = g.temperature(tn)
            TgK = Tg + 273.15
            if g.faces is not None:
                fc, Mf = g.faces
                TKf = TK[fc].mean(1)
                hb = g.h_convection + g.emissivity * SIGMA * (TgK ** 2 + TKf ** 2) * (TgK + TKf)
                Hf = (hb[:, None, None] * Mf).ravel()
                H = sp.csr_matrix((Hf, pairs[i]), shape=(n_nodes, n_nodes))
                extra = H if extra is None else extra + H
                b = b + Tg * np.asarray(H.sum(axis=1)).ravel()
            else:
                hb = g.h_convection + g.emissivity * SIGMA * (TgK ** 2 + TK ** 2) * (TgK + TK)
                diagonal = diagonal + hb * g.node_area
                b = b + hb * g.node_area * Tg
        data[diag] += diagonal
        A = sp.csr_matrix((data, indices, indptr), shape=(n_nodes, n_nodes))
        if extra is not None:
            A = (A + extra).tocsr()
        T = _cg(A, b, T)
        t = tn
        while pend and pend[0] <= t + 1e-9:
            fields[pend.pop(0)] = T.copy()
        if t >= nxt - 1e-9 or t >= t_final - 1e-9:
            record(t)
            nxt += record_every
        if verbose and step % 50 == 0:
            print(f"  t = {t / 60:6.1f} min")
        if stop_if is not None and stop_if(t, T):
            stopped = t
            if hist_t[-1] != t:
                record(t)
            break
    return ThermalResult(t=np.array(hist_t), T_gas=np.array(hist_g),
                            probes={k: np.array(v) for k, v in hist_s.items()},
                            T_final=T, fields=fields, stopped_at=stopped)


# ----------------------------------------------------------------------------
# Boundaries and shortcuts for the meshes of the package
# ----------------------------------------------------------------------------
def face_selector(m: Mesh3D, faces, tol=1e-9):
    """Function that marks the element faces (nf, k, 3) lying on the requested
    planes: x0, x1, y0, y1 (sides), base and top. With symmetry, the plane
    y = W/2 is never exposed."""
    m2 = m.m2d
    L, Wm = m2.length, m2.model_width
    z0, z1 = m.zs[0], m.zs[-1]

    def sel(P):
        s = np.zeros(len(P), bool)
        if "x0" in faces: s |= (np.abs(P[:, :, 0]) < tol).all(1)
        if "x1" in faces: s |= (np.abs(P[:, :, 0] - L) < tol).all(1)
        if "y0" in faces: s |= (np.abs(P[:, :, 1]) < tol).all(1)
        if "y1" in faces and not m2.symmetry: s |= (np.abs(P[:, :, 1] - Wm) < tol).all(1)
        if "base" in faces: s |= (np.abs(P[:, :, 2] - z0) < tol).all(1)
        if "top" in faces: s |= (np.abs(P[:, :, 2] - z1) < tol).all(1)
        return s
    return sel


def side_areas_3d(m: Mesh3D, faces=("x0", "x1", "y0", "y1"), tol=1e-9):
    """Nodal area of the exposed side faces (x = 0, x = L, y = 0, y = W),
    integrated over the element faces (W6, H8 or W15; HRZ weights).
    With symmetry, the face y = W/2 is never exposed."""
    return m.elem.face_areas(m.xyz, m.conn3, face_selector(m, faces, tol))


def top_base_areas_3d(m: Mesh3D, which="both"):
    faces = {"base": ("base",), "top": ("top",)}.get(which, ("base", "top"))
    return m.elem.face_areas(m.xyz, m.conn3, face_selector(m, faces))


def boundary_lengths_2d(m2: Mesh2D, faces=("x0", "x1", "y0", "y1"), tol=1e-9):
    if m2.vertical:      # vertical section: the top and bottom edges are periodic cuts, never exposed
        faces = tuple(f for f in faces if f in ("y0", "y1"))
    air = m2.boundary_edges()
    p = m2.xy[air]
    L, Wm = m2.length, m2.model_width
    sel = np.zeros(len(air), bool)
    if "x0" in faces: sel |= (np.abs(p[:, :, 0]) < tol).all(1)
    if "x1" in faces: sel |= (np.abs(p[:, :, 0] - L) < tol).all(1)
    if "y0" in faces: sel |= (np.abs(p[:, :, 1]) < tol).all(1)
    if "y1" in faces and not m2.symmetry: sel |= (np.abs(p[:, :, 1] - Wm) < tol).all(1)
    air = air[sel]
    comp = np.linalg.norm(m2.xy[air[:, 0]] - m2.xy[air[:, 1]], axis=1)
    return np.bincount(air.ravel(), weights=np.repeat(comp / 2, 2), minlength=m2.n_nodes)


def section_material(m2: Mesh2D, geo=None):
    """Materials of the 2D section (course 0): air in the holes, grout in the grouted
    holes (if geo is given) and unit elsewhere. In the vertical section (horizontal
    holes), the half joints at top and bottom are mortar (or dry joint)."""
    if m2.vertical:
        joint = MORTAR if (geo is None or geo.joint_type == "mortar") else JOINT
        return np.where(m2.mark == 0, AIR, np.where(m2.mark > 0, joint, CERAMIC)).astype(np.int8)
    if geo is not None:
        from .mesh import course_states
        hole, grt = course_states(m2, geo, 0)
        return np.where(grt, GROUT, np.where(hole, AIR, CERAMIC)).astype(np.int8)
    return np.where(m2.mark == 0, AIR, CERAMIC).astype(np.int8)


def boundary_faces(m: Mesh3D, node_area: np.ndarray):
    """Element faces whose nodes all have area > 0 in node_area (the faces of the
    boundary that generated node_area) and their mass matrices, integral of N_i N_j dA."""
    tp = m.elem
    F, MM = [], []
    for face_nodes, emap, rule in tp.faces_int:
        fc = m.conn3[:, face_nodes]
        sel = (node_area[fc] > 0).all(1)
        if not sel.any():
            continue
        ce = m.conn3[sel]
        st, ws = rule
        xi, dxi_ds, dxi_dt = emap(st)
        N = tp.N(xi)[:, face_nodes]
        dN = tp.dN(xi)
        X = m.xyz[ce]
        a = np.einsum("eia,pib,pb->epa", X, dN, dxi_ds)
        b = np.einsum("eia,pib,pb->epa", X, dN, dxi_dt)
        dA = np.linalg.norm(np.cross(a, b), axis=2) * ws[None, :]
        F.append(ce[:, face_nodes]); MM.append(np.einsum("ep,pi,pj->eij", dA, N, N))
    if not F:
        return None
    return np.vstack(F), np.vstack(MM)


def solve_3d(m: Mesh3D, M: Materials, boundaries, t_final, **kw) -> ThermalResult:
    """Conduction in the 3D prism. W6 and H8: capacity and heat exchange lumped
    at the nodes. W15: consistent capacity and heat-exchange matrices (lumping
    at the nodes of quadratic faces distorts the incoming flux)."""
    if m.kind == "W6":
        Ke1, vol = el.conduction_w6(m.xyz, m.conn3)
        return solve_(m.conn3, Ke1, vol, m.mat, m.n_nodes, M, boundaries, t_final, **kw)
    tp = m.elem
    Ke1, vol, wts = tp.conduction(m.xyz, m.conn3)
    if m.kind == "H8":
        return solve_(m.conn3, Ke1, vol, m.mat, m.n_nodes, M, boundaries, t_final, cap_weights=wts, **kw)
    xi, w = tp.points()
    _, detJ = tp.gradients(m.xyz, m.conn3, xi)
    N = tp.N(xi)
    Ce = np.einsum("ep,pi,pj->eij", detJ * w[None, :], N, N)
    ctr = []
    for g in boundaries:
        if g.faces is None:
            g = Boundary(g.node_area, g.temperature, g.h_convection, g.emissivity, g.name,
                         boundary_faces(m, g.node_area))
        ctr.append(g)
    return solve_(m.conn3, Ke1, vol, m.mat, m.n_nodes, M, ctr, t_final, cap_elem=Ce, **kw)


def solve_2d(m2: Mesh2D, M: Materials, boundaries, t_final, mat=None, **kw) -> ThermalResult:
    Ke1, area = el.conduction_t3(m2.xy, m2.tri)
    mat = section_material(m2) if mat is None else mat
    return solve_(m2.tri, Ke1, area, mat, m2.n_nodes, M, boundaries, t_final, **kw)
