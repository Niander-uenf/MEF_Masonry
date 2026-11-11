"""Linear structural analysis with thermal strain.

3D: W6 wedges, H8 hexahedra or W15 wedges (mesh.kind); compressive load applied on the top face of the top
plate; base restrained vertically only (free to expand in its plane),
plus the minimum against rigid-body motion. The air of the holes is left out.

2D: generalised plane strain in the section (free body), for quick studies.
"""
from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .mesh import AIR, CERAMIC, JOINT, CAPPING, STEEL, MORTAR, GROUT, Mesh2D, Mesh3D, nearest_node
from .thermal import Materials
from . import elements as el


def _properties(mat, Te, M: Materials):
    E = np.empty(len(mat)); nu = np.empty(len(mat))
    for m in np.unique(mat):
        s = mat == m
        if m == CERAMIC:
            E[s] = M.ceramic.E_T(Te[s]); nu[s] = M.ceramic.nu
        elif m == MORTAR:
            E[s] = M.mortar.E_T(Te[s]); nu[s] = M.mortar.nu
        elif m == GROUT:
            E[s] = M.grout.E_T(Te[s]); nu[s] = M.grout.nu
        elif m == JOINT:
            E[s] = M.joint.E; nu[s] = M.joint.nu
        elif m == CAPPING:
            E[s] = M.capping.E; nu[s] = M.capping.nu
        elif m == STEEL:
            E[s] = M.steel.E; nu[s] = M.steel.nu
    return E, nu


def _thermal_strain_nodes(mat, T_el_nodes, M: Materials):
    """Free thermal strain at the nodes of each element (ne, nen), from the material law.
    Dry joints and capping follow the unit; mortar and steel use their own law."""
    eps = np.empty_like(T_el_nodes)
    for m in np.unique(mat):
        s = mat == m
        if m == STEEL:
            eps[s] = M.steel.alpha * (T_el_nodes[s] - M.ceramic.T_ref)
        elif m == MORTAR:
            eps[s] = M.mortar.thermal_strain(T_el_nodes[s])
        elif m == GROUT:
            eps[s] = M.grout.thermal_strain(T_el_nodes[s])
        else:
            eps[s] = M.ceramic.thermal_strain(T_el_nodes[s])
    return eps


def _factorize(K):
    """Factorises K (symmetric positive definite). Uses pypardiso if it is
    installed (faster, Windows/Linux); otherwise SuperLU without pivoting,
    with a minimum-degree ordering (safe for an SPD matrix)."""
    try:
        import pypardiso
        Kc = K.tocsr()
        return lambda b: pypardiso.spsolve(Kc, b)
    except ImportError:
        lu = spla.splu(K.tocsc(), permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                       options=dict(SymmetricMode=True))
        return lu.solve


@dataclass
class StructuralResult:
    u: np.ndarray            # (n_nodes, 3) m
    sig: np.ndarray          # (ne, 6) Pa, nan in the air
    princ: np.ndarray        # (ne, 3) Pa, ascending
    active: np.ndarray       # solid elements


def solve_3d(m: Mesh3D, M: Materials, T_nodes: Optional[np.ndarray] = None,
                load_N: float = 0.0, include_thermal: bool = True,
                plates: str = "rigid", units: Optional[int] = None,
                fixed_base: bool = False, load_region=None) -> dict:
    """Solves up to two cases with the same stiffness: 'thermal' (field T_nodes, top
    free) and 'load' (compression load_N). Returns {'thermal': R, 'load': R};
    the cases add up with combine().

    plates = 'rigid' (default): the steel plates leave the mechanical model and the
             press platens are rigid frictionless supports — base with u_z = 0;
             in the 'load' case, top with uniform u_z, adjusted so that the reaction
             equals the load. Fewer unknowns and no steel/ceramic contrast.
    plates = 'steel': the plates enter as steel and the load is a uniform pressure
             on the top face of the top plate (or of the ceramic, if the mesh has
             no plates).
    fixed_base: all displacements of the base fixed (fixed support, like the
             "fixed support" of commercial programs), instead of u_z only.
    load_region: with plates = 'steel', function (x, y) -> mask that selects the
             part of the top face that receives the pressure (centre of each face);
             None = the whole face.
    load_N is the total test load; with symmetry, half is applied."""
    n = m.n_nodes
    if T_nodes is None:
        T_nodes = np.full(n, M.ceramic.T_ref)
    active = m.mat != AIR
    if plates == "rigid":
        active &= m.mat != STEEL
    conn3 = m.conn3[active]; mat = m.mat[active]
    used = np.unique(conn3)
    Tel = T_nodes[conn3]
    E, nu = _properties(mat, Tel.mean(1), M)
    eps_n = _thermal_strain_nodes(mat, Tel, M)
    tp = m.elem
    nd = 3 * tp.nn
    w6 = m.kind == "W6"
    if units is None:
        units = max(500, int(20000 * 324 / nd ** 2))

    ndof = 3 * n
    K = sp.csr_matrix((ndof, ndof))
    F_th = np.zeros(ndof)
    for i0 in range(0, len(conn3), units):
        s = slice(i0, i0 + units)
        eth = eps_n[s] if include_thermal else None
        if w6:
            Ke, fe = el.stiffness_w6(m.xyz, conn3[s], E[s], nu[s], eth)
        else:
            Ke, fe = tp.stiffness(m.xyz, conn3[s], E[s], nu[s], eth)
        dofs = (3 * conn3[s][:, :, None] + np.arange(3)).reshape(len(Ke), nd)
        r = np.repeat(dofs, nd, axis=1).ravel(); c = np.tile(dofs, (1, nd)).ravel()
        K = K + sp.csr_matrix((Ke.ravel(), (r, c)), shape=(ndof, ndof))
        if fe is not None:
            F_th += np.bincount(dofs.ravel(), weights=fe.ravel(), minlength=ndof)

    zmin, zmax = m.xyz[used, 2].min(), m.xyz[used, 2].max()
    factor = 0.5 if m.m2d.symmetry else 1.0
    base = used[np.isclose(m.xyz[used, 2], zmin)]
    top_nodes = used[np.isclose(m.xyz[used, 2], zmax)]

    # common boundary conditions
    fixed_dofs = np.zeros(ndof, bool)
    inactive = np.setdiff1d(np.arange(n), used)
    fixed_dofs[(3 * inactive[:, None] + np.arange(3)).ravel()] = True
    fixed_dofs[3 * base + 2] = True
    if fixed_base:
        fixed_dofs[(3 * base[:, None] + np.arange(3)).ravel()] = True
    ref = base[np.argmin(np.linalg.norm(m.xyz[base, :2] - [m.m2d.length / 2, 0.0], axis=1))]
    fixed_dofs[3 * ref] = True
    if m.m2d.symmetry:
        sym = used[np.isclose(m.xyz[used, 1], m.m2d.model_width)]
        fixed_dofs[3 * sym + 1] = True
    else:
        fixed_dofs[3 * ref + 1] = True
        other = base[np.argmax(np.abs(m.xyz[base, 0] - m.xyz[ref, 0]))]
        fixed_dofs[3 * other + 1] = True

    results = {}
    cases = []
    if include_thermal:
        cases.append(("thermal", F_th, None))
    if load_N:
        if plates == "rigid":
            # top with uniform vertical displacement (rigid frictionless platen)
            u_p = np.zeros(ndof); u_p[3 * top_nodes + 2] = -1e-3
            cases.append(("load", -(K @ u_p), u_p))
        else:
            def _sel(P):
                s = np.isclose(P[:, :, 2], zmax).all(1)
                if load_region is not None:
                    c = P.mean(1)
                    s &= np.asarray(load_region(c[:, 0], c[:, 1]), bool)
                return s
            a = tp.face_areas(m.xyz, conn3, _sel, consistent=True)
            F_c = np.zeros(ndof)
            F_c[2::3] = -load_N * factor * a / a.sum()
            cases.append(("load", F_c, None))

    fac_cache = {}
    for name, F, u_p in cases:
        fx = fixed_dofs.copy()
        if u_p is not None:
            fx[3 * top_nodes + 2] = True
        key = fx.tobytes()
        free_dofs = np.nonzero(~fx)[0]
        if key not in fac_cache:
            fac_cache[key] = _factorize(K[free_dofs][:, free_dofs])
        u = np.zeros(ndof) if u_p is None else u_p.copy()
        u[free_dofs] = fac_cache[key](F[free_dofs])
        if u_p is not None:
            reaction = (K @ u)[3 * top_nodes + 2].sum()        # force of the platen on the prism
            u *= -load_N * factor / reaction
        u3 = u.reshape(n, 3)
        epsn = eps_n if name == "thermal" else np.zeros_like(eps_n)
        sig = np.full((len(m.conn3), 6), np.nan)
        s_at = np.empty((len(conn3), 6))
        for i0 in range(0, len(conn3), units):
            s = slice(i0, i0 + units)
            ue = u3[conn3[s]].reshape(-1, nd)
            s_at[s] = (el.stresses_w6(m.xyz, conn3[s], E[s], nu[s], ue, epsn[s]) if w6
                       else tp.stresses(m.xyz, conn3[s], E[s], nu[s], ue, epsn[s]))
        sig[active] = s_at
        pr = np.full((len(m.conn3), 3), np.nan)
        pr[active] = el.principal_stresses(s_at)
        results[name] = StructuralResult(u=u3, sig=sig, princ=pr, active=active)
    return results


def combine(r1: StructuralResult, r2: StructuralResult) -> StructuralResult:
    sig = r1.sig + r2.sig
    pr = np.full_like(r1.princ, np.nan)
    pr[r1.active] = el.principal_stresses(sig[r1.active])
    return StructuralResult(u=r1.u + r2.u, sig=sig, princ=pr, active=r1.active)


# ----------------------------------------------------------------------------
# Quantities to compare with the test
# ----------------------------------------------------------------------------
def mean_uz_plane(m: Mesh3D, u, z, nodes=None):
    nodes = np.arange(m.n_nodes) if nodes is None else nodes
    sel = nodes[np.isclose(m.xyz[nodes, 2], z)]
    return u[sel, 2].mean()


def shortening_between_plates(m: Mesh3D, R: StructuralResult, H: float):
    """Change in the distance between the plate faces (z = 0 and z = H), m.
    Negative = shortening."""
    used = np.unique(m.conn3[R.active])
    return mean_uz_plane(m, R.u, H, used) - mean_uz_plane(m, R.u, 0.0, used)


def extensometer_strain(m: Mesh3D, R: StructuralResult, z1: float, z2: float,
                            y_face: float = 0.0):
    """Mean vertical strain between heights z1 and z2 on the outer face y = y_face
    (like an extensometer or DIC in the middle third). Uses the nearest node planes."""
    z1 = m.zs[np.argmin(np.abs(m.zs - z1))]; z2 = m.zs[np.argmin(np.abs(m.zs - z2))]
    face = np.nonzero(np.isclose(m.xyz[:, 1], y_face))[0]
    return (mean_uz_plane(m, R.u, z2, face) - mean_uz_plane(m, R.u, z1, face)) / (z2 - z1), z1, z2


def max_tension(m: Mesh3D, R: StructuralResult, materials=(CERAMIC,), zmin=None, zmax=None,
                  percentile=99.0):
    """Maximum principal tensile stress (Pa) in the ceramic, with the percentile and the position."""
    sel = np.isin(m.mat, materials)
    if zmin is not None or zmax is not None:
        zc = m.xyz[m.conn3].mean(1)[:, 2]
        if zmin is not None: sel &= zc >= zmin
        if zmax is not None: sel &= zc <= zmax
    s1 = R.princ[sel, 2]
    i = np.nanargmax(s1)
    pos = m.xyz[m.conn3[np.nonzero(sel)[0][i]]].mean(0)
    return float(s1[i]), float(np.nanpercentile(s1, percentile)), pos


# ----------------------------------------------------------------------------
# 2D: generalised plane strain
# ----------------------------------------------------------------------------
def solve_section_gps(m2: Mesh2D, M: Materials, T_nodes: np.ndarray, mat=None):
    """Thermal stresses in the section (free body, no loads), generalised plane
    strain: linear axial strain in x and y chosen to cancel the normal
    force and moments. Full mesh (no symmetry) recommended.
    Returns a dict with sx, sy, txy, sz, s1 (max. principal) per element, Pa."""
    mat = np.where(m2.mark == 0, AIR, CERAMIC) if mat is None else mat
    active = mat != AIR
    tri = m2.tri[active]
    xy = m2.xy
    nn = m2.n_nodes
    Te = T_nodes[tri].mean(1)
    E, nu = _properties(np.asarray(mat)[active], Te, M)
    eth = _thermal_strain_nodes(np.asarray(mat)[active], T_nodes[tri], M).mean(1)
    A, dx, dy = el.geometry_t3(xy, tri)
    A = np.abs(A)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
    D = np.zeros((len(tri), 3, 3))
    D[:, 0, 0] = D[:, 1, 1] = lam + 2 * mu; D[:, 0, 1] = D[:, 1, 0] = lam; D[:, 2, 2] = mu
    B = np.zeros((len(tri), 3, 6))
    B[:, 0, 0::2] = dx; B[:, 1, 1::2] = dy; B[:, 2, 0::2] = dy; B[:, 2, 1::2] = dx
    Ke = np.einsum("eji,ejk,ekl->eil", B, D, B) * A[:, None, None]
    # initial stress in plane strain (ez = 0): sigma0 = (3 lam + 2 mu) e on the normal components
    eps0 = np.c_[eth, eth, np.zeros(len(tri))]
    s0 = np.einsum("eij,ej->ei", D, eps0) + np.c_[lam * eth, lam * eth, np.zeros(len(tri))]
    fe = np.einsum("eji,ej->ei", B, s0) * A[:, None]
    dofs = np.zeros((len(tri), 6), int); dofs[:, 0::2] = 2 * tri; dofs[:, 1::2] = 2 * tri + 1
    K = sp.coo_matrix((Ke.ravel(), (np.repeat(dofs, 6, 1).ravel(), np.tile(dofs, 6).ravel())),
                      shape=(2 * nn, 2 * nn)).tocsr()
    F = np.bincount(dofs.ravel(), weights=fe.ravel(), minlength=2 * nn)
    used = np.unique(tri)
    fixed_dofs = np.ones(2 * nn, bool); fixed_dofs[(2 * used[:, None] + np.arange(2)).ravel()] = False
    c = xy[used]
    n0 = used[np.argmin(np.hypot(c[:, 0], c[:, 1] - m2.width / 2))]
    n1 = used[np.argmin(np.hypot(c[:, 0] - m2.length, c[:, 1] - m2.width / 2))]
    fixed_dofs[[2 * n0, 2 * n0 + 1, 2 * n1 + 1]] = True
    free_dofs = np.nonzero(~fixed_dofs)[0]
    u = np.zeros(2 * nn)
    u[free_dofs] = spla.spsolve(K[free_dofs][:, free_dofs].tocsc(), F[free_dofs])
    eps = np.einsum("eij,ej->ei", B, u[dofs])
    sig = np.einsum("eij,ej->ei", D, eps) - s0
    sx, sy, txy = sig.T
    xc = xy[tri].mean(1)
    base = nu * (sx + sy) - E * eth
    W = (A * E)
    Mw = np.array([[W.sum(), (W * xc[:, 0]).sum(), (W * xc[:, 1]).sum()],
                   [(W * xc[:, 0]).sum(), (W * xc[:, 0] ** 2).sum(), (W * xc[:, 0] * xc[:, 1]).sum()],
                   [(W * xc[:, 1]).sum(), (W * xc[:, 0] * xc[:, 1]).sum(), (W * xc[:, 1] ** 2).sum()]])
    rhs = -np.array([(A * base).sum(), (A * base * xc[:, 0]).sum(), (A * base * xc[:, 1]).sum()])
    abc = np.linalg.solve(Mw, rhs)
    sz = base + E * (abc[0] + abc[1] * xc[:, 0] + abc[2] * xc[:, 1])
    mid = (sx + sy) / 2; radius = np.sqrt(((sx - sy) / 2) ** 2 + txy ** 2)
    s1 = np.maximum(mid + radius, sz)
    out = {k: np.full(len(m2.tri), np.nan) for k in ("sx", "sy", "txy", "sz", "s1")}
    for k, v in zip(("sx", "sy", "txy", "sz", "s1"), (sx, sy, txy, sz, s1)):
        out[k][active] = v
    out["area"] = np.zeros(len(m2.tri)); out["area"][active] = A
    out["active"] = active
    return out
