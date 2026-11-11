"""Non-linear analysis of the prism: damage in the unit and mortar, dry joints with friction.

Ceramic — scalar damage with two variables, tensile and compressive, and
crack closure (Mazars, 1984; Mazars and Pijaudier-Cabot, 1989):
    sigma = (1 - d) D (eps - eps_t),  1 - d = (1 - r d_t)(1 - d_c)
    r = sum of the tensile principal effective stresses / sum of their magnitudes
    (r = 1 in pure tension, 0 in pure compression): a tensile crack that
    becomes compressed closes and transmits compression again;
    d_t grows with the Rankine criterion, <sigma1_ef>/E, and d_c with the von
    Mises norm of the compressive principal stresses, sqrt(3 J2(<sigma_ef>-))/(k E),
    k = fc/ft (effective stresses
    sigma_ef = D (eps - eps_t); see Jirásek and Bažant, 2002); damage starts
    at ft in uniaxial tension and at fc in uniaxial compression;
    exponential softening: 1 - d_i = (k0/kappa_i) exp(-(kappa_i - k0)/ku_i);
    crack band regularisation (Bažant and Oh, 1983): the energy
    dissipated per volume is G/h, with h = element width normal
    to the crack (Oliver, 1989), fixed when each damage starts; G = G_f in
    tension and G_c in compression.

Elements: W6, H8 or W15 (mesh.kind). The state (damage, joint opening)
is evaluated at the centre of each element and holds for the whole element; the
stiffness is integrated with the full rule of the element.

Dry joints and capping — thin layer (Desai et al., 1984), orthotropic in the
joint system: normal stiffness E_joint in compression and almost zero in tension
(open joint); elastic shear up to the Coulomb limit mu |sigma_n|,
without cohesion, treated as a secant stiffness.

Solution: temperature increments (actual heating, field by field) and
then imposed displacement of the top platen. In each increment, secant
iteration; the linear systems are solved by conjugate gradients
preconditioned with the factorisation of a reference stiffness, redone when
convergence slows.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .mesh import AIR, CERAMIC, JOINT, CAPPING, STEEL, MORTAR, GROUT, Mesh3D
from .thermal import Materials
from .parameters import CeramicDamage, FrictionJoint, mortar_damage as _default_mortar_damage
from .parameters import grout_damage as _default_grout_damage
from . import elements as el
from .elements3d import kind as element_type
from .mechanical import _thermal_strain_nodes, _factorize


def damage_law(kappa, k0, ku):
    d = np.zeros_like(kappa)
    s = kappa > k0
    d[s] = 1 - (k0[s] / kappa[s]) * np.exp(-(kappa[s] - k0[s]) / ku[s])
    return d


def _tensor(v):
    S = np.empty(v.shape[:-1] + (3, 3))
    S[..., 0, 0], S[..., 1, 1], S[..., 2, 2] = v[..., 0], v[..., 1], v[..., 2]
    S[..., 0, 1] = S[..., 1, 0] = v[..., 3]
    S[..., 1, 2] = S[..., 2, 1] = v[..., 4]
    S[..., 0, 2] = S[..., 2, 0] = v[..., 5]
    return S


@dataclass
class Step:
    phase: str
    label: float           # time (s) during heating; platen displacement (m) during loading
    load_N: float
    max_damage: float
    damaged_vol: float   # fraction of the ceramic volume with max(d_t, d_c) > 0.5 (cracked or crushed)
    open_joints: float   # fraction of open joints
    iterations: int
    converged: bool
    contact: float = float("nan")    # fraction of the bed-joint area in contact
    eps_ext: float = float("nan")    # strain of the virtual extensometer (if defined)


class NonlinearModel:
    """3D prism with damage in the ceramic and non-linear dry joints.

    support = 'platen': base with u_z = 0 and the minimum against rigid-body motion (frictionless platens)
    support = 'fixed_end': base with all three displacements fixed (for verification tests)
    """

    def __init__(self, m: Mesh3D, M: Materials, damage: CeramicDamage = None, joint: FrictionJoint = None,
                 support: str = "platen", ft_factor_elem: Optional[np.ndarray] = None, verbose: bool = True,
                 mortar_damage: Optional[CeramicDamage] = None, grout_damage: Optional[CeramicDamage] = None):
        self.m, self.M = m, M
        self.Pd = damage or CeramicDamage()
        self.Pa = mortar_damage or _default_mortar_damage()
        self.Pg = grout_damage or _default_grout_damage()
        self.Pj = joint or FrictionJoint()
        self.verbose = verbose
        is_active = (m.mat != AIR) & (m.mat != STEEL)
        self.id_glob = np.nonzero(is_active)[0]
        self.conn3 = m.conn3[is_active]
        self.mat = m.mat[is_active]
        # elements with damage: unit, mortar and grout (each with its own parameters)
        self.ic = np.nonzero(np.isin(self.mat, (CERAMIC, MORTAR, GROUT)))[0]
        self.mor = self.mat[self.ic] == MORTAR
        self.grt_ = self.mat[self.ic] == GROUT
        self.nu_c = self._by_material(M.ceramic.nu, M.mortar.nu, M.grout.nu)
        self.du_t = self._by_material(self.Pd.du_t, self.Pa.du_t, self.Pg.du_t)
        self.du_c = self._by_material(self.Pd.du_c, self.Pa.du_c, self.Pg.du_c)
        self.ij = np.nonzero(np.isin(self.mat, (JOINT, CAPPING)))[0]
        kinds = np.array([m.layers[c].kind for c in m.layer[is_active]])
        self.normal_x = kinds[self.ij] == "unit"          # head joints (normal x)
        Ec = self.Pj.E_contact
        self.Ej = np.where(self.mat[self.ij] == JOINT, Ec or M.joint.E, Ec or M.capping.E).astype(float)
        self.nuj = np.where(self.mat[self.ij] == JOINT, M.joint.nu, M.capping.nu)
        self.ft_factor = np.ones(len(self.ic)) if ft_factor_elem is None else ft_factor_elem[self.id_glob[self.ic]]
        n = m.n_nodes
        self.n, self.ndof = n, 3 * n
        tp = m.elem
        self.tp = tp
        self.nd = nd = 3 * tp.nn
        w6 = m.kind == "W6"
        self.dofs = (3 * self.conn3[:, :, None] + np.arange(3)).reshape(-1, nd)
        # estrutura esparsa
        r = np.repeat(self.dofs, nd, axis=1).ravel().astype(np.int64)
        c = np.tile(self.dofs, (1, nd)).ravel().astype(np.int64)
        key = r * self.ndof + c
        uniq, self.emap = np.unique(key, return_inverse=True)
        self.emap = self.emap.ravel()
        self.indices = (uniq % self.ndof).astype(np.int32)
        self.indptr = np.searchsorted(uniq // self.ndof, np.arange(self.ndof + 1)).astype(np.int64)
        self.nnz = len(uniq)
        del r, c, key, uniq
        # geometry
        nu_c = self.nu_c
        cer = self.conn3[self.ic]
        if w6:
            A, dz, dLdx, dLdy = el.geometry_w6(m.xyz, self.conn3)
            self.vol = np.abs(A) * dz
            self.Bc, self.Nc = el.B_centroid_w6(m.xyz, self.conn3)
            Gc, _ = el._gradients_w6(A[self.ic], dz[self.ic], dLdx[self.ic], dLdy[self.ic], np.full(3, 1 / 3), 0.0)
            self.Ke1 = el.stiffness_w6(m.xyz, cer, np.ones(len(self.ic)), nu_c)[0]
        else:
            self.vol = tp.volumes(m.xyz, self.conn3)
            self.Bc, self.Nc, _ = tp.at_centre(m.xyz, self.conn3)
            # band width (Oliver): linear gradients — the W15 wedge uses the W6 of its corners
            lin = element_type("W6") if m.kind == "W15" else tp
            Gc = lin.at_centre(m.xyz, cer[:, :lin.nn])[2]
            self.Ke1 = np.empty((len(self.ic), nd, nd))
            bl = max(500, int(20000 * 324 / nd ** 2))
            for i0 in range(0, len(self.ic), bl):
                s_ = slice(i0, i0 + bl)
                self.Ke1[s_] = tp.stiffness(m.xyz, cer[s_], np.ones(len(cer[s_])), nu_c[s_])[0]
        self.Gc_cer = Gc                                      # for the band width (Oliver)
        # joints: integrals of the gradient products, Gab[e, a, b] = sum_p w_p g_a g_b^T
        # (nj, 3, 3, nn, nn), to assemble the stiffness with diagonal D; and, per point, w (nj,) and G (nj, nn, 3)
        cj = self.conn3[self.ij]
        xi, wq = tp.points()
        G, detJ = tp.gradients(m.xyz, cj, xi)              # (nj, np, nn, 3), (nj, np)
        wv = detJ * wq[None, :]
        self.Gab = np.einsum("ep,epia,epjb->eabij", wv, G, G)
        self._wGj = (wv, G, tp.N(xi))
        self.Bint = np.einsum("ep,epia->eai", wv, G)          # integral of the gradients (nj, 3, nn)
        # boundary
        used = np.unique(self.conn3)
        self.used = used
        z = m.xyz[:, 2]
        self.zmin, self.zmax = z[used].min(), z[used].max()
        self.base = used[np.isclose(z[used], self.zmin)]
        self.top = used[np.isclose(z[used], self.zmax)]
        fixed_dofs = np.zeros(self.ndof, bool)
        inact = np.setdiff1d(np.arange(n), used)
        fixed_dofs[(3 * inact[:, None] + np.arange(3)).ravel()] = True
        if support == "fixed_end":
            fixed_dofs[(3 * self.base[:, None] + np.arange(3)).ravel()] = True
        else:
            fixed_dofs[3 * self.base + 2] = True
            ref = self.base[np.argmin(np.linalg.norm(m.xyz[self.base, :2] - [m.m2d.length / 2, 0.0], axis=1))]
            fixed_dofs[3 * ref] = True
            if m.m2d.symmetry:
                sym = used[np.isclose(m.xyz[used, 1], m.m2d.model_width)]
                fixed_dofs[3 * sym + 1] = True
            else:
                fixed_dofs[3 * ref + 1] = True
                other = self.base[np.argmax(np.abs(m.xyz[self.base, 0] - m.xyz[ref, 0]))]
                fixed_dofs[3 * other + 1] = True
        self.base_fixed = fixed_dofs
        self.sym_factor = 2.0 if m.m2d.symmetry else 1.0
        # tributary area of the top nodes (seating load)
        zt = self.zmax
        self.top_area = tp.face_areas(m.xyz, self.conn3, lambda P: np.isclose(P[:, :, 2], zt).all(1),
                                        consistent=True)
        # bed joints and initial gaps (irregularity of the faces)
        self.dzj = tp.height(m.xyz, cj)
        self.area_j = self.vol[self.ij] / self.dzj
        self.bed = (self.mat[self.ij] == JOINT) & ~self.normal_x
        self.gap = np.zeros(len(self.ij)); self.eg = np.zeros(len(self.ij))
        if self.Pj.gap_std > 0:
            self.set_gaps(self.Pj.gap_std, self.Pj.corr_length, self.Pj.initial_contact,
                                self.Pj.seed, set_state=False)
        self._en_j = np.zeros(len(self.ij))
        self._ext = None
        # state
        nc, nj = len(self.ic), len(self.ij)
        self.kappa = np.zeros((nc, 2))            # history [tension, compression]
        self.h = np.full((nc, 2), np.nan)         # band width [tension, compression]
        self.r = np.ones(nc)                      # Mazars tension weight
        self.d = np.zeros(nc); self.dt = np.zeros(nc); self.dc = np.zeros(nc)
        self.closed = self.eg <= 0
        self.En = np.where(self.closed, self.Ej, self.Pj.open_factor * self.Ej)
        self.Gs = self.Ej / (2 * (1 + self.nuj))
        self.u = np.zeros(self.ndof)
        # transient creep (load-induced thermal strain), off until set_creep
        self.creep = None
        self.eps_tr = np.zeros((nc, 6))           # Voigt, engineering shear strains
        self.phi_hist = np.zeros(nc)              # largest phi(T) reached so far (only the first heating counts)
        self.ftr1 = np.zeros((nc, nd))
        self.tol_tr = 2e-4                        # largest increment of transient strain per substep
        self._Gint = None
        self.set_temperature(np.full(n, M.ceramic.T_ref))
        self._pc = None; self._pc_free = None
        self.cg_max, self.cg_res_ok = 12, 1e-3   # inexact linear solution (see _solve_linear)
        self.tol_u, self.tol_d = 2e-5, 1e-3        # convergence: relative change of u and maximum change of damage
        self.history: List[Step] = []
        self.instantaneous = None   # list -> stores u, damage and T of each increment (post-processor)

    # ------------------------------------------------------------------
    def set_gaps(self, std, corr_length=0.04, initial_contact=0.20, seed=1,
                       set_state=True):
        """Initial gaps in the bed joints: the height of the face is a Gaussian field
        (standard deviation std, correlation exp(-r²/l²) with l = corr_length),
        independent in each joint. The fraction initial_contact of the area (the highest
        points) touches without load; elsewhere, the gap is the height difference."""
        from scipy.ndimage import gaussian_filter
        from scipy.interpolate import RegularGridInterpolator
        rng = np.random.default_rng(seed)
        m = self.m
        xy = m.xyz[self.conn3[self.ij]].mean(1)[:, :2]
        layer = m.layer[self.id_glob[self.ij]]
        step = min(corr_length / 4, 0.005)
        xs = np.arange(-2 * step, m.m2d.length + 3 * step, step)
        ys = np.arange(-2 * step, m.m2d.width + 3 * step, step)
        self.gap[:] = 0.0
        for c in np.unique(layer[self.bed]):
            s = self.bed & (layer == c)
            field = gaussian_filter(rng.standard_normal((len(xs), len(ys))), corr_length / step / 2, mode="wrap")
            field = (field - field.mean()) / field.std()
            z = std * RegularGridInterpolator((xs, ys), field)(xy[s])
            order = np.argsort(-z)
            cum = np.cumsum(self.area_j[s][order]) / self.area_j[s].sum()
            z_ref = z[order][min(np.searchsorted(cum, initial_contact), len(z) - 1)]
            self.gap[s] = np.maximum(z_ref - z, 0.0)
        self.eg = self.gap / self.dzj
        if set_state:
            self.closed = self.eg <= 0
            self.En = np.where(self.closed, self.Ej, self.Pj.open_factor * self.Ej)

    def set_extensometer(self, z1, z2, face_y=0.0):
        """Virtual extensometer: mean strain between heights z1 and z2 on the face y = face_y."""
        m = self.m
        zs = m.zs
        z1, z2 = zs[np.argmin(np.abs(zs - z1))], zs[np.argmin(np.abs(zs - z2))]
        face = np.nonzero(np.isclose(m.xyz[:, 1], face_y))[0]
        self._ext = (face[np.isclose(m.xyz[face, 2], z1)], face[np.isclose(m.xyz[face, 2], z2)], z2 - z1)

    def extensometer_strain(self):
        if self._ext is None:
            return float("nan")
        n1, n2, L = self._ext
        u = self.u.reshape(-1, 3)
        return float((u[n2, 2].mean() - u[n1, 2].mean()) / L)

    def _by_material(self, unit, mortar, grout):
        """Value per damaged element, chosen by material (scalars or vectors)."""
        return np.where(self.grt_, grout, np.where(self.mor, mortar, unit))

    def set_temperature(self, T_nodes):
        self.T = T_nodes
        Tel = T_nodes[self.conn3]
        self.eps_n = _thermal_strain_nodes(self.mat, Tel, self.M)          # (ne, 6)
        self.eps_c = self.eps_n @ self.Nc                                # at the centroid
        Tc = Tel[self.ic].mean(1)
        pm, M = self._by_material, self.M
        self.Ecer = pm(M.ceramic.E_T(Tc), M.mortar.E_T(Tc), M.grout.E_T(Tc))
        self.ft = pm(self.Pd.ft_T(Tc), self.Pa.ft_T(Tc), self.Pg.ft_T(Tc)) * self.ft_factor
        self.fc = pm(self.Pd.fc_T(Tc), self.Pa.fc_T(Tc), self.Pg.fc_T(Tc))
        ncer = len(self.ic)
        if self.m.kind == "W6":
            self.fth1 = el.thermal_force_w6(self.m.xyz, self.conn3[self.ic], np.ones(ncer),
                                             self.nu_c, self.eps_n[self.ic])
        else:
            self.fth1 = self.tp.thermal_force(self.m.xyz, self.conn3[self.ic], np.ones(ncer),
                                              self.nu_c, self.eps_n[self.ic])
        # joints: unit thermal vectors of the normal terms of D (k = 0, 1, 2):
        # Tj[e, k] = sum_p w_p (eps_t . N_p) B_k,p, with B_k only in component k
        wv, G, Nq = self._wGj
        s_p = wv * (self.eps_n[self.ij] @ Nq.T)                           # (nj, np)
        Tk = np.einsum("ep,epia->eai", s_p, G)                            # (nj, 3, nn)
        nj, nn = len(self.ij), self.tp.nn
        self.Tj = np.zeros((nj, 3, nn, 3))
        for k in range(3):
            self.Tj[:, k, :, k] = Tk[:, k]
        self.Tj = self.Tj.reshape(nj, 3, 3 * nn)

    def _D_joint(self, En, Gs):
        """Diagonal of the constitutive matrix of the joints in the global system (nj, 6)."""
        ein = self.Pj.plane_factor * self.Ej
        x = self.normal_x
        Dd = np.empty((len(self.ij), 6))
        Dd[:, 0] = np.where(x, En, ein); Dd[:, 1] = ein; Dd[:, 2] = np.where(x, ein, En)
        Dd[:, 3] = np.where(x, Gs, ein / 2); Dd[:, 4] = np.where(x, ein / 2, Gs); Dd[:, 5] = Gs
        return Dd

    def _joint_stiffness(self, Dd):
        """Stiffness of the joints (nj, nd, nd) with diagonal D (Voigt xx, yy, zz, xy, yz, xz)."""
        G = self.Gab
        D0, D1, D2, D3, D4, D5 = (Dd[:, k, None, None] for k in range(6))
        nj, nn = len(self.ij), self.tp.nn
        K = np.empty((nj, nn, 3, nn, 3))
        K[:, :, 0, :, 0] = D0 * G[:, 0, 0] + D3 * G[:, 1, 1] + D5 * G[:, 2, 2]
        K[:, :, 1, :, 1] = D1 * G[:, 1, 1] + D3 * G[:, 0, 0] + D4 * G[:, 2, 2]
        K[:, :, 2, :, 2] = D2 * G[:, 2, 2] + D4 * G[:, 1, 1] + D5 * G[:, 0, 0]
        K[:, :, 0, :, 1] = D3 * G[:, 1, 0]; K[:, :, 1, :, 0] = D3 * G[:, 0, 1]
        K[:, :, 1, :, 2] = D4 * G[:, 2, 1]; K[:, :, 2, :, 1] = D4 * G[:, 1, 2]
        K[:, :, 0, :, 2] = D5 * G[:, 2, 0]; K[:, :, 2, :, 0] = D5 * G[:, 0, 2]
        return K.reshape(nj, 3 * nn, 3 * nn)

    def _assemble(self, d, En, Gs, clo=None):
        clo = self.closed if clo is None else clo
        nd2 = self.nd ** 2
        vals = np.empty((len(self.conn3), nd2))
        vals[self.ic] = ((self.Ecer * (1 - d))[:, None, None] * self.Ke1).reshape(-1, nd2)
        Dd = self._D_joint(En, Gs)
        vals[self.ij] = self._joint_stiffness(Dd).reshape(-1, nd2)
        fj = np.einsum("ek,eki->ei", Dd[:, :3], self.Tj)
        # gap as an initial strain -eg in the normal direction of the closed bed joints:
        # sigma_n = E_j (eps_n + eg)
        fg = np.where(clo, -self.eg, 0.0) * Dd[:, 2]
        if np.any(fg):
            bz = np.zeros((len(self.ij), self.tp.nn, 3)); bz[:, :, 2] = self.Bint[:, 2]
            fj += fg[:, None] * bz.reshape(len(self.ij), -1)
        data = np.bincount(self.emap, weights=vals.ravel(), minlength=self.nnz)
        K = sp.csr_matrix((data, self.indices, self.indptr), shape=(self.ndof, self.ndof))
        F = np.zeros(self.ndof)
        fc = (self.fth1 + self.ftr1) * (self.Ecer * (1 - d))[:, None]
        F += np.bincount(self.dofs[self.ic].ravel(), weights=fc.ravel(), minlength=self.ndof)
        F += np.bincount(self.dofs[self.ij].ravel(), weights=fj.ravel(), minlength=self.ndof)
        return K, F

    def _strains(self, u):
        e = np.einsum("eij,ej->ei", self.Bc, u[self.dofs])
        e[:, :3] -= self.eps_c[:, None]
        if self.creep is not None:
            e[self.ic] -= self.eps_tr
        return e                                                   # mechanical strain

    # ------------------------------------------------------------------ transient creep
    def set_creep(self, k_unit=0.0, k_mortar=0.0, k_grout=0.0, table_unit=None):
        """Switches on transient creep (load-induced thermal strain).

        For each material, phi(T) = k ε_th(T) (Anderberg and Thelandersson, 1976)
        or, for the unit, a table [(T °C, phi)] measured in the transient test. On
        first heating under compression, the increment is

            Δε_tr = Δphi / fc20 · [(1 + ν) σ⁻ - ν tr(σ⁻) I]

        with σ⁻ the compressive part of the nominal stress, fc20 the compressive
        strength of the material at 20 °C and ν its Poisson's ratio (uniaxially,
        Δε_tr = Δphi σ / fc20 in the load direction). Irreversible: it does not recover on
        cooling or on reheating below the maximum temperature already reached."""
        self.creep = {"k": self._by_material(k_unit, k_mortar, k_grout),
                         "table": sorted(table_unit) if table_unit else None,
                         "fc20": self._by_material(self.Pd.fc, self.Pa.fc, self.Pg.fc)}
        self.phi_hist = self._phi()

    def _phi(self):
        """phi(T) per damaged element, at the current temperature state."""
        fl = self.creep
        phi = fl["k"] * np.maximum(self.eps_c[self.ic], 0.0)
        if fl["table"]:
            Ts, vs = zip(*fl["table"])
            Tc = self.T[self.conn3[self.ic]].mean(1)
            unit = ~(self.mor | self.grt_)
            phi = np.where(unit, np.interp(Tc, Ts, vs), phi)
        return phi

    def _gradient_integral(self):
        """∫ ∂N/∂x dV per damaged element (nc, 3, nn): force of a constant stress."""
        if self._Gint is None:
            xi, wq = self.tp.points()
            cer = self.conn3[self.ic]
            out = np.empty((len(cer), 3, self.tp.nn))
            for i0 in range(0, len(cer), 4000):
                s_ = slice(i0, i0 + 4000)
                G, detJ = self.tp.gradients(self.m.xyz, cer[s_], xi)
                out[s_] = np.einsum("ep,epia->eai", detJ * wq[None, :], G)
            self._Gint = out
        return self._Gint

    def _transient_force(self):
        """Unit force vector (E = 1) of the transient strain: ∫ Bᵀ D(1, ν) ε_tr dV."""
        S = _tensor(np.einsum("eij,ej->ei", el.matrix_D(np.ones(len(self.ic)), self.nu_c), self.eps_tr))
        self.ftr1 = np.einsum("eai,eab->eib", self._gradient_integral(), S).reshape(len(self.ic), -1)

    def nominal_stress(self):
        """Nominal stress (1 - d) D e in the damaged elements (nc, 6)."""
        e = self._strains(self.u)[self.ic]
        return (1 - self.d)[:, None] * np.einsum("eij,ej->ei", el.matrix_D(self.Ecer, self.nu_c), e)

    def _transient_increment(self, sig, phi_new):
        """Δε_tr (nc, 6) for the stress sig and the new phi (explicit in the step)."""
        dphi = np.maximum(phi_new - self.phi_hist, 0.0)
        if not np.any(dphi):
            return np.zeros_like(self.eps_tr)
        w, V = np.linalg.eigh(_tensor(sig))
        wn = np.minimum(w, 0.0)
        Sn = np.einsum("eik,ek,ejk->eij", V, wn, V)                       # compressive part
        nu = self.nu_c[:, None, None]
        tr = wn.sum(1)[:, None, None]
        D = (dphi / self.creep["fc20"])[:, None, None] * ((1 + nu) * Sn - nu * tr * np.eye(3)[None])
        return np.column_stack([D[:, 0, 0], D[:, 1, 1], D[:, 2, 2], 2 * D[:, 0, 1], 2 * D[:, 1, 2], 2 * D[:, 0, 2]])

    def _update_ceramic(self, e, kappa_hist):
        nu = self.nu_c
        E, ft, fc = self.Ecer, self.ft, self.fc
        k = fc / ft
        k0 = ft / E
        sig = np.einsum("eij,ej->ei", el.matrix_D(E, nu), e)          # effective stress
        w, V = np.linalg.eigh(_tensor(sig))
        wn = np.minimum(w, 0.0)                                           # compressive part
        J2 = ((wn[:, 0] - wn[:, 1]) ** 2 + (wn[:, 1] - wn[:, 2]) ** 2 + (wn[:, 2] - wn[:, 0]) ** 2) / 6
        et = np.maximum(w[:, 2], 0.0) / E                                 # Rankine (tension)
        ec = np.sqrt(3 * J2) / (k * E)                                    # von Mises of the compressive stresses
        kappa = np.maximum(kappa_hist, np.column_stack([et, ec]))
        h = self.h.copy()
        for j, col in ((0, 2), (1, 0)):          # direction normal to the crack: sigma1 (tension), sigma3 (compression)
            new = np.isnan(h[:, j]) & (kappa[:, j] > k0)
            if new.any():
                proj = np.abs(np.einsum("eik,ek->ei", self.Gc_cer[new], V[new, :, col])).sum(1)
                h[new, j] = 2.0 / np.maximum(proj, 1e-12)
        hh = np.where(np.isnan(h), (self.vol[self.ic] ** (1 / 3))[:, None], h)
        ku_t = np.maximum(self.du_t / hh[:, 0] - k0 / 2, 0.05 * k0)        # G_f/(h ft) - k0/2, G_f = du_t ft
        ku_c = np.maximum(self.du_c / (hh[:, 1] * k) - k0 / 2, 0.05 * k0)  # G_c/(h k fc) - k0/2, G_c = du_c fc
        dt = damage_law(kappa[:, 0], k0, ku_t)
        dc = damage_law(kappa[:, 1], k0, ku_c)
        if self.Pd.crack_closure:
            sum_ = np.abs(w).sum(1)
            r = np.where(sum_ > 0, np.maximum(w, 0.0).sum(1) / np.maximum(sum_, 1e-300), 1.0)
        else:
            r = np.ones(len(dt))
        d = np.minimum(1 - (1 - r * dt) * (1 - dc), self.Pd.d_max)
        return kappa, d, h, r, dt, dc

    def _update_joints(self, e, closed, force=None):
        ts = np.where(self.normal_x, e[self.ij, 0], e[self.ij, 2])
        g1 = np.where(self.normal_x, e[self.ij, 3], e[self.ij, 4])
        g2 = e[self.ij, 5]
        gam = np.hypot(g1, g2)
        tol = self.Pj.eps_tol
        eg = self.eg                                   # initial gap / thickness
        closed = np.where(ts < -eg - tol, True, np.where(ts > -eg + tol, False, closed))   # hysteresis
        if force is not None:          # a joint that oscillates between open and closed stays open
            closed = closed & ~force
        # closed joint: sigma_n = E_j (en + eg), the gap enters as an initial strain (_assemble)
        En = np.where(closed, self.Ej, self.Pj.open_factor * self.Ej)
        self._en_j = ts
        sn = np.where(closed, self.Ej * (ts + eg), En * ts)
        G_el = self.Ej / (2 * (1 + self.nuj))
        lim = self.Pj.mu * np.maximum(-sn, 0.0)
        Gs = np.where(G_el * gam <= lim, G_el, np.maximum(lim / np.maximum(gam, 1e-30), self.Pj.G_min_factor * G_el))
        return En, Gs, closed

    def _solve_linear(self, K, rhs, free_dofs, x0):
        """Conjugate gradients preconditioned with the factorisation of an earlier
        stiffness, with few iterations (inexact solution: the outer secant
        iteration continues from it). Refactorises when the relative residual
        stays above cg_res_ok."""
        Kff = K[free_dofs][:, free_dofs].tocsr()
        if self._pc is None or self._pc_free is None or not np.array_equal(self._pc_free, free_dofs):
            self._pc = _factorize(Kff); self._pc_free = free_dofs
            return self._pc(rhs), 0
        it = [0]
        Mop = spla.LinearOperator(Kff.shape, matvec=self._pc)
        try:
            x, info = spla.cg(Kff, rhs, x0=x0, rtol=1e-8, atol=0.0, M=Mop, maxiter=self.cg_max,
                              callback=lambda xk: it.__setitem__(0, it[0] + 1))
        except TypeError:
            x, info = spla.cg(Kff, rhs, x0=x0, tol=1e-8, M=Mop, maxiter=self.cg_max,
                              callback=lambda xk: it.__setitem__(0, it[0] + 1))
        if info != 0:
            res = np.linalg.norm(Kff @ x - rhs) / max(np.linalg.norm(rhs), 1e-300)
            if res > self.cg_res_ok:
                self._pc = _factorize(Kff)
                x = self._pc(rhs)
        return x, it[0]

    def _step(self, presc, values_, F_ext, max_it=60, tol_u=None, tol_d=None):
        tol_u = self.tol_u if tol_u is None else tol_u
        tol_d = self.tol_d if tol_d is None else tol_d
        fixed_dofs = self.base_fixed.copy()
        fixed_dofs[presc] = True
        free_dofs = np.nonzero(~fixed_dofs)[0]
        d, En, Gs, clo = self.d.copy(), self.En.copy(), self.Gs.copy(), self.closed.copy()
        u = self.u.copy(); u[presc] = values_
        tol_change = max(2, int(1e-3 * len(self.ij)))
        kappa_it = self.kappa.copy()
        swaps = np.zeros(len(self.ij), int)
        conv = False
        for it in range(1, max_it + 1):
            K, Fth = self._assemble(d, En, Gs, clo)
            F = Fth + F_ext
            ub = np.zeros(self.ndof); ub[presc] = values_
            rhs = F[free_dofs] - (K @ ub)[free_dofs]
            x, _ = self._solve_linear(K, rhs, free_dofs, u[free_dofs])
            u_new = ub; u_new[free_dofs] = x
            e = self._strains(u_new)
            # history monotonic also between iterations: avoids loading/unloading oscillation
            kappa, d_new, h, r_m, dt, dc = self._update_ceramic(e[self.ic], kappa_it)
            kappa_it = kappa
            if it > 6:   # damping of the crack closure (weight r)
                d_new = np.where(d_new < d, 0.5 * (d_new + d), d_new)
            En_n, Gs_n, clo_n = self._update_joints(e, clo)
            # a joint that changes state 3 times in the same increment stays open (avoids cycling; near the switch, sigma_n ~ 0)
            swaps += clo_n != clo
            lock = swaps >= 3
            if lock.any():
                En_n, Gs_n, clo_n = self._update_joints(e, clo, force=lock)
            if it > 4:   # damping of the friction to avoid oscillation
                Gs_n = 0.5 * (Gs_n + Gs)
            du = np.linalg.norm(u_new - u) / max(np.linalg.norm(u_new), 1e-30)
            dd = np.max(np.abs(d_new - d)) if len(d) else 0.0
            changed = np.count_nonzero(clo_n != clo)
            u, d, En, Gs, clo = u_new, d_new, En_n, Gs_n, clo_n
            if du < tol_u and dd < tol_d and changed <= tol_change:
                conv = True
                break
        self.convergence = (du, dd, changed)
        # commits the state
        self.u, self.d, self.kappa, self.h, self.En, self.Gs = u, d, kappa, h, En, Gs
        self.r, self.dt, self.dc, self.closed = r_m, dt, dc, clo
        K, Fth = self._assemble(d, En, Gs, clo)
        r = K @ u - Fth - F_ext
        load = -r[3 * self.top + 2].sum() * self.sym_factor
        return load, it, conv

    def _record(self, phase, label, load, it, conv):
        vc = self.vol[self.ic]
        lt = self.bed
        cont = float(self.area_j[lt][self._en_j[lt] < -self.eg[lt]].sum() / self.area_j[lt].sum()) if lt.any() else float("nan")
        dm = np.maximum(self.dt, self.dc)
        p = Step(phase, label, load, float(dm.max()) if len(dm) else 0.0,
                  float(vc[dm > 0.5].sum() / vc.sum()) if len(vc) else 0.0,
                  float(np.mean(~self.closed)) if len(self.En) else 0.0, it, conv,
                  cont, self.extensometer_strain())
        self.history.append(p)
        if self.instantaneous is not None:       # states for the post-processor
            f32 = lambda v: np.asarray(v, np.float32).copy()
            self.instantaneous.append({"phase": phase, "label": float(label), "load_N": float(load),
                                      "u": f32(self.u), "dt": f32(self.dt), "dc": f32(self.dc), "d": f32(self.d),
                                      "T": f32(self.T) if phase == "heating" else None})
        if self.verbose:
            lbl = f"t = {label / 60:6.1f} min" if phase == "heating" else f"delta = {label * 1e3:6.3f} mm"
            print(f"  {phase:<11s} {lbl}  load {load / 1e3:8.2f} kN  max damage {p.max_damage:.3f}  "
                  f"vol(d>0,5) {100 * p.damaged_vol:5.2f}%  "
                  + (f"contact {100 * cont:5.1f}%  " if cont == cont else "no bed joints  ")
                  + f"it {it}" +("" if conv else " (not converged: du %.1e, dd %.1e, joints %d)" % self.convergence),
                  flush=True)
        return p

    # ------------------------------------------------------------------
    def heat(self, fields: Dict[float, np.ndarray], pre_load_N: Optional[float] = None,
                top_fixed: bool = False):
        """Applies the temperature fields {time: T_nodes} in sequence, with the top
        free under the seating load (or, with top_fixed, with the vertical
        displacement of the top locked at its current value: restrained heating).

        With transient creep on, each interval between fields is
        crossed in substeps (interpolated temperature) with an adaptive step,
        so that the increment of transient strain does not exceed tol_tr.
        Each substep is a predictor (stress at the start) and a corrector (mean
        stress of the substep, trapezoidal rule)."""
        P0 = self.Pj.pre_load_N if pre_load_N is None else pre_load_N
        F_ext = np.zeros(self.ndof)
        if top_fixed:
            presc = 3 * self.top + 2
            values_ = self.u[presc].copy()
        else:
            presc, values_ = np.zeros(0, int), np.zeros(0)
            F_ext[3 * np.arange(self.n) + 2] = -P0 / self.sym_factor * self.top_area / self.top_area.sum()
        for t in sorted(fields):
            T_b = np.asarray(fields[t], float)
            if self.creep is None:
                self.set_temperature(T_b)
                load, it, conv = self._step(presc, values_, F_ext)
            else:
                load, it, conv = self._heat_with_creep(T_b, presc, values_, F_ext)
            self._record("heating", t, load if top_fixed else P0, it, conv)
        self.d_heating = np.maximum(self.dt, self.dc)

    def _heat_with_creep(self, T_b, presc, values_, F_ext, df_min=1 / 32):
        T_a = self.T.copy()
        tol = self.tol_tr
        # initial step from the current stress (without stress, starts with 1/8 of the interval); without an increase of
        # phi (temperature constant or cooling), a single step, without creep
        sig0 = self.nominal_stress()
        self.set_temperature(T_b)
        phi_b = self._phi()
        if not np.any(phi_b > self.phi_hist):
            return self._step(presc, values_, F_ext)
        sty = np.abs(self._transient_increment(sig0, phi_b)).max(initial=0.0)
        self.set_temperature(T_a)
        df = 1.0 / max(1, int(np.ceil(sty / tol))) if sty > 0 else 1 / 8
        f0 = 0.0
        while f0 < 1.0 - 1e-12:
            f1 = min(1.0, f0 + df)
            saved = self.state()
            T_saved = self.T.copy()
            eps0 = self.eps_tr.copy()
            sig0 = self.nominal_stress()
            self.set_temperature(T_a + (T_b - T_a) * f1)
            phi = self._phi()
            d0 = self._transient_increment(sig0, phi)              # predictor
            self.eps_tr = eps0 + d0
            self._transient_force()
            load, it, conv = self._step(presc, values_, F_ext)
            d1 = self._transient_increment(0.5 * (sig0 + self.nominal_stress()), phi)   # corrector
            m1 = np.abs(d1).max(initial=0.0)
            if m1 > 1.5 * tol and df > df_min:
                self.restore(saved)
                self.set_temperature(T_saved)
                df /= 2
                continue
            if np.abs(d1 - d0).max(initial=0.0) > 0.05 * tol:
                self.eps_tr = eps0 + d1
                self._transient_force()
                load, it, conv = self._step(presc, values_, F_ext)
            self.phi_hist = np.maximum(self.phi_hist, phi)
            f0 = f1
            if m1 < 0.5 * tol:
                df = min(1.0, 1.5 * df)
        return load, it, conv

    @property
    def mean_transient_strain(self):
        """Mean (by volume) of the vertical transient strain in the damaged elements."""
        vc = self.vol[self.ic]
        return float((self.eps_tr[:, 2] * vc).sum() / vc.sum()) if len(vc) else 0.0

    def apply_load(self, delta_max: float, n_steps: int = 30, stop_after_peak: float = 0.80,
                 max_cuts: int = 3, max_load: Optional[float] = None, max_no_conv: int = 25):
        """Moves the top platen (uniformly, without friction) up to delta_max, until the
        load falls to stop_after_peak x the peak or exceeds max_load (N).
        Returns (deltas, loads) in m and N."""
        presc = 3 * self.top + 2
        u0 = self.u[presc].copy()
        F0 = np.zeros(self.ndof)
        deltas, loads = [0.0], [self._step(presc, u0, F0)[0]]
        step0 = delta_max / n_steps
        min_step = step0 / 2 ** max_cuts
        step = step0
        delta, peak, no_conv = 0.0, loads[0], 0
        while delta < delta_max - 1e-12:
            state = self.state()
            dn = delta + step
            load, it, conv = self._step(presc, u0 - dn, F0)
            if not conv and step / 2 >= min_step * (1 - 1e-9):
                self.restore(state)
                step /= 2
                continue
            delta = dn
            no_conv = 0 if conv else no_conv + 1     # accepted without convergence only at the minimum step
            if conv and it <= 20:
                step = min(1.5 * step, step0)
            deltas.append(delta); loads.append(load)
            self._record("load", delta, load, it, conv)
            peak = max(peak, load)
            if load < stop_after_peak * peak and peak > 0:
                break
            if max_load is not None and load >= max_load:
                break
            if no_conv >= max_no_conv:
                if self.verbose:
                    print(f"  stopped: {max_no_conv} consecutive increments without convergence at the minimum step", flush=True)
                break
        return np.array(deltas), np.array(loads)

    _STATE_FIELDS = ("u", "d", "dt", "dc", "r", "kappa", "h", "En", "Gs", "closed", "eps_tr", "phi_hist", "ftr1")

    def state(self):
        return {c: getattr(self, c).copy() for c in self._STATE_FIELDS}

    def restore(self, state):
        for c, v in state.items():
            setattr(self, c, v)

    def global_damage(self, field="d"):
        """Damage per element on the full mesh (nan outside the ceramic).
        field: 'd' (effective, with closed cracks), 'dt' (tension), 'dc' (compression)
        or 'max' (the larger of dt and dc)."""
        v = np.maximum(self.dt, self.dc) if field == "max" else getattr(self, field)
        d = np.full(len(self.m.conn3), np.nan)
        d[self.id_glob[self.ic]] = v
        return d
