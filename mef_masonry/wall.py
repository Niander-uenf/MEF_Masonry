"""Dry-stacked masonry wall under vertical load and fire on one face.

Vertical strip model (1 m of wall length), suitable for walls
whose vertical edges do not restrain out-of-plane bending:

* 2D corotational beam elements (Crisfield 1991; Battini 2002) with a
  fibre section: large displacements and rotations, hence second order
  (P-delta), buckling and thermal bowing;
* fibres = ceramic triangles of the 2D mesh of the section (one unit
  length, with its holes), each with its own temperature, taken from the
  2D thermal analysis (ISO 834 on one face);
* homogenised material: stress (net area) × strain curve of the prism,
  computed with the 3D non-linear model at 30, 600, 800 and 1000 °C and interpolated
  in temperature; no tension (the dry joints open); secant unloading; post-peak
  branch regularised by the element length;
* free thermal strain of the ceramic (alpha or table).

Conventions: x vertical (height), y normal to the wall, exposed face at
y = -t/2. In the fibres, strain and stress are positive in tension; the material
law works with compression positive (e = -eps, s = -sigma).
Forces per metre of wall (N/m) and moments in N·m/m.
"""
from dataclasses import dataclass
from typing import Callable, Dict, Optional, Sequence, Tuple

import numpy as np

from .mesh import HOLE, Mesh2D


# ----------------------------------------------------------------------------
# Fibre section
# ----------------------------------------------------------------------------
@dataclass
class FibreSection:
    y: np.ndarray          # (nf,) distance from the mid-plane, m (negative on the exposed side)
    A: np.ndarray          # (nf,) area per metre of wall, m²/m
    tri: np.ndarray        # (nf,) index of the triangle in the 2D mesh
    thickness: float       # m
    length: float     # length of the modelled section, m
    grout: Optional[np.ndarray] = None     # (nf,) True for grout fibres
    mortar: Optional[np.ndarray] = None  # (nf,) True for joint fibres (vertical section)

    def __post_init__(self):
        if self.grout is None:
            self.grout = np.zeros(len(self.y), bool)
        if self.mortar is None:
            self.mortar = np.zeros(len(self.y), bool)

    @classmethod
    def from_mesh(cls, m2: Mesh2D, geo=None) -> "FibreSection":
        """Fibres: triangles of the unit and, if geo is given, of the grout (course 0).

        In the vertical section of units with horizontal holes, only the walls that
        run down continuously through the course carry the vertical compression: the fibres
        are the triangles outside the columns of holes (face shells, vertical
        webs and the joints over them), with an area per metre of wall equal to the
        triangle area divided by the course height. The horizontal webs
        enter the thermal analysis only."""
        if m2.vertical:
            cnt = m2.xy[m2.tri].mean(1)
            col = np.zeros(len(m2.tri), bool)
            for y0, y1 in (geo.holes_h_y() if geo is not None else []):
                col |= (cnt[:, 1] > y0) & (cnt[:, 1] < y1)
            sol = np.nonzero(~col & (m2.mark != HOLE))[0]
            grt = np.zeros(len(m2.tri), bool)
            mor = m2.mark > 0
        elif geo is not None:
            from .mesh import course_states
            hole, grt = course_states(m2, geo, 0)
            sol = np.nonzero(~hole | grt)[0]
            mor = np.zeros(len(m2.tri), bool)
        else:
            hole, grt = m2.mark == HOLE, np.zeros(len(m2.tri), bool)
            sol = np.nonzero(~hole)[0]
            mor = np.zeros(len(m2.tri), bool)
        cnt = m2.xy[m2.tri[sol]].mean(1)
        a = np.abs(m2.areas()[sol])
        t = m2.width
        return cls(y=cnt[:, 1] - t / 2, A=a / m2.length, tri=sol, thickness=t, length=m2.length,
                   grout=grt[sol], mortar=mor[sol])

    @classmethod
    def rectangle(cls, t: float, n: int = 200) -> "FibreSection":
        """Solid section of thickness t (verifications): n equal strips."""
        y = (np.arange(n) + 0.5) * t / n - t / 2
        return cls(y=y, A=np.full(n, t / n), tri=np.arange(n), thickness=t, length=1.0)

    def temperatures(self, m2: Mesh2D, T_nodes: np.ndarray) -> np.ndarray:
        return np.asarray(T_nodes)[m2.tri[self.tri]].mean(1)

    @property
    def area(self):
        return float(self.A.sum())

    @property
    def inertia(self):
        return float((self.A * self.y ** 2).sum())


# ----------------------------------------------------------------------------
# Material laws (compression positive)
# ----------------------------------------------------------------------------
def _interp_tangent(x, xp, fp):
    """Piecewise linear interpolation with its derivative; outside the table: end values and zero derivative."""
    j = np.clip(np.searchsorted(xp, x, side="right"), 1, len(xp) - 1)
    x0, x1, f0, f1 = xp[j - 1], xp[j], fp[j - 1], fp[j]
    d = (f1 - f0) / (x1 - x0)
    f = f0 + (x - x0) * d
    outside = x >= xp[-1]
    f = np.where(outside, fp[-1], f)
    d = np.where(outside, 0.0, d)
    return f, d


class LinearLaw:
    """Linear elastic with tension (verifications)."""

    def __init__(self, E: float):
        self.E = E
        self.E0 = E

    def stress(self, e, T, e_max):
        return self.E * e, np.full(np.shape(e), self.E)

    def regularised(self, L_e):
        return self


class WallLaw:
    """Homogenised law of dry-stacked masonry: s(e, T) from tabulated curves.

    T      : temperatures of the curves (°C), increasing
    e, s   : for each temperature, compressive strain and stress on the net
             area (Pa), starting at (0, 0) and ending at s = 0
    e_peak : strain at the peak of each curve
    """

    def __init__(self, T, e, s, H_ref=0.3):
        self.T = np.asarray(T, float)
        self.e = [np.asarray(x, float) for x in e]
        self.s = [np.asarray(x, float) for x in s]
        self.H_ref = H_ref
        ip = [int(np.argmax(si)) for si in self.s]
        self.e_peak = np.array([ei[i] for ei, i in zip(self.e, ip)])
        self.s_peak = np.array([si[i] for si, i in zip(self.s, ip)])
        self.E0 = float(self.s_peak[0] / self.e_peak[0])          # secant at the peak at 30 °C
        self._I = np.eye(len(self.T))
        self.d_min = 1e-6 * self.E0

    @classmethod
    def from_prism(cls, curves: Dict[float, Sequence[Tuple[float, float]]], H_ref: float, A_net: float):
        """curves {T: [(delta_mm, load_kN), ...]} of the prism of height H_ref and net area A_net."""
        Ts, es, ss = [], [], []
        for T in sorted(curves):
            a = np.asarray(curves[T], float)
            a = a[np.argsort(a[:, 0])]
            e = np.r_[0.0, a[:, 0] * 1e-3 / H_ref]
            s = np.r_[0.0, a[:, 1] * 1e3 / A_net]
            e, k = np.unique(e, return_index=True)
            s = s[k]
            # post-peak branch down to zero stress, with the final slope (or, if the curve stopped at the peak,
            # with a mean descending slope of 30% of the secant modulus at the peak)
            ip = int(np.argmax(s))
            if ip < len(s) - 1 and s[-1] < s[-2]:
                inc = (s[-1] - s[-2]) / (e[-1] - e[-2])
            else:
                inc = -0.3 * s[ip] / e[ip]
            e_zero = e[-1] + s[-1] / -inc
            e, s = np.r_[e, e_zero], np.r_[s, 0.0]
            Ts.append(T); es.append(e); ss.append(s)
        return cls(Ts, es, ss, H_ref)

    def scaled(self, k_s: float, k_e: float) -> "WallLaw":
        return WallLaw(self.T, [e * k_e for e in self.e], [s * k_s for s in self.s], self.H_ref)

    def with_curve(self, T: float, e, s) -> "WallLaw":
        """Replaces (or adds) the curve of one temperature."""
        Ts, es, ss = list(self.T), list(self.e), list(self.s)
        if T in Ts:
            i = Ts.index(T)
            es[i], ss[i] = np.asarray(e, float), np.asarray(s, float)
        else:
            Ts.append(T); es.append(np.asarray(e, float)); ss.append(np.asarray(s, float))
        o = np.argsort(Ts)
        return WallLaw([Ts[i] for i in o], [es[i] for i in o], [ss[i] for i in o], self.H_ref)

    def regularised(self, L_e: float) -> "WallLaw":
        """Post-peak branch with the post-peak displacement of the prism concentrated in
        one element of length L_e: e = e_p + (e_prism - e_p) H_ref / L_e."""
        k = self.H_ref / L_e
        es = []
        for e, ep in zip(self.e, self.e_peak):
            es.append(np.where(e > ep, ep + (e - ep) * k, e))
        return WallLaw(self.T, es, self.s, self.H_ref)

    def wts(self, T):
        T = np.asarray(T, float)
        return np.stack([np.interp(T, self.T, self._I[i]) for i in range(len(self.T))])

    def envelope(self, e, T, w=None):
        w = self.wts(T) if w is None else w
        s = np.zeros(np.broadcast(e, w[0]).shape)
        d = np.zeros_like(s)
        for i in range(len(self.T)):
            if not np.any(w[i]):
                continue
            si, di = _interp_tangent(e, self.e[i], self.s[i])
            s += w[i] * si
            d += w[i] * di
        return s, d

    def stress(self, e, T, e_max):
        """Compressive stress s(e) and tangent ds/de, with secant unloading
        from the largest compression reached so far (e_max)."""
        w = self.wts(T)
        s, d = self.envelope(e, T, w)
        desc = (e < e_max) & (e_max > 0)
        if np.any(desc):
            k = np.nonzero(desc)
            wk = np.broadcast_to(w, (len(self.T),) + np.shape(e))
            em = e_max[k]
            sm, _ = self.envelope(em, None, wk[(slice(None),) + k])
            sec = sm / em
            s[k] = sec * np.maximum(e[k], 0.0)
            d[k] = sec
        tr = e < 0
        s = np.where(tr, 0.0, s)
        d = np.where(tr | (d == 0.0), self.d_min, d)
        return s, d

    def strength(self, T):
        """Peak stress (net area) at temperature T."""
        return float(self.wts(np.asarray([T]))[:, 0] @ self.s_peak)


# ----------------------------------------------------------------------------
# Corotational beam with fibres
# ----------------------------------------------------------------------------
_DOF = {2: ([0.5 - 0.5 / np.sqrt(3), 0.5 + 0.5 / np.sqrt(3)], [0.5, 0.5]),
       3: ([0.5 - 0.5 * np.sqrt(0.6), 0.5, 0.5 + 0.5 * np.sqrt(0.6)], [5 / 18, 8 / 18, 5 / 18])}


class NotConverged(Exception):
    pass


class WallModel:
    """Vertical wall strip with n_el elements.

    base_support, top_support : "pinned" or "fixed" (the top is free vertically)
    imperfection              : sinusoidal initial deflection at mid-height (m), towards the exposed face
    thermal_strain            : function T (°C) -> free thermal strain
    creep                     : function T (°C) -> c(T) (1/Pa) of the transient creep per fibre:
                                on first heating, Δε_tr = Δc · σ (compression); None switches it off
    """

    def __init__(self, H: float, section: FibreSection, law, thermal_strain: Callable, n_el: int = 24,
                 base_support: str = "pinned", top_support: str = "pinned", imperfection: float = 0.0,
                 n_gauss: int = 3, T_initial: float = 20.0, regularise: bool = True,
                 creep: Optional[Callable] = None):
        self.H, self.sec, self.n_el = H, section, int(n_el)
        n = self.n_el
        self.x0 = np.linspace(0.0, H, n + 1)
        self.y0 = -imperfection * np.sin(np.pi * self.x0 / H)
        dx, dy = np.diff(self.x0), np.diff(self.y0)
        self.L0 = np.hypot(dx, dy)
        self.beta0 = np.arctan2(dy, dx)
        self.law = law.regularised(float(self.L0.mean())) if regularise else law
        xi, w = _DOF[n_gauss]
        self.xi, self.wg = np.asarray(xi), np.asarray(w)
        self.B1, self.B2 = -4 + 6 * self.xi, -2 + 6 * self.xi        # x L0: curvatura = (B1 th1 + B2 th2) / L0
        self.nd = 3 * (n + 1)
        fix = [0, 1] + ([2] if base_support == "fixed" else []) + \
              [3 * n + 1] + ([3 * n + 2] if top_support == "fixed" else [])
        self.fixed_dofs = np.array(sorted(fix))
        self.free_dofs = np.setdiff1d(np.arange(self.nd), self.fixed_dofs)
        self.supports = (base_support, top_support)
        self.thermal_strain = thermal_strain
        nf = len(section.y)
        self.emax = np.zeros((n, len(self.xi), nf))
        self.d = np.zeros(self.nd)
        self.T_initial = T_initial
        self.T_f = np.full(nf, float(T_initial))
        self.dofs = (3 * np.arange(n))[:, None] + np.arange(6)[None, :]
        self.last = None
        # largest increment of thermal strain per step (prevents the trial state
        # from jumping over the peak of the curve and Newton from falling onto a spurious crushed branch)
        self.deps_max = 0.25 * float(np.min(self.law.e_peak)) if hasattr(self.law, "e_peak") else np.inf
        self.N_ref = max(float(getattr(self.law, "s_peak", [self.law.E0 * 1e-3])[0]) * section.area, 1.0)
        # transient creep: compressive strain (positive) per Gauss point and fibre
        self.creep = creep
        self.e_tr = np.zeros((n, len(self.xi), nf))
        self.c_hist = creep(self.T_f) if creep else np.zeros(nf)
        self.tol_tr = 1e-4

    # ------------------------------------------------------------------ kinematics and forces
    def _state(self, d, T_f, eth, emax, e_tr=None):
        n, sec = self.n_el, self.sec
        u = d.reshape(-1, 3)
        x = self.x0 + u[:, 0]
        y = self.y0 + u[:, 1]
        dx, dy = np.diff(x), np.diff(y)
        Ln = np.hypot(dx, dy)
        beta = np.arctan2(dy, dx)
        c, s = np.cos(beta), np.sin(beta)
        cb0, sb0 = np.cos(self.beta0), np.sin(self.beta0)
        alpha = np.arctan2(cb0 * s - sb0 * c, cb0 * c + sb0 * s)
        L0 = self.L0
        ul = (Ln ** 2 - L0 ** 2) / (Ln + L0)
        th1 = u[:-1, 2] - alpha
        th2 = u[1:, 2] - alpha
        # fibres
        eps0 = ul / L0
        kap = (self.B1[None, :] * th1[:, None] + self.B2[None, :] * th2[:, None]) / L0[:, None]   # (n, ng)
        eps = eps0[:, None, None] - sec.y[None, None, :] * kap[:, :, None]                       # (n, ng, nf)
        e = -(eps - eth[None, None, :])                                                           # compression +
        e = e - (self.e_tr if e_tr is None else e_tr)                       # minus the transient creep
        sc, Et = self.law.stress(e, T_f[None, None, :], emax)
        sig = -sc
        A, yf = sec.A, sec.y
        N = sig @ A
        M = -(sig @ (A * yf))
        EA = Et @ A
        ES = -(Et @ (A * yf))
        EI = Et @ (A * yf * yf)
        wg, B1, B2 = self.wg, self.B1, self.B2
        q = np.c_[N @ wg, (M * B1) @ wg, (M * B2) @ wg]                                            # (n, 3)
        kl = np.empty((n, 3, 3))
        kl[:, 0, 0] = (EA @ wg) / L0
        kl[:, 0, 1] = kl[:, 1, 0] = ((ES * B1) @ wg) / L0
        kl[:, 0, 2] = kl[:, 2, 0] = ((ES * B2) @ wg) / L0
        kl[:, 1, 1] = ((EI * B1 * B1) @ wg) / L0
        kl[:, 1, 2] = kl[:, 2, 1] = ((EI * B1 * B2) @ wg) / L0
        kl[:, 2, 2] = ((EI * B2 * B2) @ wg) / L0
        z0 = np.zeros(n)
        r = np.stack([-c, -s, z0, c, s, z0], axis=1)
        zz = np.stack([s, -c, z0, -s, c, z0], axis=1)
        Bg = np.zeros((n, 3, 6))
        Bg[:, 0] = r
        Bg[:, 1] = -zz / Ln[:, None]; Bg[:, 1, 2] += 1.0
        Bg[:, 2] = -zz / Ln[:, None]; Bg[:, 2, 5] += 1.0
        fe = np.einsum("eij,ei->ej", Bg, q)
        Ke = np.einsum("eki,ekl,elj->eij", Bg, kl, Bg)
        Ke += (q[:, 0] / Ln)[:, None, None] * np.einsum("ei,ej->eij", zz, zz)
        Ke += ((q[:, 1] + q[:, 2]) / Ln ** 2)[:, None, None] * (np.einsum("ei,ej->eij", r, zz) + np.einsum("ei,ej->eij", zz, r))
        f = np.zeros(self.nd)
        np.add.at(f, self.dofs, fe)
        K = np.zeros((self.nd, self.nd))
        np.add.at(K, (self.dofs[:, :, None], self.dofs[:, None, :]), Ke)
        return f, K, {"e": e, "sig": sig, "N": N, "M": M}

    def external_force(self, P, e_top, e_base):
        """Load P (N/m, compression) applied at the top with eccentricity e_top
        and base reaction with e_base (m, positive towards the exposed face)."""
        f = np.zeros(self.nd)
        n = self.n_el
        f[3 * n] = -P
        if self.supports[1] != "fixed":
            f[3 * n + 2] += -P * e_top
        if self.supports[0] != "fixed":
            f[2] += P * e_base
        return f

    # ------------------------------------------------------------------ equilibrium
    def equilibrate(self, fext, T_f, d0=None, max_it=40, tol=1e-8, e_tr=None):
        """Newton-Raphson with line search. Returns (d, info, iterations) or
        raises NotConverged."""
        d = (self.d if d0 is None else d0).copy()
        eth = self.thermal_strain(T_f)
        e_tr = self.e_tr if e_tr is None else e_tr
        scl = np.ones(self.nd)
        scl[2::3] = 1.0 / self.sec.thickness           # moments as equivalent forces
        lv = self.free_dofs
        ref = self.N_ref
        f, K, info = self._state(d, T_f, eth, self.emax, e_tr)
        r = (fext - f)[lv]
        nr = np.linalg.norm(r * scl[lv])
        du = None
        best, since = nr, 0
        for it in range(max_it):
            if not np.isfinite(nr):
                raise NotConverged("non-finite residual")
            if nr < 0.5 * best:
                best, since = nr, 0
            else:
                since += 1
                if since >= 10:
                    raise NotConverged(f"residual stagnated at {nr / ref:.1e}")
            if du is not None and (nr < tol * ref or (nr < 1e-5 * ref and np.abs(du).max() < 1e-10)):
                return d, info, it
            try:
                du = np.linalg.solve(K[np.ix_(lv, lv)], r)
            except np.linalg.LinAlgError:
                raise NotConverged("singular matrix")
            if not np.all(np.isfinite(du)) or np.abs(du[0::3]).max() > 0.5 * self.H:
                raise NotConverged("unbounded displacements")
            step = 1.0
            for _ in range(6):                        # simple line search (residual reduction)
                dn = d.copy(); dn[lv] += step * du
                fn, Kn, infon = self._state(dn, T_f, eth, self.emax, e_tr)
                rn = (fext - fn)[lv]
                nrn = np.linalg.norm(rn * scl[lv])
                if np.isfinite(nrn) and (nrn < nr or nrn < tol * ref):
                    break
                step *= 0.5
            d, f, K, info, r, nr = dn, fn, Kn, infon, rn, nrn
        raise NotConverged(f"no convergence in {max_it} iterations (residual {nr / ref:.1e})")

    def accept(self, d, T_f, info):
        self.d = d
        self.T_f = np.asarray(T_f, float).copy()
        self.emax = np.maximum(self.emax, info["e"])
        self.last = info

    # ------------------------------------------------------------------ measures
    @property
    def deflection(self):
        """Horizontal displacement of the nodes (m); negative = towards the exposed face."""
        return self.d[1::3]

    @property
    def mid_deflection(self):
        return float(np.interp(self.H / 2, self.x0, self.deflection))

    @property
    def max_deflection(self):
        v = self.deflection
        return float(v[np.argmax(np.abs(v))])

    @property
    def top_displacement(self):
        """Vertical displacement of the top (m); positive = upwards."""
        return float(self.d[3 * self.n_el])

    def critical_compression(self):
        """(height x in m, ratio e/e_peak) of the Gauss point with the fibre most
        compressed relative to the peak strain at its temperature."""
        if self.last is None:
            return None, 0.0
        ep = self.law.wts(self.T_f).T @ self.law.e_peak
        r = (self.last["e"] / ep[None, None, :]).max(-1)
        e, g = np.unravel_index(int(np.argmax(r)), r.shape)
        x = self.x0[e] + self.xi[g] * (self.x0[e + 1] - self.x0[e])
        return float(x), float(r[e, g])

    def mid_height_fibres(self):
        """Stress (Pa) in the fibres at the Gauss point closest to H/2."""
        if self.last is None:
            return None
        xg = (self.x0[:-1, None] + self.xi[None, :] * np.diff(self.x0)[:, None]).ravel()
        k = int(np.argmin(np.abs(xg - self.H / 2)))
        e, g = divmod(k, len(self.xi))
        return self.last["sig"][e, g]

    # ------------------------------------------------------------------ analyses
    def apply_load(self, P, e_top, e_base, n_steps=10, min_frac=1 / 512, record=None):
        """Applies P (N/m) at the current temperature in steps, with subdivision. Returns True if it got there."""
        fu = self.external_force(1.0, e_top, e_base)
        Pa, dP = 0.0, P / n_steps
        while Pa < P * (1 - 1e-12):
            Pn = min(P, Pa + dP)
            try:
                d, info, _ = self.equilibrate(Pn * fu, self.T_f)
            except NotConverged:
                dP /= 2
                if dP < P * min_frac:
                    return False
                continue
            self.accept(d, self.T_f, info)
            Pa = Pn
            if record:
                record(Pa)
        return True

    def capacity(self, e_top, e_base, P_step=None, tol=1e-3, record=None, max_steps=400):
        """Maximum load (N/m) by load control with step bisection."""
        fu = self.external_force(1.0, e_top, e_base)
        Pref = self.N_ref
        dP = P_step or Pref / 20
        Pa = 0.0
        for _ in range(max_steps):
            if dP < tol * Pref:
                break
            try:
                d, info, _ = self.equilibrate((Pa + dP) * fu, self.T_f)
            except NotConverged:
                dP /= 2
                continue
            self.accept(d, self.T_f, info)
            Pa += dP
            if record:
                record(Pa)
        return Pa

    def _creep_step(self, fext, Tn, can_split):
        """Equilibrium at Tn with the transient creep of the step (predictor from the
        stress at the start, corrector from the mean stress). Returns None if the
        increment exceeds tol_tr and the step can be split."""
        c_new = self.creep(Tn)
        dc = np.maximum(c_new - self.c_hist, 0.0)[None, None, :]
        s0 = np.maximum(-self.last["sig"], 0.0) if self.last is not None else 0.0
        d0 = dc * s0
        if can_split and np.max(d0, initial=0.0) > self.tol_tr:
            return None
        d, info, _ = self.equilibrate(fext, Tn, e_tr=self.e_tr + d0)
        d1 = dc * 0.5 * (s0 + np.maximum(-info["sig"], 0.0))
        if can_split and np.max(d1, initial=0.0) > 1.5 * self.tol_tr:
            return None
        if np.max(np.abs(d1 - d0), initial=0.0) > 0.05 * self.tol_tr:
            d, info, _ = self.equilibrate(fext, Tn, d0=d, e_tr=self.e_tr + d1)
        return d, info, d1, c_new

    def fire(self, P, e_top, e_base, fields: Sequence[Tuple[float, np.ndarray]], min_sub=1 / 64,
                 record=None, stop_if=None):
        """Constant P; fields [(t s, T_fibres)] increasing in time, from the current state.
        Returns (t_collapse or None, last converged t)."""
        fext = self.external_force(P, e_top, e_base)
        t_a, T_a = fields[0]
        for t_b, T_b in fields[1:]:
            stack = [(0.0, 1.0)]
            while stack:
                f0, f1 = stack.pop()
                Tn = T_a + (T_b - T_a) * f1
                if f1 - f0 > 1 / 1024 and \
                        np.abs(self.thermal_strain(Tn) - self.thermal_strain(self.T_f)).max() > self.deps_max:
                    fm = 0.5 * (f0 + f1)            # large increment of thermal strain: subdivide first
                    stack += [(fm, f1), (f0, fm)]
                    continue
                try:
                    if self.creep is None:
                        d, info, _ = self.equilibrate(fext, Tn)
                    else:
                        res = self._creep_step(fext, Tn, f1 - f0 > 1 / 1024)
                        if res is None:                 # large creep increment: subdivide
                            fm = 0.5 * (f0 + f1)
                            stack += [(fm, f1), (f0, fm)]
                            continue
                        d, info, de_tr, c_new = res
                except NotConverged:
                    if f1 - f0 <= min_sub:
                        return t_a + (t_b - t_a) * f0, t_a + (t_b - t_a) * f0
                    fm = 0.5 * (f0 + f1)
                    stack += [(fm, f1), (f0, fm)]
                    continue
                if self.creep is not None:
                    self.e_tr = self.e_tr + de_tr
                    self.c_hist = np.maximum(self.c_hist, c_new)
                self.accept(d, Tn, info)
                tt = t_a + (t_b - t_a) * f1
                if record:
                    record(tt)
                if stop_if and stop_if(tt):
                    return None, tt
            t_a, T_a = t_b, T_b
        return None, t_a
