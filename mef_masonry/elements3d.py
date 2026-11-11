"""Generic isoparametric solid elements: W6, H8 and W15.

W6  : 6-node linear wedge (linear triangle x linear segment).
H8  : 8-node linear hexahedron (bilinear quadrilateral x linear segment).
W15 : 15-node quadratic wedge (6-node triangle x quadratic segment,
      serendipity), node order of VTK_QUADRATIC_WEDGE in our orientation:
      0-2 bottom corners, 3-5 top corners, 6-8 midpoints of the bottom
      edges (0-1, 1-2, 2-0), 9-11 midpoints of the top ones (3-4, 4-5, 5-3) and 12-14
      midpoints of the vertical edges (0-3, 1-4, 2-5).

Natural coordinates: wedges use (L1, L2, zeta), with L3 = 1 - L1 - L2;
the hexahedron uses (xi, eta, zeta) in [-1, 1]. All routines are vectorised
over the elements and use the full Jacobian at each Gauss point.
Heat capacity and boundary areas are lumped at the nodes by the
HRZ rule (diagonal of the consistent matrix, rescaled), which gives positive
weights also in the quadratic element.
"""
from functools import lru_cache
from typing import List, Tuple

import numpy as np

from .elements import matrix_D

# ---------------------------------------------------------------------------
# Integration rules
# ---------------------------------------------------------------------------
G2 = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)]); W2 = np.array([1.0, 1.0])
G3 = np.array([-np.sqrt(0.6), 0.0, np.sqrt(0.6)]); W3 = np.array([5 / 9, 8 / 9, 5 / 9])
# triangle: 3 points (degree 2) and 6 Dunavant points (degree 4); weights sum to 1/2
TRI3 = (np.array([[2 / 3, 1 / 6], [1 / 6, 2 / 3], [1 / 6, 1 / 6]]), np.full(3, 1 / 6))
_a, _b = 0.445948490915965, 0.091576213509771
_wa, _wb = 0.223381589678011 / 2, 0.109951743655322 / 2
TRI6 = (np.array([[_a, _a], [1 - 2 * _a, _a], [_a, 1 - 2 * _a], [_b, _b], [1 - 2 * _b, _b], [_b, 1 - 2 * _b]]),
        np.array([_wa] * 3 + [_wb] * 3))


def _wedge(tri, line):
    p, w = [], []
    for (l1, l2), wt in zip(*tri):
        for z, wz in zip(*line):
            p.append([l1, l2, z]); w.append(wt * wz)
    return np.array(p), np.array(w)


def _hexa(line):
    p, w = [], []
    for a, wa in zip(*line):
        for b, wb in zip(*line):
            for c, wc in zip(*line):
                p.append([a, b, c]); w.append(wa * wb * wc)
    return np.array(p), np.array(w)


# ---------------------------------------------------------------------------
# Element types
# ---------------------------------------------------------------------------
class ElementType:
    name = ""
    nn = 0
    vtk = 0
    vtk_order: List[int] = []
    base: List[int] = []            # corners of the bottom face
    top: List[int] = []            # corners of the top face
    faces_vis: List[List[int]] = []  # faces by their corners (visualisation)
    centre = np.zeros(3)

    def N(self, xi):  # (np, nn)
        raise NotImplementedError

    def dN(self, xi):  # (np, nn, 3)
        raise NotImplementedError

    # faces for boundary integration: (local nodes of the face, parametrisation (s, t) -> xi, rule)
    faces_int: List = []

    def points(self):
        raise NotImplementedError

    # ----------------------------------------------------------------- geometry
    def gradients(self, xyz, con, xi):
        """Physical gradients G (ne, np, nn, 3) and det J (ne, np) at the points xi (np, 3)."""
        X = xyz[con]                                   # (ne, nn, 3)
        dN = self.dN(np.atleast_2d(xi))                # (np, nn, 3)
        J = np.einsum("eia,pib->epab", X, dN)          # dx_a/dxi_b
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        G = np.einsum("pib,epba->epia", dN, Jinv)
        return G, np.abs(detJ)

    def volumes(self, xyz, con):
        xi, w = self.points()
        _, detJ = self.gradients(xyz, con, xi)
        return detJ @ w

    def height(self, xyz, con):
        """Vertical extent of each element (extruded layers)."""
        z = xyz[con][:, :, 2]
        return z.max(1) - z.min(1)

    def conduction(self, xyz, con):
        """Conduction matrices with k = 1 (ne, nn, nn), volumes (ne,) and HRZ
        weights of the lumped capacity (ne, nn), which sum to 1 in each element."""
        xi, w = self.points()
        G, detJ = self.gradients(xyz, con, xi)
        dv = detJ * w[None, :]
        K = np.einsum("ep,epia,epja->eij", dv, G, G)
        vol = dv.sum(1)
        N = self.N(xi)
        diag = np.einsum("ep,pi->ei", dv, N * N)
        return K, np.abs(vol), diag / diag.sum(1, keepdims=True)

    @staticmethod
    def B(G):
        """B (..., 6, 3 nn) from G (..., nn, 3); Voigt xx, yy, zz, xy, yz, xz."""
        sh = G.shape[:-2]
        nn = G.shape[-2]
        B = np.zeros(sh + (6, 3 * nn))
        gx, gy, gz = G[..., 0], G[..., 1], G[..., 2]
        B[..., 0, 0::3] = gx
        B[..., 1, 1::3] = gy
        B[..., 2, 2::3] = gz
        B[..., 3, 0::3] = gy; B[..., 3, 1::3] = gx
        B[..., 4, 1::3] = gz; B[..., 4, 2::3] = gy
        B[..., 5, 0::3] = gz; B[..., 5, 2::3] = gx
        return B

    def stiffness(self, xyz, con, E, nu, eps_th_nodes=None, D=None):
        """Stiffness (ne, 3nn, 3nn) and thermal forces (ne, 3nn) if eps_th_nodes (ne, nn) is given.
        D (ne, 6, 6) replaces the isotropic material (E, nu)."""
        xi, w = self.points()
        G, detJ = self.gradients(xyz, con, xi)
        Nn = self.N(xi)
        D = matrix_D(E, nu) if D is None else D
        ne, nd = len(con), 3 * self.nn
        K = np.zeros((ne, nd, nd))
        f = np.zeros((ne, nd)) if eps_th_nodes is not None else None
        one = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        for p in range(len(w)):
            B = self.B(G[:, p])
            dv = detJ[:, p] * w[p]
            DB = np.einsum("eij,ejk->eik", D, B)
            K += dv[:, None, None] * np.einsum("eji,ejk->eik", B, DB)
            if f is not None:
                sig0 = (D @ one) * (eps_th_nodes @ Nn[p])[:, None]
                f += dv[:, None] * np.einsum("eji,ej->ei", B, sig0)
        return K, f

    def thermal_force(self, xyz, con, E, nu, eps_th_nodes):
        xi, w = self.points()
        G, detJ = self.gradients(xyz, con, xi)
        Nn = self.N(xi)
        D = matrix_D(E, nu)
        one = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        f = np.zeros((len(con), 3 * self.nn))
        for p in range(len(w)):
            B = self.B(G[:, p])
            sig0 = (D @ one) * (eps_th_nodes @ Nn[p])[:, None]
            f += (detJ[:, p] * w[p])[:, None] * np.einsum("eji,ej->ei", B, sig0)
        return f

    def at_centre(self, xyz, con):
        """B (ne, 6, 3nn), N (nn,) and G (ne, nn, 3) at the element centre."""
        G, _ = self.gradients(xyz, con, self.centre[None, :])
        return self.B(G[:, 0]), self.N(self.centre[None, :])[0], G[:, 0]

    def points_B(self, xyz, con):
        """List [(dv (ne,), B (ne, 6, 3nn), N (nn,))] at the Gauss points."""
        xi, w = self.points()
        G, detJ = self.gradients(xyz, con, xi)
        Nn = self.N(xi)
        return [(detJ[:, p] * w[p], self.B(G[:, p]), Nn[p]) for p in range(len(w))]

    def stresses(self, xyz, con, E, nu, u_el, eps_th_nodes):
        B, N, _ = self.at_centre(xyz, con)
        eps = np.einsum("eij,ej->ei", B, u_el)
        eps[:, :3] -= (eps_th_nodes @ N)[:, None]
        return np.einsum("eij,ej->ei", matrix_D(E, nu), eps)

    # ----------------------------------------------------------------- faces
    def face_areas(self, xyz, con, selector, consistent=False):
        """Nodal areas of the boundary faces chosen by
        selector(coordinates of the face nodes (nf, k, 3)) -> mask (nf,).
        consistent=False: HRZ weights (positive; heat exchange on the boundary);
        consistent=True: integral of N_i (consistent load vector of a
        uniform pressure; negative at the corners of the quadratic face).
        For linear elements with flat faces the two coincide.
        Returns a vector (n_nodes,)."""
        area = np.zeros(len(xyz))
        for face_nodes, emap, rule in self.faces_int:
            fc = con[:, face_nodes]
            sel = selector(xyz[fc])
            if not sel.any():
                continue
            ce = con[sel]
            st, ws = rule
            xi, dxi_ds, dxi_dt = emap(st)
            N = self.N(xi)                                   # (np, nn)
            dN = self.dN(xi)                                 # (np, nn, 3)
            X = xyz[ce]                                      # (ne, nn, 3)
            dxds = np.einsum("eia,pib,pb->epa", X, dN, dxi_ds)
            dxdt = np.einsum("eia,pib,pb->epa", X, dN, dxi_dt)
            dA = np.linalg.norm(np.cross(dxds, dxdt), axis=2) * ws[None, :]     # (ne, np)
            Af = dA.sum(1)
            if consistent:
                wts = np.einsum("ep,pi->ei", dA, N[:, face_nodes])
            else:
                d = np.einsum("ep,pi->ei", dA, N[:, face_nodes] ** 2)
                wts = d / d.sum(1, keepdims=True) * Af[:, None]
            area += np.bincount(ce[:, face_nodes].ravel(), weights=wts.ravel(), minlength=len(xyz))
        return area


def _wedge_side_face(i, j, side_rule, line):
    """Side face of the wedge over the edge (i, j) of the triangle (L_i: 1 -> 0, L_j: 0 -> 1)."""
    def emap(st):
        s, t = st[:, 0], st[:, 1]
        xi = np.zeros((len(s), 3))
        L = np.zeros((len(s), 3)); L[:, i] = 1 - s; L[:, j] = s
        xi[:, 0], xi[:, 1], xi[:, 2] = L[:, 0], L[:, 1], t
        ds = np.zeros((len(s), 3)); dt = np.zeros((len(s), 3))
        dL = np.zeros(3); dL[i] = -1; dL[j] = 1
        ds[:, 0], ds[:, 1] = dL[0], dL[1]
        dt[:, 2] = 1.0
        return xi, ds, dt
    p, w = [], []
    for s, ws in zip(*side_rule):
        for t, wt in zip(*line):
            p.append([s, t]); w.append(ws * wt)
    return emap, (np.array(p), np.array(w))


def _wedge_tri_face(zeta, tri_rule):
    def emap(st):
        xi = np.c_[st[:, 0], st[:, 1], np.full(len(st), zeta)]
        ds = np.tile([1.0, 0.0, 0.0], (len(st), 1)); dt = np.tile([0.0, 1.0, 0.0], (len(st), 1))
        return xi, ds, dt
    return emap, tri_rule


def _rule_01(line):
    """Gauss rule on [-1, 1] mapped to [0, 1]."""
    return ((line[0] + 1) / 2, line[1] / 2)


class W6(ElementType):
    name, nn, vtk = "W6", 6, 13
    vtk_order = [0, 2, 1, 3, 5, 4]
    base, top = [0, 1, 2], [3, 4, 5]
    faces_vis = [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]]
    centre = np.array([1 / 3, 1 / 3, 0.0])

    def __init__(self):
        side = _rule_01((G2, W2))
        self.faces_int = [
            ([0, 2, 1], *_wedge_tri_face(-1.0, TRI3)), ([3, 4, 5], *_wedge_tri_face(1.0, TRI3)),
            ([0, 1, 4, 3], *_wedge_side_face(0, 1, side, (G2, W2))),
            ([1, 2, 5, 4], *_wedge_side_face(1, 2, side, (G2, W2))),
            ([2, 0, 3, 5], *_wedge_side_face(2, 0, side, (G2, W2)))]

    def points(self):
        return _wedge(TRI3, (G2, W2))

    def N(self, xi):
        xi = np.atleast_2d(xi)
        L = np.c_[xi[:, 0], xi[:, 1], 1 - xi[:, 0] - xi[:, 1]]
        fb, ft = (1 - xi[:, 2]) / 2, (1 + xi[:, 2]) / 2
        return np.c_[L * fb[:, None], L * ft[:, None]]

    def dN(self, xi):
        xi = np.atleast_2d(xi)
        L = np.c_[xi[:, 0], xi[:, 1], 1 - xi[:, 0] - xi[:, 1]]
        fb, ft = (1 - xi[:, 2]) / 2, (1 + xi[:, 2]) / 2
        dL = np.array([[1, 0], [0, 1], [-1, -1]], float)    # dL_i/d(L1, L2)
        d = np.zeros((len(xi), 6, 3))
        d[:, :3, :2] = dL[None] * fb[:, None, None]
        d[:, 3:, :2] = dL[None] * ft[:, None, None]
        d[:, :3, 2] = -L / 2
        d[:, 3:, 2] = L / 2
        return d


class H8(ElementType):
    name, nn, vtk = "H8", 8, 12
    vtk_order = list(range(8))
    base, top = [0, 1, 2, 3], [4, 5, 6, 7]
    faces_vis = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]
    centre = np.zeros(3)
    _s = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                   [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)

    def __init__(self):
        q = []
        for a, wa in zip(G2, W2):
            for b, wb in zip(G2, W2):
                q.append([a, b, wa * wb])
        q = np.array(q)
        rule = (q[:, :2], q[:, 2])

        def face(fixed, value, u, v):
            def emap(st):
                xi = np.zeros((len(st), 3)); xi[:, fixed] = value; xi[:, u] = st[:, 0]; xi[:, v] = st[:, 1]
                ds = np.zeros((len(st), 3)); ds[:, u] = 1
                dt = np.zeros((len(st), 3)); dt[:, v] = 1
                return xi, ds, dt
            return emap
        self.faces_int = [([0, 3, 2, 1], face(2, -1, 0, 1), rule), ([4, 5, 6, 7], face(2, 1, 0, 1), rule),
                          ([0, 1, 5, 4], face(1, -1, 0, 2), rule), ([1, 2, 6, 5], face(0, 1, 1, 2), rule),
                          ([2, 3, 7, 6], face(1, 1, 0, 2), rule), ([3, 0, 4, 7], face(0, -1, 1, 2), rule)]

    def points(self):
        return _hexa((G2, W2))

    def N(self, xi):
        xi = np.atleast_2d(xi)
        s = self._s
        return np.prod(1 + xi[:, None, :] * s[None, :, :], axis=2) / 8

    def dN(self, xi):
        xi = np.atleast_2d(xi)
        s = self._s
        f = 1 + xi[:, None, :] * s[None, :, :]                 # (np, 8, 3)
        d = np.empty((len(xi), 8, 3))
        d[..., 0] = s[None, :, 0] * f[..., 1] * f[..., 2] / 8
        d[..., 1] = s[None, :, 1] * f[..., 0] * f[..., 2] / 8
        d[..., 2] = s[None, :, 2] * f[..., 0] * f[..., 1] / 8
        return d


class W15(ElementType):
    name, nn, vtk = "W15", 15, 26
    # VTK_QUADRATIC_WEDGE: corners as VTK_WEDGE (0, 2, 1, 3, 5, 4), then edges in the same order
    vtk_order = [0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13]
    base, top = [0, 1, 2], [3, 4, 5]
    faces_vis = W6.faces_vis
    centre = np.array([1 / 3, 1 / 3, 0.0])

    def __init__(self):
        side = _rule_01((G3, W3))
        self.faces_int = [
            ([0, 2, 1, 8, 7, 6], *_wedge_tri_face(-1.0, TRI6)), ([3, 4, 5, 9, 10, 11], *_wedge_tri_face(1.0, TRI6)),
            ([0, 1, 4, 3, 6, 13, 9, 12], *_wedge_side_face(0, 1, side, (G3, W3))),
            ([1, 2, 5, 4, 7, 14, 10, 13], *_wedge_side_face(1, 2, side, (G3, W3))),
            ([2, 0, 3, 5, 8, 12, 11, 14], *_wedge_side_face(2, 0, side, (G3, W3)))]

    def points(self):
        return _wedge(TRI6, (G3, W3))

    def N(self, xi):
        xi = np.atleast_2d(xi)
        L1, L2, z = xi[:, 0], xi[:, 1], xi[:, 2]
        L = [L1, L2, 1 - L1 - L2]
        fb, ft, fm = 1 - z, 1 + z, 1 - z * z
        N = np.empty((len(xi), 15))
        for i in range(3):
            N[:, i] = 0.5 * L[i] * (2 * L[i] - 1) * fb - 0.5 * L[i] * fm
            N[:, 3 + i] = 0.5 * L[i] * (2 * L[i] - 1) * ft - 0.5 * L[i] * fm
            N[:, 12 + i] = L[i] * fm
        for k, (i, j) in enumerate([(0, 1), (1, 2), (2, 0)]):
            N[:, 6 + k] = 2 * L[i] * L[j] * fb
            N[:, 9 + k] = 2 * L[i] * L[j] * ft
        return N

    def dN(self, xi):
        xi = np.atleast_2d(xi)
        L1, L2, z = xi[:, 0], xi[:, 1], xi[:, 2]
        L = [L1, L2, 1 - L1 - L2]
        dL = [np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([-1.0, -1.0])]
        fb, ft, fm = 1 - z, 1 + z, 1 - z * z
        d = np.zeros((len(xi), 15, 3))
        for i in range(3):
            a = 0.5 * (4 * L[i] - 1)                        # d[0.5 L (2L - 1)]/dL
            for c in range(2):
                d[:, i, c] = (a * fb - 0.5 * fm) * dL[i][c]
                d[:, 3 + i, c] = (a * ft - 0.5 * fm) * dL[i][c]
                d[:, 12 + i, c] = fm * dL[i][c]
            g = 0.5 * L[i] * (2 * L[i] - 1)
            d[:, i, 2] = -g + L[i] * z
            d[:, 3 + i, 2] = g + L[i] * z
            d[:, 12 + i, 2] = -2 * L[i] * z
        for k, (i, j) in enumerate([(0, 1), (1, 2), (2, 0)]):
            for c in range(2):
                dLL = 2 * (dL[i][c] * L[j] + L[i] * dL[j][c])
                d[:, 6 + k, c] = dLL * fb
                d[:, 9 + k, c] = dLL * ft
            d[:, 6 + k, 2] = -2 * L[i] * L[j]
            d[:, 9 + k, 2] = 2 * L[i] * L[j]
        return d


@lru_cache(maxsize=None)
def kind(name: str) -> ElementType:
    return {"W6": W6, "H8": H8, "W15": W15}[name]()
