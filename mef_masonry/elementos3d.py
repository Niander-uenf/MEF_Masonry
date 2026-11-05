"""Elementos sólidos isoparamétricos genéricos: W6, H8 e W15.

W6  : cunha linear de 6 nós (triângulo linear x segmento linear).
H8  : hexaedro linear de 8 nós (quadrilátero bilinear x segmento linear).
W15 : cunha quadrática de 15 nós (triângulo de 6 nós x segmento quadrático,
      serendipity), ordem de nós do VTK_QUADRATIC_WEDGE no nosso sentido:
      0-2 cantos de baixo, 3-5 cantos de cima, 6-8 meios das arestas de
      baixo (0-1, 1-2, 2-0), 9-11 meios das de cima (3-4, 4-5, 5-3) e 12-14
      meios das arestas verticais (0-3, 1-4, 2-5).

Coordenadas naturais: cunhas usam (L1, L2, zeta), com L3 = 1 - L1 - L2;
hexaedro usa (xi, eta, zeta) em [-1, 1]. Todas as rotinas são vetorizadas
sobre os elementos e usam o jacobiano completo em cada ponto de Gauss.
Capacidade térmica e áreas de contorno são concentradas nos nós pela
regra HRZ (diagonal da matriz consistente, reescalada), que dá pesos
positivos também no elemento quadrático.
"""
from functools import lru_cache
from typing import List, Tuple

import numpy as np

from .elementos import matriz_D

# ---------------------------------------------------------------------------
# Regras de integração
# ---------------------------------------------------------------------------
G2 = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)]); W2 = np.array([1.0, 1.0])
G3 = np.array([-np.sqrt(0.6), 0.0, np.sqrt(0.6)]); W3 = np.array([5 / 9, 8 / 9, 5 / 9])
# triângulo: 3 pontos (grau 2) e 6 pontos de Dunavant (grau 4); pesos somam 1/2
TRI3 = (np.array([[2 / 3, 1 / 6], [1 / 6, 2 / 3], [1 / 6, 1 / 6]]), np.full(3, 1 / 6))
_a, _b = 0.445948490915965, 0.091576213509771
_wa, _wb = 0.223381589678011 / 2, 0.109951743655322 / 2
TRI6 = (np.array([[_a, _a], [1 - 2 * _a, _a], [_a, 1 - 2 * _a], [_b, _b], [1 - 2 * _b, _b], [_b, 1 - 2 * _b]]),
        np.array([_wa] * 3 + [_wb] * 3))


def _cunha(tri, linha):
    p, w = [], []
    for (l1, l2), wt in zip(*tri):
        for z, wz in zip(*linha):
            p.append([l1, l2, z]); w.append(wt * wz)
    return np.array(p), np.array(w)


def _hexa(linha):
    p, w = [], []
    for a, wa in zip(*linha):
        for b, wb in zip(*linha):
            for c, wc in zip(*linha):
                p.append([a, b, c]); w.append(wa * wb * wc)
    return np.array(p), np.array(w)


# ---------------------------------------------------------------------------
# Tipos de elemento
# ---------------------------------------------------------------------------
class TipoElemento:
    nome = ""
    nn = 0
    vtk = 0
    vtk_ordem: List[int] = []
    base: List[int] = []            # cantos da face inferior
    topo: List[int] = []            # cantos da face superior
    faces_vis: List[List[int]] = []  # faces pelos cantos (visualização)
    centro = np.zeros(3)

    def N(self, xi):  # (np, nn)
        raise NotImplementedError

    def dN(self, xi):  # (np, nn, 3)
        raise NotImplementedError

    # faces para integração no contorno: (nós locais da face, parametrização (s, t) -> xi, regra)
    faces_int: List = []

    def pontos(self):
        raise NotImplementedError

    # ----------------------------------------------------------------- geometria
    def gradientes(self, xyz, con, xi):
        """Gradientes físicos G (ne, np, nn, 3) e det J (ne, np) nos pontos xi (np, 3)."""
        X = xyz[con]                                   # (ne, nn, 3)
        dN = self.dN(np.atleast_2d(xi))                # (np, nn, 3)
        J = np.einsum("eia,pib->epab", X, dN)          # dx_a/dxi_b
        detJ = np.linalg.det(J)
        Jinv = np.linalg.inv(J)
        G = np.einsum("pib,epba->epia", dN, Jinv)
        return G, np.abs(detJ)

    def volumes(self, xyz, con):
        xi, w = self.pontos()
        _, detJ = self.gradientes(xyz, con, xi)
        return detJ @ w

    def altura(self, xyz, con):
        """Extensão vertical de cada elemento (camadas extrudadas)."""
        z = xyz[con][:, :, 2]
        return z.max(1) - z.min(1)

    def conducao(self, xyz, con):
        """Matrizes de condução com k = 1 (ne, nn, nn), volumes (ne,) e pesos
        HRZ da capacidade concentrada (ne, nn), que somam 1 em cada elemento."""
        xi, w = self.pontos()
        G, detJ = self.gradientes(xyz, con, xi)
        dv = detJ * w[None, :]
        K = np.einsum("ep,epia,epja->eij", dv, G, G)
        vol = dv.sum(1)
        N = self.N(xi)
        diag = np.einsum("ep,pi->ei", dv, N * N)
        return K, np.abs(vol), diag / diag.sum(1, keepdims=True)

    @staticmethod
    def B(G):
        """B (..., 6, 3 nn) a partir de G (..., nn, 3); Voigt xx, yy, zz, xy, yz, xz."""
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

    def rigidez(self, xyz, con, E, nu, eps_th_nos=None, D=None):
        """Rigidez (ne, 3nn, 3nn) e forças térmicas (ne, 3nn) se eps_th_nos (ne, nn) for dado.
        D (ne, 6, 6) substitui o material isotrópico (E, nu)."""
        xi, w = self.pontos()
        G, detJ = self.gradientes(xyz, con, xi)
        Nn = self.N(xi)
        D = matriz_D(E, nu) if D is None else D
        ne, nd = len(con), 3 * self.nn
        K = np.zeros((ne, nd, nd))
        f = np.zeros((ne, nd)) if eps_th_nos is not None else None
        um = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        for p in range(len(w)):
            B = self.B(G[:, p])
            dv = detJ[:, p] * w[p]
            DB = np.einsum("eij,ejk->eik", D, B)
            K += dv[:, None, None] * np.einsum("eji,ejk->eik", B, DB)
            if f is not None:
                sig0 = (D @ um) * (eps_th_nos @ Nn[p])[:, None]
                f += dv[:, None] * np.einsum("eji,ej->ei", B, sig0)
        return K, f

    def forca_termica(self, xyz, con, E, nu, eps_th_nos):
        xi, w = self.pontos()
        G, detJ = self.gradientes(xyz, con, xi)
        Nn = self.N(xi)
        D = matriz_D(E, nu)
        um = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        f = np.zeros((len(con), 3 * self.nn))
        for p in range(len(w)):
            B = self.B(G[:, p])
            sig0 = (D @ um) * (eps_th_nos @ Nn[p])[:, None]
            f += (detJ[:, p] * w[p])[:, None] * np.einsum("eji,ej->ei", B, sig0)
        return f

    def no_centro(self, xyz, con):
        """B (ne, 6, 3nn), N (nn,) e G (ne, nn, 3) no centro do elemento."""
        G, _ = self.gradientes(xyz, con, self.centro[None, :])
        return self.B(G[:, 0]), self.N(self.centro[None, :])[0], G[:, 0]

    def pontos_B(self, xyz, con):
        """Lista [(dv (ne,), B (ne, 6, 3nn), N (nn,))] nos pontos de Gauss."""
        xi, w = self.pontos()
        G, detJ = self.gradientes(xyz, con, xi)
        Nn = self.N(xi)
        return [(detJ[:, p] * w[p], self.B(G[:, p]), Nn[p]) for p in range(len(w))]

    def tensoes(self, xyz, con, E, nu, u_el, eps_th_nos):
        B, N, _ = self.no_centro(xyz, con)
        eps = np.einsum("eij,ej->ei", B, u_el)
        eps[:, :3] -= (eps_th_nos @ N)[:, None]
        return np.einsum("eij,ej->ei", matriz_D(E, nu), eps)

    # ----------------------------------------------------------------- faces
    def areas_faces(self, xyz, con, seletor, consistente=False):
        """Áreas nodais das faces de contorno escolhidas por
        seletor(coordenadas dos nós da face (nf, k, 3)) -> máscara (nf,).
        consistente=False: pesos HRZ (positivos; trocas de calor no contorno);
        consistente=True: integral de N_i (vetor de cargas consistente de uma
        pressão uniforme; nos cantos da face quadrática fica negativo).
        Nos elementos lineares de faces planas os dois coincidem.
        Devolve vetor (n_nos,)."""
        area = np.zeros(len(xyz))
        for nos_face, mapa, regra in self.faces_int:
            fc = con[:, nos_face]
            sel = seletor(xyz[fc])
            if not sel.any():
                continue
            ce = con[sel]
            st, ws = regra
            xi, dxi_ds, dxi_dt = mapa(st)
            N = self.N(xi)                                   # (np, nn)
            dN = self.dN(xi)                                 # (np, nn, 3)
            X = xyz[ce]                                      # (ne, nn, 3)
            dxds = np.einsum("eia,pib,pb->epa", X, dN, dxi_ds)
            dxdt = np.einsum("eia,pib,pb->epa", X, dN, dxi_dt)
            dA = np.linalg.norm(np.cross(dxds, dxdt), axis=2) * ws[None, :]     # (ne, np)
            Af = dA.sum(1)
            if consistente:
                pesos = np.einsum("ep,pi->ei", dA, N[:, nos_face])
            else:
                d = np.einsum("ep,pi->ei", dA, N[:, nos_face] ** 2)
                pesos = d / d.sum(1, keepdims=True) * Af[:, None]
            area += np.bincount(ce[:, nos_face].ravel(), weights=pesos.ravel(), minlength=len(xyz))
        return area


def _face_lado_cunha(i, j, regra_lado, linha):
    """Face lateral da cunha sobre a aresta (i, j) do triângulo (L_i: 1 -> 0, L_j: 0 -> 1)."""
    def mapa(st):
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
    for s, ws in zip(*regra_lado):
        for t, wt in zip(*linha):
            p.append([s, t]); w.append(ws * wt)
    return mapa, (np.array(p), np.array(w))


def _face_tri_cunha(zeta, regra_tri):
    def mapa(st):
        xi = np.c_[st[:, 0], st[:, 1], np.full(len(st), zeta)]
        ds = np.tile([1.0, 0.0, 0.0], (len(st), 1)); dt = np.tile([0.0, 1.0, 0.0], (len(st), 1))
        return xi, ds, dt
    return mapa, regra_tri


def _regra_01(linha):
    """Regra de Gauss de [-1, 1] levada para [0, 1]."""
    return ((linha[0] + 1) / 2, linha[1] / 2)


class W6(TipoElemento):
    nome, nn, vtk = "W6", 6, 13
    vtk_ordem = [0, 2, 1, 3, 5, 4]
    base, topo = [0, 1, 2], [3, 4, 5]
    faces_vis = [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]]
    centro = np.array([1 / 3, 1 / 3, 0.0])

    def __init__(self):
        lado = _regra_01((G2, W2))
        self.faces_int = [
            ([0, 2, 1], *_face_tri_cunha(-1.0, TRI3)), ([3, 4, 5], *_face_tri_cunha(1.0, TRI3)),
            ([0, 1, 4, 3], *_face_lado_cunha(0, 1, lado, (G2, W2))),
            ([1, 2, 5, 4], *_face_lado_cunha(1, 2, lado, (G2, W2))),
            ([2, 0, 3, 5], *_face_lado_cunha(2, 0, lado, (G2, W2)))]

    def pontos(self):
        return _cunha(TRI3, (G2, W2))

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


class H8(TipoElemento):
    nome, nn, vtk = "H8", 8, 12
    vtk_ordem = list(range(8))
    base, topo = [0, 1, 2, 3], [4, 5, 6, 7]
    faces_vis = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]
    centro = np.zeros(3)
    _s = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                   [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)

    def __init__(self):
        q = []
        for a, wa in zip(G2, W2):
            for b, wb in zip(G2, W2):
                q.append([a, b, wa * wb])
        q = np.array(q)
        regra = (q[:, :2], q[:, 2])

        def face(fixo, valor, u, v):
            def mapa(st):
                xi = np.zeros((len(st), 3)); xi[:, fixo] = valor; xi[:, u] = st[:, 0]; xi[:, v] = st[:, 1]
                ds = np.zeros((len(st), 3)); ds[:, u] = 1
                dt = np.zeros((len(st), 3)); dt[:, v] = 1
                return xi, ds, dt
            return mapa
        self.faces_int = [([0, 3, 2, 1], face(2, -1, 0, 1), regra), ([4, 5, 6, 7], face(2, 1, 0, 1), regra),
                          ([0, 1, 5, 4], face(1, -1, 0, 2), regra), ([1, 2, 6, 5], face(0, 1, 1, 2), regra),
                          ([2, 3, 7, 6], face(1, 1, 0, 2), regra), ([3, 0, 4, 7], face(0, -1, 1, 2), regra)]

    def pontos(self):
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


class W15(TipoElemento):
    nome, nn, vtk = "W15", 15, 26
    # VTK_QUADRATIC_WEDGE: cantos como VTK_WEDGE (0, 2, 1, 3, 5, 4), depois arestas na mesma ordem
    vtk_ordem = [0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13]
    base, topo = [0, 1, 2], [3, 4, 5]
    faces_vis = W6.faces_vis
    centro = np.array([1 / 3, 1 / 3, 0.0])

    def __init__(self):
        lado = _regra_01((G3, W3))
        self.faces_int = [
            ([0, 2, 1, 8, 7, 6], *_face_tri_cunha(-1.0, TRI6)), ([3, 4, 5, 9, 10, 11], *_face_tri_cunha(1.0, TRI6)),
            ([0, 1, 4, 3, 6, 13, 9, 12], *_face_lado_cunha(0, 1, lado, (G3, W3))),
            ([1, 2, 5, 4, 7, 14, 10, 13], *_face_lado_cunha(1, 2, lado, (G3, W3))),
            ([2, 0, 3, 5, 8, 12, 11, 14], *_face_lado_cunha(2, 0, lado, (G3, W3)))]

    def pontos(self):
        return _cunha(TRI6, (G3, W3))

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
def tipo(nome: str) -> TipoElemento:
    return {"W6": W6, "H8": H8, "W15": W15}[nome]()
