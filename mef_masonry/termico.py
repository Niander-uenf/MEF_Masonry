"""Condução de calor transiente (2D ou 3D) com convecção e radiação no contorno.

Formulação: Galerkin com capacidade concentrada (lumped) e Euler implícito.
Propriedades e coeficientes de radiação são avaliados com a temperatura do
passo anterior (semi-implícito); o passo dt controla a precisão (o exemplo de
verificação mostra a convergência).
"""
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .malha import AR, CERAMICA, JUNTA, CAPEAMENTO, ACO, ARGAMASSA, GRAUTE, Malha2D, Malha3D
from .parametros import Ceramica, Argamassa, Graute, JuntaSeca, Aco, ArFuros, SIGMA
from . import elementos as el


@dataclass
class Materiais:
    ceramica: Ceramica = field(default_factory=Ceramica)
    argamassa: Argamassa = field(default_factory=Argamassa)
    graute: Graute = field(default_factory=Graute)
    junta: JuntaSeca = field(default_factory=JuntaSeca)
    capeamento: JuntaSeca = field(default_factory=JuntaSeca)
    aco: Aco = field(default_factory=Aco)
    ar: ArFuros = field(default_factory=ArFuros)
    diametro_furo: float = 0.084

    def k(self, mat, Te):
        k = np.empty(len(mat))
        cer = self.ceramica
        for m, sel in _grupos(mat):
            T = Te[sel]
            if m == CERAMICA:
                k[sel] = cer.k_T(T)
            elif m == ARGAMASSA:
                k[sel] = self.argamassa.k_T(T)
            elif m == GRAUTE:
                k[sel] = self.graute.k_T(T)
            elif m == JUNTA:
                k[sel] = cer.k_T(T) * self.junta.fator_k
            elif m == CAPEAMENTO:
                k[sel] = cer.k_T(T) * self.capeamento.fator_k
            elif m == ACO:
                k[sel] = self.aco.k
            elif m == AR:
                TK = T + 273.15
                eps = self.ar.emissividade_paredes
                k[sel] = (0.0241 * (TK / 273.15) ** 0.81
                          + 4 * SIGMA * TK ** 3 * self.diametro_furo * eps / (2 - eps))
        return k

    @property
    def capacidade_variavel(self):
        return (self.ceramica.capacidade_variavel or self.argamassa.capacidade_variavel
                or self.graute.capacidade_variavel)

    def rho_c(self, mat, Te=None):
        """Capacidade térmica volumétrica por elemento (J/m³.K); com Te, usa as
        tabelas de c(T) e rho(T) (método da capacidade aparente)."""
        rc = np.empty(len(mat))
        for m, sel in _grupos(mat):
            if m in (CERAMICA, JUNTA, CAPEAMENTO, ARGAMASSA, GRAUTE):
                mt = {ARGAMASSA: self.argamassa, GRAUTE: self.graute}.get(m, self.ceramica)
                if Te is None:
                    rc[sel] = mt.rho * mt.c
                else:
                    T = Te[sel]
                    rc[sel] = mt.rho_T(T) * mt.c_T(T)
            elif m == ACO:
                rc[sel] = self.aco.rho * self.aco.c
            else:
                rc[sel] = self.ar.rho_c
        return rc


def _grupos(mat):
    for m in np.unique(mat):
        yield int(m), mat == m


@dataclass
class Contorno:
    """Grupo de faces expostas a um gás à temperatura T_gas(t) (°C).
    area_nos: área (3D, m²) ou comprimento (2D, m) associada a cada nó.
    faces: opcional, (conectividade das faces (nf, k), matrizes de massa das
    faces (nf, k, k)) para a troca de calor consistente (elementos quadráticos)."""
    area_nos: np.ndarray
    temperatura: Callable[[float], float]
    h_conveccao: float
    emissividade: float
    nome: str = ""
    faces: Optional[tuple] = None


@dataclass
class ResultadoTermico:
    t: np.ndarray
    T_gas: np.ndarray
    sondas: Dict[str, np.ndarray]
    T_final: np.ndarray
    campos: Dict[float, np.ndarray]
    parou_em: Optional[float] = None


def estrutura_csr(conect, n):
    """Padrão esparso e mapa elemento -> posição em data (para montagem rápida)."""
    nen = conect.shape[1]
    r = np.repeat(conect, nen, axis=1).ravel().astype(np.int64)
    c = np.tile(conect, (1, nen)).ravel().astype(np.int64)
    chave = r * n + c
    unicas, mapa = np.unique(chave, return_inverse=True)
    linhas = unicas // n
    indices = (unicas % n).astype(np.int32)
    indptr = np.searchsorted(linhas, np.arange(n + 1)).astype(np.int32)
    diag = np.searchsorted(unicas, np.arange(n, dtype=np.int64) * (n + 1))
    return indices, indptr, mapa.ravel(), diag, len(unicas)


def _cg(A, b, x0, rtol=1e-10):
    d = A.diagonal()
    M = spla.LinearOperator(A.shape, matvec=lambda v: v / d)
    try:
        x, info = spla.cg(A, b, x0=x0, rtol=rtol, atol=0.0, M=M, maxiter=5000)
    except TypeError:  # SciPy antigo
        x, info = spla.cg(A, b, x0=x0, tol=rtol, M=M, maxiter=5000)
    if info != 0:
        x = spla.spsolve(A.tocsc(), b)
    return x


def resolver(conect: np.ndarray, Ke1: np.ndarray, volumes: np.ndarray, mat: np.ndarray,
             n_nos: int, M: Materiais, contornos: Sequence[Contorno], t_final: float,
             dt: float = 30.0, T_inicial: float = 30.0, sondas: Optional[Dict] = None,
             registrar_cada: float = 60.0, instantes_campo: Sequence[float] = (),
             parar_se: Optional[Callable[[float, np.ndarray], bool]] = None,
             verboso: bool = False, pesos_cap: Optional[np.ndarray] = None,
             cap_elem: Optional[np.ndarray] = None) -> ResultadoTermico:
    """Integra a condução transiente de t = 0 a t_final.

    conect, Ke1, volumes : conectividade, matrizes de condução com k = 1 e
                           volumes (ou áreas, em 2D) dos elementos
    sondas    : {nome: índice de nó ou lista de nós (média)}
    instantes_campo : tempos em que o campo completo é guardado
    parar_se  : função (t, T) -> bool para encerrar antes (ex.: critério de isolamento)
    pesos_cap : frações (ne, nen) da capacidade concentrada em cada nó (padrão:
                partes iguais)
    cap_elem  : matrizes de capacidade consistentes por unidade de rho c,
                (ne, nen, nen); se dadas, substituem a capacidade concentrada
    Contornos com .faces usam a matriz de troca de calor consistente
    (h integrado com N_i N_j em cada face, h na temperatura média da face).
    """
    nen = conect.shape[1]
    indices, indptr, mapa, diag, nnz = estrutura_csr(conect, n_nos)
    if pesos_cap is None:
        pesos_cap = np.full(conect.shape, 1.0 / nen)

    def capacidade(rc):
        """Capacidade (consistente: dados CSR e matriz; concentrada: vetor nodal)."""
        if cap_elem is not None:
            dC = np.bincount(mapa, weights=(rc[:, None, None] * cap_elem).ravel(), minlength=nnz)
            return dC, sp.csr_matrix((dC, indices, indptr), shape=(n_nos, n_nos)), None
        return None, None, np.bincount(conect.ravel(), weights=((rc * volumes)[:, None] * pesos_cap).ravel(),
                                       minlength=n_nos)

    variavel = M.capacidade_variavel
    dados_C, Cmat, C = capacidade(M.rho_c(mat))
    pares = {}
    for i, g in enumerate(contornos):
        if g.faces is not None:
            fc = g.faces[0]
            k = fc.shape[1]
            pares[i] = (np.repeat(fc, k, 1).ravel(), np.tile(fc, (1, k)).ravel())
    T = np.full(n_nos, float(T_inicial))
    sondas = sondas or {}
    sondas = {k: np.atleast_1d(v) for k, v in sondas.items()}
    hist_t, hist_g, hist_s = [], [], {k: [] for k in sondas}
    campos = {}
    pend = sorted(instantes_campo)

    def registrar(t):
        hist_t.append(t)
        hist_g.append(contornos[0].temperatura(t) if contornos else np.nan)
        for k, v in sondas.items():
            hist_s[k].append(T[v].mean())

    registrar(0.0)
    t, prox = 0.0, registrar_cada
    n_passos = int(np.ceil(t_final / dt - 1e-9))
    parou = None
    for passo in range(n_passos):
        tn = min(t_final, t + dt)
        h = tn - t
        Te = T[conect].mean(axis=1)
        dados = np.bincount(mapa, weights=(M.k(mat, Te)[:, None, None] * Ke1).ravel(), minlength=nnz)
        if variavel:                     # capacidade aparente na temperatura do passo anterior
            dados_C, Cmat, C = capacidade(M.rho_c(mat, Te))
        if C is None:
            dados += dados_C / h
            diagonal = np.zeros(n_nos)
            b = (Cmat @ T) / h
        else:
            diagonal = C / h
            b = C / h * T
        TK = T + 273.15
        extra = None
        for i, g in enumerate(contornos):
            Tg = g.temperatura(tn)
            TgK = Tg + 273.15
            if g.faces is not None:
                fc, Mf = g.faces
                TKf = TK[fc].mean(1)
                hb = g.h_conveccao + g.emissividade * SIGMA * (TgK ** 2 + TKf ** 2) * (TgK + TKf)
                Hf = (hb[:, None, None] * Mf).ravel()
                H = sp.csr_matrix((Hf, pares[i]), shape=(n_nos, n_nos))
                extra = H if extra is None else extra + H
                b = b + Tg * np.asarray(H.sum(axis=1)).ravel()
            else:
                hb = g.h_conveccao + g.emissividade * SIGMA * (TgK ** 2 + TK ** 2) * (TgK + TK)
                diagonal = diagonal + hb * g.area_nos
                b = b + hb * g.area_nos * Tg
        dados[diag] += diagonal
        A = sp.csr_matrix((dados, indices, indptr), shape=(n_nos, n_nos))
        if extra is not None:
            A = (A + extra).tocsr()
        T = _cg(A, b, T)
        t = tn
        while pend and pend[0] <= t + 1e-9:
            campos[pend.pop(0)] = T.copy()
        if t >= prox - 1e-9 or t >= t_final - 1e-9:
            registrar(t)
            prox += registrar_cada
        if verboso and passo % 50 == 0:
            print(f"  t = {t / 60:6.1f} min")
        if parar_se is not None and parar_se(t, T):
            parou = t
            if hist_t[-1] != t:
                registrar(t)
            break
    return ResultadoTermico(t=np.array(hist_t), T_gas=np.array(hist_g),
                            sondas={k: np.array(v) for k, v in hist_s.items()},
                            T_final=T, campos=campos, parou_em=parou)


# ----------------------------------------------------------------------------
# Contornos e atalhos para as malhas do pacote
# ----------------------------------------------------------------------------
def seletor_faces(m: Malha3D, faces, tol=1e-9):
    """Função que marca as faces de elemento (nf, k, 3) contidas nos planos
    pedidos: x0, x1, y0, y1 (laterais), base e topo. Com simetria, o plano
    y = W/2 nunca é exposto."""
    m2 = m.m2d
    L, Wm = m2.comprimento, m2.largura_modelo
    z0, z1 = m.zs[0], m.zs[-1]

    def sel(P):
        s = np.zeros(len(P), bool)
        if "x0" in faces: s |= (np.abs(P[:, :, 0]) < tol).all(1)
        if "x1" in faces: s |= (np.abs(P[:, :, 0] - L) < tol).all(1)
        if "y0" in faces: s |= (np.abs(P[:, :, 1]) < tol).all(1)
        if "y1" in faces and not m2.simetria: s |= (np.abs(P[:, :, 1] - Wm) < tol).all(1)
        if "base" in faces: s |= (np.abs(P[:, :, 2] - z0) < tol).all(1)
        if "topo" in faces: s |= (np.abs(P[:, :, 2] - z1) < tol).all(1)
        return s
    return sel


def areas_laterais_3d(m: Malha3D, faces=("x0", "x1", "y0", "y1"), tol=1e-9):
    """Área nodal das faces laterais expostas (x = 0, x = L, y = 0, y = W),
    integrada nas faces dos elementos (W6, H8 ou W15; pesos HRZ).
    Com simetria, a face y = W/2 nunca é exposta."""
    return m.elem.areas_faces(m.xyz, m.cun, seletor_faces(m, faces, tol))


def areas_topo_base_3d(m: Malha3D, qual="ambas"):
    faces = {"base": ("base",), "topo": ("topo",)}.get(qual, ("base", "topo"))
    return m.elem.areas_faces(m.xyz, m.cun, seletor_faces(m, faces))


def comprimentos_contorno_2d(m2: Malha2D, faces=("x0", "x1", "y0", "y1"), tol=1e-9):
    if m2.vertical:      # seção vertical: as bordas de cima e de baixo são cortes periódicos, nunca expostas
        faces = tuple(f for f in faces if f in ("y0", "y1"))
    ar = m2.arestas_contorno()
    p = m2.xy[ar]
    L, Wm = m2.comprimento, m2.largura_modelo
    sel = np.zeros(len(ar), bool)
    if "x0" in faces: sel |= (np.abs(p[:, :, 0]) < tol).all(1)
    if "x1" in faces: sel |= (np.abs(p[:, :, 0] - L) < tol).all(1)
    if "y0" in faces: sel |= (np.abs(p[:, :, 1]) < tol).all(1)
    if "y1" in faces and not m2.simetria: sel |= (np.abs(p[:, :, 1] - Wm) < tol).all(1)
    ar = ar[sel]
    comp = np.linalg.norm(m2.xy[ar[:, 0]] - m2.xy[ar[:, 1]], axis=1)
    return np.bincount(ar.ravel(), weights=np.repeat(comp / 2, 2), minlength=m2.n_nos)


def material_secao(m2: Malha2D, geo=None):
    """Materiais da seção 2D (fiada 0): ar nos furos, graute nos furos
    grauteados (se geo for dado) e bloco no resto. Na seção vertical (furos
    horizontais), as meias juntas de cima e de baixo são argamassa (ou junta seca)."""
    if m2.vertical:
        junta = ARGAMASSA if (geo is None or geo.tipo_junta == "argamassa") else JUNTA
        return np.where(m2.marca == 0, AR, np.where(m2.marca > 0, junta, CERAMICA)).astype(np.int8)
    if geo is not None:
        from .malha import estados_fiada
        furo, gr = estados_fiada(m2, geo, 0)
        return np.where(gr, GRAUTE, np.where(furo, AR, CERAMICA)).astype(np.int8)
    return np.where(m2.marca == 0, AR, CERAMICA).astype(np.int8)


def faces_contorno(m: Malha3D, area_nos: np.ndarray):
    """Faces de elemento cujos nós têm todos área > 0 em area_nos (as faces do
    contorno que gerou area_nos) e suas matrizes de massa, integral de N_i N_j dA."""
    tp = m.elem
    F, MM = [], []
    for nos_face, mapa, regra in tp.faces_int:
        fc = m.cun[:, nos_face]
        sel = (area_nos[fc] > 0).all(1)
        if not sel.any():
            continue
        ce = m.cun[sel]
        st, ws = regra
        xi, dxi_ds, dxi_dt = mapa(st)
        N = tp.N(xi)[:, nos_face]
        dN = tp.dN(xi)
        X = m.xyz[ce]
        a = np.einsum("eia,pib,pb->epa", X, dN, dxi_ds)
        b = np.einsum("eia,pib,pb->epa", X, dN, dxi_dt)
        dA = np.linalg.norm(np.cross(a, b), axis=2) * ws[None, :]
        F.append(ce[:, nos_face]); MM.append(np.einsum("ep,pi,pj->eij", dA, N, N))
    if not F:
        return None
    return np.vstack(F), np.vstack(MM)


def resolver_3d(m: Malha3D, M: Materiais, contornos, t_final, **kw) -> ResultadoTermico:
    """Condução no prisma 3D. W6 e H8: capacidade e trocas de calor concentradas
    nos nós. W15: matrizes consistentes de capacidade e de troca de calor (a
    concentração nos nós das faces quadráticas desvia o fluxo de entrada)."""
    if m.tipo == "W6":
        Ke1, vol = el.conducao_w6(m.xyz, m.cun)
        return resolver(m.cun, Ke1, vol, m.mat, m.n_nos, M, contornos, t_final, **kw)
    tp = m.elem
    Ke1, vol, pesos = tp.conducao(m.xyz, m.cun)
    if m.tipo == "H8":
        return resolver(m.cun, Ke1, vol, m.mat, m.n_nos, M, contornos, t_final, pesos_cap=pesos, **kw)
    xi, w = tp.pontos()
    _, detJ = tp.gradientes(m.xyz, m.cun, xi)
    N = tp.N(xi)
    Ce = np.einsum("ep,pi,pj->eij", detJ * w[None, :], N, N)
    ctr = []
    for g in contornos:
        if g.faces is None:
            g = Contorno(g.area_nos, g.temperatura, g.h_conveccao, g.emissividade, g.nome,
                         faces_contorno(m, g.area_nos))
        ctr.append(g)
    return resolver(m.cun, Ke1, vol, m.mat, m.n_nos, M, ctr, t_final, cap_elem=Ce, **kw)


def resolver_2d(m2: Malha2D, M: Materiais, contornos, t_final, mat=None, **kw) -> ResultadoTermico:
    Ke1, area = el.conducao_t3(m2.xy, m2.tri)
    mat = material_secao(m2) if mat is None else mat
    return resolver(m2.tri, Ke1, area, mat, m2.n_nos, M, contornos, t_final, **kw)
