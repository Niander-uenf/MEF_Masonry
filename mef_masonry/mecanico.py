"""Análise estrutural linear com deformação térmica.

3D: cunhas W6, hexaedros H8 ou cunhas W15 (malha.tipo); carga de compressão aplicada na face superior da placa
de topo; base com restrição apenas vertical (livre para dilatar no plano),
mais o mínimo contra movimento de corpo rígido. O ar dos furos não entra.

2D: estado plano generalizado na seção (corpo livre), para estudos rápidos.
"""
from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .malha import AR, CERAMICA, JUNTA, CAPEAMENTO, ACO, ARGAMASSA, GRAUTE, Malha2D, Malha3D, no_mais_proximo
from .termico import Materiais
from . import elementos as el


def _propriedades(mat, Te, M: Materiais):
    E = np.empty(len(mat)); nu = np.empty(len(mat))
    for m in np.unique(mat):
        s = mat == m
        if m == CERAMICA:
            E[s] = M.ceramica.E_T(Te[s]); nu[s] = M.ceramica.nu
        elif m == ARGAMASSA:
            E[s] = M.argamassa.E_T(Te[s]); nu[s] = M.argamassa.nu
        elif m == GRAUTE:
            E[s] = M.graute.E_T(Te[s]); nu[s] = M.graute.nu
        elif m == JUNTA:
            E[s] = M.junta.E; nu[s] = M.junta.nu
        elif m == CAPEAMENTO:
            E[s] = M.capeamento.E; nu[s] = M.capeamento.nu
        elif m == ACO:
            E[s] = M.aco.E; nu[s] = M.aco.nu
    return E, nu


def _eps_termica_nos(mat, T_el_nos, M: Materiais):
    """Deformação térmica livre nos nós de cada elemento (ne, nen), pela lei do material.
    Juntas secas e capeamento seguem o bloco; argamassa e aço usam a sua lei."""
    eps = np.empty_like(T_el_nos)
    for m in np.unique(mat):
        s = mat == m
        if m == ACO:
            eps[s] = M.aco.alfa * (T_el_nos[s] - M.ceramica.T_ref)
        elif m == ARGAMASSA:
            eps[s] = M.argamassa.eps_termica(T_el_nos[s])
        elif m == GRAUTE:
            eps[s] = M.graute.eps_termica(T_el_nos[s])
        else:
            eps[s] = M.ceramica.eps_termica(T_el_nos[s])
    return eps


def _fatorar(K):
    """Fatora K (simétrica positiva definida). Usa pypardiso se estiver
    instalado (mais rápido, Windows/Linux); senão SuperLU sem pivotamento,
    com reordenação de grau mínimo (seguro para matriz SPD)."""
    try:
        import pypardiso
        Kc = K.tocsr()
        return lambda b: pypardiso.spsolve(Kc, b)
    except ImportError:
        lu = spla.splu(K.tocsc(), permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                       options=dict(SymmetricMode=True))
        return lu.solve


@dataclass
class ResultadoEstrutural:
    u: np.ndarray            # (n_nos, 3) m
    sig: np.ndarray          # (ne, 6) Pa, nan no ar
    princ: np.ndarray        # (ne, 3) Pa, crescente
    ativos: np.ndarray       # elementos sólidos


def resolver_3d(m: Malha3D, M: Materiais, T_nos: Optional[np.ndarray] = None,
                carga_N: float = 0.0, incluir_termico: bool = True,
                placas: str = "rigidas", blocos: Optional[int] = None,
                base_engastada: bool = False, regiao_carga=None) -> dict:
    """Resolve até dois casos com a mesma rigidez: 'termico' (campo T_nos, topo
    livre) e 'carga' (compressão carga_N). Devolve {'termico': R, 'carga': R};
    os casos se somam com combinar().

    placas = 'rigidas' (padrão): as placas de aço saem do modelo mecânico e os
             pratos da prensa são apoios rígidos sem atrito — base com u_z = 0;
             no caso 'carga', topo com u_z uniforme, ajustado para a reação
             igualar a carga. Menos incógnitas e sem contraste aço/cerâmica.
    placas = 'aco': as placas entram como aço e a carga é uma pressão uniforme
             na face superior da placa de topo (ou da cerâmica, se a malha não
             tiver placas).
    base_engastada: todos os deslocamentos da base presos (apoio fixo, como o
             "fixed support" dos programas comerciais), em vez de só u_z.
    regiao_carga: com placas = 'aco', função (x, y) -> máscara que escolhe a
             parte da face superior que recebe a pressão (centro de cada face);
             None = a face toda.
    carga_N é a carga total do ensaio; com simetria, metade é aplicada."""
    n = m.n_nos
    if T_nos is None:
        T_nos = np.full(n, M.ceramica.T_ref)
    ativos = m.mat != AR
    if placas == "rigidas":
        ativos &= m.mat != ACO
    cun = m.cun[ativos]; mat = m.mat[ativos]
    usados = np.unique(cun)
    Tel = T_nos[cun]
    E, nu = _propriedades(mat, Tel.mean(1), M)
    eps_n = _eps_termica_nos(mat, Tel, M)
    tp = m.elem
    nd = 3 * tp.nn
    w6 = m.tipo == "W6"
    if blocos is None:
        blocos = max(500, int(20000 * 324 / nd ** 2))

    ndof = 3 * n
    K = sp.csr_matrix((ndof, ndof))
    F_th = np.zeros(ndof)
    for i0 in range(0, len(cun), blocos):
        s = slice(i0, i0 + blocos)
        eth = eps_n[s] if incluir_termico else None
        if w6:
            Ke, fe = el.rigidez_w6(m.xyz, cun[s], E[s], nu[s], eth)
        else:
            Ke, fe = tp.rigidez(m.xyz, cun[s], E[s], nu[s], eth)
        dofs = (3 * cun[s][:, :, None] + np.arange(3)).reshape(len(Ke), nd)
        r = np.repeat(dofs, nd, axis=1).ravel(); c = np.tile(dofs, (1, nd)).ravel()
        K = K + sp.csr_matrix((Ke.ravel(), (r, c)), shape=(ndof, ndof))
        if fe is not None:
            F_th += np.bincount(dofs.ravel(), weights=fe.ravel(), minlength=ndof)

    zmin, zmax = m.xyz[usados, 2].min(), m.xyz[usados, 2].max()
    fator = 0.5 if m.m2d.simetria else 1.0
    base = usados[np.isclose(m.xyz[usados, 2], zmin)]
    topo_nos = usados[np.isclose(m.xyz[usados, 2], zmax)]

    # condições de contorno comuns
    fixos = np.zeros(ndof, bool)
    inativos = np.setdiff1d(np.arange(n), usados)
    fixos[(3 * inativos[:, None] + np.arange(3)).ravel()] = True
    fixos[3 * base + 2] = True
    if base_engastada:
        fixos[(3 * base[:, None] + np.arange(3)).ravel()] = True
    ref = base[np.argmin(np.linalg.norm(m.xyz[base, :2] - [m.m2d.comprimento / 2, 0.0], axis=1))]
    fixos[3 * ref] = True
    if m.m2d.simetria:
        sim = usados[np.isclose(m.xyz[usados, 1], m.m2d.largura_modelo)]
        fixos[3 * sim + 1] = True
    else:
        fixos[3 * ref + 1] = True
        outro = base[np.argmax(np.abs(m.xyz[base, 0] - m.xyz[ref, 0]))]
        fixos[3 * outro + 1] = True

    resultados = {}
    casos = []
    if incluir_termico:
        casos.append(("termico", F_th, None))
    if carga_N:
        if placas == "rigidas":
            # topo com deslocamento vertical uniforme (prato rígido sem atrito)
            u_p = np.zeros(ndof); u_p[3 * topo_nos + 2] = -1e-3
            casos.append(("carga", -(K @ u_p), u_p))
        else:
            def _sel(P):
                s = np.isclose(P[:, :, 2], zmax).all(1)
                if regiao_carga is not None:
                    c = P.mean(1)
                    s &= np.asarray(regiao_carga(c[:, 0], c[:, 1]), bool)
                return s
            a = tp.areas_faces(m.xyz, cun, _sel, consistente=True)
            F_c = np.zeros(ndof)
            F_c[2::3] = -carga_N * fator * a / a.sum()
            casos.append(("carga", F_c, None))

    fat = {}
    for nome, F, u_p in casos:
        fx = fixos.copy()
        if u_p is not None:
            fx[3 * topo_nos + 2] = True
        chave = fx.tobytes()
        livres = np.nonzero(~fx)[0]
        if chave not in fat:
            fat[chave] = _fatorar(K[livres][:, livres])
        u = np.zeros(ndof) if u_p is None else u_p.copy()
        u[livres] = fat[chave](F[livres])
        if u_p is not None:
            reacao = (K @ u)[3 * topo_nos + 2].sum()        # força do prato sobre o prisma
            u *= -carga_N * fator / reacao
        u3 = u.reshape(n, 3)
        epsn = eps_n if nome == "termico" else np.zeros_like(eps_n)
        sig = np.full((len(m.cun), 6), np.nan)
        s_at = np.empty((len(cun), 6))
        for i0 in range(0, len(cun), blocos):
            s = slice(i0, i0 + blocos)
            ue = u3[cun[s]].reshape(-1, nd)
            s_at[s] = (el.tensoes_w6(m.xyz, cun[s], E[s], nu[s], ue, epsn[s]) if w6
                       else tp.tensoes(m.xyz, cun[s], E[s], nu[s], ue, epsn[s]))
        sig[ativos] = s_at
        pr = np.full((len(m.cun), 3), np.nan)
        pr[ativos] = el.principais(s_at)
        resultados[nome] = ResultadoEstrutural(u=u3, sig=sig, princ=pr, ativos=ativos)
    return resultados


def combinar(r1: ResultadoEstrutural, r2: ResultadoEstrutural) -> ResultadoEstrutural:
    sig = r1.sig + r2.sig
    pr = np.full_like(r1.princ, np.nan)
    pr[r1.ativos] = el.principais(sig[r1.ativos])
    return ResultadoEstrutural(u=r1.u + r2.u, sig=sig, princ=pr, ativos=r1.ativos)


# ----------------------------------------------------------------------------
# Grandezas para comparar com o ensaio
# ----------------------------------------------------------------------------
def media_uz_plano(m: Malha3D, u, z, nos=None):
    nos = np.arange(m.n_nos) if nos is None else nos
    sel = nos[np.isclose(m.xyz[nos, 2], z)]
    return u[sel, 2].mean()


def encurtamento_entre_placas(m: Malha3D, R: ResultadoEstrutural, H: float):
    """Variação da distância entre as faces das placas (z = 0 e z = H), m.
    Negativo = encurtamento."""
    usados = np.unique(m.cun[R.ativos])
    return media_uz_plano(m, R.u, H, usados) - media_uz_plano(m, R.u, 0.0, usados)


def deformacao_extensometro(m: Malha3D, R: ResultadoEstrutural, z1: float, z2: float,
                            y_face: float = 0.0):
    """Deformação média vertical entre as cotas z1 e z2 na face externa y = y_face
    (como um extensômetro ou DIC no terço central). Usa os planos de nós mais próximos."""
    z1 = m.zs[np.argmin(np.abs(m.zs - z1))]; z2 = m.zs[np.argmin(np.abs(m.zs - z2))]
    face = np.nonzero(np.isclose(m.xyz[:, 1], y_face))[0]
    return (media_uz_plano(m, R.u, z2, face) - media_uz_plano(m, R.u, z1, face)) / (z2 - z1), z1, z2


def tracao_maxima(m: Malha3D, R: ResultadoEstrutural, materiais=(CERAMICA,), zmin=None, zmax=None,
                  percentil=99.0):
    """Máxima tensão principal de tração (Pa) na cerâmica, com o percentil e a posição."""
    sel = np.isin(m.mat, materiais)
    if zmin is not None or zmax is not None:
        zc = m.xyz[m.cun].mean(1)[:, 2]
        if zmin is not None: sel &= zc >= zmin
        if zmax is not None: sel &= zc <= zmax
    s1 = R.princ[sel, 2]
    i = np.nanargmax(s1)
    pos = m.xyz[m.cun[np.nonzero(sel)[0][i]]].mean(0)
    return float(s1[i]), float(np.nanpercentile(s1, percentil)), pos


# ----------------------------------------------------------------------------
# 2D: estado plano generalizado
# ----------------------------------------------------------------------------
def resolver_secao_gps(m2: Malha2D, M: Materiais, T_nos: np.ndarray, mat=None):
    """Tensões térmicas na seção (corpo livre, sem cargas), estado plano
    generalizado: deformação axial linear em x e y escolhida para anular força
    normal e momentos. Malha completa (sem simetria) recomendada.
    Devolve dict com sx, sy, txy, sz, s1 (máx. principal) por elemento, Pa."""
    mat = np.where(m2.marca == 0, AR, CERAMICA) if mat is None else mat
    ativos = mat != AR
    tri = m2.tri[ativos]
    xy = m2.xy
    nn = m2.n_nos
    Te = T_nos[tri].mean(1)
    E, nu = _propriedades(np.asarray(mat)[ativos], Te, M)
    eth = _eps_termica_nos(np.asarray(mat)[ativos], T_nos[tri], M).mean(1)
    A, dx, dy = el.geometria_t3(xy, tri)
    A = np.abs(A)
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
    D = np.zeros((len(tri), 3, 3))
    D[:, 0, 0] = D[:, 1, 1] = lam + 2 * mu; D[:, 0, 1] = D[:, 1, 0] = lam; D[:, 2, 2] = mu
    B = np.zeros((len(tri), 3, 6))
    B[:, 0, 0::2] = dx; B[:, 1, 1::2] = dy; B[:, 2, 0::2] = dy; B[:, 2, 1::2] = dx
    Ke = np.einsum("eji,ejk,ekl->eil", B, D, B) * A[:, None, None]
    # tensão inicial em deformação plana (ez = 0): sigma0 = (3 lam + 2 mu) e nas normais
    eps0 = np.c_[eth, eth, np.zeros(len(tri))]
    s0 = np.einsum("eij,ej->ei", D, eps0) + np.c_[lam * eth, lam * eth, np.zeros(len(tri))]
    fe = np.einsum("eji,ej->ei", B, s0) * A[:, None]
    dofs = np.zeros((len(tri), 6), int); dofs[:, 0::2] = 2 * tri; dofs[:, 1::2] = 2 * tri + 1
    K = sp.coo_matrix((Ke.ravel(), (np.repeat(dofs, 6, 1).ravel(), np.tile(dofs, 6).ravel())),
                      shape=(2 * nn, 2 * nn)).tocsr()
    F = np.bincount(dofs.ravel(), weights=fe.ravel(), minlength=2 * nn)
    usados = np.unique(tri)
    fixos = np.ones(2 * nn, bool); fixos[(2 * usados[:, None] + np.arange(2)).ravel()] = False
    c = xy[usados]
    n0 = usados[np.argmin(np.hypot(c[:, 0], c[:, 1] - m2.largura / 2))]
    n1 = usados[np.argmin(np.hypot(c[:, 0] - m2.comprimento, c[:, 1] - m2.largura / 2))]
    fixos[[2 * n0, 2 * n0 + 1, 2 * n1 + 1]] = True
    livres = np.nonzero(~fixos)[0]
    u = np.zeros(2 * nn)
    u[livres] = spla.spsolve(K[livres][:, livres].tocsc(), F[livres])
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
    meio = (sx + sy) / 2; raio = np.sqrt(((sx - sy) / 2) ** 2 + txy ** 2)
    s1 = np.maximum(meio + raio, sz)
    out = {k: np.full(len(m2.tri), np.nan) for k in ("sx", "sy", "txy", "sz", "s1")}
    for k, v in zip(("sx", "sy", "txy", "sz", "s1"), (sx, sy, txy, sz, s1)):
        out[k][ativos] = v
    out["area"] = np.zeros(len(m2.tri)); out["area"][ativos] = A
    out["ativos"] = ativos
    return out
