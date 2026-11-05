"""Análise não linear do prisma: dano no bloco e na argamassa, juntas secas com atrito.

Cerâmica — dano escalar com duas variáveis, de tração e de compressão, e
fechamento de fissuras (Mazars, 1984; Mazars e Pijaudier-Cabot, 1989):
    sigma = (1 - d) D (eps - eps_t),  1 - d = (1 - r d_t)(1 - d_c)
    r = soma das tensões efetivas principais de tração / soma dos módulos
    (r = 1 em tração pura, 0 em compressão pura): uma fissura de tração que
    passa a ser comprimida fecha e volta a transmitir compressão;
    d_t cresce com o critério de Rankine, <sigma1_ef>/E, e d_c com o de von
    Mises das tensões principais de compressão, sqrt(3 J2(<sigma_ef>-))/(k E),
    k = fc/ft (tensões efetivas
    sigma_ef = D (eps - eps_t); ver Jirásek e Bažant, 2002); o dano começa
    com ft em tração uniaxial e fc em compressão uniaxial;
    amolecimento exponencial: 1 - d_i = (k0/kappa_i) exp(-(kappa_i - k0)/ku_i);
    regularização pela banda de fissuração (Bažant e Oh, 1983): a energia
    dissipada por volume é G/h, com h = largura do elemento na direção normal
    à fissura (Oliver, 1989), fixada quando cada dano começa; G = G_f em
    tração e G_c em compressão.

Elementos: W6, H8 ou W15 (malha.tipo). O estado (dano, abertura das juntas)
é avaliado no centro de cada elemento e vale para o elemento inteiro; a
rigidez é integrada com a regra completa do elemento.

Juntas secas e capeamento — camada fina (Desai et al., 1984), ortotrópica no
sistema da junta: rigidez normal E_junta em compressão e quase nula em tração
(junta aberta); cisalhamento elástico até o limite de Coulomb mu |sigma_n|,
sem coesão, tratado como rigidez secante.

Solução: incrementos de temperatura (aquecimento real, campo a campo) e
depois deslocamento imposto no prato de topo. Em cada incremento, iteração
secante; os sistemas lineares são resolvidos por gradientes conjugados
precondicionados pela fatoração de uma rigidez de referência, refeita quando
a convergência piora.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .malha import AR, CERAMICA, JUNTA, CAPEAMENTO, ACO, ARGAMASSA, GRAUTE, Malha3D
from .termico import Materiais
from .parametros import DanoCeramica, JuntaAtrito, dano_argamassa as _dano_argamassa_padrao
from .parametros import dano_graute as _dano_graute_padrao
from . import elementos as el
from .elementos3d import tipo as tipo_elemento
from .mecanico import _eps_termica_nos, _fatorar


def lei_dano(kappa, k0, ku):
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
class Passo:
    fase: str
    rotulo: float           # tempo (s) no aquecimento; deslocamento do prato (m) na carga
    carga_N: float
    dano_max: float
    vol_danificado: float   # fração do volume de cerâmica com max(d_t, d_c) > 0,5 (fissurada ou esmagada)
    juntas_abertas: float   # fração das juntas abertas
    iteracoes: int
    convergiu: bool
    contato: float = float("nan")    # fração da área das juntas de leito em contato
    eps_ext: float = float("nan")    # deformação no extensômetro virtual (se definido)


class ModeloNaoLinear:
    """Prisma 3D com dano na cerâmica e juntas secas não lineares.

    apoio = 'prato': base com u_z = 0 e o mínimo contra corpo rígido (pratos sem atrito)
    apoio = 'engaste': base com os três deslocamentos nulos (para ensaios de verificação)
    """

    def __init__(self, m: Malha3D, M: Materiais, dano: DanoCeramica = None, junta: JuntaAtrito = None,
                 apoio: str = "prato", fator_ft_elem: Optional[np.ndarray] = None, verboso: bool = True,
                 dano_argamassa: Optional[DanoCeramica] = None, dano_graute: Optional[DanoCeramica] = None):
        self.m, self.M = m, M
        self.Pd = dano or DanoCeramica()
        self.Pa = dano_argamassa or _dano_argamassa_padrao()
        self.Pg = dano_graute or _dano_graute_padrao()
        self.Pj = junta or JuntaAtrito()
        self.verboso = verboso
        ativo = (m.mat != AR) & (m.mat != ACO)
        self.id_glob = np.nonzero(ativo)[0]
        self.cun = m.cun[ativo]
        self.mat = m.mat[ativo]
        # elementos com dano: bloco, argamassa e graute (cada um com seus parâmetros)
        self.ic = np.nonzero(np.isin(self.mat, (CERAMICA, ARGAMASSA, GRAUTE)))[0]
        self.arg = self.mat[self.ic] == ARGAMASSA
        self.gra = self.mat[self.ic] == GRAUTE
        self.nu_c = self._por_material(M.ceramica.nu, M.argamassa.nu, M.graute.nu)
        self.du_t = self._por_material(self.Pd.du_t, self.Pa.du_t, self.Pg.du_t)
        self.du_c = self._por_material(self.Pd.du_c, self.Pa.du_c, self.Pg.du_c)
        self.ij = np.nonzero(np.isin(self.mat, (JUNTA, CAPEAMENTO)))[0]
        tipos = np.array([m.camadas[c].tipo for c in m.camada[ativo]])
        self.normal_x = tipos[self.ij] == "unidade"          # juntas verticais (normal x)
        Ec = self.Pj.E_contato
        self.Ej = np.where(self.mat[self.ij] == JUNTA, Ec or M.junta.E, Ec or M.capeamento.E).astype(float)
        self.nuj = np.where(self.mat[self.ij] == JUNTA, M.junta.nu, M.capeamento.nu)
        self.fator_ft = np.ones(len(self.ic)) if fator_ft_elem is None else fator_ft_elem[self.id_glob[self.ic]]
        n = m.n_nos
        self.n, self.ndof = n, 3 * n
        tp = m.elem
        self.tp = tp
        self.nd = nd = 3 * tp.nn
        w6 = m.tipo == "W6"
        self.dofs = (3 * self.cun[:, :, None] + np.arange(3)).reshape(-1, nd)
        # estrutura esparsa
        r = np.repeat(self.dofs, nd, axis=1).ravel().astype(np.int64)
        c = np.tile(self.dofs, (1, nd)).ravel().astype(np.int64)
        chave = r * self.ndof + c
        unicas, self.mapa = np.unique(chave, return_inverse=True)
        self.mapa = self.mapa.ravel()
        self.indices = (unicas % self.ndof).astype(np.int32)
        self.indptr = np.searchsorted(unicas // self.ndof, np.arange(self.ndof + 1)).astype(np.int64)
        self.nnz = len(unicas)
        del r, c, chave, unicas
        # geometria
        nu_c = self.nu_c
        cer = self.cun[self.ic]
        if w6:
            A, dz, dLdx, dLdy = el.geometria_w6(m.xyz, self.cun)
            self.vol = np.abs(A) * dz
            self.Bc, self.Nc = el.B_centroide_w6(m.xyz, self.cun)
            Gc, _ = el._gradientes_w6(A[self.ic], dz[self.ic], dLdx[self.ic], dLdy[self.ic], np.full(3, 1 / 3), 0.0)
            self.Ke1 = el.rigidez_w6(m.xyz, cer, np.ones(len(self.ic)), nu_c)[0]
        else:
            self.vol = tp.volumes(m.xyz, self.cun)
            self.Bc, self.Nc, _ = tp.no_centro(m.xyz, self.cun)
            # largura da banda (Oliver): gradientes lineares — a cunha W15 usa a W6 dos cantos
            lin = tipo_elemento("W6") if m.tipo == "W15" else tp
            Gc = lin.no_centro(m.xyz, cer[:, :lin.nn])[2]
            self.Ke1 = np.empty((len(self.ic), nd, nd))
            bl = max(500, int(20000 * 324 / nd ** 2))
            for i0 in range(0, len(self.ic), bl):
                s_ = slice(i0, i0 + bl)
                self.Ke1[s_] = tp.rigidez(m.xyz, cer[s_], np.ones(len(cer[s_])), nu_c[s_])[0]
        self.Gc_cer = Gc                                      # para a largura da banda (Oliver)
        # juntas: integrais dos produtos de gradientes, Gab[e, a, b] = soma_p w_p g_a g_b^T
        # (nj, 3, 3, nn, nn), para montar a rigidez com D diagonal; e, por ponto, w (nj,) e G (nj, nn, 3)
        cj = self.cun[self.ij]
        xi, wq = tp.pontos()
        G, detJ = tp.gradientes(m.xyz, cj, xi)              # (nj, np, nn, 3), (nj, np)
        wv = detJ * wq[None, :]
        self.Gab = np.einsum("ep,epia,epjb->eabij", wv, G, G)
        self._wGj = (wv, G, tp.N(xi))
        self.Bint = np.einsum("ep,epia->eai", wv, G)          # integral dos gradientes (nj, 3, nn)
        # contorno
        usados = np.unique(self.cun)
        self.usados = usados
        z = m.xyz[:, 2]
        self.zmin, self.zmax = z[usados].min(), z[usados].max()
        self.base = usados[np.isclose(z[usados], self.zmin)]
        self.topo = usados[np.isclose(z[usados], self.zmax)]
        fixos = np.zeros(self.ndof, bool)
        inat = np.setdiff1d(np.arange(n), usados)
        fixos[(3 * inat[:, None] + np.arange(3)).ravel()] = True
        if apoio == "engaste":
            fixos[(3 * self.base[:, None] + np.arange(3)).ravel()] = True
        else:
            fixos[3 * self.base + 2] = True
            ref = self.base[np.argmin(np.linalg.norm(m.xyz[self.base, :2] - [m.m2d.comprimento / 2, 0.0], axis=1))]
            fixos[3 * ref] = True
            if m.m2d.simetria:
                sim = usados[np.isclose(m.xyz[usados, 1], m.m2d.largura_modelo)]
                fixos[3 * sim + 1] = True
            else:
                fixos[3 * ref + 1] = True
                outro = self.base[np.argmax(np.abs(m.xyz[self.base, 0] - m.xyz[ref, 0]))]
                fixos[3 * outro + 1] = True
        self.fixos_base = fixos
        self.fator_sim = 2.0 if m.m2d.simetria else 1.0
        # área de influência dos nós do topo (carga de assentamento)
        zt = self.zmax
        self.area_topo = tp.areas_faces(m.xyz, self.cun, lambda P: np.isclose(P[:, :, 2], zt).all(1),
                                        consistente=True)
        # juntas de leito e folgas iniciais (irregularidade das faces)
        self.dzj = tp.altura(m.xyz, cj)
        self.area_j = self.vol[self.ij] / self.dzj
        self.leito = (self.mat[self.ij] == JUNTA) & ~self.normal_x
        self.folga = np.zeros(len(self.ij)); self.eg = np.zeros(len(self.ij))
        if self.Pj.desvio_folga > 0:
            self.definir_folgas(self.Pj.desvio_folga, self.Pj.comprimento_corr, self.Pj.contato_inicial,
                                self.Pj.semente, ajustar_estado=False)
        self._en_j = np.zeros(len(self.ij))
        self._ext = None
        # estado
        nc, nj = len(self.ic), len(self.ij)
        self.kappa = np.zeros((nc, 2))            # histórico [tração, compressão]
        self.h = np.full((nc, 2), np.nan)         # largura da banda [tração, compressão]
        self.r = np.ones(nc)                      # peso de tração de Mazars
        self.d = np.zeros(nc); self.dt = np.zeros(nc); self.dc = np.zeros(nc)
        self.fechada = self.eg <= 0
        self.En = np.where(self.fechada, self.Ej, self.Pj.fator_aberta * self.Ej)
        self.Gs = self.Ej / (2 * (1 + self.nuj))
        self.u = np.zeros(self.ndof)
        # fluência transitória (deformação térmica induzida pela carga), desligada até definir_fluencia
        self.fluencia = None
        self.eps_tr = np.zeros((nc, 6))           # Voigt, distorções de engenharia
        self.phi_hist = np.zeros(nc)              # maior phi(T) já atingido (só o primeiro aquecimento conta)
        self.ftr1 = np.zeros((nc, nd))
        self.tol_tr = 2e-4                        # maior incremento de deformação transitória por subpasso
        self._Gint = None
        self.definir_temperatura(np.full(n, M.ceramica.T_ref))
        self._pc = None; self._pc_livres = None
        self.cg_max, self.cg_res_ok = 12, 1e-3   # solução linear inexata (ver _resolver_linear)
        self.tol_u, self.tol_d = 2e-5, 1e-3        # convergência: variação relativa de u e máxima do dano
        self.historico: List[Passo] = []
        self.instantaneos = None   # lista -> guarda u, dano e T de cada incremento (pós-processador)

    # ------------------------------------------------------------------
    def definir_folgas(self, desvio, comprimento_corr=0.04, contato_inicial=0.20, semente=1,
                       ajustar_estado=True):
        """Folgas iniciais nas juntas de leito: a cota da face é um campo gaussiano
        (desvio-padrão desvio, correlação exp(-r²/l²) com l = comprimento_corr),
        independente em cada junta. A fração contato_inicial da área (os pontos
        mais altos) toca sem carga; no resto, a folga é a diferença de cota."""
        from scipy.ndimage import gaussian_filter
        from scipy.interpolate import RegularGridInterpolator
        rng = np.random.default_rng(semente)
        m = self.m
        xy = m.xyz[self.cun[self.ij]].mean(1)[:, :2]
        camada = m.camada[self.id_glob[self.ij]]
        passo = min(comprimento_corr / 4, 0.005)
        xs = np.arange(-2 * passo, m.m2d.comprimento + 3 * passo, passo)
        ys = np.arange(-2 * passo, m.m2d.largura + 3 * passo, passo)
        self.folga[:] = 0.0
        for c in np.unique(camada[self.leito]):
            s = self.leito & (camada == c)
            campo = gaussian_filter(rng.standard_normal((len(xs), len(ys))), comprimento_corr / passo / 2, mode="wrap")
            campo = (campo - campo.mean()) / campo.std()
            z = desvio * RegularGridInterpolator((xs, ys), campo)(xy[s])
            ordem = np.argsort(-z)
            acum = np.cumsum(self.area_j[s][ordem]) / self.area_j[s].sum()
            z_ref = z[ordem][min(np.searchsorted(acum, contato_inicial), len(z) - 1)]
            self.folga[s] = np.maximum(z_ref - z, 0.0)
        self.eg = self.folga / self.dzj
        if ajustar_estado:
            self.fechada = self.eg <= 0
            self.En = np.where(self.fechada, self.Ej, self.Pj.fator_aberta * self.Ej)

    def definir_extensometro(self, z1, z2, face_y=0.0):
        """Extensômetro virtual: deformação média entre as cotas z1 e z2 na face y = face_y."""
        m = self.m
        zs = m.zs
        z1, z2 = zs[np.argmin(np.abs(zs - z1))], zs[np.argmin(np.abs(zs - z2))]
        face = np.nonzero(np.isclose(m.xyz[:, 1], face_y))[0]
        self._ext = (face[np.isclose(m.xyz[face, 2], z1)], face[np.isclose(m.xyz[face, 2], z2)], z2 - z1)

    def deformacao_extensometro(self):
        if self._ext is None:
            return float("nan")
        n1, n2, L = self._ext
        u = self.u.reshape(-1, 3)
        return float((u[n2, 2].mean() - u[n1, 2].mean()) / L)

    def _por_material(self, bloco, argamassa, graute):
        """Valor por elemento com dano, escolhido pelo material (escalares ou vetores)."""
        return np.where(self.gra, graute, np.where(self.arg, argamassa, bloco))

    def definir_temperatura(self, T_nos):
        self.T = T_nos
        Tel = T_nos[self.cun]
        self.eps_n = _eps_termica_nos(self.mat, Tel, self.M)          # (ne, 6)
        self.eps_c = self.eps_n @ self.Nc                                # no centróide
        Tc = Tel[self.ic].mean(1)
        pm, M = self._por_material, self.M
        self.Ecer = pm(M.ceramica.E_T(Tc), M.argamassa.E_T(Tc), M.graute.E_T(Tc))
        self.ft = pm(self.Pd.ft_T(Tc), self.Pa.ft_T(Tc), self.Pg.ft_T(Tc)) * self.fator_ft
        self.fc = pm(self.Pd.fc_T(Tc), self.Pa.fc_T(Tc), self.Pg.fc_T(Tc))
        ncer = len(self.ic)
        if self.m.tipo == "W6":
            self.fth1 = el.forca_termica_w6(self.m.xyz, self.cun[self.ic], np.ones(ncer),
                                             self.nu_c, self.eps_n[self.ic])
        else:
            self.fth1 = self.tp.forca_termica(self.m.xyz, self.cun[self.ic], np.ones(ncer),
                                              self.nu_c, self.eps_n[self.ic])
        # juntas: vetores térmicos unitários dos termos normais de D (k = 0, 1, 2):
        # Tj[e, k] = soma_p w_p (eps_t . N_p) B_k,p, com B_k só na componente k
        wv, G, Nq = self._wGj
        s_p = wv * (self.eps_n[self.ij] @ Nq.T)                           # (nj, np)
        Tk = np.einsum("ep,epia->eai", s_p, G)                            # (nj, 3, nn)
        nj, nn = len(self.ij), self.tp.nn
        self.Tj = np.zeros((nj, 3, nn, 3))
        for k in range(3):
            self.Tj[:, k, :, k] = Tk[:, k]
        self.Tj = self.Tj.reshape(nj, 3, 3 * nn)

    def _D_junta(self, En, Gs):
        """Diagonal da matriz constitutiva das juntas no sistema global (nj, 6)."""
        ein = self.Pj.fator_plano * self.Ej
        x = self.normal_x
        Dd = np.empty((len(self.ij), 6))
        Dd[:, 0] = np.where(x, En, ein); Dd[:, 1] = ein; Dd[:, 2] = np.where(x, ein, En)
        Dd[:, 3] = np.where(x, Gs, ein / 2); Dd[:, 4] = np.where(x, ein / 2, Gs); Dd[:, 5] = Gs
        return Dd

    def _rigidez_juntas(self, Dd):
        """Rigidez das juntas (nj, nd, nd) com D diagonal (Voigt xx, yy, zz, xy, yz, xz)."""
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

    def _montar(self, d, En, Gs, fech=None):
        fech = self.fechada if fech is None else fech
        nd2 = self.nd ** 2
        vals = np.empty((len(self.cun), nd2))
        vals[self.ic] = ((self.Ecer * (1 - d))[:, None, None] * self.Ke1).reshape(-1, nd2)
        Dd = self._D_junta(En, Gs)
        vals[self.ij] = self._rigidez_juntas(Dd).reshape(-1, nd2)
        fj = np.einsum("ek,eki->ei", Dd[:, :3], self.Tj)
        # folga como deformação inicial -eg na direção normal das juntas de leito fechadas:
        # sigma_n = E_j (eps_n + eg)
        fg = np.where(fech, -self.eg, 0.0) * Dd[:, 2]
        if np.any(fg):
            bz = np.zeros((len(self.ij), self.tp.nn, 3)); bz[:, :, 2] = self.Bint[:, 2]
            fj += fg[:, None] * bz.reshape(len(self.ij), -1)
        dados = np.bincount(self.mapa, weights=vals.ravel(), minlength=self.nnz)
        K = sp.csr_matrix((dados, self.indices, self.indptr), shape=(self.ndof, self.ndof))
        F = np.zeros(self.ndof)
        fc = (self.fth1 + self.ftr1) * (self.Ecer * (1 - d))[:, None]
        F += np.bincount(self.dofs[self.ic].ravel(), weights=fc.ravel(), minlength=self.ndof)
        F += np.bincount(self.dofs[self.ij].ravel(), weights=fj.ravel(), minlength=self.ndof)
        return K, F

    def _deformacoes(self, u):
        e = np.einsum("eij,ej->ei", self.Bc, u[self.dofs])
        e[:, :3] -= self.eps_c[:, None]
        if self.fluencia is not None:
            e[self.ic] -= self.eps_tr
        return e                                                   # deformação mecânica

    # ------------------------------------------------------------------ fluência transitória
    def definir_fluencia(self, k_bloco=0.0, k_argamassa=0.0, k_graute=0.0, tabela_bloco=None):
        """Liga a fluência transitória (deformação térmica induzida pela carga).

        Para cada material, phi(T) = k ε_th(T) (Anderberg e Thelandersson, 1976)
        ou, no bloco, uma tabela [(T °C, phi)] medida no ensaio transiente. No
        primeiro aquecimento sob compressão, o incremento é

            Δε_tr = Δphi / fc20 · [(1 + ν) σ⁻ - ν tr(σ⁻) I]

        com σ⁻ a parte compressiva da tensão nominal, fc20 a resistência à
        compressão a 20 °C do material e ν o seu Poisson (no uniaxial,
        Δε_tr = Δphi σ / fc20 na direção da carga). Irreversível: não volta no
        resfriamento nem ao reaquecer abaixo da temperatura máxima já atingida."""
        self.fluencia = {"k": self._por_material(k_bloco, k_argamassa, k_graute),
                         "tabela": sorted(tabela_bloco) if tabela_bloco else None,
                         "fc20": self._por_material(self.Pd.fc, self.Pa.fc, self.Pg.fc)}
        self.phi_hist = self._phi()

    def _phi(self):
        """phi(T) por elemento com dano, no estado de temperatura atual."""
        fl = self.fluencia
        phi = fl["k"] * np.maximum(self.eps_c[self.ic], 0.0)
        if fl["tabela"]:
            Ts, vs = zip(*fl["tabela"])
            Tc = self.T[self.cun[self.ic]].mean(1)
            bloco = ~(self.arg | self.gra)
            phi = np.where(bloco, np.interp(Tc, Ts, vs), phi)
        return phi

    def _integral_gradientes(self):
        """∫ ∂N/∂x dV por elemento com dano (nc, 3, nn): força de uma tensão constante."""
        if self._Gint is None:
            xi, wq = self.tp.pontos()
            cer = self.cun[self.ic]
            out = np.empty((len(cer), 3, self.tp.nn))
            for i0 in range(0, len(cer), 4000):
                s_ = slice(i0, i0 + 4000)
                G, detJ = self.tp.gradientes(self.m.xyz, cer[s_], xi)
                out[s_] = np.einsum("ep,epia->eai", detJ * wq[None, :], G)
            self._Gint = out
        return self._Gint

    def _forca_transitoria(self):
        """Vetor de força unitário (E = 1) da deformação transitória: ∫ Bᵀ D(1, ν) ε_tr dV."""
        S = _tensor(np.einsum("eij,ej->ei", el.matriz_D(np.ones(len(self.ic)), self.nu_c), self.eps_tr))
        self.ftr1 = np.einsum("eai,eab->eib", self._integral_gradientes(), S).reshape(len(self.ic), -1)

    def tensao_nominal(self):
        """Tensão nominal (1 - d) D e nos elementos com dano (nc, 6)."""
        e = self._deformacoes(self.u)[self.ic]
        return (1 - self.d)[:, None] * np.einsum("eij,ej->ei", el.matriz_D(self.Ecer, self.nu_c), e)

    def _incremento_transitorio(self, sig, phi_novo):
        """Δε_tr (nc, 6) para a tensão sig e o novo phi (explícito no passo)."""
        dphi = np.maximum(phi_novo - self.phi_hist, 0.0)
        if not np.any(dphi):
            return np.zeros_like(self.eps_tr)
        w, V = np.linalg.eigh(_tensor(sig))
        wn = np.minimum(w, 0.0)
        Sn = np.einsum("eik,ek,ejk->eij", V, wn, V)                       # parte compressiva
        nu = self.nu_c[:, None, None]
        tr = wn.sum(1)[:, None, None]
        D = (dphi / self.fluencia["fc20"])[:, None, None] * ((1 + nu) * Sn - nu * tr * np.eye(3)[None])
        return np.column_stack([D[:, 0, 0], D[:, 1, 1], D[:, 2, 2], 2 * D[:, 0, 1], 2 * D[:, 1, 2], 2 * D[:, 0, 2]])

    def _atualizar_ceramica(self, e, kappa_hist):
        nu = self.nu_c
        E, ft, fc = self.Ecer, self.ft, self.fc
        k = fc / ft
        k0 = ft / E
        sig = np.einsum("eij,ej->ei", el.matriz_D(E, nu), e)          # tensão efetiva
        w, V = np.linalg.eigh(_tensor(sig))
        wn = np.minimum(w, 0.0)                                           # parte compressiva
        J2 = ((wn[:, 0] - wn[:, 1]) ** 2 + (wn[:, 1] - wn[:, 2]) ** 2 + (wn[:, 2] - wn[:, 0]) ** 2) / 6
        et = np.maximum(w[:, 2], 0.0) / E                                 # Rankine (tração)
        ec = np.sqrt(3 * J2) / (k * E)                                    # von Mises das tensões de compressão
        kappa = np.maximum(kappa_hist, np.column_stack([et, ec]))
        h = self.h.copy()
        for j, col in ((0, 2), (1, 0)):          # direção normal à fissura: sigma1 (tração), sigma3 (compressão)
            novo = np.isnan(h[:, j]) & (kappa[:, j] > k0)
            if novo.any():
                proj = np.abs(np.einsum("eik,ek->ei", self.Gc_cer[novo], V[novo, :, col])).sum(1)
                h[novo, j] = 2.0 / np.maximum(proj, 1e-12)
        hh = np.where(np.isnan(h), (self.vol[self.ic] ** (1 / 3))[:, None], h)
        ku_t = np.maximum(self.du_t / hh[:, 0] - k0 / 2, 0.05 * k0)        # G_f/(h ft) - k0/2, G_f = du_t ft
        ku_c = np.maximum(self.du_c / (hh[:, 1] * k) - k0 / 2, 0.05 * k0)  # G_c/(h k fc) - k0/2, G_c = du_c fc
        dt = lei_dano(kappa[:, 0], k0, ku_t)
        dc = lei_dano(kappa[:, 1], k0, ku_c)
        if self.Pd.fechamento_fissuras:
            soma = np.abs(w).sum(1)
            r = np.where(soma > 0, np.maximum(w, 0.0).sum(1) / np.maximum(soma, 1e-300), 1.0)
        else:
            r = np.ones(len(dt))
        d = np.minimum(1 - (1 - r * dt) * (1 - dc), self.Pd.d_max)
        return kappa, d, h, r, dt, dc

    def _atualizar_juntas(self, e, fechada, forcar=None):
        en = np.where(self.normal_x, e[self.ij, 0], e[self.ij, 2])
        g1 = np.where(self.normal_x, e[self.ij, 3], e[self.ij, 4])
        g2 = e[self.ij, 5]
        gam = np.hypot(g1, g2)
        tol = self.Pj.eps_tol
        eg = self.eg                                   # folga inicial / espessura
        fechada = np.where(en < -eg - tol, True, np.where(en > -eg + tol, False, fechada))   # histerese
        if forcar is not None:          # junta que oscila entre aberta e fechada fica aberta
            fechada = fechada & ~forcar
        # junta fechada: sigma_n = E_j (en + eg), a folga entra como deformação inicial (_montar)
        En = np.where(fechada, self.Ej, self.Pj.fator_aberta * self.Ej)
        self._en_j = en
        sn = np.where(fechada, self.Ej * (en + eg), En * en)
        G_el = self.Ej / (2 * (1 + self.nuj))
        lim = self.Pj.mu * np.maximum(-sn, 0.0)
        Gs = np.where(G_el * gam <= lim, G_el, np.maximum(lim / np.maximum(gam, 1e-30), self.Pj.fator_G_min * G_el))
        return En, Gs, fechada

    def _resolver_linear(self, K, rhs, livres, x0):
        """Gradientes conjugados precondicionados pela fatoração de uma rigidez
        anterior, com poucas iterações (solução inexata: a iteração secante
        externa continua a partir dela). Refatora quando o resíduo relativo
        fica acima de cg_res_ok."""
        Kff = K[livres][:, livres].tocsr()
        if self._pc is None or self._pc_livres is None or not np.array_equal(self._pc_livres, livres):
            self._pc = _fatorar(Kff); self._pc_livres = livres
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
                self._pc = _fatorar(Kff)
                x = self._pc(rhs)
        return x, it[0]

    def _passo(self, presc, valores, F_ext, max_it=60, tol_u=None, tol_d=None):
        tol_u = self.tol_u if tol_u is None else tol_u
        tol_d = self.tol_d if tol_d is None else tol_d
        fixos = self.fixos_base.copy()
        fixos[presc] = True
        livres = np.nonzero(~fixos)[0]
        d, En, Gs, fech = self.d.copy(), self.En.copy(), self.Gs.copy(), self.fechada.copy()
        u = self.u.copy(); u[presc] = valores
        tol_mud = max(2, int(1e-3 * len(self.ij)))
        kappa_it = self.kappa.copy()
        trocas = np.zeros(len(self.ij), int)
        conv = False
        for it in range(1, max_it + 1):
            K, Fth = self._montar(d, En, Gs, fech)
            F = Fth + F_ext
            ub = np.zeros(self.ndof); ub[presc] = valores
            rhs = F[livres] - (K @ ub)[livres]
            x, _ = self._resolver_linear(K, rhs, livres, u[livres])
            u_novo = ub; u_novo[livres] = x
            e = self._deformacoes(u_novo)
            # histórico monotônico também entre iterações: evita a oscilação carga/descarga
            kappa, d_novo, h, r_m, dt, dc = self._atualizar_ceramica(e[self.ic], kappa_it)
            kappa_it = kappa
            if it > 6:   # amortecimento do fechamento de fissuras (peso r)
                d_novo = np.where(d_novo < d, 0.5 * (d_novo + d), d_novo)
            En_n, Gs_n, fech_n = self._atualizar_juntas(e, fech)
            # junta que troca de estado 3 vezes no mesmo incremento fica aberta (evita ciclo; perto da troca, sigma_n ~ 0)
            trocas += fech_n != fech
            travar = trocas >= 3
            if travar.any():
                En_n, Gs_n, fech_n = self._atualizar_juntas(e, fech, forcar=travar)
            if it > 4:   # amortecimento do atrito para evitar oscilação
                Gs_n = 0.5 * (Gs_n + Gs)
            du = np.linalg.norm(u_novo - u) / max(np.linalg.norm(u_novo), 1e-30)
            dd = np.max(np.abs(d_novo - d)) if len(d) else 0.0
            mudou = np.count_nonzero(fech_n != fech)
            u, d, En, Gs, fech = u_novo, d_novo, En_n, Gs_n, fech_n
            if du < tol_u and dd < tol_d and mudou <= tol_mud:
                conv = True
                break
        self.convergencia = (du, dd, mudou)
        # consolida o estado
        self.u, self.d, self.kappa, self.h, self.En, self.Gs = u, d, kappa, h, En, Gs
        self.r, self.dt, self.dc, self.fechada = r_m, dt, dc, fech
        K, Fth = self._montar(d, En, Gs, fech)
        r = K @ u - Fth - F_ext
        carga = -r[3 * self.topo + 2].sum() * self.fator_sim
        return carga, it, conv

    def _registrar(self, fase, rotulo, carga, it, conv):
        vc = self.vol[self.ic]
        lt = self.leito
        cont = float(self.area_j[lt][self._en_j[lt] < -self.eg[lt]].sum() / self.area_j[lt].sum()) if lt.any() else float("nan")
        dm = np.maximum(self.dt, self.dc)
        p = Passo(fase, rotulo, carga, float(dm.max()) if len(dm) else 0.0,
                  float(vc[dm > 0.5].sum() / vc.sum()) if len(vc) else 0.0,
                  float(np.mean(~self.fechada)) if len(self.En) else 0.0, it, conv,
                  cont, self.deformacao_extensometro())
        self.historico.append(p)
        if self.instantaneos is not None:       # estados para o pós-processador
            f32 = lambda v: np.asarray(v, np.float32).copy()
            self.instantaneos.append({"fase": fase, "rotulo": float(rotulo), "carga_N": float(carga),
                                      "u": f32(self.u), "dt": f32(self.dt), "dc": f32(self.dc), "d": f32(self.d),
                                      "T": f32(self.T) if fase == "aquecimento" else None})
        if self.verboso:
            rot = f"t = {rotulo / 60:6.1f} min" if fase == "aquecimento" else f"delta = {rotulo * 1e3:6.3f} mm"
            print(f"  {fase:<11s} {rot}  carga {carga / 1e3:8.2f} kN  dano máx {p.dano_max:.3f}  "
                  f"vol(d>0,5) {100 * p.vol_danificado:5.2f}%  "
                  + (f"contato {100 * cont:5.1f}%  " if cont == cont else "sem juntas de leito  ")
                  + f"it {it}" +("" if conv else " (não convergiu: du %.1e, dd %.1e, juntas %d)" % self.convergencia),
                  flush=True)
        return p

    # ------------------------------------------------------------------
    def aquecer(self, campos: Dict[float, np.ndarray], pre_carga_N: Optional[float] = None,
                topo_fixo: bool = False):
        """Aplica em sequência os campos de temperatura {tempo: T_nos}, com o topo
        livre sob a carga de assentamento (ou, com topo_fixo, com o deslocamento
        vertical do topo travado no valor atual: aquecimento restringido).

        Com a fluência transitória ligada, cada intervalo entre campos é
        percorrido em subpassos (temperatura interpolada) com passo adaptativo,
        de modo que o incremento de deformação transitória não passe de tol_tr.
        Cada subpasso é um preditor (tensão do início) e um corretor (tensão
        média do subpasso, regra do trapézio)."""
        P0 = self.Pj.pre_carga_N if pre_carga_N is None else pre_carga_N
        F_ext = np.zeros(self.ndof)
        if topo_fixo:
            presc = 3 * self.topo + 2
            valores = self.u[presc].copy()
        else:
            presc, valores = np.zeros(0, int), np.zeros(0)
            F_ext[3 * np.arange(self.n) + 2] = -P0 / self.fator_sim * self.area_topo / self.area_topo.sum()
        for t in sorted(campos):
            T_b = np.asarray(campos[t], float)
            if self.fluencia is None:
                self.definir_temperatura(T_b)
                carga, it, conv = self._passo(presc, valores, F_ext)
            else:
                carga, it, conv = self._aquecer_com_fluencia(T_b, presc, valores, F_ext)
            self._registrar("aquecimento", t, carga if topo_fixo else P0, it, conv)
        self.d_aquecimento = np.maximum(self.dt, self.dc)

    def _aquecer_com_fluencia(self, T_b, presc, valores, F_ext, df_min=1 / 32):
        T_a = self.T.copy()
        tol = self.tol_tr
        # passo inicial pela tensão atual (sem tensão, começa com 1/8 do intervalo); sem aumento de
        # phi (temperatura parada ou resfriamento), um passo só, sem fluência
        sig0 = self.tensao_nominal()
        self.definir_temperatura(T_b)
        phi_b = self._phi()
        if not np.any(phi_b > self.phi_hist):
            return self._passo(presc, valores, F_ext)
        est = np.abs(self._incremento_transitorio(sig0, phi_b)).max(initial=0.0)
        self.definir_temperatura(T_a)
        df = 1.0 / max(1, int(np.ceil(est / tol))) if est > 0 else 1 / 8
        f0 = 0.0
        while f0 < 1.0 - 1e-12:
            f1 = min(1.0, f0 + df)
            salvo = self.estado()
            T_salvo = self.T.copy()
            eps0 = self.eps_tr.copy()
            sig0 = self.tensao_nominal()
            self.definir_temperatura(T_a + (T_b - T_a) * f1)
            phi = self._phi()
            d0 = self._incremento_transitorio(sig0, phi)              # preditor
            self.eps_tr = eps0 + d0
            self._forca_transitoria()
            carga, it, conv = self._passo(presc, valores, F_ext)
            d1 = self._incremento_transitorio(0.5 * (sig0 + self.tensao_nominal()), phi)   # corretor
            m1 = np.abs(d1).max(initial=0.0)
            if m1 > 1.5 * tol and df > df_min:
                self.restaurar(salvo)
                self.definir_temperatura(T_salvo)
                df /= 2
                continue
            if np.abs(d1 - d0).max(initial=0.0) > 0.05 * tol:
                self.eps_tr = eps0 + d1
                self._forca_transitoria()
                carga, it, conv = self._passo(presc, valores, F_ext)
            self.phi_hist = np.maximum(self.phi_hist, phi)
            f0 = f1
            if m1 < 0.5 * tol:
                df = min(1.0, 1.5 * df)
        return carga, it, conv

    @property
    def deformacao_transitoria_media(self):
        """Média (pelo volume) da deformação transitória vertical nos elementos com dano."""
        vc = self.vol[self.ic]
        return float((self.eps_tr[:, 2] * vc).sum() / vc.sum()) if len(vc) else 0.0

    def carregar(self, delta_max: float, n_passos: int = 30, parar_apos_pico: float = 0.80,
                 max_cortes: int = 3, carga_max: Optional[float] = None, max_sem_conv: int = 25):
        """Desloca o prato de topo (uniforme, sem atrito) até delta_max, até a
        carga cair a parar_apos_pico x o pico ou até passar de carga_max (N).
        Devolve (deltas, cargas) em m e N."""
        presc = 3 * self.topo + 2
        u0 = self.u[presc].copy()
        F0 = np.zeros(self.ndof)
        deltas, cargas = [0.0], [self._passo(presc, u0, F0)[0]]
        passo0 = delta_max / n_passos
        passo_min = passo0 / 2 ** max_cortes
        passo = passo0
        delta, pico, sem_conv = 0.0, cargas[0], 0
        while delta < delta_max - 1e-12:
            estado = self.estado()
            dn = delta + passo
            carga, it, conv = self._passo(presc, u0 - dn, F0)
            if not conv and passo / 2 >= passo_min * (1 - 1e-9):
                self.restaurar(estado)
                passo /= 2
                continue
            delta = dn
            sem_conv = 0 if conv else sem_conv + 1     # aceito sem convergir só com o passo mínimo
            if conv and it <= 20:
                passo = min(1.5 * passo, passo0)
            deltas.append(delta); cargas.append(carga)
            self._registrar("carga", delta, carga, it, conv)
            pico = max(pico, carga)
            if carga < parar_apos_pico * pico and pico > 0:
                break
            if carga_max is not None and carga >= carga_max:
                break
            if sem_conv >= max_sem_conv:
                if self.verboso:
                    print(f"  parado: {max_sem_conv} incrementos seguidos sem convergir com o passo mínimo", flush=True)
                break
        return np.array(deltas), np.array(cargas)

    _CAMPOS_ESTADO = ("u", "d", "dt", "dc", "r", "kappa", "h", "En", "Gs", "fechada", "eps_tr", "phi_hist", "ftr1")

    def estado(self):
        return {c: getattr(self, c).copy() for c in self._CAMPOS_ESTADO}

    def restaurar(self, estado):
        for c, v in estado.items():
            setattr(self, c, v)

    def dano_global(self, campo="d"):
        """Dano por elemento na malha completa (nan fora da cerâmica).
        campo: 'd' (efetivo, com fissuras fechadas), 'dt' (tração), 'dc' (compressão)
        ou 'max' (o maior entre dt e dc)."""
        v = np.maximum(self.dt, self.dc) if campo == "max" else getattr(self, campo)
        d = np.full(len(self.m.cun), np.nan)
        d[self.id_glob[self.ic]] = v
        return d
