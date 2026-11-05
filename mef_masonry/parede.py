"""Parede de alvenaria seca sob carga vertical e incêndio numa face.

Modelo de faixa vertical (1 m de comprimento de parede), adequado a paredes
em que as bordas verticais não restringem a flexão fora do plano:

* elementos de viga 2D corrotacionais (Crisfield 1991; Battini 2002) com
  seção de fibras: grandes deslocamentos e rotações, portanto segunda ordem
  (P-delta), flambagem e arqueamento térmico;
* fibras = triângulos de cerâmica da malha 2D da seção (um bloco de
  comprimento, com os furos), cada um com a sua temperatura, tirada da
  análise térmica 2D (ISO 834 numa face);
* material homogeneizado: curva tensão (área líquida) × deformação do prisma,
  calculada pelo modelo não linear 3D a 30, 600, 800 e 1000 °C e interpolada
  na temperatura; sem tração (as juntas secas abrem); descarga secante; ramo
  pós-pico regularizado pelo comprimento do elemento;
* deformação térmica livre da cerâmica (alfa ou tabela).

Convenções: x vertical (altura), y perpendicular à parede, face exposta em
y = -t/2. Nas fibras, deformação e tensão positivas em tração; a lei do
material trabalha com compressão positiva (e = -eps, s = -sigma).
Forças por metro de parede (N/m) e momentos em N·m/m.
"""
from dataclasses import dataclass
from typing import Callable, Dict, Optional, Sequence, Tuple

import numpy as np

from .malha import FURO, Malha2D


# ----------------------------------------------------------------------------
# Seção de fibras
# ----------------------------------------------------------------------------
@dataclass
class SecaoFibras:
    y: np.ndarray          # (nf,) distância ao plano médio, m (negativa do lado exposto)
    A: np.ndarray          # (nf,) área por metro de parede, m²/m
    tri: np.ndarray        # (nf,) índice do triângulo na malha 2D
    espessura: float       # m
    comprimento: float     # comprimento da seção modelada, m
    graute: Optional[np.ndarray] = None     # (nf,) True nas fibras de graute
    argamassa: Optional[np.ndarray] = None  # (nf,) True nas fibras de junta (seção vertical)

    def __post_init__(self):
        if self.graute is None:
            self.graute = np.zeros(len(self.y), bool)
        if self.argamassa is None:
            self.argamassa = np.zeros(len(self.y), bool)

    @classmethod
    def da_malha(cls, m2: Malha2D, geo=None) -> "SecaoFibras":
        """Fibras: triângulos do bloco e, se geo for dado, do graute (fiada 0).

        Na seção vertical de blocos de furos horizontais, só as paredes que
        descem contínuas pela fiada carregam a compressão vertical: as fibras
        são os triângulos fora das colunas de furos (paredes externas, septos
        verticais e as juntas sobre eles), com área por metro de parede igual à
        área do triângulo dividida pela altura da fiada. Os septos horizontais
        entram só na térmica."""
        if m2.vertical:
            cen = m2.xy[m2.tri].mean(1)
            col = np.zeros(len(m2.tri), bool)
            for y0, y1 in (geo.furos_h_y() if geo is not None else []):
                col |= (cen[:, 1] > y0) & (cen[:, 1] < y1)
            sol = np.nonzero(~col & (m2.marca != FURO))[0]
            gr = np.zeros(len(m2.tri), bool)
            arg = m2.marca > 0
        elif geo is not None:
            from .malha import estados_fiada
            furo, gr = estados_fiada(m2, geo, 0)
            sol = np.nonzero(~furo | gr)[0]
            arg = np.zeros(len(m2.tri), bool)
        else:
            furo, gr = m2.marca == FURO, np.zeros(len(m2.tri), bool)
            sol = np.nonzero(~furo)[0]
            arg = np.zeros(len(m2.tri), bool)
        cen = m2.xy[m2.tri[sol]].mean(1)
        a = np.abs(m2.areas()[sol])
        t = m2.largura
        return cls(y=cen[:, 1] - t / 2, A=a / m2.comprimento, tri=sol, espessura=t, comprimento=m2.comprimento,
                   graute=gr[sol], argamassa=arg[sol])

    @classmethod
    def retangulo(cls, t: float, n: int = 200) -> "SecaoFibras":
        """Seção maciça de espessura t (verificações): n faixas iguais."""
        y = (np.arange(n) + 0.5) * t / n - t / 2
        return cls(y=y, A=np.full(n, t / n), tri=np.arange(n), espessura=t, comprimento=1.0)

    def temperaturas(self, m2: Malha2D, T_nos: np.ndarray) -> np.ndarray:
        return np.asarray(T_nos)[m2.tri[self.tri]].mean(1)

    @property
    def area(self):
        return float(self.A.sum())

    @property
    def inercia(self):
        return float((self.A * self.y ** 2).sum())


# ----------------------------------------------------------------------------
# Leis do material (compressão positiva)
# ----------------------------------------------------------------------------
def _interp_tangente(x, xp, fp):
    """Interpolação linear por partes com a derivada; fora da tabela: extremos e derivada nula."""
    j = np.clip(np.searchsorted(xp, x, side="right"), 1, len(xp) - 1)
    x0, x1, f0, f1 = xp[j - 1], xp[j], fp[j - 1], fp[j]
    d = (f1 - f0) / (x1 - x0)
    f = f0 + (x - x0) * d
    fora = x >= xp[-1]
    f = np.where(fora, fp[-1], f)
    d = np.where(fora, 0.0, d)
    return f, d


class LeiLinear:
    """Elástica linear com tração (verificações)."""

    def __init__(self, E: float):
        self.E = E
        self.E0 = E

    def tensao(self, e, T, e_max):
        return self.E * e, np.full(np.shape(e), self.E)

    def regularizada(self, L_e):
        return self


class LeiParede:
    """Lei homogeneizada da alvenaria seca: s(e, T) por curvas tabeladas.

    T      : temperaturas das curvas (°C), crescentes
    e, s   : para cada temperatura, deformação de compressão e tensão na área
             líquida (Pa), começando em (0, 0) e terminando em s = 0
    e_pico : deformação no pico de cada curva
    """

    def __init__(self, T, e, s, H_ref=0.3):
        self.T = np.asarray(T, float)
        self.e = [np.asarray(x, float) for x in e]
        self.s = [np.asarray(x, float) for x in s]
        self.H_ref = H_ref
        ip = [int(np.argmax(si)) for si in self.s]
        self.e_pico = np.array([ei[i] for ei, i in zip(self.e, ip)])
        self.s_pico = np.array([si[i] for si, i in zip(self.s, ip)])
        self.E0 = float(self.s_pico[0] / self.e_pico[0])          # secante no pico a 30 °C
        self._I = np.eye(len(self.T))
        self.d_min = 1e-6 * self.E0

    @classmethod
    def de_prisma(cls, curvas: Dict[float, Sequence[Tuple[float, float]]], H_ref: float, A_liq: float):
        """curvas {T: [(delta_mm, carga_kN), ...]} do prisma de altura H_ref e área líquida A_liq."""
        Ts, es, ss = [], [], []
        for T in sorted(curvas):
            a = np.asarray(curvas[T], float)
            a = a[np.argsort(a[:, 0])]
            e = np.r_[0.0, a[:, 0] * 1e-3 / H_ref]
            s = np.r_[0.0, a[:, 1] * 1e3 / A_liq]
            e, k = np.unique(e, return_index=True)
            s = s[k]
            # ramo pós-pico até tensão nula, pela inclinação final (ou, se a curva parou no pico,
            # pela inclinação média de descida 30% do módulo secante no pico)
            ip = int(np.argmax(s))
            if ip < len(s) - 1 and s[-1] < s[-2]:
                inc = (s[-1] - s[-2]) / (e[-1] - e[-2])
            else:
                inc = -0.3 * s[ip] / e[ip]
            e_zero = e[-1] + s[-1] / -inc
            e, s = np.r_[e, e_zero], np.r_[s, 0.0]
            Ts.append(T); es.append(e); ss.append(s)
        return cls(Ts, es, ss, H_ref)

    def escalada(self, k_s: float, k_e: float) -> "LeiParede":
        return LeiParede(self.T, [e * k_e for e in self.e], [s * k_s for s in self.s], self.H_ref)

    def com_curva(self, T: float, e, s) -> "LeiParede":
        """Troca (ou acrescenta) a curva de uma temperatura."""
        Ts, es, ss = list(self.T), list(self.e), list(self.s)
        if T in Ts:
            i = Ts.index(T)
            es[i], ss[i] = np.asarray(e, float), np.asarray(s, float)
        else:
            Ts.append(T); es.append(np.asarray(e, float)); ss.append(np.asarray(s, float))
        o = np.argsort(Ts)
        return LeiParede([Ts[i] for i in o], [es[i] for i in o], [ss[i] for i in o], self.H_ref)

    def regularizada(self, L_e: float) -> "LeiParede":
        """Ramo pós-pico com o deslocamento pós-pico do prisma concentrado num
        elemento de comprimento L_e: e = e_p + (e_prisma - e_p) H_ref / L_e."""
        k = self.H_ref / L_e
        es = []
        for e, ep in zip(self.e, self.e_pico):
            es.append(np.where(e > ep, ep + (e - ep) * k, e))
        return LeiParede(self.T, es, self.s, self.H_ref)

    def pesos(self, T):
        T = np.asarray(T, float)
        return np.stack([np.interp(T, self.T, self._I[i]) for i in range(len(self.T))])

    def envoltoria(self, e, T, w=None):
        w = self.pesos(T) if w is None else w
        s = np.zeros(np.broadcast(e, w[0]).shape)
        d = np.zeros_like(s)
        for i in range(len(self.T)):
            if not np.any(w[i]):
                continue
            si, di = _interp_tangente(e, self.e[i], self.s[i])
            s += w[i] * si
            d += w[i] * di
        return s, d

    def tensao(self, e, T, e_max):
        """Tensão de compressão s(e) e tangente ds/de, com descarga secante a
        partir da maior compressão já atingida (e_max)."""
        w = self.pesos(T)
        s, d = self.envoltoria(e, T, w)
        desc = (e < e_max) & (e_max > 0)
        if np.any(desc):
            k = np.nonzero(desc)
            wk = np.broadcast_to(w, (len(self.T),) + np.shape(e))
            em = e_max[k]
            sm, _ = self.envoltoria(em, None, wk[(slice(None),) + k])
            sec = sm / em
            s[k] = sec * np.maximum(e[k], 0.0)
            d[k] = sec
        tr = e < 0
        s = np.where(tr, 0.0, s)
        d = np.where(tr | (d == 0.0), self.d_min, d)
        return s, d

    def resistencia(self, T):
        """Tensão de pico (área líquida) na temperatura T."""
        return float(self.pesos(np.asarray([T]))[:, 0] @ self.s_pico)


# ----------------------------------------------------------------------------
# Viga corrotacional com fibras
# ----------------------------------------------------------------------------
_GL = {2: ([0.5 - 0.5 / np.sqrt(3), 0.5 + 0.5 / np.sqrt(3)], [0.5, 0.5]),
       3: ([0.5 - 0.5 * np.sqrt(0.6), 0.5, 0.5 + 0.5 * np.sqrt(0.6)], [5 / 18, 8 / 18, 5 / 18])}


class NaoConvergiu(Exception):
    pass


class ModeloParede:
    """Faixa vertical de parede com n_el elementos.

    apoio_base, apoio_topo : "articulado" ou "engastado" (o topo é livre na vertical)
    imperfeicao            : flecha inicial senoidal no meio da altura (m), para o lado exposto
    eps_termica            : função T (°C) -> deformação térmica livre
    fluencia               : função T (°C) -> c(T) (1/Pa) da fluência transitória por fibra:
                             no primeiro aquecimento, Δε_tr = Δc · σ (compressão); None desliga
    """

    def __init__(self, H: float, secao: SecaoFibras, lei, eps_termica: Callable, n_el: int = 24,
                 apoio_base: str = "articulado", apoio_topo: str = "articulado", imperfeicao: float = 0.0,
                 n_gauss: int = 3, T_inicial: float = 20.0, regularizar: bool = True,
                 fluencia: Optional[Callable] = None):
        self.H, self.sec, self.n_el = H, secao, int(n_el)
        n = self.n_el
        self.x0 = np.linspace(0.0, H, n + 1)
        self.y0 = -imperfeicao * np.sin(np.pi * self.x0 / H)
        dx, dy = np.diff(self.x0), np.diff(self.y0)
        self.L0 = np.hypot(dx, dy)
        self.beta0 = np.arctan2(dy, dx)
        self.lei = lei.regularizada(float(self.L0.mean())) if regularizar else lei
        xi, w = _GL[n_gauss]
        self.xi, self.wg = np.asarray(xi), np.asarray(w)
        self.B1, self.B2 = -4 + 6 * self.xi, -2 + 6 * self.xi        # x L0: curvatura = (B1 th1 + B2 th2) / L0
        self.nd = 3 * (n + 1)
        fix = [0, 1] + ([2] if apoio_base == "engastado" else []) + \
              [3 * n + 1] + ([3 * n + 2] if apoio_topo == "engastado" else [])
        self.fixos = np.array(sorted(fix))
        self.livres = np.setdiff1d(np.arange(self.nd), self.fixos)
        self.apoios = (apoio_base, apoio_topo)
        self.eps_termica = eps_termica
        nf = len(secao.y)
        self.emax = np.zeros((n, len(self.xi), nf))
        self.d = np.zeros(self.nd)
        self.T_inicial = T_inicial
        self.T_f = np.full(nf, float(T_inicial))
        self.dofs = (3 * np.arange(n))[:, None] + np.arange(6)[None, :]
        self.ultimo = None
        # maior incremento de deformação térmica por passo (evita que o estado de tentativa
        # salte o pico da curva e o Newton caia num ramo esmagado espúrio)
        self.deps_max = 0.25 * float(np.min(self.lei.e_pico)) if hasattr(self.lei, "e_pico") else np.inf
        self.N_ref = max(float(getattr(self.lei, "s_pico", [self.lei.E0 * 1e-3])[0]) * secao.area, 1.0)
        # fluência transitória: deformação de compressão (positiva) por ponto de Gauss e fibra
        self.fluencia = fluencia
        self.e_tr = np.zeros((n, len(self.xi), nf))
        self.c_hist = fluencia(self.T_f) if fluencia else np.zeros(nf)
        self.tol_tr = 1e-4

    # ------------------------------------------------------------------ cinemática e forças
    def _estado(self, d, T_f, eth, emax, e_tr=None):
        n, sec = self.n_el, self.sec
        u = d.reshape(-1, 3)
        x = self.x0 + u[:, 0]
        y = self.y0 + u[:, 1]
        dx, dy = np.diff(x), np.diff(y)
        Ln = np.hypot(dx, dy)
        beta = np.arctan2(dy, dx)
        c, s = np.cos(beta), np.sin(beta)
        cb0, sb0 = np.cos(self.beta0), np.sin(self.beta0)
        alfa = np.arctan2(cb0 * s - sb0 * c, cb0 * c + sb0 * s)
        L0 = self.L0
        ul = (Ln ** 2 - L0 ** 2) / (Ln + L0)
        th1 = u[:-1, 2] - alfa
        th2 = u[1:, 2] - alfa
        # fibras
        eps0 = ul / L0
        kap = (self.B1[None, :] * th1[:, None] + self.B2[None, :] * th2[:, None]) / L0[:, None]   # (n, ng)
        eps = eps0[:, None, None] - sec.y[None, None, :] * kap[:, :, None]                       # (n, ng, nf)
        e = -(eps - eth[None, None, :])                                                           # compressão +
        e = e - (self.e_tr if e_tr is None else e_tr)                       # menos a fluência transitória
        sc, Et = self.lei.tensao(e, T_f[None, None, :], emax)
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

    def forca_externa(self, P, e_topo, e_base):
        """Carga P (N/m, compressão) aplicada no topo com excentricidade e_topo
        e reação da base com e_base (m, positivas para o lado exposto)."""
        f = np.zeros(self.nd)
        n = self.n_el
        f[3 * n] = -P
        if self.apoios[1] != "engastado":
            f[3 * n + 2] += -P * e_topo
        if self.apoios[0] != "engastado":
            f[2] += P * e_base
        return f

    # ------------------------------------------------------------------ equilíbrio
    def equilibrar(self, fext, T_f, d0=None, max_it=40, tol=1e-8, e_tr=None):
        """Newton-Raphson com busca linear. Devolve (d, info, iterações) ou
        levanta NaoConvergiu."""
        d = (self.d if d0 is None else d0).copy()
        eth = self.eps_termica(T_f)
        e_tr = self.e_tr if e_tr is None else e_tr
        esc = np.ones(self.nd)
        esc[2::3] = 1.0 / self.sec.espessura           # momentos em força equivalente
        lv = self.livres
        ref = self.N_ref
        f, K, info = self._estado(d, T_f, eth, self.emax, e_tr)
        r = (fext - f)[lv]
        nr = np.linalg.norm(r * esc[lv])
        du = None
        melhor, desde = nr, 0
        for it in range(max_it):
            if not np.isfinite(nr):
                raise NaoConvergiu("resíduo não finito")
            if nr < 0.5 * melhor:
                melhor, desde = nr, 0
            else:
                desde += 1
                if desde >= 10:
                    raise NaoConvergiu(f"resíduo estagnado em {nr / ref:.1e}")
            if du is not None and (nr < tol * ref or (nr < 1e-5 * ref and np.abs(du).max() < 1e-10)):
                return d, info, it
            try:
                du = np.linalg.solve(K[np.ix_(lv, lv)], r)
            except np.linalg.LinAlgError:
                raise NaoConvergiu("matriz singular")
            if not np.all(np.isfinite(du)) or np.abs(du[0::3]).max() > 0.5 * self.H:
                raise NaoConvergiu("deslocamentos sem limite")
            passo = 1.0
            for _ in range(6):                        # busca linear simples (redução do resíduo)
                dn = d.copy(); dn[lv] += passo * du
                fn, Kn, infon = self._estado(dn, T_f, eth, self.emax, e_tr)
                rn = (fext - fn)[lv]
                nrn = np.linalg.norm(rn * esc[lv])
                if np.isfinite(nrn) and (nrn < nr or nrn < tol * ref):
                    break
                passo *= 0.5
            d, f, K, info, r, nr = dn, fn, Kn, infon, rn, nrn
        raise NaoConvergiu(f"sem convergência em {max_it} iterações (resíduo {nr / ref:.1e})")

    def aceitar(self, d, T_f, info):
        self.d = d
        self.T_f = np.asarray(T_f, float).copy()
        self.emax = np.maximum(self.emax, info["e"])
        self.ultimo = info

    # ------------------------------------------------------------------ medidas
    @property
    def flecha(self):
        """Deslocamento horizontal dos nós (m); negativo = para o lado exposto."""
        return self.d[1::3]

    @property
    def flecha_meio(self):
        return float(np.interp(self.H / 2, self.x0, self.flecha))

    @property
    def flecha_max(self):
        v = self.flecha
        return float(v[np.argmax(np.abs(v))])

    @property
    def deslocamento_topo(self):
        """Deslocamento vertical do topo (m); positivo = para cima."""
        return float(self.d[3 * self.n_el])

    def compressao_critica(self):
        """(altura x em m, razão e/e_pico) do ponto de Gauss com a fibra mais
        comprimida em relação à deformação de pico na sua temperatura."""
        if self.ultimo is None:
            return None, 0.0
        ep = self.lei.pesos(self.T_f).T @ self.lei.e_pico
        r = (self.ultimo["e"] / ep[None, None, :]).max(-1)
        e, g = np.unravel_index(int(np.argmax(r)), r.shape)
        x = self.x0[e] + self.xi[g] * (self.x0[e + 1] - self.x0[e])
        return float(x), float(r[e, g])

    def fibras_meia_altura(self):
        """Tensão (Pa) nas fibras no ponto de Gauss mais próximo de H/2."""
        if self.ultimo is None:
            return None
        xg = (self.x0[:-1, None] + self.xi[None, :] * np.diff(self.x0)[:, None]).ravel()
        k = int(np.argmin(np.abs(xg - self.H / 2)))
        e, g = divmod(k, len(self.xi))
        return self.ultimo["sig"][e, g]

    # ------------------------------------------------------------------ análises
    def carregar(self, P, e_topo, e_base, n_passos=10, min_frac=1 / 512, registrar=None):
        """Aplica P (N/m) à temperatura atual em passos, com subdivisão. Devolve True se chegou."""
        fu = self.forca_externa(1.0, e_topo, e_base)
        Pa, dP = 0.0, P / n_passos
        while Pa < P * (1 - 1e-12):
            Pn = min(P, Pa + dP)
            try:
                d, info, _ = self.equilibrar(Pn * fu, self.T_f)
            except NaoConvergiu:
                dP /= 2
                if dP < P * min_frac:
                    return False
                continue
            self.aceitar(d, self.T_f, info)
            Pa = Pn
            if registrar:
                registrar(Pa)
        return True

    def capacidade(self, e_topo, e_base, P_passo=None, tol=1e-3, registrar=None, max_passos=400):
        """Carga máxima (N/m) por controle de carga com bisseção do passo."""
        fu = self.forca_externa(1.0, e_topo, e_base)
        Pref = self.N_ref
        dP = P_passo or Pref / 20
        Pa = 0.0
        for _ in range(max_passos):
            if dP < tol * Pref:
                break
            try:
                d, info, _ = self.equilibrar((Pa + dP) * fu, self.T_f)
            except NaoConvergiu:
                dP /= 2
                continue
            self.aceitar(d, self.T_f, info)
            Pa += dP
            if registrar:
                registrar(Pa)
        return Pa

    def _passo_fluencia(self, fext, Tn, pode_dividir):
        """Equilíbrio em Tn com a fluência transitória do passo (preditor pela
        tensão do início, corretor pela tensão média). Devolve None se o
        incremento passar de tol_tr e o passo puder ser dividido."""
        c_novo = self.fluencia(Tn)
        dc = np.maximum(c_novo - self.c_hist, 0.0)[None, None, :]
        s0 = np.maximum(-self.ultimo["sig"], 0.0) if self.ultimo is not None else 0.0
        d0 = dc * s0
        if pode_dividir and np.max(d0, initial=0.0) > self.tol_tr:
            return None
        d, info, _ = self.equilibrar(fext, Tn, e_tr=self.e_tr + d0)
        d1 = dc * 0.5 * (s0 + np.maximum(-info["sig"], 0.0))
        if pode_dividir and np.max(d1, initial=0.0) > 1.5 * self.tol_tr:
            return None
        if np.max(np.abs(d1 - d0), initial=0.0) > 0.05 * self.tol_tr:
            d, info, _ = self.equilibrar(fext, Tn, d0=d, e_tr=self.e_tr + d1)
        return d, info, d1, c_novo

    def incendio(self, P, e_topo, e_base, campos: Sequence[Tuple[float, np.ndarray]], min_sub=1 / 64,
                 registrar=None, parar_se=None):
        """P constante; campos [(t s, T_fibras)] crescentes no tempo, a partir do estado atual.
        Devolve (t_colapso ou None, último t convergido)."""
        fext = self.forca_externa(P, e_topo, e_base)
        t_a, T_a = campos[0]
        for t_b, T_b in campos[1:]:
            pilha = [(0.0, 1.0)]
            while pilha:
                f0, f1 = pilha.pop()
                Tn = T_a + (T_b - T_a) * f1
                if f1 - f0 > 1 / 1024 and \
                        np.abs(self.eps_termica(Tn) - self.eps_termica(self.T_f)).max() > self.deps_max:
                    fm = 0.5 * (f0 + f1)            # incremento de deformação térmica grande: subdivide antes
                    pilha += [(fm, f1), (f0, fm)]
                    continue
                try:
                    if self.fluencia is None:
                        d, info, _ = self.equilibrar(fext, Tn)
                    else:
                        res = self._passo_fluencia(fext, Tn, f1 - f0 > 1 / 1024)
                        if res is None:                 # incremento de fluência grande: subdivide
                            fm = 0.5 * (f0 + f1)
                            pilha += [(fm, f1), (f0, fm)]
                            continue
                        d, info, de_tr, c_novo = res
                except NaoConvergiu:
                    if f1 - f0 <= min_sub:
                        return t_a + (t_b - t_a) * f0, t_a + (t_b - t_a) * f0
                    fm = 0.5 * (f0 + f1)
                    pilha += [(fm, f1), (f0, fm)]
                    continue
                if self.fluencia is not None:
                    self.e_tr = self.e_tr + de_tr
                    self.c_hist = np.maximum(self.c_hist, c_novo)
                self.aceitar(d, Tn, info)
                tt = t_a + (t_b - t_a) * f1
                if registrar:
                    registrar(tt)
                if parar_se and parar_se(tt):
                    return None, tt
            t_a, T_a = t_b, T_b
        return None, t_a
