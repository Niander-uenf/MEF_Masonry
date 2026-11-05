"""Funções de alto nível: calibração das juntas, simulação do ensaio da UFV em 3D
e análises rápidas da seção 2D (forno e parede em incêndio-padrão)."""
import copy
import os
from dataclasses import replace
from typing import Optional

import numpy as np

from . import malha as ml, termico as tm, mecanico as mc, curvas
from .malha import AR, CERAMICA, JUNTA, CAPEAMENTO, ACO, FURO
from .parametros import Geometria, Ensaio
from .termico import Materiais, Contorno


# ----------------------------------------------------------------------------
# Rigidez do prisma e calibração das juntas secas
# ----------------------------------------------------------------------------
def modulo_prisma(m: ml.Malha3D, M: Materiais, carga_N=1e5, geo: Geometria = Geometria(),
                  terco_central=True):
    """Módulo secante do prisma em área bruta (MPa): pelo extensômetro no terço
    central da face externa e pelo encurtamento entre placas."""
    R = mc.resolver_3d(m, M, carga_N=carga_N, incluir_termico=False)["carga"]
    A_bruta = geo.comprimento * geo.largura
    sig = carga_N / A_bruta
    H = geo.altura_prisma
    eps, z1, z2 = mc.deformacao_extensometro(m, R, H / 3, 2 * H / 3)
    d = mc.encurtamento_entre_placas(m, R, H)
    return {"E_extensometro": sig / -eps / 1e6, "E_placas": sig / (-d / H) / 1e6,
            "z1": z1, "z2": z2, "resultado": R}


def estimativa_serie(M: Materiais, E_alvo_MPa, base_mm, n_juntas):
    """E da junta (material, Pa) pela associação em série na área bruta:
    base/E_alvo = (base - n t)/E_bloco + n t/(E_j r), r = área líquida/bruta."""
    t = M.junta.espessura
    r = M.ceramica.razao_area_liquida
    base = base_mm / 1000
    flex_j = base / (E_alvo_MPa * 1e6) - (base - n_juntas * t) / M.ceramica.E_bloco_bruto
    if flex_j <= 0:
        raise ValueError("Alvo mais rígido que o próprio bloco; reveja os dados.")
    return n_juntas * t / flex_j / r


def calibrar_junta(m: ml.Malha3D, M: Materiais, E_alvo_MPa=1100.0, geo: Geometria = Geometria(),
                   tol=0.002, max_iter=8, verboso=True):
    """Ajusta M.junta.E (e o capeamento, se igual) para o módulo no terço central
    igualar E_alvo_MPa. Usa secante em 1/E_junta. Devolve (E_junta, histórico)."""
    base = modulo_prisma(m, M, geo=geo)
    n_j = int(sum(1 for c in m.camadas if c.tipo == "junta"
                  and base["z1"] <= 0.5 * (c.z0 + c.z1) <= base["z2"]))
    Ej = estimativa_serie(M, E_alvo_MPa, (base["z2"] - base["z1"]) * 1000, max(n_j, 1))
    hist = []
    x_prev = f_prev = None
    for it in range(max_iter):
        M.junta.E = Ej; M.capeamento.E = Ej
        Ep = modulo_prisma(m, M, geo=geo)["E_extensometro"]
        hist.append((Ej / 1e6, Ep))
        if verboso:
            print(f"  iteração {it}: E_junta = {Ej / 1e6:8.2f} MPa -> E_prisma = {Ep:8.1f} MPa")
        f = 1 / Ep - 1 / E_alvo_MPa
        if abs(Ep - E_alvo_MPa) / E_alvo_MPa < tol:
            break
        x = 1 / Ej
        if x_prev is None:
            x_new = x * (1 + 0.3 * np.sign(f))
        else:
            x_new = x - f * (x - x_prev) / (f - f_prev)
        x_prev, f_prev = x, f
        Ej = 1 / x_new
    return Ej, hist


# ----------------------------------------------------------------------------
# Ensaio da UFV em 3D
# ----------------------------------------------------------------------------
def simular_ensaio_3d(m: ml.Malha3D, M: Materiais, ensaio: Ensaio, geo: Geometria = Geometria(),
                      dt=30.0, h_extremos=0.0, pasta_vtk: Optional[str] = None, verboso=False):
    """Aquece o prisma pelo protocolo do ensaio e aplica a carga máxima ao fim
    do patamar. Devolve um dicionário com históricos e indicadores."""
    f, t_fim, t_rampa = curvas.rampa_patamar(ensaio.T_inicial, ensaio.taxa_C_min,
                                             ensaio.T_alvo, ensaio.patamar_s)
    contornos = [Contorno(tm.areas_laterais_3d(m), f, ensaio.h_conveccao, ensaio.emissividade, "laterais")]
    if h_extremos > 0:
        contornos.append(Contorno(tm.areas_topo_base_3d(m), curvas.constante(ensaio.T_inicial),
                                  h_extremos, 0.0, "topo e base"))
    H = geo.altura_prisma
    zm = m.zs[np.argmin(np.abs(m.zs - H / 2))]
    L, Wm = m.m2d.comprimento, m.m2d.largura_modelo
    sondas = {
        "superfície (meia altura, x = 140 mm)": ml.no_mais_proximo(m.xyz, (0.140, 0.0, zm)),
        "centro do septo (x = 280 mm, meia altura)": ml.no_mais_proximo(m.xyz, (0.280, Wm, zm)),
    }
    cer_meio = np.unique(m.cun[(m.mat == CERAMICA) & (np.abs(m.xyz[m.cun[:, 0], 2] - H / 2) < H / 6)])
    sondas["média da cerâmica no terço central"] = cer_meio
    rt = tm.resolver_3d(m, M, contornos, t_fim, dt=dt, T_inicial=ensaio.T_inicial, sondas=sondas,
                        registrar_cada=60.0, instantes_campo=[t_rampa], verboso=verboso)
    T = rt.T_final
    res = mc.resolver_3d(m, M, T_nos=T, carga_N=ensaio.carga_N)
    Rt, Rc = res["termico"], res["carga"]
    Rtot = mc.combinar(Rt, Rc)
    A_bruta = geo.comprimento * geo.largura
    eps_c, _, _ = mc.deformacao_extensometro(m, Rc, H / 3, 2 * H / 3)
    ind = {
        "T_alvo": ensaio.T_alvo, "tempo_ate_carga_min": t_fim / 60,
        "T_superficie": rt.sondas[list(sondas)[0]][-1],
        "T_septo": rt.sondas[list(sondas)[1]][-1],
        "T_media_ceramica_terco": rt.sondas[list(sondas)[2]][-1],
        "T_min_ceramica": T[np.unique(m.cun[m.mat == CERAMICA])].min(),
        "dilatacao_entre_placas_mm": mc.encurtamento_entre_placas(m, Rt, H) * 1e3,
        "encurtamento_carga_mm": -mc.encurtamento_entre_placas(m, Rc, H) * 1e3,
        "E_terco_central_MPa": ensaio.carga_N / A_bruta / -eps_c / 1e6,
    }
    for nome, R in (("termico", Rt), ("carga", Rc), ("total", Rtot)):
        s, p99, pos = mc.tracao_maxima(m, R)
        sm, p99m, posm = mc.tracao_maxima(m, R, zmin=geo.altura_fiada, zmax=H - geo.altura_fiada)
        ind[f"tracao_max_{nome}_MPa"] = s / 1e6
        ind[f"tracao_p99_{nome}_MPa"] = p99 / 1e6
        ind[f"tracao_max_{nome}_fiadas_centrais_MPa"] = sm / 1e6
        ind[f"posicao_tracao_max_{nome}_mm"] = tuple(np.round(pos * 1000, 1))
    s1 = Rtot.princ[:, 2]
    vol = m.elem.volumes(m.xyz, m.cun)
    cer = m.mat == CERAMICA
    ind["fracao_ceramica_tracao_acima_1MPa"] = float(vol[cer & (s1 > 1e6)].sum() / vol[cer].sum())
    if pasta_vtk:
        from .saida import salvar_vtk
        os.makedirs(pasta_vtk, exist_ok=True)
        salvar_vtk(os.path.join(pasta_vtk, f"prisma_{int(ensaio.T_alvo)}C.vtk"), m,
                   dados_nos={"T_C": T, "u_total_mm": Rtot.u * 1e3, "u_termico_mm": Rt.u * 1e3},
                   dados_elem={"material": m.mat.astype(float),
                               "s1_total_MPa": Rtot.princ[:, 2] / 1e6,
                               "s1_termico_MPa": Rt.princ[:, 2] / 1e6,
                               "s3_total_MPa": Rtot.princ[:, 0] / 1e6,
                               "szz_total_MPa": Rtot.sig[:, 2] / 1e6})
    return {"termico": rt, "mecanico": res, "total": Rtot, "indicadores": ind, "sondas": list(sondas),
            "t_fim_rampa": t_rampa}


# ----------------------------------------------------------------------------
# Seção 2D: ensaio no forno e parede em incêndio-padrão
# ----------------------------------------------------------------------------
def _geo_furo(geo: Geometria, diametro):
    g = copy.copy(geo)
    g.diametro_furo = diametro
    return g


def secao_forno(M: Materiais, ensaio: Ensaio, geo: Geometria = Geometria(), h=0.005, dt=30.0,
                diametro=None, y_furos=None):
    """Seção do prisma no forno (todas as faces expostas) até o fim do patamar;
    tensões térmicas em estado plano generalizado."""
    g = _geo_furo(geo, geo.diametro_furo if diametro is None else diametro)
    furos = None if y_furos is None else [(x, y_furos) for x in g.centros_furos]
    m2 = ml.malha_secao(g, h=h, simetria=False, faixas=(), furos=furos)
    M2 = copy.deepcopy(M); M2.diametro_furo = g.diametro_furo
    f, t_fim, _ = curvas.rampa_patamar(ensaio.T_inicial, ensaio.taxa_C_min, ensaio.T_alvo, ensaio.patamar_s)
    ctr = [Contorno(tm.comprimentos_contorno_2d(m2), f, ensaio.h_conveccao, ensaio.emissividade)]
    r = tm.resolver_2d(m2, M2, ctr, t_fim, dt=dt, T_inicial=ensaio.T_inicial)
    s = mc.resolver_secao_gps(m2, M2, r.T_final)
    sol = m2.marca != FURO
    nos_sol = np.unique(m2.tri[sol])
    a = s["area"]
    return {"T_min": r.T_final[nos_sol].min(), "T_media": float((r.T_final[m2.tri].mean(1)[sol] * a[sol]).sum() / a[sol].sum()),
            "s1_max_MPa": np.nanmax(s["s1"]) / 1e6, "s1_p99_MPa": np.nanpercentile(s["s1"][sol], 99) / 1e6,
            "fracao_s1_acima_1MPa": float(a[sol & (s["s1"] > 1e6)].sum() / a[sol].sum()),
            "razao_area_liquida": float(a[sol].sum() / (g.comprimento * g.largura)),
            "malha": m2, "T": r.T_final, "tensoes": s}


def secao_parede_iso834(M: Materiais, geo: Geometria = Geometria(), h=0.005, dt=30.0, t_max_min=240.0,
                        diametro=None, y_furos=None, h_exposta=25.0, eps_exposta=0.8,
                        h_nao_exposta=4.0, eps_nao_exposta=0.8, T_amb=20.0, t_tensao_min=30.0):
    """Parede de BCEPQ (espessura 140 mm) com incêndio-padrão ISO 834 na face
    y = 0 e ambiente na face y = W; faces x = 0 e x = L adiabáticas
    (repetição ao longo da parede). Isolamento: elevação média de 140 °C ou
    máxima de 180 °C na face não exposta (NBR 5628 / ISO 834-1).
    Furos tratados como cavidades fechadas (sem convecção vertical)."""
    g = _geo_furo(geo, geo.diametro_furo if diametro is None else diametro)
    furos = None if y_furos is None else [(x, y_furos) for x in g.centros_furos]
    m2 = ml.malha_secao(g, h=h, simetria=False, faixas=(), furos=furos)
    M2 = copy.deepcopy(M); M2.diametro_furo = g.diametro_furo
    M2.ceramica.T_ref = T_amb
    exp = tm.comprimentos_contorno_2d(m2, faces=("y0",))
    nexp = tm.comprimentos_contorno_2d(m2, faces=("y1",))
    ctr = [Contorno(exp, curvas.iso834(T_amb), h_exposta, eps_exposta, "exposta"),
           Contorno(nexp, curvas.constante(T_amb), h_nao_exposta, eps_nao_exposta, "não exposta")]
    face = np.nonzero(nexp > 0)[0]
    w = nexp[face] / nexp[face].sum()

    def criterio(t, T):
        return (T[face] @ w - T_amb >= 140.0) or (T[face].max() - T_amb >= 180.0)
    r = tm.resolver_2d(m2, M2, ctr, t_max_min * 60, dt=dt, T_inicial=T_amb, parar_se=criterio,
                       sondas={"face não exposta (média)": face}, registrar_cada=60.0,
                       instantes_campo=[t_tensao_min * 60])
    out = {"t_isolamento_min": (r.parou_em / 60 if r.parou_em else np.inf),
           "razao_area_liquida": float(np.abs(m2.areas())[m2.marca != FURO].sum() / (g.comprimento * g.largura)),
           "malha": m2, "resultado": r}
    if t_tensao_min * 60 in r.campos:
        s = mc.resolver_secao_gps(m2, M2, r.campos[t_tensao_min * 60])
        out["s1_max_MPa"] = np.nanmax(s["s1"]) / 1e6
        out["s1_p99_MPa"] = np.nanpercentile(s["s1"][m2.marca != FURO], 99) / 1e6
    return out


# ----------------------------------------------------------------------------
# Análise não linear do prisma
# ----------------------------------------------------------------------------
def modulo_nao_linear(mod, geo: Geometria = Geometria(), delta=2e-6):
    """Módulo do prisma (MPa, área bruta, terço central da face externa) no
    regime elástico do modelo não linear, com as juntas fechadas."""
    from .naolinear import ModeloNaoLinear  # noqa
    m = mod.m
    estado = mod.estado()
    presc = 3 * mod.topo + 2
    P, _, _ = mod._passo(presc, mod.u[presc] - delta, np.zeros(mod.ndof))
    u = mod.u.reshape(-1, 3)
    H = geo.altura_prisma
    z1 = m.zs[np.argmin(np.abs(m.zs - H / 3))]; z2 = m.zs[np.argmin(np.abs(m.zs - 2 * H / 3))]
    face = np.nonzero(np.isclose(m.xyz[:, 1], 0.0))[0]
    uz = lambda z: u[face[np.isclose(m.xyz[face, 2], z)], 2].mean()
    eps = (uz(z2) - uz(z1)) / (z2 - z1)
    mod.restaurar(estado)
    return P / (geo.comprimento * geo.largura) / -eps / 1e6


def calibrar_junta_nao_linear(mod, E_alvo_MPa=1100.0, geo: Geometria = Geometria(), tol=0.002, verboso=True):
    """Ajusta o E das juntas da camada fina (modelo não linear) para o prisma
    reproduzir E_alvo no terço central. Secante em 1/E_junta."""
    def aplicar(Ej):
        mod.Ej[:] = Ej
        mod.En = mod.Ej.copy(); mod.Gs = mod.Ej / (2 * (1 + mod.nuj)); mod.fechada[:] = True
        mod._pc = None
    Ej = float(mod.Ej[0])
    x_prev = f_prev = None
    for it in range(8):
        aplicar(Ej)
        Ep = modulo_nao_linear(mod, geo)
        if verboso:
            print(f"  iteração {it}: E_junta = {Ej / 1e6:7.2f} MPa -> E_prisma = {Ep:7.1f} MPa")
        f = 1 / Ep - 1 / E_alvo_MPa
        if abs(Ep - E_alvo_MPa) / E_alvo_MPa < tol:
            break
        x = 1 / Ej
        x_new = x * (1 + 0.3 * np.sign(f)) if x_prev is None else x - f * (x - x_prev) / (f - f_prev)
        x_prev, f_prev = x, f
        Ej = 1 / x_new
    aplicar(Ej)
    return Ej


def campos_aquecimento(m, M: Materiais, ensaio: Ensaio, n_rampa=6, n_patamar=2, dt=30.0):
    """Roda a térmica 3D e devolve {tempo: campo} em instantes da rampa e do patamar."""
    f, t_fim, t_rampa = curvas.rampa_patamar(ensaio.T_inicial, ensaio.taxa_C_min, ensaio.T_alvo, ensaio.patamar_s)
    tempos = list(np.linspace(0, t_rampa, n_rampa + 1)[1:]) + list(np.linspace(t_rampa, t_fim, n_patamar + 1)[1:])
    tempos = [round(t / dt) * dt for t in tempos]
    ctr = [Contorno(tm.areas_laterais_3d(m), f, ensaio.h_conveccao, ensaio.emissividade)]
    r = tm.resolver_3d(m, M, ctr, t_fim, dt=dt, T_inicial=ensaio.T_inicial, instantes_campo=tempos)
    return r.campos


def modulo_corda(mod, P_ref_N, geo: Geometria = Geometria(), f1=0.05, f2=0.30):
    """Módulo de corda do prisma (MPa, área bruta) entre f1 e f2 de P_ref, com a
    deformação do extensômetro virtual (mod.definir_extensometro)."""
    h = [p for p in mod.historico if p.fase == "carga"]
    P = np.r_[0.0, [p.carga_N for p in h]]; e = np.r_[0.0, [-p.eps_ext for p in h]]
    if len(P) < 3 or P.max() < f2 * P_ref_N:
        return float("nan")
    i = int(np.argmax(P >= f2 * P_ref_N)) + 1
    e1, e2 = np.interp([f1 * P_ref_N, f2 * P_ref_N], P[:i], e[:i])
    return (f2 - f1) * P_ref_N / (geo.comprimento * geo.largura) / (e2 - e1) / 1e6


def ensaio_nao_linear(m, M: Materiais, Pd, Pj, campos=None, delta_max=0.9e-3, n_passos=30,
                      carga_max=None, geo: Geometria = Geometria(), verboso=True, tol_u=1e-4, tol_d=5e-3):
    """Um ensaio do prisma com o modelo não linear: aquecimento (campos
    {tempo: T_nos}, ou temperatura uniforme T_ref) sob a carga de assentamento
    e depois carga em deslocamento controlado até depois do pico. Devolve o
    modelo (com o histórico e o extensômetro no terço central) e (deltas, cargas)."""
    from .naolinear import ModeloNaoLinear
    mod = ModeloNaoLinear(m, M, Pd, Pj, verboso=verboso)
    mod.tol_u, mod.tol_d = tol_u, tol_d
    H = geo.altura_prisma
    mod.definir_extensometro(H / 3, 2 * H / 3)
    if campos is None:
        campos = {0.0: np.full(m.n_nos, M.ceramica.T_ref)}
    mod.aquecer(campos)
    deltas, cargas = mod.carregar(delta_max, n_passos, carga_max=carga_max)
    return mod, (deltas, cargas)


def calibrar_folga(m, M: Materiais, Pd, Pj, P_ref_N, E_alvo_MPa=1100.0, folgas=(0.024e-3, 0.032e-3),
                   sementes=(1, 2, 3), geo: Geometria = Geometria(), verboso=True):
    """Desvio-padrão das folgas das juntas de leito que reproduz, na média das
    sementes do campo aleatório, o módulo de corda 5-30% de P_ref (E_alvo).
    Cada caso carrega o prisma até 0,32 P_ref. A rigidez varia cerca de 20%
    entre sementes, por isso a média. Devolve (desvio, [(desvio, E_médio, [E por semente])])."""
    from dataclasses import replace
    pts = []
    for sg in folgas:
        Es = []
        for sem in sementes:
            mod, _ = ensaio_nao_linear(m, M, Pd, replace(Pj, desvio_folga=sg, semente=sem), delta_max=0.45e-3,
                                       n_passos=30, carga_max=0.32 * P_ref_N, geo=geo, verboso=False)
            Es.append(modulo_corda(mod, P_ref_N, geo))
        pts.append((sg, float(np.mean(Es)), Es))
        if verboso:
            print(f"    folga {sg * 1e3:.3f} mm -> E_corda = {np.mean(Es):7.1f} MPa (sementes: "
                  + ", ".join(f"{E:.0f}" for E in Es) + ")", flush=True)
    a = np.array([(p[0], p[1]) for p in pts])
    if len(a) == 1:
        return float(a[0, 0]), pts
    c = np.polyfit(a[:, 1], a[:, 0], 1)                 # desvio linear em E na faixa testada
    if not (a[:, 1].min() <= E_alvo_MPa <= a[:, 1].max()) and verboso:
        print("    aviso: E_alvo fora da faixa testada; valor extrapolado", flush=True)
    return float(np.polyval(c, E_alvo_MPa)), pts


def calibrar_fc(m, M: Materiais, Pd, Pj, P_alvo_N, geo: Geometria = Geometria(), tol=0.03, max_it=3, verboso=True):
    """Resistência à compressão efetiva da cerâmica (fc) que reproduz a carga
    máxima P_alvo a 30 °C. A carga máxima do modelo é quase proporcional a fc;
    a iteração usa essa proporção. Devolve (fc, pico, [(fc, pico)])."""
    from dataclasses import replace
    fc = Pd.fc
    pts = []
    for it in range(max_it):
        mod, (d, P) = ensaio_nao_linear(m, M, replace(Pd, fc=fc), Pj, geo=geo, verboso=False)
        pts.append((fc, P.max()))
        if verboso:
            print(f"    fc = {fc / 1e6:5.2f} MPa -> carga máxima {P.max() / 1e3:6.1f} kN", flush=True)
        if abs(P.max() / P_alvo_N - 1) < tol:
            break
        fc = fc * P_alvo_N / P.max()
    return pts[-1][0], pts[-1][1], pts
