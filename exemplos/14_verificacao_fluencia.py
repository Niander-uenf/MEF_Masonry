"""Verificação da fase 4: fluência transitória (deformação térmica induzida pela carga).

Modelo: no primeiro aquecimento sob compressão, Δε_tr = Δφ / fc20 · [(1 + ν) σ⁻ - ν tr(σ⁻) I],
com φ(T) = k_tr ε_th(T) (Anderberg e Thelandersson, 1976) ou tabelado.

V1  Prisma homogêneo sob tensão constante, aquecido uniformemente (W6, H8 e W15):
    ε_z = ε_th + σ/E + k σ/fc ε_th e ε_x = ε_th - ν σ/E - ν k σ/fc ε_th, exatos.
V2  Aquecimento com o topo travado (relaxação): σ = fc/k [exp(-E k ε_th / fc) - 1], a
    solução da equação dσ/dε_th = -E (1 + k σ/fc); o passo adaptativo controla o erro.
V3  Irreversibilidade: aquecer a 520 °C, resfriar a 270 °C e reaquecer a 520 °C não muda
    ε_tr; seguir até 620 °C soma só a parte de 520 a 620 °C.
V4  Sem compressão não há fluência: sob tração, ε_tr = 0.
V5  φ(T) tabelado substitui k ε_th: ε_tr = φ(T) σ / fc.
V6  Parede (fibras): coluna curta sob carga constante, aquecida por igual; encurtamento
    exato. O coeficiente da parede com argamassa é a média em série de bloco e junta.
V7  Parede BCEPQ no ISO 834 (k_tr = 0 por padrão): a fluência ligada e desligada dão o
    mesmo tempo; com k_tr hipotético, o tempo muda (registrado à parte, sem tolerância).
V8  Padrões e biblioteca: argamassa 2,0 (Prakash et al., 2020), bloco e graute 0; os
    projetos prontos passam na validação.
V9  (--completo) Ensaio transiente do prisma de tijolo maciço: a fluência da argamassa
    alivia as tensões térmicas entre tijolo e junta, e o dano no aquecimento cai.

Uso: python exemplos/14_verificacao_fluencia.py [--completo]
Grava resultados/verificacao_fluencia.csv.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, projeto as pj, biblioteca as bib, curvas
from mef_masonry.parametros import Geometria, DanoCeramica
from mef_masonry.termico import Materiais
from mef_masonry.naolinear import ModeloNaoLinear
from mef_masonry.parede import ModeloParede, SecaoFibras, LeiLinear
from mef_masonry.saida import salvar_csv

COMPLETO = "--completo" in sys.argv
linhas = []


def registrar(caso, modelo, ref, tol_rel, unidade="", tol_abs=0.0):
    erro = abs(modelo - ref) / max(abs(ref), 1e-30)
    ok = erro <= tol_rel or abs(modelo - ref) <= tol_abs
    tol = f"{100 * tol_rel:.2g}%" if tol_abs <= 0 else f"{tol_abs:g} {unidade}"
    linhas.append([caso, round(modelo, 9), round(ref, 9), unidade, f"{100 * erro:.3g}%", tol,
                   "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: modelo {modelo:.6g} {unidade}, referência {ref:.6g} ({100 * erro:.3g}%)")


def condicao(caso, ok, a, b, unidade, criterio):
    linhas.append([caso, round(a, 6), round(b, 6), unidade, "—", criterio, "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: {a:.5g} e {b:.5g} {unidade} ({criterio})")


def informar(caso, a, unidade, nota):
    linhas.append([caso, round(a, 4), "", unidade, "—", nota, "informativo"])
    print(f"--  {caso}: {a:.4g} {unidade} ({nota})")


# ----------------------------------------------------------------------------- prisma de teste
E, NU, FC, ALFA, K = 10e9, 0.2, 20e6, 8e-6, 2.0
geo = Geometria(comprimento_bloco=0.10, largura=0.10, altura_fiada=0.20, diametro_furo=0, furos_por_bloco=0,
                blocos_por_fiada=1, n_fiadas=1, amarracao="prumo", tipo_junta="argamassa")
A, H = geo.comprimento * geo.largura, geo.altura_prisma


def materiais():
    M = Materiais()
    c = M.ceramica
    c.E_bloco_bruto, c.razao_area_liquida, c.nu, c.alfa, c.T_ref = E, 1.0, NU, ALFA, 20.0
    return M


def modelo(m, k=K, tabela=None, ft=50e6):
    mod = ModeloNaoLinear(m, materiais(), DanoCeramica(ft=ft, fc=FC), None, verboso=False)
    mod.definir_fluencia(k_bloco=k, tabela_bloco=tabela)
    return mod


def uniforme(m, T):
    return np.full(m.n_nos, float(T))


def aquecer_ate(mod, m, T0, T1, n, **kw):
    mod.aquecer({i: uniforme(m, T0 + (T1 - T0) * i / n) for i in range(1, n + 1)}, **kw)


def lateral(mod, m):
    ux, x = mod.u[0::3], m.xyz[:, 0]
    return (ux[np.isclose(x, geo.comprimento)].mean() - ux[np.isclose(x, 0.0)].mean()) / geo.comprimento


# ---------------------------------------------------------------------------- V1
sig = -6e6
eth = ALFA * 500
ref_z = eth + sig / E + K * sig / FC * eth
ref_x = eth - NU * sig / E - NU * K * sig / FC * eth
for el in ("W6", "H8", "W15"):
    m = ml.malha_prisma(geo, h=0.025, dz=0.05, t_junta=0.01, simetria=True, placas=False, capeamento=False, elemento=el)
    mod = modelo(m)
    mod.aquecer({0.0: uniforme(m, 20)}, pre_carga_N=-sig * A)
    aquecer_ate(mod, m, 20, 520, 10, pre_carga_N=-sig * A)
    registrar(f"V1 deformação axial sob carga constante, {el}", mod.u[3 * mod.topo + 2].mean() / H * 1e3,
              ref_z * 1e3, 1e-9, "‰")
    registrar(f"V1 deformação lateral (Poisson transitório = ν), {el}", lateral(mod, m) * 1e3, ref_x * 1e3, 1e-9, "‰")

m = ml.malha_prisma(geo, h=0.025, dz=0.05, t_junta=0.01, simetria=True, placas=False, capeamento=False)

# ---------------------------------------------------------------------------- V2
ref2 = FC / K * (np.exp(-E * K * ALFA * 500 / FC) - 1)
for tol_tr, tol in ((None, 1e-2), (1e-5, 1e-3)):
    mod = modelo(m)
    if tol_tr:
        mod.tol_tr = tol_tr
    mod.aquecer({0.0: uniforme(m, 20)}, pre_carga_N=0.0)
    aquecer_ate(mod, m, 20, 520, 10, topo_fixo=True)
    s2 = float(mod.tensao_nominal()[:, 2].mean())
    registrar(f"V2 tensão no aquecimento restringido a 520 °C (relaxação), tol_tr = {mod.tol_tr:g}", s2 / 1e6,
              ref2 / 1e6, tol, "MPa")
informar("V2 tensão sem fluência no mesmo aquecimento (-E α ΔT)", -E * ALFA * 500 / 1e6, "MPa", "para comparação")

# ---------------------------------------------------------------------------- V3
mod = modelo(m)
P = 6e6 * A
mod.aquecer({0.0: uniforme(m, 20)}, pre_carga_N=P)
aquecer_ate(mod, m, 20, 520, 5, pre_carga_N=P)
e520 = float(mod.eps_tr[:, 2].mean())
aquecer_ate(mod, m, 520, 270, 5, pre_carga_N=P)
e_resf = float(mod.eps_tr[:, 2].mean())
aquecer_ate(mod, m, 270, 520, 5, pre_carga_N=P)
e_reaq = float(mod.eps_tr[:, 2].mean())
aquecer_ate(mod, m, 520, 620, 2, pre_carga_N=P)
e620 = float(mod.eps_tr[:, 2].mean())
registrar("V3 ε_tr a 520 °C = k σ/fc ε_th", e520 * 1e3, -K * 0.3 * ALFA * 500 * 1e3, 1e-9, "‰")
registrar("V3 ε_tr depois de resfriar a 270 °C (não volta)", e_resf * 1e3, e520 * 1e3, 1e-12, "‰")
registrar("V3 ε_tr depois de reaquecer a 520 °C (não soma)", e_reaq * 1e3, e520 * 1e3, 1e-12, "‰")
registrar("V3 ε_tr a 620 °C (soma só de 520 a 620 °C)", e620 * 1e3, -K * 0.3 * ALFA * 600 * 1e3, 1e-9, "‰")

# ---------------------------------------------------------------------------- V4
mod = modelo(m)
mod.aquecer({0.0: uniforme(m, 20)}, pre_carga_N=-1e6 * A)          # tração de 1 MPa
aquecer_ate(mod, m, 20, 520, 5, pre_carga_N=-1e6 * A)
condicao("V4 sob tração não há fluência transitória", float(np.abs(mod.eps_tr).max()) < 1e-12,
         float(np.abs(mod.eps_tr).max()), 0.0, "", "< 1e-12 (arredondamento)")

# ---------------------------------------------------------------------------- V5
tab = [(20.0, 0.0), (300.0, 0.002), (520.0, 0.010)]
mod = modelo(m, k=0.0, tabela=tab)
mod.aquecer({0.0: uniforme(m, 20)}, pre_carga_N=P)
aquecer_ate(mod, m, 20, 410, 6, pre_carga_N=P)
phi410 = float(np.interp(410, *zip(*tab)))
registrar("V5 ε_tr com φ(T) tabelado, a 410 °C", float(mod.eps_tr[:, 2].mean()) * 1e3, -phi410 * 0.3 * 1e3, 1e-9, "‰")

# ---------------------------------------------------------------------------- V6
t = 0.14
sec = SecaoFibras.retangulo(t, 20)
Ew, a_w, fc_w = 5e9, 6e-6, 10e6
c = lambda T: K * a_w * np.maximum(np.asarray(T, float) - 20, 0.0) / fc_w   # noqa: E731
mw = ModeloParede(1.0, sec, LeiLinear(Ew), lambda T: a_w * (np.asarray(T, float) - 20), n_el=4, T_inicial=20,
                  regularizar=False, fluencia=c)
Pw = 3e6 * t
mw.carregar(Pw, 0.0, 0.0)
mw.incendio(Pw, 0.0, 0.0, [(60.0 * i, np.full(len(sec.y), 20 + 50.0 * i)) for i in range(11)])
ref6 = a_w * 500 - 3e6 / Ew - K * a_w * 500 * 3e6 / fc_w
registrar("V6 parede: deslocamento do topo sob carga constante e aquecimento", mw.deslocamento_topo * 1e3,
          ref6 * 1e3, 1e-9, "mm")
p6 = pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", "macico_10x6x22"))
p6["modelo"].update({"analise": "parede", "k_tr_bloco": 1.5})
g6 = pj.geometria(p6)
M6 = pj.materiais(p6, g6)
_, _, sec6 = pj._secao_parede(p6, g6)
coef = pj._fluencia_parede(p6, g6, M6, sec6)
T6 = 400.0
fa = g6.espessura_junta / g6.altura_fiada
ref_c = ((1 - fa) * 1.5 * float(M6.ceramica.eps_termica(np.array([T6]))[0]) / (p6["modelo"]["fc"] * 1e6)
         + fa * 2.0 * float(M6.argamassa.eps_termica(np.array([T6]))[0]) / (p6["ceramica"]["arg_fc"] * 1e6))
registrar("V6 coeficiente da parede: bloco e junta em série, a 400 °C", float(coef(np.full(len(sec6.y), T6))[0]) * 1e12,
          ref_c * 1e12, 1e-12, "×10⁻¹² /Pa")

# ---------------------------------------------------------------------------- V7
base = pj.completar(pj.aplicar(pj.projeto_padrao(), pj.PREDEFINIDOS["Parede BCEPQ 2,94 × 3,00 m — ISO 834 sob carga"]))
base["malha"]["n_elementos"] = 12
base["ensaio"]["duracao_iso"] = 120


def parede(**mods):
    """(tempo de resistência, flecha máxima no fim em mm) da parede BCEPQ no ISO 834."""
    import tempfile
    q = {**base, "modelo": {**base["modelo"], **mods}}
    with tempfile.TemporaryDirectory() as d:
        r = pj.executar(q, d, log=lambda s: None)
    v = {nome: x for nome, x, _ in r["resumo"]}
    return v["Tempo de resistência (perda de equilíbrio)"], float(v["Flecha máxima no fim"])


t_off, f_off = parede(fluencia=False)
t_on, f_on = parede(fluencia=True)
condicao(f"V7 parede BCEPQ, k_tr = 0: flecha no fim igual com a fluência ligada (tempo {t_on} min)",
         t_on == t_off and f_on == f_off, f_on, f_off, "mm", "iguais")
if COMPLETO:
    base["ensaio"]["duracao_iso"] = 240
    t0, f0 = parede(fluencia=False)
    for k in (1.0, 2.0):
        tk, fk = parede(fluencia=True, k_tr_bloco=k)
        informar(f"V7 parede BCEPQ, 240 min, k_tr = {k:g} (hipotético): tempo {tk} min; flecha no fim", fk, "mm",
                 f"sem fluência: {t0} min e {f0:.1f} mm")

# ---------------------------------------------------------------------------- V8
p8 = pj.completar(pj.projeto_padrao())
condicao("V8 padrões: k_tr da argamassa 2,0; bloco e graute 0",
         (p8["modelo"]["k_tr_argamassa"], p8["modelo"]["k_tr_bloco"], p8["modelo"]["k_tr_graute"]) == (2.0, 0.0, 0.0),
         p8["modelo"]["k_tr_argamassa"], p8["modelo"]["k_tr_bloco"], "", "2,0 / 0 / 0")
ruins = [k for k, v in pj.PREDEFINIDOS.items() if pj.validar(pj.completar(pj.aplicar(pj.projeto_padrao(), v)))]
condicao("V8 projetos prontos passam na validação", not ruins, len(pj.PREDEFINIDOS) - len(ruins),
         len(pj.PREDEFINIDOS), "projetos", "todos")

# ---------------------------------------------------------------------------- V9
if COMPLETO:
    import tempfile

    def dano_transiente(fluencia):
        q = pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", "macico_10x6x22"))
        q["modelo"].update({"analise": "nao_linear", "fluencia": fluencia})
        q["ensaio"].update({"tipo": "transiente", "carga": 46, "taxa": 5, "T_alvo": 800, "patamar": 30,
                            "T_inicial": 20})
        with tempfile.TemporaryDirectory() as d:
            pj.executar(q, d, log=lambda s: None)
            a = np.loadtxt(os.path.join(d, "ensaio_transiente.csv"), delimiter=",", skiprows=1)
        return float(a[-1, -1])
    d_off, d_on = dano_transiente(False), dano_transiente(True)
    condicao("V9 tijolo maciço no ensaio transiente: a fluência da argamassa reduz o dano (%)", d_on < d_off,
             d_on, d_off, "%", "com < sem")

saida = os.path.join(os.path.dirname(__file__), "..", "resultados", "verificacao_fluencia.csv")
salvar_csv(saida, linhas, ["caso", "modelo", "referência", "unidade", "erro", "tolerância", "situação"])
nf = sum(1 for l in linhas if l[-1] == "FALHOU")
nv = sum(1 for l in linhas if l[-1] != "informativo")
print(f"\n{nv - nf}/{nv} casos dentro da tolerância"
      + (f" (mais {len(linhas) - nv} resultado(s) informativo(s))" if len(linhas) > nv else "")
      + f". Gravado em {os.path.abspath(saida)}")
sys.exit(1 if nf else 0)
