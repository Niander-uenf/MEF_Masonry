"""Verificação da fase 2: células retangulares e graute.

V1  Área líquida das células retangulares: a área dos elementos de bloco na
    malha 2D, fiada a fiada (com contrafiada as células mudam de lugar), é igual
    à área bruta menos as células calculadas pela geometria. Inclui um bloco de
    paredes vazadas (células pequenas nas paredes externas).
V2  Volume de graute no prisma 3D (junta a prumo, argamassa, todas as células
    grauteadas): área grauteada × altura entre os capeamentos.
V3  Rigidez em paralelo: uma fiada de blocos com células alternadas grauteadas,
    com o mesmo Poisson nos dois materiais (a deformação lateral é a mesma, sem
    interação): encurtamento = P H / (Eb Ab + Eg Ag), exato.
V4  Propriedades do graute chegam aos solvers: condutividade, capacidade,
    módulo e deformação térmica, inclusive as tabelas da EN 1992-1-2.
V5  Seção 2D e parede: o graute entra como material na térmica e como fibra
    na parede; a área das fibras é a área resistente (bloco + graute).
V6  Graute aumenta a resistência do prisma não linear (bloco de concreto B1 de
    Camacho et al., 2015) e comparação, informativa, com os ensaios: prisma oco
    5,63 MPa e grauteado com G1 11,29 MPa (área bruta). Só com --completo
    (cerca de 10 minutos).
V7  Biblioteca: os blocos já liberados passam na validação; os ainda não liberados
    (tijolo de solo-cimento, sem os furos medidos) são recusados.

Uso: python exemplos/12_verificacao_celulas_graute.py [--completo]
Grava resultados/verificacao_celulas_graute.csv.
"""
import os
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, termico as tm, mecanico as mc, projeto as pj, biblioteca as bib
from mef_masonry.malha import CERAMICA, GRAUTE, AR
from mef_masonry.parametros import Geometria, EN1992_E_FATOR, EN1992_EPS_TERMICA
from mef_masonry.termico import Materiais
from mef_masonry.parede import SecaoFibras
from mef_masonry.saida import salvar_csv

COMPLETO = "--completo" in sys.argv
linhas = []


def registrar(caso, modelo, ref, tol_rel, unidade=""):
    erro = abs(modelo - ref) / max(abs(ref), 1e-30)
    ok = erro <= tol_rel
    linhas.append([caso, round(modelo, 6), round(ref, 6), unidade, f"{100 * erro:.3f}%", f"{100 * tol_rel:.1f}%",
                   "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: modelo {modelo:.6g} {unidade}, referência {ref:.6g} ({100 * erro:.3f}%)")


def comparar(caso, modelo, ref, unidade=""):
    """Comparação com ensaio (informativa: não conta como falha)."""
    erro = (modelo - ref) / ref
    linhas.append([caso, round(modelo, 4), round(ref, 4), unidade, f"{100 * erro:+.1f}%", "—", "comparação"])
    print(f"--  {caso}: modelo {modelo:.4g} {unidade}, ensaio {ref:.4g} ({100 * erro:+.1f}%)")


def condicao(caso, ok, a, b, unidade, criterio):
    linhas.append([caso, round(a, 6), round(b, 6), unidade, "—", criterio, "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: {a:.4g} e {b:.4g} {unidade} ({criterio})")


# ---------------------------------------------------------------------------- V1
tj = 0.010
casos = {
    "concreto 14×19×39, bloco e meio, contrafiada": Geometria(
        forma_furo="retangular", comprimento_bloco=0.40, largura=0.14, altura_fiada=0.20, furos_por_bloco=2,
        blocos_por_fiada=1.5, n_fiadas=2, amarracao="contrafiada", espessura_junta=tj, tipo_junta="argamassa",
        parede_externa=0.032, septo=0.030, parede_topo=0.030),
    "cerâmico de paredes vazadas": Geometria(
        forma_furo="retangular", comprimento_bloco=0.30, largura=0.14, altura_fiada=0.20, furos_por_bloco=2,
        blocos_por_fiada=1, n_fiadas=1, amarracao="prumo", espessura_junta=tj, tipo_junta="argamassa",
        parede_externa=0.030, septo=0.008, parede_topo=0.008, celulas_parede=6, parede_minima=0.008),
}
for nome, g in casos.items():
    assert not g.validar(), g.validar()
    m2 = ml.malha_secao_retangular(g, h=0.008, simetria=False, faixas=g.posicoes_juntas, t_junta=tj)
    a = np.abs(m2.areas())
    for i in range(g.n_fiadas):
        furo, _ = ml.estados_fiada(m2, g, i)
        ref = g.comprimento * g.largura - g.area_furos(i, t_junta=tj)
        registrar(f"V1 área do bloco na fiada {i + 1}, {nome}", a[~furo].sum() * 1e4, ref * 1e4, 1e-9, "cm²")

# ---------------------------------------------------------------------------- V2
gv = Geometria(forma_furo="retangular", comprimento_bloco=0.30, largura=0.14, altura_fiada=0.20, furos_por_bloco=2,
               blocos_por_fiada=1, n_fiadas=3, amarracao="prumo", espessura_junta=tj, tipo_junta="argamassa",
               parede_externa=0.028, septo=0.028, parede_topo=0.028, grauteamento="todos")
t_cap = 0.003
mv = ml.malha_prisma(gv, h=0.01, dz=0.02, t_junta=tj, simetria=False, placas=False, capeamento=True)
xyz = mv.xyz[mv.cun]
p = xyz[:, :3, :2]
ar = 0.5 * np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
vol = ar * (xyz[:, 3, 2] - xyz[:, 0, 2])
A_g = gv.area_furos(0, t_junta=tj) - gv.area_furos(0, graute=True, t_junta=tj)
registrar("V2 volume de graute no prisma de 3 fiadas", vol[mv.mat == GRAUTE].sum() * 1e6,
          A_g * (gv.altura_prisma - 2 * t_cap) * 1e6, 1e-9, "cm³")
v_ar = float(vol[mv.mat == AR].sum() * 1e6)
condicao("V2 nenhum vazio no prisma todo grauteado com assentamento total", v_ar == 0.0, v_ar, 0.0, "cm³", "zero")

# ---------------------------------------------------------------------------- V3
g3 = replace(gv, n_fiadas=1, grauteamento="alternados", furos_por_bloco=4, blocos_por_fiada=1)
m3 = ml.malha_prisma(g3, h=0.01, dz=0.02, t_junta=tj, simetria=True, placas=False, capeamento=False)
M3 = Materiais()
M3.ceramica.E_bloco_bruto, M3.ceramica.razao_area_liquida, M3.ceramica.nu = 9000e6, 1.0, 0.2
M3.graute.E_bloco_bruto, M3.graute.nu = 26000e6, 0.2
P_N = 100e3
R3 = mc.resolver_3d(m3, M3, carga_N=P_N, incluir_termico=False)["carga"]
H3 = g3.altura_prisma
d3 = -mc.encurtamento_entre_placas(m3, R3, H3)
A_b = g3.comprimento * g3.largura - g3.area_furos(0, t_junta=tj)
A_g3 = g3.area_furos(0, t_junta=tj) - g3.area_furos(0, graute=True, t_junta=tj)
registrar("V3 encurtamento em paralelo bloco + graute (ν iguais)", d3 * 1e6,
          P_N * H3 / (9000e6 * A_b + 26000e6 * A_g3) * 1e6, 1e-6, "µm")

# ---------------------------------------------------------------------------- V4
M4 = Materiais()
g = M4.graute
g.k, g.rho, g.c, g.T_ref = 1.75, 2400.0, 1000.0, 20.0
g.E_bloco_bruto, g.razao_area_liquida = 26.1e9, 1.0
mat4 = np.array([GRAUTE], dtype=np.int8)
registrar("V4 condutividade do graute", float(M4.k(mat4, np.array([300.0]))[0]), 1.75, 1e-12, "W/m·K")
registrar("V4 capacidade volumétrica do graute", float(M4.rho_c(mat4)[0]) / 1e6, 2.4, 1e-12, "MJ/m³·K")
g.tabela_E_fator, g.tabela_eps_termica = list(EN1992_E_FATOR), list(EN1992_EPS_TERMICA)
E4, _ = mc._propriedades(mat4, np.array([600.0]), M4)
registrar("V4 E do graute a 600 °C (EN 1992-1-2: 0,45 × 0,0025/0,025)", E4[0] / 1e9, 26.1 * 0.045, 1e-12, "GPa")
eps4 = mc._eps_termica_nos(mat4, np.array([[400.0] * 6]), M4)[0, 0]
registrar("V4 deformação térmica do graute a 400 °C (silicoso)", eps4 * 1e3,
          (-1.8e-4 + 9e-6 * 400 + 2.3e-11 * 400 ** 3) * 1e3, 1e-3, "‰")

# ---------------------------------------------------------------------------- V5
p5 = pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", "concreto_m15_29"))
geo5 = pj.geometria(p5)
geo1 = replace(geo5, blocos_por_fiada=1.0)
m25 = ml.malha_secao(geo1, h=0.005, simetria=False, faixas=())
mat5 = tm.material_secao(m25, geo1)
a5 = np.abs(m25.areas())
A_g5 = geo1.area_furos(0) - geo1.area_furos(0, graute=True)
registrar("V5 área de graute na seção 2D (material da térmica)", a5[mat5 == GRAUTE].sum() * 1e4, A_g5 * 1e4, 1e-9,
          "cm²")
sec = SecaoFibras.da_malha(m25, geo1)
registrar("V5 área das fibras da parede = área resistente (bloco + graute)", sec.area * sec.comprimento * 1e4,
          geo1.comprimento * geo1.largura * pj.razao_area_solida(geo1) * 1e4, 1e-9, "cm²")
registrar("V5 fibras de graute na parede", float((sec.A[sec.graute] * sec.comprimento).sum() * 1e4), A_g5 * 1e4,
          1e-9, "cm²")

# ---------------------------------------------------------------------------- V6
if COMPLETO:
    from mef_masonry.naolinear import ModeloNaoLinear
    res = {}
    for gr in ("nenhum", "todos"):
        q = pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", "concreto_m15_29"))
        q["geometria"]["grauteamento"] = gr
        q["modelo"]["analise"] = "nao_linear"
        q["ensaio"].update({"tipo": "ambiente", "T_inicial": 20})
        q["malha"].update({"h": 12, "dz": 20})
        geo = pj.geometria(q)
        M = pj.materiais(q, geo)
        m = pj.malha_3d(q, geo)
        mod = ModeloNaoLinear(m, M, pj.dano(q), pj.junta_atrito(q), verboso=False, dano_argamassa=pj.dano_arg(q),
                              dano_graute=pj.dano_gra(q))
        mod.aquecer({0.0: np.full(m.n_nos, M.ceramica.T_ref)})
        d, P = mod.carregar(0.006 * geo.altura_prisma, 60)
        res[gr] = float(P.max()) / (geo.comprimento * geo.largura) / 1e6
        print(f"    prisma {gr}: {res[gr]:.2f} MPa (área bruta)")
    condicao("V6 graute aumenta a resistência do prisma (MPa, área bruta)", res["todos"] > 1.5 * res["nenhum"],
             res["todos"], res["nenhum"], "MPa", "grauteado > 1,5 × oco")
    comparar("V6 prisma oco B1 contra o ensaio (Camacho et al., 2015)", res["nenhum"], 5.63, "MPa")
    comparar("V6 prisma grauteado B1-G1 contra o ensaio", res["todos"], 11.29, "MPa")

# ---------------------------------------------------------------------------- V7
validos = recusados = 0
for k, item in bib.BLOCOS.items():
    q = pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", k))
    erros = pj.validar(q)
    if item.get("fase", 0) <= bib.FASE_ATUAL:
        validos += not erros
        if erros:
            print("   ", k, erros)
    else:
        recusados += bool(erros)
n_ok = sum(1 for v in bib.BLOCOS.values() if v.get("fase", 0) <= bib.FASE_ATUAL)
registrar("V7 blocos da biblioteca que rodam e passam na validação", validos, n_ok, 0.0, "itens")
registrar("V7 blocos de fases futuras recusados", recusados, len(bib.BLOCOS) - n_ok, 0.0, "itens")

saida = os.path.join(os.path.dirname(__file__), "..", "resultados", "verificacao_celulas_graute.csv")
salvar_csv(saida, linhas, ["caso", "modelo", "referência", "unidade", "erro", "tolerância", "situação"])
nf = sum(1 for l in linhas if l[-1] == "FALHOU")
nv = sum(1 for l in linhas if l[-1] != "comparação")
print(f"\n{nv - nf}/{nv} casos dentro da tolerância"
      + (f" (mais {len(linhas) - nv} comparações com ensaios)" if len(linhas) > nv else "")
      + f". Gravado em {os.path.abspath(saida)}")
sys.exit(1 if nf else 0)
