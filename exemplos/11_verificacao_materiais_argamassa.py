"""Verificação da fase 1: calor específico com a temperatura, juntas de
argamassa e biblioteca.

V1  Capacidade térmica variável c(T) (pico entre 100 e 200 °C): corpo de
    condutividade muito alta (Biot << 1) aquecido por convecção; compara com a
    EDO  rho c(T) A dT/dt = h P (Tg - T), integrada com tolerância 1e-10.
V2  Coluna de tijolos e argamassa (linear 3D) com nu/E igual nos dois materiais
    (mesma deformação lateral, logo tensão uniaxial uniforme): encurtamento =
    soma de h_i / E_i (associação em série, exata).
V3  Mesma coluna com nu/E maior na argamassa (0,30/1 GPa contra 0,15/3 GPa): a
    argamassa tenta expandir mais que o tijolo e o traciona lateralmente
    (mecanismo de Hilsdorf): tensão lateral de tração no tijolo junto à junta.
V4  Coluna não linear com nu/E igual: sem confinamento, a ruptura vem da
    argamassa (mais fraca), com carga máxima = fa × área. Com Poisson iguais e
    E diferentes, a argamassa fica confinada e a carga sobe (registrado à parte).
V5  Propriedades por material: condutividade, capacidade e deformação térmica
    da argamassa chegam aos solvers.
V6  Biblioteca: os itens que rodam passam na validação; os da fase 2 em diante
    são recusados.

Uso: python exemplos/11_verificacao_materiais_argamassa.py
Grava resultados/verificacao_materiais_argamassa.csv.
"""
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, termico as tm, mecanico as mc, projeto as pj, biblioteca as bib, curvas
from mef_masonry.malha import CERAMICA, ARGAMASSA
from mef_masonry.parametros import Geometria, DanoCeramica
from mef_masonry.termico import Materiais, Contorno
from mef_masonry.naolinear import ModeloNaoLinear
from mef_masonry.saida import salvar_csv

linhas = []


def registrar(caso, modelo, ref, tol_rel, unidade=""):
    erro = abs(modelo - ref) / max(abs(ref), 1e-30)
    ok = erro <= tol_rel
    linhas.append([caso, round(modelo, 6), round(ref, 6), unidade, f"{100 * erro:.3f}%", f"{100 * tol_rel:.1f}%",
                   "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: modelo {modelo:.6g} {unidade}, referência {ref:.6g} ({100 * erro:.3f}%)")


# ---------------------------------------------------------------------------- V1
a = 0.05
m2 = ml.malha_retangulo(a, a, a / 10)
M = Materiais()
M.ceramica.k, M.ceramica.rho, M.ceramica.c = 2000.0, 1000.0, 1000.0
tab_c = [(20, 1000.0), (100, 1000.0), (115, 3000.0), (200, 1000.0), (1200, 1000.0)]
M.ceramica.tabela_c = tab_c
M.ceramica.T_ref = 20.0
h, Tg = 25.0, 500.0
comp = tm.comprimentos_contorno_2d(m2)
ctr = [Contorno(comp, curvas.constante(Tg), h, 0.0)]
inst = [600.0, 1200.0, 2400.0]
r = tm.resolver_2d(m2, M, ctr, 2400.0, dt=1.0, T_inicial=20.0, instantes_campo=inst,
                   mat=np.full(len(m2.tri), CERAMICA, dtype=np.int8))
A, P = a * a, 4 * a
cT = lambda T: np.interp(T, [x[0] for x in tab_c], [x[1] for x in tab_c])  # noqa: E731
sol = solve_ivp(lambda t, T: h * P * (Tg - T) / (1000.0 * cT(T) * A), (0, 2400), [20.0], rtol=1e-10, atol=1e-10,
                dense_output=True, max_step=1.0)
for t in inst:
    registrar(f"V1 temperatura com c(T) aos {t / 60:.0f} min", float(r.campos[t].mean()), float(sol.sol(t)[0]), 0.005,
              "°C")

# ---------------------------------------------------------------------------- V2 e V3
geo = Geometria(comprimento_bloco=0.23, largura=0.10, altura_fiada=0.07, diametro_furo=0.0, furos_por_bloco=0,
                blocos_por_fiada=1, n_fiadas=3, amarracao="prumo", tipo_junta="argamassa")
mm = ml.malha_prisma(geo, h=0.01, dz=0.005, t_junta=0.01, simetria=True, placas=False, capeamento=False)
Mm = Materiais()
Mm.ceramica.E_bloco_bruto, Mm.ceramica.razao_area_liquida, Mm.ceramica.nu = 3000e6, 1.0, 0.15
Mm.argamassa.E_bloco_bruto, Mm.argamassa.nu = 1000e6, 0.05
P_N = 50e3
R = mc.resolver_3d(mm, Mm, carga_N=P_N, incluir_termico=False)["carga"]
H = geo.altura_prisma
d_mef = -mc.encurtamento_entre_placas(mm, R, H)
hs = {CERAMICA: 0.0, ARGAMASSA: 0.0}
for c in mm.camadas:
    hs[ARGAMASSA if c.tipo == "junta" else CERAMICA] += c.z1 - c.z0
Aarea = geo.comprimento * geo.largura
d_ref = P_N / Aarea * (hs[CERAMICA] / 3000e6 + hs[ARGAMASSA] / 1000e6)
registrar("V2 encurtamento da coluna tijolo + argamassa (série)", d_mef * 1e6, d_ref * 1e6, 0.002, "µm")

Mm.argamassa.nu, Mm.ceramica.nu = 0.30, 0.15
R3 = mc.resolver_3d(mm, Mm, carga_N=P_N, incluir_termico=False)["carga"]
ativos = np.nonzero(R3.ativos)[0]
mat_at = mm.mat[ativos]
zc = mm.xyz[mm.cun[ativos]].mean(1)[:, 2]
zj = [0.5 * (c.z0 + c.z1) for c in mm.camadas if c.tipo == "junta"][0]
perto = (mat_at == CERAMICA) & (np.abs(zc - zj) < 0.012)
sx_tijolo = float(np.median(R3.sig[perto, 0]))
sx_arg = float(np.median(R3.sig[mat_at == ARGAMASSA, 0]))
ok3 = sx_tijolo > 0 and sx_arg < 0
linhas.append(["V3 tijolo tracionado e argamassa comprimida lateralmente (Hilsdorf)",
               round(sx_tijolo / 1e6, 4), round(sx_arg / 1e6, 4), "MPa (σx tijolo / argamassa)", "—", "sinais",
               "ok" if ok3 else "FALHOU"])
print(f"{'ok ' if ok3 else 'XX '} V3 σx no tijolo junto à junta {sx_tijolo / 1e6:+.4f} MPa; na argamassa {sx_arg / 1e6:+.4f} MPa")

# ---------------------------------------------------------------------------- V4
Mn = Materiais()
mn = ml.malha_prisma(geo, h=0.02, dz=0.01, t_junta=0.01, simetria=True, placas=False, capeamento=False)
Pd = DanoCeramica(ft=5e6, fc=20e6)
Pa = DanoCeramica(ft=2e6, fc=4e6, du_t=0.065e-3, du_c=1.6e-3)
P4max = {}
for nome, nu_b, nu_m in (("nu/E igual", 0.15, 0.05), ("Poisson iguais", 0.18, 0.18)):
    Mn = Materiais()
    Mn.ceramica.E_bloco_bruto, Mn.ceramica.razao_area_liquida, Mn.ceramica.nu = 3000e6, 1.0, nu_b
    Mn.argamassa.E_bloco_bruto, Mn.argamassa.nu = 1000e6, nu_m
    mod = ModeloNaoLinear(mn, Mn, Pd, None, verboso=False, dano_argamassa=Pa)
    mod.aquecer({0.0: np.full(mn.n_nos, Mn.ceramica.T_ref)}, pre_carga_N=0.0)
    d4, P4 = mod.carregar(1.0e-3, 50)
    P4max[nome] = float(P4.max())
registrar("V4 carga máxima da coluna = fa × área (nu/E igual)", P4max["nu/E igual"] / 1e3, 4e6 * Aarea / 1e3, 0.02,
          "kN")
ok4 = P4max["Poisson iguais"] > P4max["nu/E igual"]
linhas.append(["V4 argamassa confinada (Poisson iguais, E menor) resiste mais", round(P4max["Poisson iguais"] / 1e3, 3),
               round(P4max["nu/E igual"] / 1e3, 3), "kN", "—", "maior", "ok" if ok4 else "FALHOU"])
print(f"{'ok ' if ok4 else 'XX '} V4 com Poisson iguais: {P4max['Poisson iguais'] / 1e3:.1f} kN "
      f"(> {P4max['nu/E igual'] / 1e3:.1f} kN sem confinamento)")

# ---------------------------------------------------------------------------- V5
M5 = Materiais()
M5.argamassa.k, M5.argamassa.rho, M5.argamassa.c, M5.argamassa.alfa = 1.15, 2000.0, 1000.0, 10e-6
M5.argamassa.T_ref = 20.0
mat5 = np.array([ARGAMASSA], dtype=np.int8)
registrar("V5 condutividade da argamassa no solver térmico", float(M5.k(mat5, np.array([300.0]))[0]), 1.15, 1e-12,
          "W/m·K")
registrar("V5 capacidade volumétrica da argamassa", float(M5.rho_c(mat5)[0]) / 1e6, 2.0, 1e-12, "MJ/m³·K")
eps5 = mc._eps_termica_nos(mat5, np.array([[520.0] * 6]), M5)[0, 0]
registrar("V5 deformação térmica da argamassa a 520 °C", eps5 * 1e3, 10e-6 * 500 * 1e3, 1e-12, "‰")

# ---------------------------------------------------------------------------- V6
validos = recusados = 0
for k, item in bib.BLOCOS.items():
    p = pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", k))
    erros = pj.validar(p)
    if item.get("fase", 0) <= bib.FASE_ATUAL:
        validos += not erros
        if erros:
            print("   ", k, erros)
    else:
        recusados += bool(erros)
n_ok = sum(1 for v in bib.BLOCOS.values() if v.get("fase", 0) <= bib.FASE_ATUAL)
n_fut = len(bib.BLOCOS) - n_ok
registrar("V6 itens da biblioteca que rodam e passam na validação", validos, n_ok, 0.0, "itens")
registrar("V6 itens de fases futuras recusados pela validação", recusados, n_fut, 0.0, "itens")

saida = os.path.join(os.path.dirname(__file__), "..", "resultados", "verificacao_materiais_argamassa.csv")
salvar_csv(saida, linhas, ["caso", "modelo", "referência", "unidade", "erro", "tolerância", "situação"])
nf = sum(1 for l in linhas if l[-1] != "ok")
print(f"\n{len(linhas) - nf}/{len(linhas)} casos dentro da tolerância. Gravado em {os.path.abspath(saida)}")
sys.exit(1 if nf else 0)
