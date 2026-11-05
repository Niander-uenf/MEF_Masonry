"""Verificação do modelo de parede (parede.py) contra soluções conhecidas.

V1  Amplificação da imperfeição senoidal (coluna biarticulada elástica):
    flecha = a0 (P/Pcr') / (1 - P/Pcr'), com Pcr' = pi² EI / H'² e H' encurtado por P/EA
V2  Arqueamento térmico livre (gradiente linear): arco de raio t / (alfa dT)
V3  Carga excêntrica nas duas pontas (fórmula da secante): e [sec(pi/2 sqrt(P/Pcr)) - 1]
V4  Seção de fibras da malha 2D: razão área líquida/bruta e inércia contra a geometria
V5  Parede curta a 20 °C, carga centrada: capacidade = resistência do prisma (área bruta)
V6  Parede curta aquecida por igual de 20 a 600 °C sob 1% da carga: dilatação livre e
    capacidade = resistência do prisma a 600 °C
V7  Convergência com o número de elementos (parede de 3 m): capacidade a frio e tempo de
    resistência no ISO 834 com 40% da capacidade

Uso: python exemplos/10_verificacao_parede.py [--rapido]   (--rapido pula a V7)
Grava resultados/verificacao_parede.csv.
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, termico as tm, curvas
from mef_masonry.parametros import Geometria
from mef_masonry.parede import SecaoFibras, LeiLinear, LeiParede, ModeloParede
from mef_masonry.dados_parede import CURVAS_PRISMA_BCEPQ, H_REF, A_BRUTA_REF
from mef_masonry.termico import Materiais, Contorno
from mef_masonry.saida import salvar_csv

ap = argparse.ArgumentParser()
ap.add_argument("--rapido", action="store_true")
a = ap.parse_args()

linhas = []


def registrar(caso, modelo, ref, tol_rel, unidade=""):
    erro = abs(modelo - ref) / abs(ref)
    ok = erro <= tol_rel
    linhas.append([caso, round(modelo, 6), round(ref, 6), unidade, f"{100 * erro:.3f}%", f"{100 * tol_rel:.1f}%",
                   "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: modelo {modelo:.6g} {unidade}, referência {ref:.6g} ({100 * erro:.3f}%)")


zero = lambda T: 0.0 * np.asarray(T, float)  # noqa: E731
t, E, H = 0.14, 1e9, 3.0
ret = SecaoFibras.retangulo(t, 200)
EI, EA = E * ret.inercia, E * ret.area
Pcr = math.pi ** 2 * EI / H ** 2

# V1
a0 = H / 1000
for r in (0.5, 0.9):
    m = ModeloParede(H, ret, LeiLinear(E), zero, n_el=24, imperfeicao=a0)
    m.carregar(r * Pcr, 0, 0, n_passos=20)
    Pc = Pcr / (1 - r * Pcr / EA) ** 2
    ref = a0 * (r * Pcr / Pc) / (1 - r * Pcr / Pc)
    registrar(f"V1 flecha adicional com P = {r:.1f} Pcr", -m.flecha_meio * 1e3, ref * 1e3, 0.005, "mm")

# V2
alfa, dT = 6e-6, 800.0
m = ModeloParede(H, ret, LeiLinear(E), lambda T: alfa * (np.asarray(T) - 20.0), n_el=24, T_inicial=20.0)
Tf = 20.0 + dT * (0.5 - ret.y / t)
d, info, _ = m.equilibrar(np.zeros(m.nd), Tf)
m.aceitar(d, Tf, info)
R = t / (alfa * dT)
registrar("V2 flecha térmica livre", -m.flecha_meio * 1e3, (R - math.sqrt(R * R - H * H / 4)) * 1e3, 0.005, "mm")
Lm = H * (1 + alfa * dT / 2)                     # eixo médio alongado, curvado em arco de raio R
corda = 2 * R * math.sin(Lm / (2 * R))
registrar("V2 deslocamento do topo (corda do arco)", m.deslocamento_topo * 1e3, (corda - H) * 1e3, 0.005, "mm")

# V3
e = 0.010
m = ModeloParede(H, ret, LeiLinear(E), zero, n_el=24)
m.carregar(0.5 * Pcr, e, e, n_passos=20)
Pc = Pcr / (1 - 0.5 * Pcr / EA) ** 2
registrar("V3 flecha com excentricidade de 10 mm, P = 0,5 Pcr", m.flecha_meio * 1e3,
          e * (1 / math.cos(math.pi / 2 * math.sqrt(0.5 * Pcr / Pc)) - 1) * 1e3, 0.005, "mm")

# V4
geo1 = Geometria(blocos_por_fiada=1.0)
m2 = ml.malha_secao(geo1, h=0.005, simetria=False, faixas=())
sec = SecaoFibras.da_malha(m2)
Rf = geo1.diametro_furo / 2
nfu = len(geo1.centros_furos)
razao = 1 - nfu * math.pi * Rf ** 2 / (geo1.comprimento * t)
I_an = t ** 3 / 12 - nfu * math.pi * Rf ** 4 / 4 / geo1.comprimento
registrar("V4 razão área líquida/bruta das fibras", sec.area / t, razao, 0.002)
registrar("V4 inércia por metro das fibras", sec.inercia * 1e6, I_an * 1e6, 0.002, "×10⁻⁶ m⁴/m")

# V5 e V6
A_liq_ref = A_BRUTA_REF * (1 - 3 * math.pi * 0.042 ** 2 / A_BRUTA_REF)
lei = LeiParede.de_prisma(CURVAS_PRISMA_BCEPQ, H_REF, A_liq_ref)
M = Materiais(); M.ceramica.T_ref = 20.0
f_bruta = {T: max(c[1] for c in v) * 1e3 / A_BRUTA_REF / 1e6 for T, v in CURVAS_PRISMA_BCEPQ.items()}
m = ModeloParede(0.3, sec, lei, M.ceramica.eps_termica, n_el=4, T_inicial=20.0)
P = m.capacidade(0, 0, tol=1e-4)
registrar("V5 capacidade da parede curta a 20 °C (área bruta)", P / t / 1e6, f_bruta[30], 0.003, "MPa")
m = ModeloParede(0.3, sec, lei, M.ceramica.eps_termica, n_el=4, T_inicial=20.0)
P0 = 0.01 * f_bruta[30] * 1e6 * t                 # 1% da capacidade, para o equilíbrio ser único
m.carregar(P0, 0, 0)
u0 = m.deslocamento_topo
campos = [(60.0 * k, np.full(len(sec.y), 20.0 + 29.0 * k)) for k in range(21)]      # 20 -> 600 °C
tc, _ = m.incendio(P0, 0, 0, campos)
registrar("V6 dilatação livre de 20 a 600 °C sob 1% da carga", (m.deslocamento_topo - u0) * 1e3,
          0.3 * float(M.ceramica.eps_termica(600.0)) * 1e3, 0.01, "mm")
P = m.capacidade(0, 0, tol=1e-4)
registrar("V6 capacidade a 600 °C (área bruta)", P / t / 1e6, f_bruta[600], 0.003, "MPa")

# V7
if not a.rapido:
    Ti = 20.0
    exp = tm.comprimentos_contorno_2d(m2, faces=("y0",)); nexp = tm.comprimentos_contorno_2d(m2, faces=("y1",))
    ctr = [Contorno(exp, curvas.iso834(Ti), 25, 0.8), Contorno(nexp, curvas.constante(Ti), 4.0, 0.8)]
    r = tm.resolver_2d(m2, M, ctr, 120 * 60, dt=30, T_inicial=Ti, instantes_campo=list(np.arange(60.0, 7201, 60.0)))
    campos = [(0.0, np.full(len(sec.y), Ti))] + [(tt, sec.temperaturas(m2, r.campos[tt])) for tt in sorted(r.campos)]
    res = {}
    for n in (16, 32, 64):
        mo = ModeloParede(H, sec, lei, M.ceramica.eps_termica, n_el=n, imperfeicao=H / 450, T_inicial=Ti)
        Pm = mo.capacidade(0, 0)
        mo = ModeloParede(H, sec, lei, M.ceramica.eps_termica, n_el=n, imperfeicao=H / 450, T_inicial=Ti)
        mo.carregar(0.4 * Pm, 0, 0)
        tc, _ = mo.incendio(0.4 * Pm, 0, 0, campos)
        res[n] = (Pm, tc / 60 if tc else np.inf)
        print(f"    n = {n}: capacidade {Pm / 1e3:.2f} kN/m, tempo de resistência com 40%: {res[n][1]:.1f} min")
    for n in (16, 32):
        registrar(f"V7 capacidade a frio com {n} elementos (ref. 64)", res[n][0] / 1e3, res[64][0] / 1e3, 0.01, "kN/m")
        registrar(f"V7 tempo de resistência com {n} elementos (ref. 64)", res[n][1], res[64][1],
                  0.02 if n == 32 else 0.20, "min")

saida = os.path.join(os.path.dirname(__file__), "..", "resultados", "verificacao_parede.csv")
salvar_csv(saida, linhas, ["caso", "modelo", "referência", "unidade", "erro", "tolerância", "situação"])
nf = sum(1 for l in linhas if l[-1] != "ok")
print(f"\n{len(linhas) - nf}/{len(linhas)} casos dentro da tolerância. Gravado em {os.path.abspath(saida)}")
sys.exit(1 if nf else 0)
