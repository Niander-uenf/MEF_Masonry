"""Estudo paramétrico na seção 2D (rápido): diâmetro dos furos x condutividade.

Para cada combinação calcula:
  - ensaio da UFV (forno a 5 °C/min até 600 °C + 30 min): temperatura mínima e
    média na carga, tração térmica máxima e percentil 99;
  - parede em incêndio-padrão ISO 834 numa face: tempo de isolamento térmico
    (critério de 140/180 °C na face não exposta) e tração térmica aos 30 min.

Uso: python exemplos/04_parametrico_secao.py [--rapido]
"""
import argparse
import itertools
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import estudos
from mef_masonry.parametros import Geometria, Ensaio
from mef_masonry.termico import Materiais
from mef_masonry.saida import salvar_csv

ap = argparse.ArgumentParser()
ap.add_argument("--rapido", action="store_true", help="malha e passo mais grossos")
ap.add_argument("--diametros", type=float, nargs="+", default=[0, 60, 70, 84, 95], help="mm")
ap.add_argument("--k", type=float, nargs="+", default=[0.40, 0.52, 0.70], help="W/m.K")
ap.add_argument("--T_alvo", type=float, default=600.0)
a = ap.parse_args()
h, dt = (0.008, 60.0) if a.rapido else (0.005, 30.0)

SAIDA = os.path.join(os.path.dirname(__file__), "..", "resultados")
os.makedirs(SAIDA, exist_ok=True)
geo = Geometria()
linhas, res = [], {}
t0 = time.time()
for D, k in itertools.product(a.diametros, a.k):
    M = Materiais(); M.ceramica.k = k
    f = estudos.secao_forno(M, Ensaio(T_alvo=a.T_alvo), geo, h=h, dt=dt, diametro=D / 1000)
    p = estudos.secao_parede_iso834(M, geo, h=h, dt=dt, diametro=D / 1000, t_max_min=480)
    res[(D, k)] = (f, p)
    linhas.append([D, k, f"{f['razao_area_liquida']:.3f}", f"{f['T_min']:.1f}", f"{f['T_media']:.1f}",
                   f"{f['s1_max_MPa']:.2f}", f"{f['s1_p99_MPa']:.2f}", f"{p['t_isolamento_min']:.1f}",
                   f"{p.get('s1_p99_MPa', np.nan):.2f}"])
    print(f"D = {D:3.0f} mm, k = {k:.2f}: T_min {f['T_min']:5.0f} °C, tração p99 {f['s1_p99_MPa']:.2f} MPa, "
          f"isolamento ISO 834 {p['t_isolamento_min']:.0f} min  ({time.time() - t0:.0f} s)")

cab = ["D_furo_mm", "k_W_mK", "razao_area_liquida", "forno_T_min_C", "forno_T_media_C",
       "forno_tracao_max_MPa", "forno_tracao_p99_MPa", "ISO834_isolamento_min", "ISO834_tracao_p99_30min_MPa"]
salvar_csv(os.path.join(SAIDA, "parametrico_secao.csv"), linhas, cab)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 2, figsize=(10, 7))
paineis = [(axs[0, 0], lambda f, p: f["T_min"], f"Ponto mais frio na carga (ensaio {a.T_alvo:.0f} °C), °C"),
           (axs[0, 1], lambda f, p: f["s1_p99_MPa"], "Tração térmica, percentil 99 (ensaio), MPa"),
           (axs[1, 0], lambda f, p: p["t_isolamento_min"], "Isolamento em ISO 834 (parede), min"),
           (axs[1, 1], lambda f, p: f["razao_area_liquida"], "Área líquida / bruta")]
cores = plt.cm.Blues(np.linspace(0.45, 0.95, len(a.k)))
for ax, fn, tit in paineis:
    for c, k in zip(cores, a.k):
        ax.plot(a.diametros, [fn(*res[(D, k)]) for D in a.diametros], "o-", color=c, label=f"k = {k:.2f}")
    ax.set_title(tit, fontsize=9); ax.set_xlabel("diâmetro do furo (mm)"); ax.grid(alpha=0.3)
axs[0, 0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(SAIDA, "parametrico_secao.png"), dpi=130)
print("Resultados em", os.path.abspath(SAIDA), f"({time.time() - t0:.0f} s)")
