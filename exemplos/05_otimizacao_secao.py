"""Exemplo de otimização da geometria da seção (modelo 2D).

Variáveis: diâmetro dos furos D e deslocamento dy dos furos em direção à face
não exposta. Objetivo: bloco mais leve (maior área de furos) que ainda atenda
a um tempo mínimo de isolamento em ISO 834 e a uma parede mínima.
Também traça a fronteira de Pareto (área líquida x isolamento) numa grade.

É um modelo de partida: troque o objetivo e as restrições pelos do estudo
(ex.: tração térmica máxima, critério de fabricação, NBR 15270).

Uso: python exemplos/05_otimizacao_secao.py [--meta 90] [--completo]
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import estudos
from mef_masonry.parametros import Geometria
from mef_masonry.termico import Materiais
from mef_masonry.saida import salvar_csv

ap = argparse.ArgumentParser()
ap.add_argument("--meta", type=float, default=90.0, help="isolamento mínimo, min")
ap.add_argument("--parede_min", type=float, default=15.0, help="parede mínima, mm")
ap.add_argument("--completo", action="store_true", help="malha fina (lento)")
a = ap.parse_args()
h, dt = (0.005, 30.0) if a.completo else (0.008, 60.0)

SAIDA = os.path.join(os.path.dirname(__file__), "..", "resultados")
os.makedirs(SAIDA, exist_ok=True)
geo = Geometria()
W = geo.largura
M = Materiais()
cache = {}


def avaliar(D_mm, dy_mm):
    chave = (round(D_mm, 2), round(dy_mm, 2))
    if chave not in cache:
        p = estudos.secao_parede_iso834(M, geo, h=h, dt=dt, diametro=D_mm / 1000,
                                        y_furos=W / 2 + dy_mm / 1000, t_max_min=300, t_tensao_min=1e9)
        cache[chave] = (p["t_isolamento_min"], p["razao_area_liquida"])
    return cache[chave]


def parede(D_mm, dy_mm):
    return W * 1000 / 2 - abs(dy_mm) - D_mm / 2


def objetivo(x):
    D, dy = x
    pen = 0.0
    if parede(D, dy) < a.parede_min:
        return 1e3 + (a.parede_min - parede(D, dy)) * 10
    t, ral = avaliar(D, dy)
    if t < a.meta:
        pen += 10 * (a.meta - t)
    return ral + pen / 100    # minimizar a área líquida (massa) com isolamento >= meta


t0 = time.time()
print("Grade para a fronteira de Pareto...")
grade = []
for D in np.arange(50, 96, 5):
    for dy in np.arange(-15, 16, 5):
        if parede(D, dy) >= a.parede_min:
            t, ral = avaliar(D, dy)
            grade.append((D, dy, t, ral))
grade = np.array(grade)
# Pareto: menor área líquida e maior isolamento
pareto = [i for i, (D, dy, t, r) in enumerate(grade)
          if not any((g[3] <= r and g[2] >= t) and (g[3] < r or g[2] > t) for g in grade)]
print(f"  {len(grade)} geometrias em {time.time() - t0:.0f} s; {len(pareto)} na fronteira de Pareto")

print(f"Otimização (evolução diferencial): isolamento >= {a.meta:.0f} min, parede >= {a.parede_min:.0f} mm")
res = differential_evolution(objetivo, bounds=[(50, 100), (-20, 20)], seed=1, maxiter=15, popsize=8,
                             tol=1e-3, polish=False)
D, dy = res.x
t, ral = avaliar(D, dy)
t_atual, ral_atual = avaliar(geo.diametro_furo * 1000, 0.0)
print(f"\nGeometria atual: D = 84 mm, dy = 0 -> isolamento {t_atual:.0f} min, área líquida {ral_atual:.3f}")
print(f"Ótimo:           D = {D:.1f} mm, dy = {dy:+.1f} mm -> isolamento {t:.0f} min, área líquida {ral:.3f}"
      f" (parede mínima {parede(D, dy):.1f} mm)")
print(f"Avaliações do modelo: {len(cache)}; tempo total {time.time() - t0:.0f} s")

salvar_csv(os.path.join(SAIDA, "otimizacao_grade.csv"),
           [[f"{g[0]:.0f}", f"{g[1]:.0f}", f"{g[2]:.1f}", f"{g[3]:.3f}", "sim" if i in pareto else ""]
            for i, g in enumerate(grade)],
           ["D_mm", "dy_mm", "isolamento_min", "razao_area_liquida", "pareto"])

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.scatter(grade[:, 3], grade[:, 2], c="0.75", s=18, label="geometrias da grade")
P = grade[pareto][np.argsort(grade[pareto][:, 3])]
ax.plot(P[:, 3], P[:, 2], "o-", color="tab:blue", ms=4, label="fronteira de Pareto")
ax.plot(ral_atual, t_atual, "s", color="k", label="BCEPQ atual (D = 84 mm)")
ax.plot(ral, t, "*", color="tab:orange", ms=12, label=f"ótimo (meta {a.meta:.0f} min)")
ax.axhline(a.meta, color="0.5", lw=0.8, ls="--")
ax.set_xlabel("área líquida / bruta (proporcional à massa)")
ax.set_ylabel("isolamento em ISO 834 (min)")
ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(SAIDA, "otimizacao_pareto.png"), dpi=130)
print("Resultados em", os.path.abspath(SAIDA))
