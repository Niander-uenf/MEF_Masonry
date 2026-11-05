"""Simula os ensaios da UFV no prisma 3D: aquecimento a 5 °C/min, 30 min de
patamar e carga máxima do ensaio aplicada ao fim do patamar.

Uso: python exemplos/03_ensaio_ufv_3d.py [--temperaturas 600 800 1000] [--h 8] [--dz 12]
Se existirem arquivos dados/termopares_<T>.csv (colunas: tempo_min, T_C),
eles são sobrepostos ao gráfico para validação.
"""
import argparse
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, estudos
from mef_masonry.parametros import Geometria, Ensaio, CARGAS_UFV_KGF
from mef_masonry.termico import Materiais
from mef_masonry.saida import salvar_csv

ap = argparse.ArgumentParser()
ap.add_argument("--temperaturas", type=float, nargs="+", default=[600, 800, 1000])
ap.add_argument("--h", type=float, default=8.0, help="tamanho de elemento na seção, mm")
ap.add_argument("--dz", type=float, default=12.0, help="passo vertical, mm")
ap.add_argument("--dt", type=float, default=30.0, help="passo de tempo, s")
ap.add_argument("--hc", type=float, default=10.0, help="convecção no forno, W/m².K")
ap.add_argument("--emissividade", type=float, default=0.8)
ap.add_argument("--sem_vtk", action="store_true")
a = ap.parse_args()

BASE = os.path.join(os.path.dirname(__file__), "..")
SAIDA = os.path.join(BASE, "resultados")
os.makedirs(SAIDA, exist_ok=True)
geo = Geometria()
M = Materiais()
m = ml.malha_prisma(geo, h=a.h / 1000, dz=a.dz / 1000, t_junta=M.junta.espessura)
print("Malha:", ml.resumo(m))
print(f"E da junta = {M.junta.E / 1e6:.1f} MPa (recalibre com o exemplo 02 se mudar a malha)\n")

todos, historicos = [], {}
for T_alvo in a.temperaturas:
    t0 = time.time()
    ens = Ensaio(T_alvo=T_alvo, carga_kgf=CARGAS_UFV_KGF.get(int(T_alvo), 10000.0),
                 h_conveccao=a.hc, emissividade=a.emissividade)
    r = estudos.simular_ensaio_3d(m, M, ens, geo, dt=a.dt,
                                  pasta_vtk=None if a.sem_vtk else os.path.join(SAIDA, "vtk"))
    ind = r["indicadores"]
    todos.append(ind)
    historicos[T_alvo] = r
    print(f"Ensaio de {T_alvo:.0f} °C (carga {ens.carga_kgf:.1f} kgf), {time.time() - t0:.0f} s")
    print(f"  temperatura na carga: superfície {ind['T_superficie']:.0f} °C, "
          f"centro do septo {ind['T_septo']:.0f} °C, "
          f"média do terço central {ind['T_media_ceramica_terco']:.0f} °C")
    print(f"  dilatação térmica entre placas {ind['dilatacao_entre_placas_mm']:.3f} mm; "
          f"encurtamento pela carga {ind['encurtamento_carga_mm']:.3f} mm")
    print(f"  tração máxima na cerâmica: térmica {ind['tracao_max_termico_MPa']:.2f} MPa "
          f"(p99 {ind['tracao_p99_termico_MPa']:.2f}); total {ind['tracao_max_total_MPa']:.2f} MPa "
          f"(p99 {ind['tracao_p99_total_MPa']:.2f}); fiadas centrais {ind['tracao_max_total_fiadas_centrais_MPa']:.2f} MPa")
    print(f"  cerâmica com tração > 1 MPa: {100 * ind['fracao_ceramica_tracao_acima_1MPa']:.0f}% do volume\n")

    rt = r["termico"]
    linhas = [[f"{t / 60:.2f}", f"{g:.1f}"] + [f"{rt.sondas[k][i]:.1f}" for k in r["sondas"]]
              for i, (t, g) in enumerate(zip(rt.t, rt.T_gas))]
    salvar_csv(os.path.join(SAIDA, f"historico_{int(T_alvo)}C.csv"), linhas,
               ["tempo_min", "forno_C"] + r["sondas"])

chaves = list(todos[0].keys())
salvar_csv(os.path.join(SAIDA, "ensaio_ufv_3d.csv"), [[str(i[k]) for k in chaves] for i in todos], chaves)

# figura: temperaturas ao longo do ensaio, com termopares se houver
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, len(a.temperaturas), figsize=(5 * len(a.temperaturas), 3.8), squeeze=False)
for ax, T_alvo in zip(axs[0], a.temperaturas):
    r = historicos[T_alvo]; rt = r["termico"]
    ax.plot(rt.t / 60, rt.T_gas, "k--", lw=1, label="forno")
    for k, est in zip(r["sondas"], ("-", "-", ":")):
        ax.plot(rt.t / 60, rt.sondas[k], est, lw=1.5, label=k)
    arq = os.path.join(BASE, "dados", f"termopares_{int(T_alvo)}.csv")
    if os.path.exists(arq):
        d = np.genfromtxt(arq, delimiter=",", names=True)
        ax.plot(d[d.dtype.names[0]], d[d.dtype.names[1]], "o", ms=3, label="termopar (ensaio)")
    ax.axvline(r["t_fim_rampa"] / 60, color="0.6", lw=0.8)
    ax.set_title(f"Ensaio de {T_alvo:.0f} °C"); ax.set_xlabel("tempo (min)"); ax.set_ylabel("°C")
    ax.grid(alpha=0.3)
axs[0][0].legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(SAIDA, "ensaio_ufv_temperaturas.png"), dpi=130)
print("Resultados em", os.path.abspath(SAIDA))
