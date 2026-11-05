"""Calibração do modelo não linear pelo prisma a 30 °C, em duas etapas
independentes:

  1. Rigidez: com o contato tão rígido quanto a cerâmica, procura o
     desvio-padrão das folgas das juntas de leito que reproduz o módulo de
     corda do prisma entre 5% e 30% da carga máxima (E = 1.100 MPa, BPQ), na
     média de três sementes do campo aleatório (a rigidez varia cerca de 20%
     entre sementes). Usa só o trecho inicial da curva.
  2. Resistência: com essa folga, procura a resistência à compressão efetiva
     da cerâmica (fc) que reproduz a carga máxima do ensaio a 30 °C.

Os valores calibrados vão para mef_masonry/parametros.py (JuntaAtrito.desvio_folga
e CALIBRACAO_30C) e servem, sem mudança, para prever os ensaios aquecidos
(07_ensaio_ufv_nao_linear.py --calibrado).

Uso: python exemplos/08_calibrar_imperfeicao.py [--E_prisma 1100] [--P30 98.7]
     [--folgas 0.024 0.032] [--sementes 1 2 3] [--h 12] [--dz 20]
Tempo: cerca de 15 min com a malha padrão (2 núcleos).
"""
import argparse
import os
import sys
import time
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, estudos
from mef_masonry.parametros import Geometria, DanoCeramica, JuntaAtrito, CARGAS_UFV_KGF
from mef_masonry.termico import Materiais
from mef_masonry.saida import salvar_csv

ap = argparse.ArgumentParser()
ap.add_argument("--E_prisma", type=float, default=1100.0, help="módulo de corda 5-30%% do prisma, MPa")
ap.add_argument("--P30", type=float, default=CARGAS_UFV_KGF[30] * 9.80665 / 1e3, help="carga máxima a 30 °C, kN")
ap.add_argument("--folgas", type=float, nargs="+", default=[0.024, 0.032], help="desvios a testar, mm")
ap.add_argument("--sementes", type=int, nargs="+", default=[1, 2, 3], help="sementes do campo de folgas")
ap.add_argument("--h", type=float, default=12.0)
ap.add_argument("--dz", type=float, default=20.0)
a = ap.parse_args()

SAIDA = os.path.join(os.path.dirname(__file__), "..", "resultados")
os.makedirs(SAIDA, exist_ok=True)
geo, M, Pd = Geometria(), Materiais(), DanoCeramica()
m = ml.malha_prisma(geo, h=a.h / 1000, dz=a.dz / 1000, t_junta=M.junta.espessura)
P30 = a.P30 * 1e3
Pj = JuntaAtrito(semente=a.sementes[0])
print(f"Malha: {ml.resumo(m)}")
print(f"Alvos: E de corda (5-30%) = {a.E_prisma:.0f} MPa; carga máxima = {a.P30:.1f} kN; "
      f"E_contato = {Pj.E_contato / 1e6:.0f} MPa\n")

t0 = time.time()
print("1. Folga das juntas pelo módulo do prisma")
sg, pts_E = estudos.calibrar_folga(m, M, Pd, Pj, P30, a.E_prisma, folgas=[f * 1e-3 for f in a.folgas],
                                   sementes=a.sementes, geo=geo)
print(f"   desvio_folga = {sg * 1e3:.4f} mm ({time.time() - t0:.0f} s)\n")
sg = round(sg * 1e6) / 1e6          # 0,001 mm

t0 = time.time()
print("2. Resistência à compressão efetiva da cerâmica pela carga máxima a 30 °C")
Pj = replace(Pj, desvio_folga=sg)
fc, Pmax, pts_fc = estudos.calibrar_fc(m, M, Pd, Pj, P30, geo=geo)
print(f"   fc efetivo = {fc / 1e6:.2f} MPa (literatura: {Pd.fc / 1e6:.2f} MPa); "
      f"carga máxima {Pmax / 1e3:.1f} kN ({time.time() - t0:.0f} s)\n")

linhas = [["etapa 1", f"{s * 1e3:.4f}", "", f"{E:.1f}", ""] for s, E, _ in pts_E]
linhas += [["etapa 2", f"{sg * 1e3:.4f}", f"{f / 1e6:.3f}", "", f"{P / 1e3:.2f}"] for f, P in pts_fc]
linhas.append(["calibrado", f"{sg * 1e3:.4f}", f"{fc / 1e6:.3f}", f"{a.E_prisma:.0f}", f"{Pmax / 1e3:.2f}"])
salvar_csv(os.path.join(SAIDA, "calibracao_30C.csv"), linhas,
           ["etapa", "desvio_folga_mm", "fc_MPa", "E_corda_MPa", "carga_maxima_kN"])
print("Coloque em mef_masonry/parametros.py:")
print(f"  JuntaAtrito.desvio_folga = {sg:.3e}")
print(f"  CALIBRACAO_30C = {{'fc': {fc:.3e}}}")
