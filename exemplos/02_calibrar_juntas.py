"""Calibra a rigidez das juntas secas para o prisma reproduzir o módulo medido,
e compara experimental x numérico x fórmulas (série e NBR 16868).

Uso: python exemplos/02_calibrar_juntas.py [--E_alvo 1100] [--h 8] [--dz 12]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, estudos
from mef_masonry.parametros import Geometria, E_PRISMA_SECO_MPA, CARGAS_UFV_KGF
from mef_masonry.termico import Materiais
from mef_masonry.saida import salvar_csv

ap = argparse.ArgumentParser()
ap.add_argument("--E_alvo", type=float, default=1100.0, help="módulo do prisma seco, MPa (BPQ: 1.001 a 1.216)")
ap.add_argument("--h", type=float, default=8.0, help="tamanho de elemento na seção, mm")
ap.add_argument("--dz", type=float, default=12.0, help="passo vertical, mm")
a = ap.parse_args()

SAIDA = os.path.join(os.path.dirname(__file__), "..", "resultados")
os.makedirs(SAIDA, exist_ok=True)
geo = Geometria()
M = Materiais()
t0 = time.time()
m = ml.malha_prisma(geo, h=a.h / 1000, dz=a.dz / 1000, t_junta=M.junta.espessura)
print("Malha:", ml.resumo(m))

print(f"\nCalibrando E da junta para E_prisma = {a.E_alvo:.0f} MPa (terço central, área bruta)")
Ej, hist = estudos.calibrar_junta(m, M, a.E_alvo, geo)
r = estudos.modulo_prisma(m, M, geo=geo)
base = (r["z2"] - r["z1"]) * 1000
Ej_serie = estudos.estimativa_serie(M, a.E_alvo, base, 2)

fpk = CARGAS_UFV_KGF[30] * 9.80665 / (geo.comprimento * geo.largura) / 1e6
linhas = [
    ["E do bloco (área bruta)", "experimental (BPQ, n = 6)", f"{M.ceramica.E_bloco_bruto / 1e6:.0f}"],
    ["E do material cerâmico (área líquida)", "derivado: E_bloco / 0,717", f"{M.ceramica.E_material / 1e6:.0f}"],
]
for k, v in E_PRISMA_SECO_MPA.items():
    linhas.append([f"E do prisma seco {k}", "experimental (BPQ, n = 4)", f"{v:.0f}"])
linhas += [
    ["E do prisma = 600 fpk", f"NBR 16868-1, fpk = {fpk:.2f} MPa (UFV, 30 °C)", f"{600 * fpk:.0f}"],
    ["E da junta: associação em série", f"teórico, base {base:.0f} mm com 2 juntas", f"{Ej_serie / 1e6:.1f}"],
    ["E da junta: calibrado no MEF", "numérico 3D", f"{Ej / 1e6:.1f}"],
    ["E do prisma no terço central", "numérico 3D, junta calibrada", f"{r['E_extensometro']:.0f}"],
    ["E do prisma entre placas", "numérico 3D, junta calibrada", f"{r['E_placas']:.0f}"],
    ["Rigidez normal da junta k_n = E_j/t", "numérico 3D (N/mm³)", f"{Ej / M.junta.espessura / 1e9:.1f}"],
]
print("\nResumo (MPa, exceto onde indicado):")
for l in linhas:
    print(f"  {l[0]:<40s} {l[1]:<45s} {l[2]:>8s}")
salvar_csv(os.path.join(SAIDA, "calibracao_juntas.csv"), linhas, ["grandeza", "origem", "valor"])
salvar_csv(os.path.join(SAIDA, "calibracao_iteracoes.csv"), [[f"{x:.3f}", f"{y:.2f}"] for x, y in hist],
           ["E_junta_MPa", "E_prisma_MPa"])
print(f"\nUse E_junta = {Ej / 1e6:.1f} MPa nos demais exemplos (parametros.JuntaSeca.E). "
      f"Tempo: {time.time() - t0:.0f} s")
