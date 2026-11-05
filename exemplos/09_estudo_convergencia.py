"""Estudos de convergência do MEF (Tabela S2 do artigo).

  malha_nl : prisma a 30 °C (ou --T) com o modelo não linear calibrado, em
             várias malhas; carga máxima, deslocamento no pico e E de corda.
             Os parâmetros calibrados (folga, fc efetiva) ficam fixos.
  passo_nl : mesma análise na malha padrão com incrementos de carga diferentes.
  tol_nl   : mesma análise com tolerâncias de convergência diferentes.
  termico  : ensaio aquecido (térmica + termoelástico), em várias malhas e
             passos de tempo; temperaturas no momento da carga e tração térmica.
             A rigidez da junta é recalibrada em cada malha (E_prisma = 1.100 MPa).

Uso:
  python exemplos/09_estudo_convergencia.py malha_nl --malhas 12x20 10x15 8x12
  python exemplos/09_estudo_convergencia.py passo_nl --passos 20 40 80
  python exemplos/09_estudo_convergencia.py tol_nl
  python exemplos/09_estudo_convergencia.py termico --malhas 12x20 8x12 --dts 15 30 60 --T 600
Cada estudo acrescenta linhas a resultados/convergencia_<estudo>.csv.
"""
import argparse
import csv
import os
import sys
import time
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, estudos
from mef_masonry.parametros import Geometria, Ensaio, DanoCeramica, JuntaAtrito, CARGAS_UFV_KGF, CALIBRACAO_30C
from mef_masonry.termico import Materiais

ap = argparse.ArgumentParser()
ap.add_argument("estudo", choices=["malha_nl", "passo_nl", "tol_nl", "termico"])
ap.add_argument("--malhas", nargs="+", default=["12x20"], help="h x dz em mm, por exemplo 12x20 8x12")
ap.add_argument("--passos", type=int, nargs="+", default=[40])
ap.add_argument("--dts", type=float, nargs="+", default=[30.0], help="passos de tempo da térmica, s")
ap.add_argument("--T", type=float, default=30.0, help="temperatura-alvo do ensaio, °C")
ap.add_argument("--delta_max", type=float, default=0.8, help="mm")
a = ap.parse_args()

BASE = os.path.join(os.path.dirname(__file__), "..")
SAIDA = os.path.join(BASE, "resultados")
os.makedirs(SAIDA, exist_ok=True)
arq = os.path.join(SAIDA, f"convergencia_{a.estudo}.csv")
geo = Geometria()
P30 = CARGAS_UFV_KGF[30] * 9.80665


def gravar(linha: dict):
    novo = not os.path.exists(arq)
    with open(arq, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(linha))
        if novo:
            w.writeheader()
        w.writerow(linha)
    print("  ->", linha, flush=True)


def malha(txt):
    h, dz = (float(v) / 1000 for v in txt.lower().split("x"))
    return ml.malha_prisma(geo, h=h, dz=dz, t_junta=Materiais().junta.espessura), h, dz


def ensaio_nl(m, n_passos=40, tol_u=1e-4, tol_d=5e-3):
    M = Materiais()
    Pd = DanoCeramica(fc=CALIBRACAO_30C["fc"])
    Pj = JuntaAtrito()
    campos = estudos.campos_aquecimento(m, M, Ensaio(T_alvo=a.T)) if a.T > M.ceramica.T_ref + 1 else None
    t0 = time.time()
    mod, (d, P) = estudos.ensaio_nao_linear(m, M, Pd, Pj, campos=campos, delta_max=a.delta_max / 1000,
                                            n_passos=n_passos, geo=geo, verboso=True, tol_u=tol_u, tol_d=tol_d)
    i = int(np.argmax(P))
    hist = [p for p in mod.historico if p.fase == "carga"]
    return {"T_C": a.T, "nos": m.n_nos, "gl_livres": int(np.count_nonzero(~mod.fixos_base)),
            "elementos_ceramica": len(mod.ic), "pico_kN": round(P[i] / 1e3, 2), "delta_pico_mm": round(d[i] * 1e3, 3),
            "E_corda_MPa": round(estudos.modulo_corda(mod, P30, geo), 0),
            "contato_no_pico_pct": round(100 * hist[i - 1].contato, 1) if i >= 1 else "",
            "incrementos": len(hist), "sem_convergir": sum(1 for p in hist if not p.convergiu),
            "tempo_min": round((time.time() - t0) / 60, 1)}


if a.estudo == "malha_nl":
    for txt in a.malhas:
        m, h, dz = malha(txt)
        print(f"Malha {txt}: {ml.resumo(m)}", flush=True)
        gravar({"malha_mm": txt, **ensaio_nl(m)})

elif a.estudo == "passo_nl":
    m, h, dz = malha(a.malhas[0])
    for n in a.passos:
        print(f"{n} incrementos até {a.delta_max} mm", flush=True)
        gravar({"malha_mm": a.malhas[0], "incremento_mm": round(a.delta_max / n, 4), **ensaio_nl(m, n_passos=n)})

elif a.estudo == "tol_nl":
    m, h, dz = malha(a.malhas[0])
    for tu, td in ((1e-4, 5e-3), (2e-5, 1e-3)):
        print(f"tol_u = {tu}, tol_d = {td}", flush=True)
        gravar({"malha_mm": a.malhas[0], "tol_u": tu, "tol_d": td, **ensaio_nl(m, tol_u=tu, tol_d=td)})

elif a.estudo == "termico":
    for txt in a.malhas:
        m, h, dz = malha(txt)
        M = Materiais()
        Ej, _ = estudos.calibrar_junta(m, M, 1100.0, geo, verboso=False)
        for dt in a.dts:
            t0 = time.time()
            ind = estudos.simular_ensaio_3d(m, M, Ensaio(T_alvo=a.T), geo, dt=dt)["indicadores"]
            gravar({"malha_mm": txt, "dt_s": dt, "nos": m.n_nos, "E_junta_MPa": round(Ej / 1e6, 1),
                    "T_superficie_C": round(ind["T_superficie"], 1), "T_media_terco_C": round(ind["T_media_ceramica_terco"], 1),
                    "T_septo_C": round(ind["T_septo"], 1),
                    "tracao_max_termica_MPa": round(ind["tracao_max_termico_MPa"], 3),
                    "tracao_p99_termica_MPa": round(ind["tracao_p99_termico_MPa"], 3),
                    "tempo_min": round((time.time() - t0) / 60, 1)})
