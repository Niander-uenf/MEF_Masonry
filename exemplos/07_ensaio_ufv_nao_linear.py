"""Ensaio da UFV com o modelo não linear: dano na cerâmica e juntas secas com
atrito. Aquece o prisma campo a campo (rampa de 5 °C/min e patamar) e depois
comprime com o prato de topo em deslocamento controlado, até depois do pico.

As juntas de leito têm folgas iniciais (irregularidade das faces), com o
desvio-padrão calibrado pelo módulo do prisma a 30 °C (script 08; padrão em
JuntaAtrito). Com --folga 0, o contato é perfeito e a rigidez das juntas é
recalibrada pelo E do prisma.

Resistências: as da literatura (DanoCeramica) por padrão; com --calibrado, a
resistência à compressão efetiva que reproduz a carga máxima a 30 °C
(CALIBRACAO_30C, script 08).

Uso: python exemplos/07_ensaio_ufv_nao_linear.py [--temperaturas 30 600 800 1000]
     [--calibrado] [--fc 5.76] [--ft 0.66] [--folga 0.028] [--semente 1]
     [--h 12] [--dz 20] [--delta_max 0.8] [--passos 40]
Tempo: cerca de 10 min por temperatura com a malha padrão.
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, estudos
from mef_masonry.parametros import Geometria, Ensaio, DanoCeramica, JuntaAtrito, CARGAS_UFV_KGF, CALIBRACAO_30C
from mef_masonry.termico import Materiais
from mef_masonry.naolinear import ModeloNaoLinear
from mef_masonry.saida import salvar_csv, salvar_vtk

ap = argparse.ArgumentParser()
ap.add_argument("--temperaturas", type=float, nargs="+", default=[30, 600, 800, 1000])
ap.add_argument("--h", type=float, default=12.0)
ap.add_argument("--dz", type=float, default=20.0)
ap.add_argument("--ft", type=float, default=DanoCeramica().ft / 1e6, help="resistência à tração da cerâmica, MPa")
ap.add_argument("--fc", type=float, default=None, help="resistência à compressão da cerâmica, MPa")
ap.add_argument("--calibrado", action="store_true", help="usa a fc efetiva calibrada a 30 °C")
ap.add_argument("--E_prisma", type=float, default=1100.0, help="módulo-alvo do prisma a 30 °C, MPa")
ap.add_argument("--delta_max", type=float, default=0.8, help="deslocamento máximo do prato, mm")
ap.add_argument("--passos", type=int, default=40)
ap.add_argument("--folga", type=float, default=JuntaAtrito().desvio_folga * 1e3,
                help="desvio-padrão das folgas das juntas de leito, mm (0 = contato perfeito)")
ap.add_argument("--E_contato", type=float, default=(JuntaAtrito().E_contato or 0) / 1e6,
                help="E das juntas onde há contato, MPa")
ap.add_argument("--semente", type=int, default=JuntaAtrito().semente, help="semente do campo de folgas")
ap.add_argument("--rotulo", default="")
a = ap.parse_args()

BASE = os.path.join(os.path.dirname(__file__), "..")
SAIDA = os.path.join(BASE, "resultados")
os.makedirs(os.path.join(SAIDA, "vtk"), exist_ok=True)
geo = Geometria()
M = Materiais()
fc = a.fc * 1e6 if a.fc else (CALIBRACAO_30C["fc"] if a.calibrado else DanoCeramica().fc)
Pd = DanoCeramica(ft=a.ft * 1e6, fc=fc)
Pj = JuntaAtrito(desvio_folga=a.folga * 1e-3, semente=a.semente,
                 E_contato=(a.E_contato * 1e6 or None) if a.folga > 0 else None)
P30 = CARGAS_UFV_KGF[30] * 9.80665
m = ml.malha_prisma(geo, h=a.h / 1000, dz=a.dz / 1000, t_junta=M.junta.espessura)
print("Malha:", ml.resumo(m))
print(f"Cerâmica: ft = {Pd.ft / 1e6:.2f} MPa, fc = {Pd.fc / 1e6:.2f} MPa, G_f = {Pd.G_f:.1f} N/m, "
      f"G_c = {Pd.G_c:.0f} N/m; juntas: mu = {Pj.mu}")

t0 = time.time()
if Pj.desvio_folga <= 0:
    mod0 = ModeloNaoLinear(m, M, Pd, Pj, verboso=False)
    print("Contato perfeito. Calibrando a rigidez das juntas para E_prisma =", a.E_prisma, "MPa")
    Ej = estudos.calibrar_junta_nao_linear(mod0, a.E_prisma, geo)
    M.junta.E = Ej; M.capeamento.E = Ej
    del mod0
    print(f"  E_junta = {Ej / 1e6:.1f} MPa ({time.time() - t0:.0f} s)\n")
else:
    print(f"Juntas com folga: desvio {Pj.desvio_folga * 1e3:.3f} mm, contato inicial {100 * Pj.contato_inicial:.0f}%, "
          f"correlação {Pj.comprimento_corr * 1e3:.0f} mm, E_contato {Pj.E_contato / 1e6:.0f} MPa, semente {Pj.semente}\n")

resumo, curvas_pd = [], {}
for T_alvo in a.temperaturas:
    t0 = time.time()
    print(f"Ensaio de {T_alvo:.0f} °C")
    campos = estudos.campos_aquecimento(m, M, Ensaio(T_alvo=T_alvo)) if T_alvo > M.ceramica.T_ref + 1 else None
    mod, (deltas, cargas) = estudos.ensaio_nao_linear(m, M, Pd, Pj, campos=campos, delta_max=a.delta_max / 1000,
                                                      n_passos=a.passos, geo=geo, verboso=True)
    vc = mod.vol[mod.ic]
    frac_aq = float(vc[mod.d_aquecimento > 0.5].sum() / vc.sum())
    d_aq = np.full(len(m.cun), np.nan); d_aq[mod.id_glob[mod.ic]] = mod.d_aquecimento
    i = int(np.argmax(cargas))
    P_exp = CARGAS_UFV_KGF.get(int(T_alvo), np.nan) * 9.80665
    resumo.append([T_alvo, f"{cargas[i] / 1e3:.1f}", f"{P_exp / 1e3:.1f}", f"{(cargas[i] / P_exp - 1) * 100:+.0f}%",
                   f"{deltas[i] * 1e3:.3f}", f"{100 * frac_aq:.2f}", f"{estudos.modulo_corda(mod, P30, geo):.0f}",
                   sum(1 for p in mod.historico if p.fase == "carga" and not p.convergiu), f"{time.time() - t0:.0f}"])
    curvas_pd[T_alvo] = (deltas, cargas)
    hist = [p for p in mod.historico if p.fase == "carga"]
    salvar_csv(os.path.join(SAIDA, f"nl_carga_deslocamento_{int(T_alvo)}C{a.rotulo}.csv"),
               [[f"{p.rotulo * 1e3:.4f}", f"{p.carga_N / 1e3:.3f}", f"{-p.eps_ext * 1e3:.4f}", f"{100 * p.contato:.1f}",
                 f"{100 * p.vol_danificado:.2f}"] for p in hist],
               ["delta_mm", "carga_kN", "deformacao_terco_central_por_mil", "contato_juntas_pct", "ceramica_d_maior_0.5_pct"])
    folga_el = np.full(len(m.cun), np.nan); folga_el[mod.id_glob[mod.ij]] = mod.folga * 1e3
    salvar_vtk(os.path.join(SAIDA, "vtk", f"nl_prisma_{int(T_alvo)}C{a.rotulo}.vtk"), m,
               dados_nos={"T_C": mod.T, "u_mm": mod.u.reshape(-1, 3) * 1e3},
               dados_elem={"material": m.mat.astype(float), "dano_final": mod.dano_global(),
                           "dano_apos_aquecimento": d_aq, "folga_inicial_mm": folga_el})
    print(f"  pico {cargas[i] / 1e3:.1f} kN com delta = {deltas[i] * 1e3:.3f} mm; ensaio {P_exp / 1e3:.1f} kN; "
          f"cerâmica com d > 0,5 após o aquecimento: {100 * frac_aq:.2f}%  ({time.time() - t0:.0f} s)\n")

cab = ["T_alvo_C", "pico_modelo_kN", "pico_ensaio_kN", "diferenca", "delta_no_pico_mm",
       "ceramica_danificada_apos_aquecimento_pct", "E_corda_5_30_MPa", "incrementos_sem_convergir", "tempo_s"]
salvar_csv(os.path.join(SAIDA, f"nl_resumo{a.rotulo}.csv"), resumo, cab)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.5, 4.3))
cores = plt.cm.Blues(np.linspace(0.35, 0.95, len(a.temperaturas)))
for c, T in zip(cores, a.temperaturas):
    d, P = curvas_pd[T]
    ax.plot(d * 1e3, P / 1e3, "-", color=c, lw=1.8, label=f"modelo, {T:.0f} °C")
    Pe = CARGAS_UFV_KGF.get(int(T))
    if Pe:
        ax.axhline(Pe * 9.80665 / 1e3, color=c, ls=":", lw=1)
ax.set_xlabel("deslocamento do prato após o aquecimento (mm)"); ax.set_ylabel("carga (kN)")
ax.set_title("Prisma BCEPQ: modelo não linear (linhas) e carga máxima do ensaio (pontilhado)", fontsize=9)
ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(SAIDA, f"nl_carga_deslocamento{a.rotulo}.png"), dpi=130)
print("Resumo:")
for l in resumo:
    print("  ", dict(zip(cab, l)))
