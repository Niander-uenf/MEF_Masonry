"""Verificação do modelo não linear (dano na cerâmica e juntas secas).

  N1. Tração uniaxial com uma camada enfraquecida: resistência de pico e
      energia dissipada = G_f x área, para dois tamanhos de elemento
      (objetividade da malha).
  N2. Compressão uniaxial homogênea: tensão de pico = fc.
  N3. Junta seca: não transmite tração; em compressão, rigidez igual à do
      modelo linear; ao cisalhar sob compressão, a força satura em mu x N.
  N4. Junta com folga inicial g: quase sem força até fechar a folga; depois,
      força = rigidez da série x (deslocamento - g).
  N5. Fechamento de fissura: barra fissurada em tração e depois comprimida
      recupera a rigidez intacta em compressão.

Uso: python exemplos/06_verificacao_nao_linear.py [--elemento H8|W15]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, mecanico as mc
from mef_masonry.malha import Camada, CERAMICA
from mef_masonry.termico import Materiais
from mef_masonry.naolinear import ModeloNaoLinear
from mef_masonry.parametros import DanoCeramica, JuntaAtrito
from mef_masonry.saida import salvar_csv

ap = argparse.ArgumentParser()
ap.add_argument("--elemento", default="W6", choices=["W6", "H8", "W15"])
EL = ap.parse_args().elemento
_extrudar = ml.extrudar


def extrudar(*a, **k):
    return _extrudar(*a, tipo=EL, **k)


ml.extrudar = extrudar
print(f"Elemento 3D: {EL}\n")
SAIDA = os.path.join(os.path.dirname(__file__), "..", "resultados")
os.makedirs(SAIDA, exist_ok=True)
linhas = []


def registrar(caso, grandeza, mef, ref, erro, tol, un):
    ok = abs(erro) <= tol
    linhas.append([caso, grandeza, f"{mef:.5g}", f"{ref:.5g}", f"{erro:.3g}", f"{tol:.3g}", un, "OK" if ok else "FALHOU"])
    print(f"  {grandeza:<46s} MEF {mef:11.5g}  referência {ref:11.5g}  erro {erro:8.3g} {un}  {'OK' if ok else 'FALHOU'}")


M = Materiais()
Pd = DanoCeramica()
Pj0 = JuntaAtrito(desvio_folga=0.0, E_contato=None)   # contato perfeito, E da junta de Materiais
B = 0.02
A = B * B

print("N1. Tração uniaxial: pico e energia dissipada (objetividade da malha)")
for dz in (0.010, 0.005):
    m2 = ml.malha_retangulo(B, B, 0.01)
    L = 0.10
    m3 = ml.extrudar(m2, [Camada(0.0, L, "uniforme", material=CERAMICA)], dz=dz)
    zc = m3.xyz[m3.cun].mean(1)[:, 2]
    fraco = np.abs(zc - (L / 2 - dz / 2)) < dz / 4      # uma camada no meio
    fator = np.where(fraco, 0.95, 1.0)
    mod = ModeloNaoLinear(m3, M, Pd, Pj0, apoio="prato", fator_ft_elem=fator, verboso=False)
    presc = 3 * mod.topo + 2
    ds = np.linspace(0, 0.12e-3, 121)
    Ps = []
    for dlt in ds[1:]:
        P, it, conv = mod._passo(presc, np.full(len(presc), dlt), np.zeros(mod.ndof))
        Ps.append(-P)
    Ps = np.r_[0.0, Ps]
    W = np.trapezoid(Ps, ds) if hasattr(np, "trapezoid") else np.trapz(Ps, ds)
    Gf = Pd.du_t * Pd.ft * 0.95
    registrar(f"N1 dz={dz * 1e3:.0f}mm", "tensão de pico (MPa)", Ps.max() / A / 1e6, 0.95 * Pd.ft / 1e6,
              (Ps.max() / A - 0.95 * Pd.ft) / (0.95 * Pd.ft) * 100, 2.0, "%")
    # tolerância de 10%: na banda, a camada danificada tem a contração lateral
    # restringida pelas vizinhas, o que eleva a tensão efetiva principal em até ~9%
    registrar(f"N1 dz={dz * 1e3:.0f}mm", "energia dissipada / área (N/m)", W / A, Gf, (W / A - Gf) / Gf * 100, 10.0, "%")

print("\nN5. Fechamento de fissura: a barra do N1, fissurada em tração, é comprimida")
E0 = float(M.ceramica.E_T(np.array([M.ceramica.T_ref]))[0])
k_el = E0 * A / L
dlt_c = 2e-6
Pc5 = mod._passo(presc, np.full(len(presc), -dlt_c), np.zeros(mod.ndof))[0]
registrar("N5", "rigidez em compressão após fissurar / intacta", Pc5 / (k_el * dlt_c), 1.0,
          (Pc5 / (k_el * dlt_c) - 1) * 100, 2.0, "%")

print("\nN2. Compressão uniaxial homogênea: pico = fc")
m2 = ml.malha_retangulo(B, B, 0.01)
m3 = ml.extrudar(m2, [Camada(0.0, 0.04, "uniforme", material=CERAMICA)], dz=0.01)
mod = ModeloNaoLinear(m3, M, Pd, Pj0, apoio="prato", verboso=False)
presc = 3 * mod.topo + 2
Ps = []
for dlt in np.linspace(0, 0.12e-3, 61)[1:]:
    Ps.append(mod._passo(presc, np.full(len(presc), -dlt), np.zeros(mod.ndof))[0])
registrar("N2", "tensão de pico (MPa)", max(Ps) / A / 1e6, Pd.fc / 1e6, (max(Ps) / A - Pd.fc) / Pd.fc * 100, 2.0, "%")

print("\nN3. Junta seca entre dois blocos (corpo largo e baixo, como num ensaio de cisalhamento direto)")
Bj, tb, tj = 0.10, 0.01, 0.002
cam = [Camada(0.0, tb, "uniforme", material=CERAMICA), Camada(tb, tb + tj, "junta"),
       Camada(tb + tj, 2 * tb + tj, "uniforme", material=CERAMICA)]
m3 = ml.extrudar(ml.malha_retangulo(Bj, Bj, 0.02), cam, dz=0.005)
Aj = Bj * Bj
mod = ModeloNaoLinear(m3, M, Pd, Pj0, apoio="prato", verboso=False)
presc = 3 * mod.topo + 2
d0 = 2e-7
Pc = mod._passo(presc, np.full(len(presc), -d0), np.zeros(mod.ndof))[0]
mod = ModeloNaoLinear(m3, M, Pd, Pj0, apoio="prato", verboso=False)
Pt = mod._passo(presc, np.full(len(presc), +d0), np.zeros(mod.ndof))[0]
k_an = Aj / (2 * tb / M.ceramica.E_material + tj / M.junta.E)      # série: cerâmica + camada fina
registrar("N3", "rigidez em compressão (kN/mm)", Pc / d0 / 1e6, k_an / 1e6, (Pc / d0 - k_an) / k_an * 100, 1.0, "%")
registrar("N3", "força em tração / força em compressão", abs(Pt) / Pc, 0.0, abs(Pt) / Pc, 1e-3, "-")

mod = ModeloNaoLinear(m3, M, Pd, Pj0, apoio="engaste", verboso=False)
N = 0.5e6 * Aj                     # compressão de 0,5 MPa
F = np.zeros(mod.ndof)
F[3 * mod.topo + 2] = -N * mod.area_topo[mod.topo] / mod.area_topo[mod.topo].sum()
presc = 3 * mod.topo
for s in (1e-6, 3e-6, 1e-5, 3e-5, 1e-4):
    mod._passo(presc, np.full(len(presc), s), F)
K, Fth = mod._montar(mod.d, mod.En, mod.Gs)
V = (K @ mod.u - Fth - F)[presc].sum()
mu = JuntaAtrito().mu
registrar("N3", "cisalhamento final / (mu N)", V / (mu * N), 1.0, (V / (mu * N) - 1) * 100, 2.0, "%")

print("\nN4. Junta com folga uniforme g: força quase nula até fechar; depois, rigidez da série")
g = 4e-6
mod = ModeloNaoLinear(m3, M, Pd, Pj0, apoio="prato", verboso=False)
mod.folga[mod.leito] = g; mod.eg = mod.folga / mod.dzj
mod.fechada = mod.eg <= 0; mod.En = np.where(mod.fechada, mod.Ej, mod.Pj.fator_aberta * mod.Ej)
presc = 3 * mod.topo + 2
P1 = mod._passo(presc, np.full(len(presc), -0.5 * g), np.zeros(mod.ndof))[0]
d1 = 1e-6
P2 = mod._passo(presc, np.full(len(presc), -(g + d1)), np.zeros(mod.ndof))[0]
registrar("N4", "força com 0,5 g / (k x g)", P1 / (k_an * g), 0.0, P1 / (k_an * g), 1e-3, "-")
registrar("N4", "força com g + d / (k x d)", P2 / (k_an * d1), 1.0, (P2 / (k_an * d1) - 1) * 100, 1.0, "%")

salvar_csv(os.path.join(SAIDA, "verificacao_nao_linear.csv" if EL == "W6" else f"verificacao_nao_linear_{EL}.csv"), linhas,
           ["caso", "grandeza", "MEF", "referência", "erro", "tolerância", "unidade", "situação"])
falhas = sum(1 for l in linhas if l[-1] != "OK")
print(f"\n{len(linhas) - falhas} de {len(linhas)} verificações OK.")
sys.exit(1 if falhas else 0)
