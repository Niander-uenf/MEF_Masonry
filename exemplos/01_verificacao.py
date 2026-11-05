"""Verificação do código contra soluções analíticas.

Rode antes de qualquer estudo:  python exemplos/01_verificacao.py
Gera resultados/verificacao.csv. Todos os casos devem aparecer como OK.
Outros elementos 3D:  python exemplos/01_verificacao.py --elemento H8   (ou W15)
-> resultados/verificacao_H8.csv
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, termico as tm, mecanico as mc
from mef_masonry.malha import Camada, CERAMICA, JUNTA
from mef_masonry.termico import Materiais, Contorno
from mef_masonry.parametros import Geometria
from mef_masonry import curvas
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


def registrar(caso, grandeza, mef, analitico, erro, tolerancia, unidade):
    ok = abs(erro) <= tolerancia
    linhas.append([caso, grandeza, f"{mef:.6g}", f"{analitico:.6g}", f"{erro:.3g}", f"{tolerancia:.3g}",
                   unidade, "OK" if ok else "FALHOU"])
    print(f"  {grandeza:<42s} MEF {mef:12.5g}  analítico {analitico:12.5g}  erro {erro:9.3g} {unidade}"
          f"  {'OK' if ok else 'FALHOU'}")


def placa_1d_centro(Fo, termos=200):
    n = np.arange(termos)
    return np.sum(4 * (-1) ** n / ((2 * n + 1) * np.pi) * np.exp(-((2 * n + 1) * np.pi / 2) ** 2 * Fo))


# ---------------------------------------------------------------------------
print("V1. Condução transiente numa placa com temperatura superficial imposta")
M = Materiais()
k, rc = M.ceramica.k, M.ceramica.rho * M.ceramica.c
esp = 0.14; l = esp / 2; a = k / rc
Ti, Ts = 30.0, 600.0
Fos = [0.05, 0.1, 0.2, 0.5, 1.0]
tempos = [F * l ** 2 / a for F in Fos]
for nome, tipo in ((f"3D ({EL})", "3d"), ("2D (triângulos)", "2d")):
    for dt in (60.0, 15.0):
        m2 = ml.malha_retangulo(esp, 0.01, 0.0035, hy=0.005)
        if tipo == "3d":
            m3 = ml.extrudar(m2, [Camada(0.0, 0.01, "uniforme", material=CERAMICA)], dz=0.01)
            area = tm.areas_laterais_3d(m3, faces=("x0", "x1"))
            ctr = [Contorno(area, curvas.constante(Ts), 1e8, 0.0)]
            centro = np.nonzero(np.isclose(m3.xyz[:, 0], l))[0]
            r = tm.resolver_3d(m3, M, ctr, max(tempos), dt=dt, T_inicial=Ti, sondas={"c": centro},
                               registrar_cada=dt, instantes_campo=tempos)
            T_c = [r.campos[t][centro].mean() for t in tempos]
        else:
            area = tm.comprimentos_contorno_2d(m2, faces=("x0", "x1"))
            ctr = [Contorno(area, curvas.constante(Ts), 1e8, 0.0)]
            centro = np.nonzero(np.isclose(m2.xy[:, 0], l))[0]
            mat = np.full(len(m2.tri), CERAMICA, dtype=np.int8)
            r = tm.resolver_2d(m2, M, ctr, max(tempos), mat=mat, dt=dt, T_inicial=Ti,
                               registrar_cada=dt, instantes_campo=tempos)
            T_c = [r.campos[t][centro].mean() for t in tempos]
        print(f" {nome}, dt = {dt:.0f} s")
        for F, Tc in zip(Fos, T_c):
            Tan = Ts + (Ti - Ts) * placa_1d_centro(F)
            tol = 6.0 if dt == 60.0 else 2.0
            registrar(f"V1 {nome} dt={dt:.0f}s", f"T no centro, Fo = {F}", Tc, Tan, Tc - Tan, tol, "°C")

# ---------------------------------------------------------------------------
print("\nV2. Compressão uniaxial de bloco maciço (3D)")
M = Materiais()
E = M.ceramica.E_material; nu = M.ceramica.nu
B, H, P = 0.10, 0.20, 10e3
m2 = ml.malha_retangulo(B, B, 0.02)
m3 = ml.extrudar(m2, [Camada(0.0, H, "uniforme", material=CERAMICA)], dz=0.04)
sig = P / B ** 2
for modo, rotulo in (("rigidas", "prato rígido"), ("aco", "pressão no topo")):
    R = mc.resolver_3d(m3, M, carga_N=P, incluir_termico=False, placas=modo)["carga"]
    uz_top = mc.media_uz_plano(m3, R.u, H)
    registrar(f"V2 {rotulo}", "deslocamento do topo (mm)", uz_top * 1e3, -sig * H / E * 1e3,
              (uz_top + sig * H / E) / (sig * H / E) * 100, 0.1, "%")
    dux = R.u[np.isclose(m3.xyz[:, 0], B), 0].mean() - R.u[np.isclose(m3.xyz[:, 0], 0), 0].mean()
    registrar(f"V2 {rotulo}", "expansão lateral (mm)", dux * 1e3, nu * sig / E * B * 1e3,
              (dux - nu * sig / E * B) / (nu * sig / E * B) * 100, 0.1, "%")
    registrar(f"V2 {rotulo}", "tensão szz média (MPa)", np.nanmean(R.sig[:, 2]) / 1e6, -sig / 1e6,
              (np.nanmean(R.sig[:, 2]) + sig) / sig * 100, 0.1, "%")

# ---------------------------------------------------------------------------
print("\nV3. Dilatação térmica livre uniforme (3D)")
dT = 500.0
T = np.full(m3.n_nos, M.ceramica.T_ref + dT)
R = mc.resolver_3d(m3, M, T_nos=T)["termico"]
uz_top = mc.media_uz_plano(m3, R.u, H)
esperado = M.ceramica.alfa * dT * H
registrar("V3", "alongamento livre (mm)", uz_top * 1e3, esperado * 1e3,
          (uz_top - esperado) / esperado * 100, 0.1, "%")
registrar("V3", "maior |tensão| (MPa)", np.nanmax(np.abs(R.sig)) / 1e6, 0.0,
          np.nanmax(np.abs(R.sig)) / 1e6, 1e-3, "MPa")

# ---------------------------------------------------------------------------
print("\nV4. Coluna em camadas: cerâmica + juntas finas, nu = 0 (associação em série)")
M = Materiais()
M.ceramica.nu = 0.0; M.junta.nu = 0.0; M.junta.E = 100e6
cam, z = [], 0.0
for i in range(5):
    cam.append(Camada(z, z + 0.058, "uniforme", material=CERAMICA)); z += 0.058
    if i < 4:
        cam.append(Camada(z, z + 0.002, "uniforme", material=JUNTA)); z += 0.002
m3 = ml.extrudar(m2, cam, dz=0.02)
R = mc.resolver_3d(m3, M, carga_N=P, incluir_termico=False)["carga"]
flex = 5 * 0.058 / M.ceramica.E_material + 4 * 0.002 / M.junta.E
uz_top = mc.media_uz_plano(m3, R.u, z)
registrar("V4", "encurtamento total (mm)", -uz_top * 1e3, sig * flex * 1e3,
          (-uz_top - sig * flex) / (sig * flex) * 100, 0.1, "%")
Eeq_mef = sig * z / -uz_top; Eeq_an = z / flex
registrar("V4", "módulo equivalente (MPa)", Eeq_mef / 1e6, Eeq_an / 1e6,
          (Eeq_mef - Eeq_an) / Eeq_an * 100, 0.1, "%")

# ---------------------------------------------------------------------------
print("\nV5. Placa livre com perfil parabólico de temperatura na espessura (3D)")
M = Materiais()
E = M.ceramica.E_material; nu = M.ceramica.nu; alfa = M.ceramica.alfa
Lx, t, Hz, dT0 = 0.24, 0.04, 0.24, 300.0
m2 = ml.malha_retangulo(Lx, t, 0.01, hy=0.0025)
# a cunha quadrática usa camadas de 2 cm (mesmo número de nós na vertical; a
# fatoração direta da W15 com camadas de 1 cm passa da memória sem o pypardiso)
m3 = ml.extrudar(m2, [Camada(0.0, Hz, "uniforme", material=CERAMICA)], dz=0.02 if EL == "W15" else 0.01)
y = m3.xyz[:, 1]
T = M.ceramica.T_ref + dT0 * (2 * y / t - 1) ** 2
t0 = time.time()
R = mc.resolver_3d(m3, M, T_nos=T)["termico"]
cen = m3.xyz[m3.cun].mean(1)
perto = (np.abs(cen[:, 0] - Lx / 2) < 0.011) & (np.abs(cen[:, 2] - Hz / 2) < 0.011)
yc = cen[perto, 1]
an = E * alfa * (dT0 / 3 - dT0 * (2 * yc / t - 1) ** 2) / (1 - nu)
ref = E * alfa * dT0
for comp, nome in ((0, "sxx"), (2, "szz")):
    err = np.max(np.abs(R.sig[perto, comp] - an)) / ref * 100
    registrar("V5", f"{nome}: maior erro / (E alfa dT0)", np.max(np.abs(R.sig[perto, comp])) / 1e6,
              np.max(np.abs(an)) / 1e6, err, 3.0, "%")
print(f"  ({len(m3.cun)} elementos {EL}, {time.time() - t0:.1f} s)")

# ---------------------------------------------------------------------------
print("\nV6. O mesmo perfil parabólico na seção 2D (estado plano generalizado)")
m2 = ml.malha_retangulo(Lx, t, 0.0025)
T = M.ceramica.T_ref + dT0 * (2 * m2.xy[:, 1] / t - 1) ** 2
mat = np.full(len(m2.tri), CERAMICA, dtype=np.int8)
s = mc.resolver_secao_gps(m2, M, T, mat=mat)
cen = m2.xy[m2.tri].mean(1)
perto = np.abs(cen[:, 0] - Lx / 2) < 0.006
an = E * alfa * (dT0 / 3 - dT0 * (2 * cen[perto, 1] / t - 1) ** 2) / (1 - nu)
for comp in ("sx", "sz"):
    err = np.max(np.abs(s[comp][perto] - an)) / ref * 100
    registrar("V6", f"{comp}: maior erro / (E alfa dT0)", np.max(np.abs(s[comp][perto])) / 1e6,
              np.max(np.abs(an)) / 1e6, err, 3.0, "%")

# ---------------------------------------------------------------------------
print("\nV7. Seção com furos: térmica 2D x 3D extrudado (devem coincidir)")
M = Materiais()
geo = Geometria()
f, t_fim, _ = curvas.rampa_patamar(T_alvo=600.0)
m2 = ml.malha_secao(geo, h=0.006, simetria=True)
ctr2 = [Contorno(tm.comprimentos_contorno_2d(m2), f, 10.0, 0.8)]
r2 = tm.resolver_2d(m2, M, ctr2, t_fim, dt=30.0)
m3 = ml.extrudar(m2, [Camada(0.0, 0.02, "unidade", fiada=0)], dz=0.01, geo=geo)
m3.mat[m3.mat == JUNTA] = CERAMICA
ctr3 = [Contorno(tm.areas_laterais_3d(m3), f, 10.0, 0.8)]
r3 = tm.resolver_3d(m3, M, ctr3, t_fim, dt=30.0)
T3 = r3.T_final[:m2.n_nos]          # nós de canto do plano z = 0 (vêm primeiro em todas as malhas)
for nome, v2, v3 in (("mínima da seção", r2.T_final.min(), T3.min()),
                     ("média nodal", r2.T_final.mean(), T3.mean())):
    # H8 e W15: a seção 3D (quadriláteros, ou triângulos quadráticos) não é a malha T3
    # do modelo 2D; as duas convergem para a mesma solução, mas não coincidem nó a nó
    registrar("V7", f"T {nome} no fim (°C)", v3, v2, v3 - v2, 0.5 if EL == "W6" else 1.0, "°C")

salvar_csv(os.path.join(SAIDA, "verificacao.csv" if EL == "W6" else f"verificacao_{EL}.csv"), linhas,
           ["caso", "grandeza", "MEF", "analítico/referência", "erro", "tolerância", "unidade", "situação"])
falhas = sum(1 for l in linhas if l[-1] != "OK")
print(f"\n{len(linhas) - falhas} de {len(linhas)} verificações OK.")
sys.exit(1 if falhas else 0)
