"""Seção 7.3 do artigo: teste barato da hipótese de transporte de umidade nas paredes de
Oliveira et al. (2021). O modelo térmico não transporta água; aqui a água é redistribuída
à mão antes do ensaio, como se já tivesse migrado para a casca fria: o pico de calor
específico da EN 1992-1-2 (100 a 115 °C) sai da casca exposta e vai para a casca não exposta,
com a massa de água conservada.

  A: a água da casca exposta passa para a casca não exposta (exposta 0 %, septos 1,5 %, fria 3 %)
  B: toda a água da seção vai para a casca não exposta (exposta e septos 0 %, fria ~5 %)

Pico de c(T): EN 1992-1-2 (0 % -> 900, 1,5 % -> 1470, 3 % -> 2020 J/kg.K), extrapolado
linearmente acima de 3 %. Base: projeto "oliveira" (k limite inferior, 1,5 %, emissividade 0,7)."""
import os, sys, json, copy
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oliveira import projeto, AQUI
from mef_masonry import projeto as pj, termico as tm, malha as ml, curvas
from mef_masonry.malha import AR, CERAMICA, ARGAMASSA, GRAUTE
from mef_masonry.termico import Contorno


FORMA = "en"   # "en": pico da EN 1992-1-2; "latente": todo o calor de vaporização entre 100 e 110 °C


def tabela_c(u):
    global FORMA
    if FORMA == "latente":
        P = 2.257e6 * u / 100 / 10.0          # J/kg.K a mais entre 100 e 110 °C
        return [(20, 900), (100, 900), (100.01, 900 + P), (110, 910 + P), (110.01, 910), (200, 1000), (400, 1100), (1200, 1100)]
    pico = 900 + u / 1.5 * 570 if u <= 1.5 else 1470 + (u - 1.5) / 1.5 * 550
    return [(20, 900), (100, 900), (100.01, pico), (115, pico), (200, 1000), (400, 1100), (1200, 1100)]


def roda(u_exp, u_septo, u_fria, nome):
    p = projeto("oliveira")
    geo = pj.geometria(p)
    M = pj.materiais(p, geo)
    en, ma = p["ensaio"], p["malha"]
    m2 = ml.malha_secao(geo, h=ma["h"] / 1e3, simetria=False, faixas=())
    mat = tm.material_secao(m2, geo).copy()
    yc = m2.xy[m2.tri][:, :, 1].mean(1)
    W, tc = geo.largura, p["geometria"]["parede_externa"] / 1e3
    sol = mat == CERAMICA
    exp, fria = sol & (yc < tc), sol & (yc > W - tc)
    mat[exp], mat[fria] = ARGAMASSA, GRAUTE          # só rótulos: as três zonas têm o mesmo concreto
    for alvo in (M.argamassa, M.graute):
        alvo.k, alvo.rho, alvo.c = M.ceramica.k, M.ceramica.rho, M.ceramica.c
        alvo.tabela_k, alvo.tabela_rho = M.ceramica.tabela_k, M.ceramica.tabela_rho
    M.argamassa.tabela_c = tabela_c(u_exp)
    M.ceramica.tabela_c = tabela_c(u_septo)
    M.graute.tabela_c = tabela_c(u_fria)
    d1, d2 = m2.xy[m2.tri[:, 1]] - m2.xy[m2.tri[:, 0]], m2.xy[m2.tri[:, 2]] - m2.xy[m2.tri[:, 0]]
    area = np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]) / 2
    agua = sum(u * area[s].sum() for u, s in ((u_exp, exp), (u_septo, sol & ~exp & ~fria), (u_fria, fria)))
    f, tf, _ = pj.programa_termico(p)
    nexp = tm.comprimentos_contorno_2d(m2, faces=("y1",))
    ctr = [Contorno(tm.comprimentos_contorno_2d(m2, faces=("y0",)), f, en["h_conveccao"], en["emissividade"], "exposta"),
           Contorno(nexp, curvas.constante(en["T_inicial"]), 4.0, en["emissividade"], "não exposta")]
    face = np.nonzero(nexp > 0)[0]; w = nexp[face] / nexp[face].sum(); Ti = en["T_inicial"]
    r = tm.resolver_2d(m2, M, ctr, tf, dt=ma["dt"], T_inicial=Ti, mat=mat, registrar_cada=60.0,
                       sondas={"face": face}, parar_se=lambda t, T: (T[face] @ w - Ti >= 140) or (T[face].max() - Ti >= 180))
    I = r.parou_em / 60 if r.parou_em else None
    print(f"{nome:40s} u = {u_exp:.2f} / {u_septo:.2f} / {u_fria:.2f} %   água relativa = {agua:.4e}   isolamento = {I} min", flush=True)
    return {"u_exposta": u_exp, "u_septos": u_septo, "u_fria": u_fria, "agua_rel": agua, "I_min": I,
            "areas_mm2": {"exposta": area[exp].sum() * 1e6, "septos": area[sol & ~exp & ~fria].sum() * 1e6, "fria": area[fria].sum() * 1e6}}


if __name__ == "__main__":
    out = {}
    for forma in ("en", "latente"):
        FORMA = forma
        print("pico de c(T):", forma)
        o = out[forma] = {"referencia": roda(1.5, 1.5, 1.5, "referência (1,5 % em toda a seção)")}
        a = o["referencia"]["areas_mm2"]
        o["A"] = roda(0.0, 1.5, 1.5 * (a["exposta"] + a["fria"]) / a["fria"], "A: água da casca exposta na fria")
        o["B"] = roda(0.0, 0.0, 1.5 * (a["exposta"] + a["septos"] + a["fria"]) / a["fria"], "B: toda a água na casca fria")
        for u in (10, 20):
            o[f"fria_{u}"] = roda(0.0, 0.0, u, f"casca fria com {u} % (água a mais)")
    json.dump(out, open(os.path.join(AQUI, "oliveira_umidade_migrada.json"), "w"), indent=1)
