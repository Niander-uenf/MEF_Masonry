"""Seção 7.3 do artigo: de onde vem o calor a mais da seção explícita do bloco de concreto
de Oliveira et al. (2021)? Varia o modelo das cavidades e a condução pelos septos e mede o
tempo de isolamento. Base: projeto "oliveira" (k limite inferior, 1,5 %, emissividade 0,7).

A cavidade é um material com condutividade equivalente k_ar + 4 sigma T^3 D eps/(2 - eps),
a troca entre placas paralelas infinitas a uma distância D. Numa célula de dimensões finitas
só uma parte da radiação da casca quente chega direto à casca fria: o fator de forma entre as
faces opostas da célula. Escalar o termo radiativo por esse fator equivale a usar uma
emissividade efetiva eps' com eps'/(2 - eps') = F eps/(2 - eps)."""
import os, sys, json
from math import atan, log, sqrt, pi
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oliveira import projeto, AQUI
from mef_masonry import projeto as pj, termico as tm, malha as ml, curvas
from mef_masonry.malha import CERAMICA, ARGAMASSA
from mef_masonry.termico import Contorno


def F_retangulos(a, b, c):
    """Fator de forma entre dois retângulos paralelos a x b, alinhados, a uma distância c."""
    X, Y = a / c, b / c
    return 2 / (pi * X * Y) * (log(sqrt((1 + X * X) * (1 + Y * Y) / (1 + X * X + Y * Y)))
        + X * sqrt(1 + Y * Y) * atan(X / sqrt(1 + Y * Y)) + Y * sqrt(1 + X * X) * atan(Y / sqrt(1 + X * X))
        - X * atan(X) - Y * atan(Y))


def F_faixas(w, c):
    """Fator de forma 2D entre duas faixas paralelas de largura w a uma distância c."""
    return sqrt(1 + (c / w) ** 2) - c / w


def eps_efetiva(F, eps=0.9):
    r = F * eps / (2 - eps)
    return 2 * r / (1 + r)


def roda(nome, eps_furos=0.9, fator_k_septos=1.0):
    p = projeto("oliveira")
    p["ensaio"]["emissividade_furos"] = eps_furos
    geo = pj.geometria(p)
    M = pj.materiais(p, geo)
    en, ma = p["ensaio"], p["malha"]
    m2 = ml.malha_secao(geo, h=ma["h"] / 1e3, simetria=False, faixas=())
    mat = tm.material_secao(m2, geo).copy()
    if fator_k_septos != 1.0:
        yc = m2.xy[m2.tri][:, :, 1].mean(1)
        W, tc = geo.largura, p["geometria"]["parede_externa"] / 1e3
        septo = (mat == CERAMICA) & (yc > tc) & (yc < W - tc)
        mat[septo] = ARGAMASSA
        c = M.ceramica
        M.argamassa.k, M.argamassa.rho, M.argamassa.c = c.k * fator_k_septos, c.rho, c.c
        M.argamassa.tabela_k = [(T, v * fator_k_septos) for T, v in c.tabela_k]
        M.argamassa.tabela_rho, M.argamassa.tabela_c = c.tabela_rho, c.tabela_c
    f, tf, _ = pj.programa_termico(p)
    nexp = tm.comprimentos_contorno_2d(m2, faces=("y1",))
    ctr = [Contorno(tm.comprimentos_contorno_2d(m2, faces=("y0",)), f, en["h_conveccao"], en["emissividade"], "exposta"),
           Contorno(nexp, curvas.constante(en["T_inicial"]), 4.0, en["emissividade"], "não exposta")]
    face = np.nonzero(nexp > 0)[0]; w = nexp[face] / nexp[face].sum(); Ti = en["T_inicial"]
    r = tm.resolver_2d(m2, M, ctr, tf, dt=ma["dt"], T_inicial=Ti, mat=mat, registrar_cada=60.0,
                       sondas={"face": face}, parar_se=lambda t, T: (T[face] @ w - Ti >= 140) or (T[face].max() - Ti >= 180))
    I = r.parou_em / 60 if r.parou_em else None
    print(f"{nome:58s} eps_furos = {eps_furos:.3f}  k septos x {fator_k_septos:.2f}  ->  isolamento = {I} min", flush=True)
    return {"eps_furos": eps_furos, "fator_k_septos": fator_k_septos, "I_min": I}


if __name__ == "__main__":
    p = projeto("oliveira"); g = p["geometria"]
    D = g["largura"] - 2 * g["parede_externa"]                                   # 68 mm, espessura da célula
    a = (g["comprimento_bloco"] - 2 * g["parede_topo"] - (g["furos_por_bloco"] - 1) * g["septo"]) / g["furos_por_bloco"]
    b = g["altura_fiada"] - 2 * g["parede_leito"] if g.get("parede_leito") else g["altura_fiada"]
    F2, F3 = F_faixas(a, D), F_retangulos(a, b, D)
    print(f"célula: {a:.1f} x {b:.1f} mm, faces opostas a {D:.0f} mm; fator de forma 2D = {F2:.3f}, 3D = {F3:.3f}")
    out = {"celula_mm": [a, b, D], "fator_forma_2D": F2, "fator_forma_3D": F3,
           "eps_efetiva_2D": eps_efetiva(F2), "eps_efetiva_3D": eps_efetiva(F3)}
    out["referencia"] = roda("referência: placas paralelas infinitas, eps = 0,9")
    out["F_2D"] = roda("radiação direta entre as faces, fator de forma 2D", eps_efetiva(F2))
    out["F_3D"] = roda("radiação direta entre as faces, fator de forma 3D", eps_efetiva(F3))
    out["sem_radiacao"] = roda("sem radiação nas cavidades (só condução do ar)", 0.0)
    out["septos_metade"] = roda("septos com metade da condutividade", 0.9, 0.5)
    out["F_3D_septos_metade"] = roda("fator de forma 3D e septos com metade da condutividade", eps_efetiva(F3), 0.5)
    out["sem_radiacao_septos_metade"] = roda("sem radiação e septos com metade da condutividade", 0.0, 0.5)
    json.dump(out, open(os.path.join(AQUI, "oliveira_secao_explicita.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
