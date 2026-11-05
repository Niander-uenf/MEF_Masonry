"""Verificação da análise de sensibilidade (método de Morris) e das rodadas da parede.

V1  Plano de Morris: cada trajetória tem k + 1 pontos na grade de p níveis; entre pontos
    vizinhos muda um parâmetro só, de ±Δ = p / [2 (p - 1)]; cada parâmetro muda uma vez.
V2  Função linear y = Σ a_i x_i: EE_i = a_i em toda trajetória (μ = a_i, μ* = |a_i|, σ = 0).
V3  Função aditiva y = Σ b_i x_i²: μ* igual à média, trajetória a trajetória, das diferenças
    nos pares de níveis (0, 2/3) e (1/3, 1), calculada à parte; com r = 2000, μ* tende ao
    valor esperado b_i (os dois pares valem 2/3 b_i e 4/3 b_i).
V4  Interação y = x1 x2 (+ x3 sem efeito): σ > 0 nos dois que interagem; μ* = σ = 0 no x3.
V5  Função g de Sobol (a = 0, 1, 4,5, 9, 99, 99, 99, 99): a ordem dos μ* segue a ordem de
    importância conhecida, e os quatro últimos ficam abaixo de 5% do maior.
V6  Valores ausentes (rodada que falhou): só os dois efeitos vizinhos saem da conta.
V7  Reprodutibilidade: a mesma semente repete o plano; o intervalo de 95% contém μ*.
V8  Faixas físicas: extremos, imperfeição com amplitude linear, cópia do projeto sem mexer no
    original, tabelas escaladas junto com k e α, excentricidade no topo e na base.
V9  Fatores da lei da parede: a resistência de pico multiplica pelo fator de resistência e o
    módulo inicial, pelo fator de rigidez.
V10 Parede BCEPQ no ISO 834: a rodada rápida do estudo dá os mesmos R, I e Pmax que a análise
    Parede com a mesma malha; 10% a mais de resistência aumenta a capacidade a frio.
V11 (--completo) Estudo pequeno da parede (2 parâmetros, 2 trajetórias): em paralelo e em
    sequência, resultados idênticos.
V12 (--completo) Malha rápida × malha recomendada na parede BCEPQ nominal (informativo).

Uso: python exemplos/15_verificacao_sensibilidade.py [--completo]
Grava resultados/verificacao_sensibilidade.csv.
"""
import copy
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import projeto as pj, sensibilidade as sb  # noqa: E402
from mef_masonry.saida import salvar_csv  # noqa: E402

COMPLETO = "--completo" in sys.argv
linhas = []


def registrar(caso, modelo, ref, tol_rel, unidade="", tol_abs=0.0):
    erro = abs(modelo - ref) / max(abs(ref), 1e-30)
    ok = erro <= tol_rel or abs(modelo - ref) <= tol_abs
    tol = f"{100 * tol_rel:.2g}%" if tol_abs <= 0 else f"{tol_abs:g} {unidade}"
    txt = f"{100 * erro:.3g}%" if abs(ref) > 1e-30 else f"{abs(modelo - ref):.2g} (absoluto)"
    linhas.append([caso, round(modelo, 9), round(ref, 9), unidade, txt, tol, "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: modelo {modelo:.6g} {unidade}, referência {ref:.6g} ({txt})")


def condicao(caso, ok, a, b, unidade, criterio):
    linhas.append([caso, round(float(a), 6), round(float(b), 6), unidade, "—", criterio, "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: {float(a):.5g} e {float(b):.5g} {unidade} ({criterio})")


def informar(caso, a, unidade, nota):
    linhas.append([caso, round(a, 4), "", unidade, "—", nota, "informativo"])
    print(f"--  {caso}: {a:.4g} {unidade} ({nota})")


def avaliar(f, X):
    r, n, k = X.shape
    return np.array([[f(X[j, s]) for s in range(n)] for j in range(r)])


def principal():
    # ---------------------------------------------------------------------------- V1
    k, r, p = 6, 30, 4
    X = sb.trajetorias(k, r, p, semente=3)
    D = p / (2 * (p - 1))
    dif = np.diff(X, axis=1)                                    # (r, k, k)
    mudou = np.abs(dif) > 1e-12
    grade = np.abs(X * (p - 1) - np.round(X * (p - 1))).max()
    condicao("V1 um parâmetro por passo e cada parâmetro uma vez por trajetória",
             bool(np.all(mudou.sum(2) == 1) and np.all(mudou.sum(1) == 1)), mudou.sum(), r * k, "mudanças", "r × k")
    registrar("V1 tamanho do passo Δ", float(np.abs(dif[mudou]).max()), D, 1e-12)
    condicao("V1 pontos na grade de p níveis dentro de [0, 1]", grade < 1e-12 and X.min() >= 0 and X.max() <= 1,
             X.min(), X.max(), "", "0 ≤ x ≤ 1, na grade")

    # ---------------------------------------------------------------------------- V2
    a = np.array([3.0, -2.0, 0.5, 0.0, 7.0, -0.1])
    e = sb.efeitos(X, avaliar(lambda x: a @ x, X))
    registrar("V2 linear: maior erro de μ em relação a a_i", float(np.abs(e["mu"] - a).max()), 0.0, 0, tol_abs=1e-12)
    registrar("V2 linear: maior erro de μ* em relação a |a_i|", float(np.abs(e["mu_estrela"] - np.abs(a)).max()), 0.0, 0,
              tol_abs=1e-12)
    registrar("V2 linear: maior σ", float(np.nanmax(e["sigma"])), 0.0, 0, tol_abs=1e-12)

    # ---------------------------------------------------------------------------- V3
    b = np.array([1.0, 2.0, -3.0, 0.5, 0.0, 4.0])
    e = sb.efeitos(X, avaliar(lambda x: b @ x ** 2, X))
    x0 = X[:, 0, :]                                              # nível de partida de cada parâmetro
    par_baixo = x0 < 0.5                                         # 0 ou 1/3 sobe; 2/3 ou 1 desce
    lo_n = np.where(np.isclose(x0, 0) | np.isclose(x0, 2 / 3), 0.0, 1 / 3)
    ee_ref = b * ((lo_n + D) ** 2 - lo_n ** 2) / D
    registrar("V3 aditiva: μ* contra a média pelos pares de níveis", float(np.abs(e["mu_estrela"] - np.abs(ee_ref).mean(0)).max()),
              0.0, 0, tol_abs=1e-12)
    X2 = sb.trajetorias(2, 2000, 4, semente=11)
    e2 = sb.efeitos(X2, avaliar(lambda x: x[0] ** 2 + 2 * x[1] ** 2, X2))
    registrar("V3 aditiva, r = 2000: μ* de x1² tende a 1", float(e2["mu_estrela"][0]), 1.0, 0.02)
    registrar("V3 aditiva, r = 2000: μ* de 2 x2² tende a 2", float(e2["mu_estrela"][1]), 2.0, 0.02)
    condicao("V3 aditiva: σ > 0 (efeito não linear)", bool(np.all(e2["sigma"] > 0.2 * e2["mu_estrela"])), e2["sigma"][0],
             e2["mu_estrela"][0], "", "σ/μ* > 0,2")

    # ---------------------------------------------------------------------------- V4
    X3 = sb.trajetorias(3, 40, 4, semente=5)
    e = sb.efeitos(X3, avaliar(lambda x: x[0] * x[1], X3))
    condicao("V4 interação x1 x2: σ/μ* nos dois", bool(e["sigma"][0] / e["mu_estrela"][0] > 0.3 and
                                                        e["sigma"][1] / e["mu_estrela"][1] > 0.3),
             e["sigma"][0] / e["mu_estrela"][0], e["sigma"][1] / e["mu_estrela"][1], "", "> 0,3")
    registrar("V4 parâmetro sem efeito: μ* + σ", float(e["mu_estrela"][2] + e["sigma"][2]), 0.0, 0, tol_abs=1e-14)

    # ---------------------------------------------------------------------------- V5
    ag = np.array([0, 1, 4.5, 9, 99, 99, 99, 99])


    def g_sobol(x):
        return float(np.prod((np.abs(4 * x - 2) + ag) / (1 + ag)))


    X5 = sb.trajetorias(8, 50, 4, semente=2)
    e = sb.efeitos(X5, avaliar(g_sobol, X5))
    ms = e["mu_estrela"]
    condicao("V5 função g de Sobol: ordem dos quatro primeiros", bool(np.all(np.diff(ms[:4]) < 0)), ms[0], ms[3], "",
             "μ*1 > μ*2 > μ*3 > μ*4")
    condicao("V5 função g de Sobol: os quatro últimos < 5% do maior", bool(ms[4:].max() < 0.05 * ms.max()),
             ms[4:].max(), 0.05 * ms.max(), "", "desprezíveis")

    # ---------------------------------------------------------------------------- V6
    Y = avaliar(lambda x: a @ x, X)
    Y[4, 3] = np.nan
    e = sb.efeitos(X, Y)
    fora = np.isnan(e["EE"]).sum()
    condicao("V6 rodada ausente: dois efeitos fora da conta, os demais exatos",
             fora == 2 and float(np.abs(e["mu"] - a).max()) < 1e-12, fora, 2, "efeitos", "2 e μ = a_i")

    # ---------------------------------------------------------------------------- V7
    igual = np.array_equal(sb.trajetorias(5, 8, 4, 9), sb.trajetorias(5, 8, 4, 9))
    difer = not np.array_equal(sb.trajetorias(5, 8, 4, 9), sb.trajetorias(5, 8, 4, 10))
    condicao("V7 mesma semente, mesmo plano; outra semente, outro plano", igual and difer, int(igual), int(difer), "", "1 e 1")
    e = sb.efeitos(X5, avaliar(g_sobol, X5))
    dentro = bool(np.all((e["ic_lo"] <= e["mu_estrela"] + 1e-12) & (e["mu_estrela"] <= e["ic_hi"] + 1e-12)))
    condicao("V7 intervalo de 95% (bootstrap) contém μ*", dentro, e["ic_lo"][0], e["ic_hi"][0], "", "ic_lo ≤ μ* ≤ ic_hi")

    # ---------------------------------------------------------------------------- V8
    nome = [k for k in pj.PREDEFINIDOS if "sensibilidade no ISO" in k][0]
    proj = pj.aplicar(pj.projeto_padrao(), pj.PREDEFINIDOS[nome])
    condicao("V8 projeto pronto da sensibilidade passa na validação", not pj.validar(proj), 0, len(pj.validar(proj)),
             "erros", "nenhum")
    pars = sb.ativos(proj["sensibilidade"])
    P = {d["chave"]: d for d in pars}
    registrar("V8 imperfeição no meio da faixa: 1/n é a média de 1/1000 e 1/250", 1 / sb.valor(P["imp"], 0.5),
              0.5 * (1 / 1000 + 1 / 250), 1e-12)
    registrar("V8 emissividade dos furos no extremo superior", sb.valor(P["eps_furos"], 1.0), 0.95, 1e-12)
    base = pj.completar(copy.deepcopy(proj))
    base["ceramica"]["tabela_k"] = [[20, 0.5], [800, 0.6]]
    base["ceramica"]["tabela_eps_termica"] = [[20, 0.0], [820, 4.8]]
    copia = copy.deepcopy(base)
    x = np.ones(len(pars))
    q = sb.aplicar(base, pars, x)
    condicao("V8 o projeto original não muda", base == copia, 1, 1, "", "igual")
    registrar("V8 k e a tabela k(T) escalados juntos (+20%)", q["ceramica"]["tabela_k"][1][1] / 0.6,
              q["ceramica"]["k"] / base["ceramica"]["k"], 1e-12)
    registrar("V8 α = 8: tabela de deformação térmica escalada (6 → 8)", q["ceramica"]["tabela_eps_termica"][1][1], 4.8 * 8 / 6,
              1e-12, "‰")
    condicao("V8 excentricidade no topo e na base", q["ensaio"]["excentricidade_topo"] == q["ensaio"]["excentricidade_base"] == 10,
             q["ensaio"]["excentricidade_topo"], q["ensaio"]["excentricidade_base"], "mm", "10 e 10")
    condicao("V8 k_tr liga a fluência", q["modelo"]["fluencia"] and q["modelo"]["k_tr_bloco"] == 2.35,
             q["modelo"]["k_tr_bloco"], 2.35, "", "fluência ligada")

    # ---------------------------------------------------------------------------- V9
    from mef_masonry.parede import LeiParede  # noqa: E402
    from mef_masonry.dados_parede import CURVAS_PRISMA_BCEPQ, H_REF, A_BRUTA_REF  # noqa: E402
    from mef_masonry.parametros import Geometria  # noqa: E402

    lei = LeiParede.de_prisma(CURVAS_PRISMA_BCEPQ, H_REF, A_BRUTA_REF * pj.razao_area_liquida(Geometria()))
    fr, fk = 1.15, 0.75
    l2 = lei.escalada(fr, fr / fk)
    registrar("V9 resistência de pico × fator de resistência", float(l2.s_pico[0] / lei.s_pico[0]), fr, 1e-12)
    registrar("V9 módulo inicial × fator de rigidez", float(l2.E0 / lei.E0), fk, 1e-12)

    # ---------------------------------------------------------------------------- V10
    b0 = pj._base_sensibilidade(proj)
    t0 = time.time()
    r0 = pj.avaliar_parede(b0)
    t_rap = time.time() - t0
    import tempfile  # noqa: E402
    with tempfile.TemporaryDirectory() as d:
        rp = pj.executar(copy.deepcopy(b0), d, log=lambda s: None)
    res = {it[0]: it[1] for it in rp["resumo"]}
    # o resumo da análise Parede guarda 4 algarismos significativos
    registrar("V10 rodada rápida = análise Parede: tempo de resistência", r0["R"],
              float(res["Tempo de resistência (perda de equilíbrio)"]), 0, "min", tol_abs=0.05)
    registrar("V10 rodada rápida = análise Parede: isolamento", r0["I"],
              float(res["Tempo de isolamento (140/180 °C na face não exposta)"]), 0, "min", tol_abs=0.05)
    registrar("V10 rodada rápida = análise Parede: capacidade a frio", r0["Pmax"],
              float([v for k_, v in res.items() if k_.startswith("Capacidade da parede")][0]), 1e-3, "kN")
    b1 = copy.deepcopy(b0)
    b1["modelo"]["fator_resistencia_parede"] = 1.1
    r1 = pj.avaliar_parede(b1)
    condicao("V10 10% a mais de resistência: a capacidade a frio cresce", r1["Pmax"] > r0["Pmax"], r1["Pmax"], r0["Pmax"],
             "kN", "maior")
    informar("V10 tempo de uma rodada rápida da parede", t_rap, "s", "um processo")

    # ---------------------------------------------------------------------------- V11 e V12
    if COMPLETO:
        se = {f"var_{d['chave']}": d["chave"] in ("h", "eps") for d in sb.PARAMETROS}
        pp = sb.ativos(se)
        Xw = sb.trajetorias(len(pp), 2, 4, semente=4)
        proj_w = [sb.aplicar(b0, pp, x) for x in Xw.reshape(-1, len(pp))]
        seq = sb.rodar(proj_w, None, processos=1, log=lambda s: None)
        par = sb.rodar(proj_w, None, processos=2, log=lambda s: None)
        dif = max(abs(a_[c] - b_[c]) for a_, b_ in zip(seq, par) for c in ("R", "I", "Pmax"))
        registrar("V11 estudo pequeno da parede: paralelo = sequencial (maior diferença)", dif, 0.0, 0, "", tol_abs=0.0)
        rec = copy.deepcopy(b0)
        rec["malha"].update({"h": 5, "n_elementos": 32})
        t0 = time.time()
        rr = pj.avaliar_parede(rec)
        informar("V12 malha recomendada (5 mm, 32 elementos): R − R da malha rápida", rr["R"] - r0["R"], "min",
                 f"R = {rr['R']:.1f} e {r0['R']:.1f} min")
        informar("V12 malha recomendada: I − I da malha rápida", rr["I"] - r0["I"], "min",
                 f"I = {rr['I']:.1f} e {r0['I']:.1f} min")
        informar("V12 tempo da rodada com a malha recomendada / rápida", (time.time() - t0) / t_rap, "",
                 "razão de tempos")

    saida = os.path.join(os.path.dirname(__file__), "..", "resultados", "verificacao_sensibilidade.csv")
    salvar_csv(saida, linhas, ["caso", "modelo", "referência", "unidade", "erro", "tolerância", "situação"])
    nf = sum(1 for l in linhas if l[-1] == "FALHOU")
    nv = sum(1 for l in linhas if l[-1] != "informativo")
    print(f"\n{nv - nf}/{nv} casos dentro da tolerância"
          + (f" (mais {len(linhas) - nv} resultado(s) informativo(s))" if len(linhas) > nv else "")
          + f". Gravado em {os.path.abspath(saida)}")
    sys.exit(1 if nf else 0)


if __name__ == "__main__":      # os processos do paralelo (spawn) importam este arquivo
    principal()
