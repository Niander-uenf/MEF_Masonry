"""Validação térmica: paredes de blocos de concreto de 3 células de Oliveira et al. (2021), ISO 834."""
import sys, os, json, csv, copy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from mef_masonry import projeto as pj, biblioteca as bib
AQUI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "resultados", "artigo")

def base():
    p = pj.projeto_padrao()
    p = bib.aplicar(p, "ceramica", "material_biblioteca", "concreto_nbr15220") if hasattr(bib, "aplicar") else p
    return p

def projeto(var):
    p = pj.projeto_padrao()
    m = copy.deepcopy(bib.MATERIAIS["concreto_nbr15220"]["mudancas"])
    for g, d in m.items(): p[g].update(d)
    p["ceramica"]["material_biblioteca"] = "concreto_nbr15220"
    p["geometria"].update({"bloco_biblioteca": "personalizado", "forma_furo": "retangular", "comprimento_bloco": 201,
        "largura": 100, "altura_fiada": 100, "furos_por_bloco": 3, "parede_externa": 16, "septo": 16, "parede_topo": 16,
        "celulas_parede": 0, "grauteamento": "nenhum", "blocos_por_fiada": 1, "n_fiadas": 10, "tipo_junta": "seca",
        "espessura_junta": 2})
    p["modelo"].update({"analise": "secao_2d", "elemento": "T3"})
    p["ensaio"].update({"tipo": "iso834", "T_inicial": 20, "duracao_iso": 150, "h_conveccao": 25})
    p["malha"].update({"h": 4, "dt": 30})
    c = p["ceramica"]
    if var in ("en_umidade", "oliveira"):
        # EN 1992-1-2: k limite inferior; c com o pico de umidade u = 1,5 % (1470 J/kg.K entre 100 e 115 °C); rho(T)
        c["tabela_k"] = [[T, 1.36 - 0.136 * (T / 100) + 0.0057 * (T / 100) ** 2] for T in range(20, 1201, 100)]
        c["tabela_c"] = [[20, 900], [100, 900], [100.01, 1470], [115, 1470], [200, 1000], [400, 1100], [1200, 1100]]
        c["rho"] = 2300
        c["tabela_rho"] = [[20, 2300], [115, 2300], [200, 2300 * 0.98], [400, 2300 * 0.95], [1200, 2300 * 0.88]]
    if var == "oliveira":
        p["ensaio"]["emissividade"] = 0.7
    return p

if __name__ == "__main__":
    out = {}
    for var in ("biblioteca", "en_umidade", "oliveira"):
        p = projeto(var)
        e = pj.validar(p)
        if e: print(var, e); continue
        r = pj.executar(p, os.path.join(AQUI, "oliveira_" + var), log=lambda s: None)
        res = {it[0]: it[1] for it in r["resumo"]}
        rows = list(csv.reader(open(os.path.join(AQUI, "oliveira_" + var, "temperaturas_2d.csv"), encoding="utf-8")))
        out[var] = {"I": res.get("Tempo de isolamento (140/180 °C na face não exposta)"), "razao": res.get("Razão área líquida/bruta"),
                    "T": [[float(a), float(c)] for a, b, c in rows[1:]]}
        print(var, out[var]["I"], out[var]["razao"], flush=True)
    json.dump(out, open(os.path.join(AQUI, "oliveira.json"), "w"))
