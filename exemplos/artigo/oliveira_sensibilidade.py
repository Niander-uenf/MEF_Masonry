"""Validação térmica (seção 7.3 do artigo): sensibilidade da previsão de isolamento das paredes de
Oliveira et al. (2021) e seção maciça de 100 mm (comparação com a EI 90 da EN 1992-1-2)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oliveira import projeto, AQUI
from mef_masonry import projeto as pj

def roda(nome, p):
    r = pj.executar(p, os.path.join(AQUI, "oliveira_" + nome), log=lambda s: None)
    res = {it[0]: it[1] for it in r["resumo"]}
    print(f"{nome:28s} isolamento = {res.get('Tempo de isolamento (140/180 °C na face não exposta)')} min", flush=True)

if __name__ == "__main__":
    for var in ("biblioteca", "en_umidade"):
        p = projeto(var); p["geometria"].update(furos_por_bloco=0, forma_furo="circular", diametro_furo=0)
        roda("macica_" + var, p)
    p = projeto("oliveira"); p["ceramica"]["tabela_c"] = [[20, 900], [100, 900], [100.01, 2020], [115, 2020], [200, 1000], [400, 1100], [1200, 1100]]
    roda("umidade_3pct", p)
    p = projeto("oliveira"); p["ensaio"]["emissividade_furos"] = 0.7; roda("eps_furos_07", p)
    p = projeto("oliveira"); p["malha"]["h"] = 2; roda("malha_2mm", p)
    p = projeto("oliveira"); p["geometria"].update(furos_por_bloco=0, forma_furo="circular", diametro_furo=0)
    p["ceramica"].update(rho=2300 * 0.5365, tabela_rho=[], k=1.36 * 0.5365, tabela_k=[]); roda("homogenea_equivalente", p)
