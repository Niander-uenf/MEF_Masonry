"""Seção 8.1 do artigo: a triagem de Morris da parede BCEPQ repetida com a faixa de α
estreitada de 4 a 8 para 5 a 7 ×10⁻⁶/°C (as outras dez entradas, as trajetórias e a
semente iguais às da triagem do artigo). Mostra quanto o primeiro lugar de α depende da
faixa que lhe foi atribuída. ~12 min em dois núcleos."""
import os, sys, json
RAIZ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, RAIZ)
from mef_masonry import projeto as pj

p = json.load(open(os.path.join(RAIZ, "resultados", "sensibilidade_parede_bcepq", "projeto.json")))
p["sensibilidade"].update(min_alfa=5.0, max_alfa=7.0)
p["nome"] = "Parede BCEPQ — Morris com α de 5 a 7 ×10⁻⁶/°C"
pj.executar(p, os.path.join(RAIZ, "resultados", "artigo", "morris_alfa_5_7"))
