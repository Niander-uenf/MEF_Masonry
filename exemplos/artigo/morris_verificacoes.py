"""Seção 8.1 do artigo: duas verificações da triagem de Morris da parede BCEPQ.

(a) Rodadas truncadas. Uma rodada que chega ao fim do ensaio (300 min) sem ruína entra na
    triagem com R = 300 min, um limite inferior. Aqui os efeitos elementares de R são
    recalculados (i) como no artigo e (ii) descartando todo efeito que envolve uma rodada
    truncada, para ver se a ordem dos parâmetros muda.
(b) k_tr. A rodada nominal da triagem usa k_tr = 0 (nenhum valor publicado para a cerâmica
    queimada); aqui a mesma parede roda com k_tr = 2,35 (o limite superior da faixa, valor de
    concreto do EN 1992-1-2 / Anderberg e Thelandersson), para dar o sinal e o tamanho do efeito
    na rodada nominal.

Lê resultados/sensibilidade_parede_bcepq/ e grava resultados/artigo/morris_verificacoes.json.
A triagem com a faixa de α estreitada está em morris_alfa_estreito.py."""
import csv
import json
import os
import sys

import numpy as np

RAIZ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, RAIZ)
from mef_masonry import projeto as pj, sensibilidade as sb  # noqa: E402

PASTA = os.path.join(RAIZ, "resultados", "sensibilidade_parede_bcepq")
p = json.load(open(os.path.join(PASTA, "projeto.json")))
se = p["sensibilidade"]
pars = sb.ativos(se)
k, r, nv = len(pars), int(se["trajetorias"]), int(se["niveis"])
X = sb.trajetorias(k, r, nv, int(se["semente"]))

# ---------------------------------------------------------------- (a) rodadas truncadas
with open(os.path.join(PASTA, "sensibilidade_rodadas.csv"), encoding="utf-8") as f:
    linhas = list(csv.DictReader(f))
corpo = linhas[1:]
assert len(corpo) == r * (k + 1)
R = np.array([float(x["R (min)"]) for x in corpo]).reshape(r, k + 1)
trunc = np.array(["sem ruína" in x["observação"] for x in corpo]).reshape(r, k + 1)
R_sem = np.where(trunc, np.nan, R)

ef_todos = sb.efeitos(X, R)
ef_sem = sb.efeitos(X, R_sem)


def tabela(ef):
    o = np.argsort(-np.nan_to_num(ef["mu_estrela"], nan=-1))
    num = lambda v: round(float(v), 2) if np.isfinite(v) else None  # sem efeito válido: null
    return [{"parametro": pars[i]["curto"], "mu_estrela": num(ef["mu_estrela"][i]), "mu": num(ef["mu"][i]),
             "sigma": num(ef["sigma"][i])} for i in o]


# número de efeitos elementares válidos por parâmetro, sem as rodadas truncadas
n_ef = np.zeros(k, int)
for j in range(r):
    for s in range(k):
        i = int(np.argmax(np.abs(X[j, s + 1] - X[j, s]) > 1e-9))
        if np.isfinite(R_sem[j, s]) and np.isfinite(R_sem[j, s + 1]):
            n_ef[i] += 1

saida = {
    "rodadas": int(R.size), "rodadas_truncadas": int(trunc.sum()),
    "R_como_no_artigo": tabela(ef_todos),
    "R_sem_rodadas_truncadas": tabela(ef_sem),
    "efeitos_validos_por_parametro_sem_truncadas": {pars[i]["curto"]: int(n_ef[i]) for i in range(k)},
}

# ---------------------------------------------------------------- (b) k_tr na rodada nominal
class _S:  # resumo mínimo pedido por _lei_parede
    def item(self, *a, **kw):
        pass


base = pj._base_sensibilidade(p)
lei = pj._lei_parede(base, pj.geometria(base), print, _S())
par_ktr = sb.POR_CHAVE["ktr"]
res = {}
for v in (0.0, 2.35):
    q = sb.aplicar(base, [{**par_ktr, "lo": v, "hi": v}], [0.0])
    out = sb.rodar([q], lei, 1, print)[0]
    res[f"k_tr = {v:g}"] = {"R": out.get("R"), "R_censurado": out.get("R_censurado"), "I": out.get("I")}
    print(v, res[f"k_tr = {v:g}"])
saida["k_tr_na_rodada_nominal"] = res

os.makedirs(os.path.join(RAIZ, "resultados", "artigo"), exist_ok=True)
with open(os.path.join(RAIZ, "resultados", "artigo", "morris_verificacoes.json"), "w", encoding="utf-8") as f:
    json.dump(saida, f, ensure_ascii=False, indent=1)
print(json.dumps({kk: saida[kk] for kk in ("rodadas_truncadas", "R_sem_rodadas_truncadas")}, ensure_ascii=False,
                 indent=1))
