import os
import sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from mef_masonry import projeto as pj, sensibilidade as sb
if __name__ == "__main__":
    nome = [k for k in pj.PREDEFINIDOS if "sensibilidade no ISO" in k][0]
    p = pj.aplicar(pj.projeto_padrao(), pj.PREDEFINIDOS[nome]); p["ensaio"]["duracao_iso"] = 240
    b = pj._base_sensibilidade(p)
    pars = [d for d in sb.ativos(p["sensibilidade"]) if d["chave"] in ("k", "c", "eps")]
    X = sb.trajetorias(len(pars), 2, 4, 3).reshape(-1, len(pars))
    projs = [sb.aplicar(b, pars, x) for x in X]
    for procs in (1, 2):
        t = time.time(); r = sb.rodar(projs, None, procs, log=lambda s: None); dt = time.time() - t
        print(procs, len(projs), f"{dt:.1f} s", [round(x["R"], 1) for x in r], flush=True)
