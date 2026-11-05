"""Ábacos de projeto da parede BCEPQ: tempo de resistência × nível de carga, por altura e α."""
import sys, os, time, json, copy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from mef_masonry import projeto as pj
AQUI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "resultados", "artigo")

def um(args):
    nf, alfa = args
    import contextlib, io
    p = pj.aplicar(pj.projeto_padrao(), pj.PREDEFINIDOS["Parede BCEPQ 2,94 × 3,00 m — ISO 834 sob carga"])
    p["geometria"]["n_fiadas"] = nf
    p["ceramica"]["alfa"] = alfa
    p["ensaio"].update({"duracao_iso": 240, "niveis_carga": 9, "carga": 100})
    p["malha"].update({"h": 10, "n_elementos": 16})
    pasta = os.path.join(AQUI, "abacos", f"H{nf}_a{alfa:g}")
    t = time.time()
    with contextlib.redirect_stdout(io.StringIO()):
        r = pj.executar(p, pasta, log=lambda s: None)
    return nf, alfa, time.time() - t, r["resumo"]

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    casos = [(nf, a) for a in (5.0, 6.0, 7.0) for nf in (40, 47, 54, 60, 67)]
    with ProcessPoolExecutor(2, mp_context=mp.get_context("spawn")) as ex:
        for nf, a, dt, res in ex.map(um, casos):
            print(nf, a, f"{dt:.0f} s", flush=True)
    print("fim", flush=True)
