"""Convergência da parede BCEPQ no ISO 834: malha da seção, elementos de viga e passo de tempo."""
import sys, os, time, json, copy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from mef_masonry import projeto as pj
AQUI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "resultados", "artigo")

def um(caso):
    h, n, dt = caso
    p = pj.aplicar(pj.projeto_padrao(), pj.PREDEFINIDOS["Parede BCEPQ 2,94 × 3,00 m — ISO 834 sob carga"])
    p["ensaio"].update({"duracao_iso": 240, "niveis_carga": 0})
    p["malha"].update({"h": h, "n_elementos": n, "dt": dt})
    q = pj.completar(p)
    geo = pj.geometria(q); _, m2, sec = pj._secao_parede(q, geo)
    t = time.time()
    r = pj.avaliar_parede(q)
    return {"h": h, "n_el": n, "dt": dt, "fibras": len(sec.y), "nos_secao": m2.n_nos, "tempo_s": time.time() - t,
            **{k: (float(v) if isinstance(v, (int, float)) else v) for k, v in r.items()}}

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    casos = sorted({(h, 32, 30) for h in (10, 7.5, 5, 3.5)} | {(5, n, 30) for n in (8, 16, 32, 64)} | {(5, 32, d) for d in (60, 30, 15)})
    out = []
    with ProcessPoolExecutor(2, mp_context=mp.get_context("spawn")) as ex:
        for r in ex.map(um, casos):
            print(json.dumps(r), flush=True); out.append(r)
    json.dump(out, open(os.path.join(AQUI, "conv_parede.json"), "w"), indent=1)
    print("fim", flush=True)
