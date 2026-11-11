"""Convergence of the BCEPQ wall under ISO 834: section mesh, beam elements and time step."""
import sys, os, time, json, copy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from mef_masonry import project as pj
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results", "paper")

def one(case):
    h, n, dt = case
    p = pj.apply(pj.default_project(), pj.PRESETS["BCEPQ wall 2.94 × 3.00 m — ISO 834 under load"])
    p["test"].update({"iso_duration": 240, "load_levels": 0})
    p["mesh"].update({"h": h, "n_elements": n, "dt": dt})
    q = pj.complete(p)
    geo = pj.geometry(q); _, m2, sec = pj._wall_section(q, geo)
    t = time.time()
    r = pj.evaluate_wall(q)
    return {"h": h, "n_el": n, "dt": dt, "fibres": len(sec.y), "section_nodes": m2.n_nodes, "time_s": time.time() - t,
            **{k: (float(v) if isinstance(v, (int, float)) else v) for k, v in r.items()}}

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    cases = sorted({(h, 32, 30) for h in (10, 7.5, 5, 3.5)} | {(5, n, 30) for n in (8, 16, 32, 64)} | {(5, 32, d) for d in (60, 30, 15)})
    out = []
    with ProcessPoolExecutor(2, mp_context=mp.get_context("spawn")) as ex:
        for r in ex.map(one, cases):
            print(json.dumps(r), flush=True); out.append(r)
    json.dump(out, open(os.path.join(HERE, "wall_convergence.json"), "w"), indent=1)
    print("end", flush=True)
