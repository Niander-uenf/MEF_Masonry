"""Design charts of the BCEPQ wall: fire resistance time × load level, by height and α."""
import sys, os, time, json, copy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from mef_masonry import project as pj
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results", "paper")

def one(args):
    nf, alpha = args
    import contextlib, io
    p = pj.apply(pj.default_project(), pj.PRESETS["BCEPQ wall 2.94 × 3.00 m — ISO 834 under load"])
    p["geometry"]["n_courses"] = nf
    p["ceramic"]["alpha"] = alpha
    p["test"].update({"iso_duration": 240, "load_levels": 9, "load": 100})
    p["mesh"].update({"h": 10, "n_elements": 16})
    folder = os.path.join(HERE, "charts", f"H{nf}_a{alpha:g}")
    t = time.time()
    with contextlib.redirect_stdout(io.StringIO()):
        r = pj.run(p, folder, log=lambda s: None)
    return nf, alpha, time.time() - t, r["summary"]

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    cases = [(nf, a) for a in (5.0, 6.0, 7.0) for nf in (40, 47, 54, 60, 67)]
    with ProcessPoolExecutor(2, mp_context=mp.get_context("spawn")) as ex:
        for nf, a, dt, res in ex.map(one, cases):
            print(nf, a, f"{dt:.0f} s", flush=True)
    print("end", flush=True)
