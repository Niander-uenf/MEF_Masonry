import os
import sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from mef_masonry import project as pj, sensitivity as sens
if __name__ == "__main__":
    name = [k for k in pj.PRESETS if "sensitivity under ISO" in k][0]
    p = pj.apply(pj.default_project(), pj.PRESETS[name]); p["test"]["iso_duration"] = 240
    b = pj._sensitivity_base(p)
    pars = [d for d in sens.active(p["sensitivity"]) if d["key"] in ("k", "c", "eps")]
    X = sens.trajectories(len(pars), 2, 4, 3).reshape(-1, len(pars))
    projs = [sens.apply(b, pars, x) for x in X]
    for procs in (1, 2):
        t = time.time(); r = sens.run_many(projs, None, procs, log=lambda s: None); dt = time.time() - t
        print(procs, len(projs), f"{dt:.1f} s", [round(x["R"], 1) for x in r], flush=True)
