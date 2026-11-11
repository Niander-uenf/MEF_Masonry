"""Section 8.1 of the paper: two checks of the Morris screening of the BCEPQ wall.

(a) Truncated runs. A run that reaches the end of the test (300 min) without failure enters the
    screening with R = 300 min, a lower bound. Here the elementary effects of R are
    recomputed (i) as in the paper and (ii) discarding every effect that involves a truncated
    run, to see whether the order of the parameters changes.
(b) k_tr. The nominal run of the screening uses k_tr = 0 (no published value for fired
    clay); here the same wall runs with k_tr = 2.35 (the upper limit of the range, a concrete
    value from EN 1992-1-2 / Anderberg and Thelandersson), to give the sign and size of the effect
    on the nominal run.

Reads results/sensitivity_wall_bcepq/ and writes results/paper/morris_checks.json.
The screening with the narrowed α range is in morris_narrow_alpha.py."""
import csv
import json
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
from mef_masonry import project as pj, sensitivity as sens  # noqa: E402

FOLDER = os.path.join(ROOT, "results", "sensitivity_wall_bcepq")
p = json.load(open(os.path.join(FOLDER, "project.json")))
sset = p["sensitivity"]
pars = sens.active(sset)
k, r, nv = len(pars), int(sset["trajectories"]), int(sset["levels"])
X = sens.trajectories(k, r, nv, int(sset["seed"]))

# ---------------------------------------------------------------- (a) truncated runs
with open(os.path.join(FOLDER, "sensitivity_runs.csv"), encoding="utf-8") as f:
    lines = list(csv.DictReader(f))
body = lines[1:]
assert len(body) == r * (k + 1)
R = np.array([float(x["R (min)"]) for x in body]).reshape(r, k + 1)
trunc = np.array(["no failure" in x["remark"] for x in body]).reshape(r, k + 1)
R_excl = np.where(trunc, np.nan, R)

ef_all = sens.effects(X, R)
ef_excl = sens.effects(X, R_excl)


def table(ef):
    o = np.argsort(-np.nan_to_num(ef["mu_star"], nan=-1))
    num = lambda v: round(float(v), 2) if np.isfinite(v) else None  # no valid effect: null
    return [{"parameter": pars[i]["short"], "mu_star": num(ef["mu_star"][i]), "mu": num(ef["mu"][i]),
             "sigma": num(ef["sigma"][i])} for i in o]


# number of valid elementary effects per parameter, without the truncated runs
n_ef = np.zeros(k, int)
for j in range(r):
    for s in range(k):
        i = int(np.argmax(np.abs(X[j, s + 1] - X[j, s]) > 1e-9))
        if np.isfinite(R_excl[j, s]) and np.isfinite(R_excl[j, s + 1]):
            n_ef[i] += 1

output = {
    "runs": int(R.size), "truncated_runs": int(trunc.sum()),
    "R_as_in_paper": table(ef_all),
    "R_without_truncated_runs": table(ef_excl),
    "valid_effects_per_parameter_without_truncated": {pars[i]["short"]: int(n_ef[i]) for i in range(k)},
}

# ---------------------------------------------------------------- (b) k_tr on the nominal run
class _S:  # minimal summary required by _wall_law
    def item(self, *a, **kw):
        pass


base = pj._sensitivity_base(p)
law = pj._wall_law(base, pj.geometry(base), print, _S())
par_ktr = sens.BY_KEY["ktr"]
res = {}
for v in (0.0, 2.35):
    q = sens.apply(base, [{**par_ktr, "lo": v, "hi": v}], [0.0])
    out = sens.run_many([q], law, 1, print)[0]
    res[f"k_tr = {v:g}"] = {"R": out.get("R"), "R_censored": out.get("R_censored"), "I": out.get("I")}
    print(v, res[f"k_tr = {v:g}"])
output["k_tr_in_nominal_run"] = res

os.makedirs(os.path.join(ROOT, "results", "paper"), exist_ok=True)
with open(os.path.join(ROOT, "results", "paper", "morris_checks.json"), "w", encoding="utf-8") as f:
    json.dump(output, f, ensure_ascii=False, indent=1)
print(json.dumps({kk: output[kk] for kk in ("truncated_runs", "R_without_truncated_runs")}, ensure_ascii=False,
                 indent=1))
