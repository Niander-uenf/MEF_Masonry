"""Verification of the sensitivity analysis (Morris method) and of the wall runs.

V1  Morris plan: each trajectory has k + 1 points on the p-level grid; between neighbouring
    points only one parameter changes, by ±Δ = p / [2 (p - 1)]; each parameter changes once.
V2  Linear function y = Σ a_i x_i: EE_i = a_i on every trajectory (μ = a_i, μ* = |a_i|, σ = 0).
V3  Additive function y = Σ b_i x_i²: μ* equal to the mean, trajectory by trajectory, of the differences
    over the level pairs (0, 2/3) and (1/3, 1), computed separately; with r = 2000, μ* tends to the
    expected value b_i (the two pairs give 2/3 b_i and 4/3 b_i).
V4  Interaction y = x1 x2 (+ x3 without effect): σ > 0 for the two that interact; μ* = σ = 0 for x3.
V5  Sobol' g-function (a = 0, 1, 4.5, 9, 99, 99, 99, 99): the order of μ* follows the known order
    of importance, and the last four stay below 5% of the largest.
V6  Missing values (a run that failed): only the two neighbouring effects leave the count.
V7  Reproducibility: the same seed repeats the plan; the 95% interval contains μ*.
V8  Physical ranges: limits, imperfection with linear amplitude, copy of the project without touching the
    original, tables scaled together with k and α, eccentricity at top and base.
V9  Factors of the wall law: the peak strength is multiplied by the strength factor and the
    initial modulus by the stiffness factor.
V10 BCEPQ wall under ISO 834: the fast run of the study gives the same R, I and Pmax as the Wall
    analysis with the same mesh; 10% more strength increases the cold capacity.
V11 (--full) Small wall study (2 parameters, 2 trajectories): in parallel and
    sequentially, identical results.
V12 (--full) Fast mesh × recommended mesh for the nominal BCEPQ wall (informative).

Usage: python examples/15_verification_sensitivity.py [--full]
Writes results/verification_sensitivity.csv.
"""
import copy
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import project as pj, sensitivity as sens  # noqa: E402
from mef_masonry.output import save_csv  # noqa: E402

FULL = "--full" in sys.argv
lines = []


def record(case, model, ref, tol_rel, unit_str="", tol_abs=0.0):
    error = abs(model - ref) / max(abs(ref), 1e-30)
    ok = error <= tol_rel or abs(model - ref) <= tol_abs
    tol = f"{100 * tol_rel:.2g}%" if tol_abs <= 0 else f"{tol_abs:g} {unit_str}"
    txt = f"{100 * error:.3g}%" if abs(ref) > 1e-30 else f"{abs(model - ref):.2g} (absolute)"
    lines.append([case, round(model, 9), round(ref, 9), unit_str, txt, tol, "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: model {model:.6g} {unit_str}, reference {ref:.6g} ({txt})")


def condition(case, ok, a, b, unit_str, criterion):
    lines.append([case, round(float(a), 6), round(float(b), 6), unit_str, "—", criterion, "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: {float(a):.5g} and {float(b):.5g} {unit_str} ({criterion})")


def report(case, a, unit_str, note):
    lines.append([case, round(a, 4), "", unit_str, "—", note, "informative"])
    print(f"--  {case}: {a:.4g} {unit_str} ({note})")


def evaluate(f, X):
    r, n, k = X.shape
    return np.array([[f(X[j, s]) for s in range(n)] for j in range(r)])


def main():
    # ---------------------------------------------------------------------------- V1
    k, r, p = 6, 30, 4
    X = sens.trajectories(k, r, p, seed=3)
    D = p / (2 * (p - 1))
    diff = np.diff(X, axis=1)                                    # (r, k, k)
    changed = np.abs(diff) > 1e-12
    grid = np.abs(X * (p - 1) - np.round(X * (p - 1))).max()
    condition("V1 one parameter per step and each parameter once per trajectory",
             bool(np.all(changed.sum(2) == 1) and np.all(changed.sum(1) == 1)), changed.sum(), r * k, "changes", "r × k")
    record("V1 step size Δ", float(np.abs(diff[changed]).max()), D, 1e-12)
    condition("V1 points on the p-level grid within [0, 1]", grid < 1e-12 and X.min() >= 0 and X.max() <= 1,
             X.min(), X.max(), "", "0 ≤ x ≤ 1, on the grid")

    # ---------------------------------------------------------------------------- V2
    a = np.array([3.0, -2.0, 0.5, 0.0, 7.0, -0.1])
    e = sens.effects(X, evaluate(lambda x: a @ x, X))
    record("V2 linear: largest error of μ relative to a_i", float(np.abs(e["mu"] - a).max()), 0.0, 0, tol_abs=1e-12)
    record("V2 linear: largest error of μ* relative to |a_i|", float(np.abs(e["mu_star"] - np.abs(a)).max()), 0.0, 0,
              tol_abs=1e-12)
    record("V2 linear: largest σ", float(np.nanmax(e["sigma"])), 0.0, 0, tol_abs=1e-12)

    # ---------------------------------------------------------------------------- V3
    b = np.array([1.0, 2.0, -3.0, 0.5, 0.0, 4.0])
    e = sens.effects(X, evaluate(lambda x: b @ x ** 2, X))
    x0 = X[:, 0, :]                                              # starting level of each parameter
    par_low = x0 < 0.5                                         # 0 or 1/3 goes up; 2/3 or 1 goes down
    lo_n = np.where(np.isclose(x0, 0) | np.isclose(x0, 2 / 3), 0.0, 1 / 3)
    ee_ref = b * ((lo_n + D) ** 2 - lo_n ** 2) / D
    record("V3 additive: μ* against the mean over the level pairs", float(np.abs(e["mu_star"] - np.abs(ee_ref).mean(0)).max()),
              0.0, 0, tol_abs=1e-12)
    X2 = sens.trajectories(2, 2000, 4, seed=11)
    e2 = sens.effects(X2, evaluate(lambda x: x[0] ** 2 + 2 * x[1] ** 2, X2))
    record("V3 additive, r = 2000: μ* of x1² tends to 1", float(e2["mu_star"][0]), 1.0, 0.02)
    record("V3 additive, r = 2000: μ* of 2 x2² tends to 2", float(e2["mu_star"][1]), 2.0, 0.02)
    condition("V3 additive: σ > 0 (non-linear effect)", bool(np.all(e2["sigma"] > 0.2 * e2["mu_star"])), e2["sigma"][0],
             e2["mu_star"][0], "", "σ/μ* > 0.2")

    # ---------------------------------------------------------------------------- V4
    X3 = sens.trajectories(3, 40, 4, seed=5)
    e = sens.effects(X3, evaluate(lambda x: x[0] * x[1], X3))
    condition("V4 interaction x1 x2: σ/μ* for both", bool(e["sigma"][0] / e["mu_star"][0] > 0.3 and
                                                        e["sigma"][1] / e["mu_star"][1] > 0.3),
             e["sigma"][0] / e["mu_star"][0], e["sigma"][1] / e["mu_star"][1], "", "> 0.3")
    record("V4 parameter without effect: μ* + σ", float(e["mu_star"][2] + e["sigma"][2]), 0.0, 0, tol_abs=1e-14)

    # ---------------------------------------------------------------------------- V5
    ag = np.array([0, 1, 4.5, 9, 99, 99, 99, 99])


    def g_sobol(x):
        return float(np.prod((np.abs(4 * x - 2) + ag) / (1 + ag)))


    X5 = sens.trajectories(8, 50, 4, seed=2)
    e = sens.effects(X5, evaluate(g_sobol, X5))
    ms = e["mu_star"]
    condition("V5 Sobol' g-function: order of the first four", bool(np.all(np.diff(ms[:4]) < 0)), ms[0], ms[3], "",
             "μ*1 > μ*2 > μ*3 > μ*4")
    condition("V5 Sobol' g-function: the last four < 5% of the largest", bool(ms[4:].max() < 0.05 * ms.max()),
             ms[4:].max(), 0.05 * ms.max(), "", "negligible")

    # ---------------------------------------------------------------------------- V6
    Y = evaluate(lambda x: a @ x, X)
    Y[4, 3] = np.nan
    e = sens.effects(X, Y)
    outside = np.isnan(e["EE"]).sum()
    condition("V6 missing run: two effects out of the count, the others exact",
             outside == 2 and float(np.abs(e["mu"] - a).max()) < 1e-12, outside, 2, "effects", "2 and μ = a_i")

    # ---------------------------------------------------------------------------- V7
    same = np.array_equal(sens.trajectories(5, 8, 4, 9), sens.trajectories(5, 8, 4, 9))
    differ = not np.array_equal(sens.trajectories(5, 8, 4, 9), sens.trajectories(5, 8, 4, 10))
    condition("V7 same seed, same plan; other seed, other plan", same and differ, int(same), int(differ), "", "1 and 1")
    e = sens.effects(X5, evaluate(g_sobol, X5))
    inside = bool(np.all((e["ic_lo"] <= e["mu_star"] + 1e-12) & (e["mu_star"] <= e["ic_hi"] + 1e-12)))
    condition("V7 95% interval (bootstrap) contains μ*", inside, e["ic_lo"][0], e["ic_hi"][0], "", "ic_lo ≤ μ* ≤ ic_hi")

    # ---------------------------------------------------------------------------- V8
    name = [k for k in pj.PRESETS if "sensitivity under ISO" in k][0]
    proj = pj.apply(pj.default_project(), pj.PRESETS[name])
    condition("V8 ready-made sensitivity project passes validation", not pj.validate(proj), 0, len(pj.validate(proj)),
             "errors", "none")
    pars = sens.active(proj["sensitivity"])
    P = {d["key"]: d for d in pars}
    record("V8 imperfection at mid-range: 1/n is the mean of 1/1000 and 1/250", 1 / sens.value(P["imp"], 0.5),
              0.5 * (1 / 1000 + 1 / 250), 1e-12)
    record("V8 hole emissivity at the upper limit", sens.value(P["eps_holes"], 1.0), 0.95, 1e-12)
    base = pj.complete(copy.deepcopy(proj))
    base["ceramic"]["table_k"] = [[20, 0.5], [800, 0.6]]
    base["ceramic"]["table_thermal_strain"] = [[20, 0.0], [820, 4.8]]
    dup = copy.deepcopy(base)
    x = np.ones(len(pars))
    q = sens.apply(base, pars, x)
    condition("V8 the original project does not change", base == dup, 1, 1, "", "same")
    record("V8 k and the k(T) table scaled together (+20%)", q["ceramic"]["table_k"][1][1] / 0.6,
              q["ceramic"]["k"] / base["ceramic"]["k"], 1e-12)
    record("V8 α = 8: thermal strain table scaled (6 → 8)", q["ceramic"]["table_thermal_strain"][1][1], 4.8 * 8 / 6,
              1e-12, "‰")
    condition("V8 eccentricity at top and base", q["test"]["top_eccentricity"] == q["test"]["base_eccentricity"] == 10,
             q["test"]["top_eccentricity"], q["test"]["base_eccentricity"], "mm", "10 and 10")
    condition("V8 k_tr switches creep on", q["model"]["creep"] and q["model"]["k_tr_unit"] == 2.35,
             q["model"]["k_tr_unit"], 2.35, "", "creep on")

    # ---------------------------------------------------------------------------- V9
    from mef_masonry.wall import WallLaw  # noqa: E402
    from mef_masonry.wall_data import PRISM_CURVES_BCEPQ, H_REF, A_GROSS_REF  # noqa: E402
    from mef_masonry.parameters import Geometry  # noqa: E402

    law = WallLaw.from_prism(PRISM_CURVES_BCEPQ, H_REF, A_GROSS_REF * pj.net_area_ratio(Geometry()))
    fr, fk = 1.15, 0.75
    l2 = law.scaled(fr, fr / fk)
    record("V9 peak strength × strength factor", float(l2.s_peak[0] / law.s_peak[0]), fr, 1e-12)
    record("V9 initial modulus × stiffness factor", float(l2.E0 / law.E0), fk, 1e-12)

    # ---------------------------------------------------------------------------- V10
    b0 = pj._sensitivity_base(proj)
    t0 = time.time()
    r0 = pj.evaluate_wall(b0)
    t_rap = time.time() - t0
    import tempfile  # noqa: E402
    with tempfile.TemporaryDirectory() as d:
        rp = pj.run(copy.deepcopy(b0), d, log=lambda s: None)
    res = {it[0]: it[1] for it in rp["summary"]}
    # the summary of the Wall analysis keeps 4 significant digits
    record("V10 fast run = Wall analysis: fire resistance", r0["R"],
              float(res["Fire resistance time (loss of equilibrium)"]), 0, "min", tol_abs=0.05)
    record("V10 fast run = Wall analysis: insulation", r0["I"],
              float(res["Insulation time (140/180 °C on the unexposed face)"]), 0, "min", tol_abs=0.05)
    record("V10 fast run = Wall analysis: cold capacity", r0["Pmax"],
              float([v for k_, v in res.items() if k_.startswith("Wall capacity")][0]), 1e-3, "kN")
    b1 = copy.deepcopy(b0)
    b1["model"]["wall_strength_factor"] = 1.1
    r1 = pj.evaluate_wall(b1)
    condition("V10 10% more strength: the cold capacity grows", r1["Pmax"] > r0["Pmax"], r1["Pmax"], r0["Pmax"],
             "kN", "largest")
    report("V10 time of one fast wall run", t_rap, "s", "one process")

    # ---------------------------------------------------------------------------- V11 and V12
    if FULL:
        sset = {f"var_{d['key']}": d["key"] in ("h", "eps") for d in sens.PARAMETERS}
        pp = sens.active(sset)
        Xw = sens.trajectories(len(pp), 2, 4, seed=4)
        proj_w = [sens.apply(b0, pp, x) for x in Xw.reshape(-1, len(pp))]
        seq = sens.run_many(proj_w, None, processes=1, log=lambda s: None)
        par = sens.run_many(proj_w, None, processes=2, log=lambda s: None)
        diff = max(abs(a_[c] - b_[c]) for a_, b_ in zip(seq, par) for c in ("R", "I", "Pmax"))
        record("V11 small wall study: parallel = sequential (largest difference)", diff, 0.0, 0, "", tol_abs=0.0)
        rec = copy.deepcopy(b0)
        rec["mesh"].update({"h": 5, "n_elements": 32})
        t0 = time.time()
        rr = pj.evaluate_wall(rec)
        report("V12 recommended mesh (5 mm, 32 elements): R − R of the fast mesh", rr["R"] - r0["R"], "min",
                 f"R = {rr['R']:.1f} and {r0['R']:.1f} min")
        report("V12 recommended mesh: I − I of the fast mesh", rr["I"] - r0["I"], "min",
                 f"I = {rr['I']:.1f} and {r0['I']:.1f} min")
        report("V12 run time with the recommended / fast mesh", (time.time() - t0) / t_rap, "",
                 "time ratio")

    output = os.path.join(os.path.dirname(__file__), "..", "results", "verification_sensitivity.csv")
    save_csv(output, lines, ["case", "model", "reference", "unit", "error", "tolerance", "status"])
    nf = sum(1 for l in lines if l[-1] == "FAILED")
    nv = sum(1 for l in lines if l[-1] != "informative")
    print(f"\n{nv - nf}/{nv} cases within tolerance"
          + (f" (plus {len(lines) - nv} informative result(s))" if len(lines) > nv else "")
          + f". Written to {os.path.abspath(output)}")
    sys.exit(1 if nf else 0)


if __name__ == "__main__":      # the parallel processes (spawn) import this file
    main()
