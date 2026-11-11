"""Global sensitivity analysis by the Morris method (elementary effects).

Morris (1991), with the μ* measure of Campolongo, Cariboni and Saltelli (2007):
each trajectory starts at a point of the p-level grid in the hypercube [0, 1]^k and
changes one parameter at a time, in random order, by Δ = p / [2 (p - 1)].
The elementary effect of parameter i is EE_i = [y(x + Δ e_i) - y(x)] / (±Δ), on the
scale of the whole range of the parameter: μ* ≈ how much the output changes, on average, when
the parameter spans its range; σ measures non-linearity and interactions.

The study applied here is that of the wall under ISO 834 (project.evaluate_wall): for each
point, the project is built with the parameter values and run without figures.
"""
import copy
import os
import time
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Uncertain inputs of the wall (default ranges and sources in the interface)
# ---------------------------------------------------------------------------
# kind "rel": range in ± % of the project value; "abs": absolute minimum and maximum
PARAMETERS = [
    {"key": "k", "name": "Thermal conductivity of the unit", "short": "k", "kind": "rel", "unit": "%",
     "default": (20.0,)},
    {"key": "c", "name": "Specific heat of the unit", "short": "c", "kind": "rel", "unit": "%",
     "default": (15.0,)},
    {"key": "eps_holes", "name": "Emissivity of the hole walls", "short": "ε holes", "kind": "abs",
     "unit": "", "default": (0.70, 0.95)},
    {"key": "h", "name": "Convection on the exposed face", "short": "h", "kind": "abs", "unit": "W/m²·K",
     "default": (15.0, 35.0)},
    {"key": "eps", "name": "Surface emissivity", "short": "ε face", "kind": "abs", "unit": "",
     "default": (0.60, 0.95)},
    {"key": "alpha", "name": "Thermal expansion of the unit", "short": "α", "kind": "abs", "unit": "×10⁻⁶/°C",
     "default": (4.0, 8.0)},
    {"key": "ktr", "name": "Transient creep of the unit (k_tr)", "short": "k_tr", "kind": "abs", "unit": "",
     "default": (0.0, 2.35)},
    {"key": "strength", "name": "Masonry strength (factor on the law)", "short": "f_str", "kind": "rel",
     "unit": "%", "default": (15.0,)},
    {"key": "stiffness", "name": "Masonry stiffness (factor on the law)", "short": "f_stf", "kind": "rel", "unit": "%",
     "default": (25.0,)},
    {"key": "imp", "name": "Initial imperfection (H divided by)", "short": "imperf.", "kind": "abs", "unit": "",
     "default": (1000.0, 250.0)},
    {"key": "exc", "name": "Load eccentricity (top and base)", "short": "excentr.", "kind": "abs",
     "unit": "mm", "default": (-10.0, 10.0)},
]
BY_KEY = {d["key"]: d for d in PARAMETERS}
OUTPUTS = [("R", "fire resistance", "min"), ("I", "insulation time", "min"),
          ("Pmax", "cold capacity", "kN")]


def active(sens: Dict) -> List[Dict]:
    """Parameters selected on the Sensitivity page, with the given range (lo, hi)."""
    out = []
    for d in PARAMETERS:
        ch = d["key"]
        if not sens.get(f"var_{ch}", True):
            continue
        if d["kind"] == "rel":
            f = sens.get(f"range_{ch}", d["default"][0])
            out.append({**d, "lo": -f, "hi": f})
        else:
            out.append({**d, "lo": sens.get(f"min_{ch}", d["default"][0]), "hi": sens.get(f"max_{ch}", d["default"][1])})
    return out


def value(par: Dict, x: float) -> float:
    """Physical value of the parameter at position x ∈ [0, 1] of the range. For the imperfection, the
    amplitude H/n varies linearly (n goes from lo to hi)."""
    lo, hi = par["lo"], par["hi"]
    if par["key"] == "imp":
        a = (1 - x) / lo + x / hi
        return 1.0 / a
    return lo + (hi - lo) * x


def apply(base: Dict, pars: Sequence[Dict], x: Sequence[float]) -> Dict:
    """Project with the parameters at position x (copy)."""
    q = copy.deepcopy(base)
    c, ts, mdl = q["ceramic"], q["test"], q["model"]
    for par, xi in zip(pars, x):
        v = value(par, float(xi))
        ch = par["key"]
        if ch in ("k", "c"):
            f = 1 + v / 100
            c[ch] = c[ch] * f
            tab = f"table_{ch}"
            if c.get(tab):
                c[tab] = [[T, val * f] for T, val in c[tab]]
        elif ch == "eps_holes":
            ts["hole_emissivity"] = v
        elif ch == "h":
            ts["h_convection"] = v
        elif ch == "eps":
            ts["emissivity"] = v
        elif ch == "alpha":
            if c.get("table_thermal_strain") and c.get("alpha", 0) > 0:
                f = v / c["alpha"]
                c["table_thermal_strain"] = [[T, e * f] for T, e in c["table_thermal_strain"]]
            c["alpha"] = v
        elif ch == "ktr":
            mdl["creep"] = True
            mdl["k_tr_unit"] = v
            mdl["table_creep_unit"] = []
        elif ch == "strength":
            mdl["wall_strength_factor"] = mdl.get("wall_strength_factor", 1.0) * (1 + v / 100)
        elif ch == "stiffness":
            mdl["wall_stiffness_factor"] = mdl.get("wall_stiffness_factor", 1.0) * (1 + v / 100)
        elif ch == "imp":
            ts["imperfection"] = v
        elif ch == "exc":
            ts["top_eccentricity"] = ts["base_eccentricity"] = v
    return q


# ---------------------------------------------------------------------------
# Morris plan and elementary effects
# ---------------------------------------------------------------------------
def trajectories(k: int, r: int, p: int = 4, seed: int = 1) -> np.ndarray:
    """r trajectories of k + 1 points in the hypercube [0, 1]^k (grid of p levels)."""
    rng = np.random.default_rng(seed)
    levels = np.arange(p) / (p - 1)
    delta = p / (2 * (p - 1))
    X = np.empty((r, k + 1, k))
    for j in range(r):
        x = rng.choice(levels, size=k)
        X[j, 0] = x
        for s, i in enumerate(rng.permutation(k)):
            x = x.copy()
            x[i] = x[i] + delta if x[i] + delta <= 1 + 1e-12 else x[i] - delta
            X[j, s + 1] = x
    return X


def effects(X: np.ndarray, Y: np.ndarray, n_boot: int = 1000, seed: int = 7) -> Dict[str, np.ndarray]:
    """μ*, μ, σ and the 95% interval of μ* (bootstrap over trajectories) for the
    output Y (r, k + 1). Pairs with NaN are left out."""
    r, n, k = X.shape
    EE = np.full((r, k), np.nan)
    for j in range(r):
        for s in range(k):
            dx = X[j, s + 1] - X[j, s]
            i = int(np.argmax(np.abs(dx)))
            EE[j, i] = (Y[j, s + 1] - Y[j, s]) / dx[i]
    a = np.abs(EE)
    mu_s = np.nanmean(a, axis=0)
    mu = np.nanmean(EE, axis=0)
    sig = np.array([np.nanstd(EE[:, i], ddof=1) if np.sum(np.isfinite(EE[:, i])) > 1 else np.nan
                    for i in range(k)])
    rng = np.random.default_rng(seed)
    bs = np.array([np.nanmean(a[rng.integers(0, r, r)], axis=0) for _ in range(n_boot)])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
    return {"EE": EE, "mu_star": mu_s, "mu": mu, "sigma": sig, "ic_lo": lo, "ic_hi": hi,
            "valid": np.sum(np.isfinite(EE), axis=0)}


# ---------------------------------------------------------------------------
# Execution of the runs (in parallel when possible)
# ---------------------------------------------------------------------------
def _evaluate(args):
    import contextlib
    import io
    from .project import evaluate_wall
    q, law = args
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            return evaluate_wall(q, law)
    except Exception as e:  # noqa: BLE001 - a run that fails becomes NaN
        return {"error": f"{type(e).__name__}: {e}"}


def default_processes():
    """Cores minus one (the interface keeps responding); with 2 cores, both."""
    n = os.cpu_count() or 1
    return n if n <= 2 else n - 1


_SINGLE_THREAD = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def run_many(projects: Sequence[Dict], law, processes: int = 0, log: Callable[[str], None] = print) -> List[Dict]:
    """Evaluates the projects (parallel with processes > 1; 0 = automatic)."""
    n = len(projects)
    procs = processes or default_processes()
    res: List[Optional[Dict]] = [None] * n
    t0 = time.time()

    def progress(i_done):
        dt = time.time() - t0
        remaining = dt / i_done * (n - i_done)
        log(f"  {i_done}/{n} runs ({dt / 60:.1f} min; about {remaining / 60:.0f} min left)")

    pool = None
    if procs > 1 and n > 1:
        # new processes ("spawn", as on Windows): the server has other threads. Each
        # process uses a single BLAS thread, so as not to compete for the cores; the variable holds until the
        # processes are spawned, in submit
        before = {k: os.environ.get(k) for k in _SINGLE_THREAD}
        os.environ.update({k: "1" for k in _SINGLE_THREAD})
        try:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor, as_completed
            pool = ProcessPoolExecutor(max_workers=min(procs, n), mp_context=mp.get_context("spawn"))
            futs = {pool.submit(_evaluate, (q, law)): i for i, q in enumerate(projects)}
        except Exception as e:  # noqa: BLE001
            log(f"  No parallel processes ({e}); running sequentially.")
            if pool is not None:
                pool.shutdown(wait=False, cancel_futures=True)
            pool = None
        finally:
            for k, v in before.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
    if pool is not None:
        try:
            done = 0
            for f in as_completed(futs):
                try:
                    res[futs[f]] = f.result()
                except Exception as e:  # noqa: BLE001 - process that crashed
                    res[futs[f]] = {"error": f"{type(e).__name__}: {e}"}
                done += 1
                if done % max(1, n // 20) == 0 or done == n:
                    progress(done)
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
    else:
        for i, q in enumerate(projects):
            res[i] = _evaluate((q, law))
            if (i + 1) % max(1, n // 20) == 0 or i + 1 == n:
                progress(i + 1)
    return res
