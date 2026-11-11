"""Example of optimising the section geometry (2D model).

Variables: hole diameter D and shift dy of the holes towards the
unexposed face. Objective: a lighter unit (larger hole area) that still meets
a minimum insulation time under ISO 834 and a minimum wall thickness.
Also traces the Pareto front (net area x insulation) on a grid.

It is a starting model: replace the objective and the constraints with those of the study
(e.g. maximum thermal tension, manufacturing criterion, NBR 15270).

Usage: python examples/05_section_optimisation.py [--target 90] [--full]
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import studies
from mef_masonry.parameters import Geometry
from mef_masonry.thermal import Materials
from mef_masonry.output import save_csv

ap = argparse.ArgumentParser()
ap.add_argument("--target", type=float, default=90.0, help="minimum insulation, min")
ap.add_argument("--min_wall", type=float, default=15.0, help="minimum wall, mm")
ap.add_argument("--full", action="store_true", help="fine mesh (slow)")
a = ap.parse_args()
h, dt = (0.005, 30.0) if a.full else (0.008, 60.0)

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT_DIR, exist_ok=True)
geo = Geometry()
W = geo.width
M = Materials()
cache = {}


def evaluate(D_mm, dy_mm):
    key = (round(D_mm, 2), round(dy_mm, 2))
    if key not in cache:
        p = studies.wall_section_iso834(M, geo, h=h, dt=dt, diameter=D_mm / 1000,
                                        y_holes=W / 2 + dy_mm / 1000, t_max_min=300, t_min_stress=1e9)
        cache[key] = (p["t_insulation_min"], p["net_area_ratio"])
    return cache[key]


def wall(D_mm, dy_mm):
    return W * 1000 / 2 - abs(dy_mm) - D_mm / 2


def objective(x):
    D, dy = x
    pen = 0.0
    if wall(D, dy) < a.min_wall:
        return 1e3 + (a.min_wall - wall(D, dy)) * 10
    t, ral = evaluate(D, dy)
    if t < a.target:
        pen += 10 * (a.target - t)
    return ral + pen / 100    # minimise the net area (mass) with insulation >= target


t0 = time.time()
print("Grid for the Pareto front...")
grid = []
for D in np.arange(50, 96, 5):
    for dy in np.arange(-15, 16, 5):
        if wall(D, dy) >= a.min_wall:
            t, ral = evaluate(D, dy)
            grid.append((D, dy, t, ral))
grid = np.array(grid)
# Pareto: smallest net area and largest insulation
pareto = [i for i, (D, dy, t, r) in enumerate(grid)
          if not any((g[3] <= r and g[2] >= t) and (g[3] < r or g[2] > t) for g in grid)]
print(f"  {len(grid)} geometries in {time.time() - t0:.0f} s; {len(pareto)} on the Pareto front")

print(f"Optimisation (differential evolution): insulation >= {a.target:.0f} min, wall >= {a.min_wall:.0f} mm")
res = differential_evolution(objective, bounds=[(50, 100), (-20, 20)], seed=1, maxiter=15, popsize=8,
                             tol=1e-3, polish=False)
D, dy = res.x
t, ral = evaluate(D, dy)
t_now, ral_current = evaluate(geo.hole_diameter * 1000, 0.0)
print(f"\nCurrent geometry: D = 84 mm, dy = 0 -> insulation {t_now:.0f} min, net area {ral_current:.3f}")
print(f"Optimum:          D = {D:.1f} mm, dy = {dy:+.1f} mm -> insulation {t:.0f} min, net area {ral:.3f}"
      f" (minimum wall {wall(D, dy):.1f} mm)")
print(f"Model evaluations: {len(cache)}; total time {time.time() - t0:.0f} s")

save_csv(os.path.join(OUT_DIR, "optimisation_grid.csv"),
           [[f"{g[0]:.0f}", f"{g[1]:.0f}", f"{g[2]:.1f}", f"{g[3]:.3f}", "yes" if i in pareto else ""]
            for i, g in enumerate(grid)],
           ["D_mm", "dy_mm", "insulation_min", "net_area_ratio", "pareto"])

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.scatter(grid[:, 3], grid[:, 2], c="0.75", s=18, label="grid geometries")
P = grid[pareto][np.argsort(grid[pareto][:, 3])]
ax.plot(P[:, 3], P[:, 2], "o-", color="tab:blue", ms=4, label="Pareto front")
ax.plot(ral_current, t_now, "s", color="k", label="current BCEPQ (D = 84 mm)")
ax.plot(ral, t, "*", color="tab:orange", ms=12, label=f"optimum (target {a.target:.0f} min)")
ax.axhline(a.target, color="0.5", lw=0.8, ls="--")
ax.set_xlabel("net / gross area (proportional to the mass)")
ax.set_ylabel("insulation under ISO 834 (min)")
ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, "optimisation_pareto.png"), dpi=130)
print("Results in", os.path.abspath(OUT_DIR))
