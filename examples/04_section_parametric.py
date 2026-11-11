"""Parametric study on the 2D section (fast): hole diameter x conductivity.

For each combination it computes:
  - UFV test (furnace at 5 °C/min up to 600 °C + 30 min): minimum and
    mean temperature at loading, maximum thermal tension and 99th percentile;
  - wall under the ISO 834 standard fire on one face: insulation time
    (140/180 °C criterion on the unexposed face) and thermal tension at 30 min.

Usage: python examples/04_section_parametric.py [--fast]
"""
import argparse
import itertools
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import studies
from mef_masonry.parameters import Geometry, Test
from mef_masonry.thermal import Materials
from mef_masonry.output import save_csv

ap = argparse.ArgumentParser()
ap.add_argument("--fast", action="store_true", help="coarser mesh and time step")
ap.add_argument("--diameters", type=float, nargs="+", default=[0, 60, 70, 84, 95], help="mm")
ap.add_argument("--k", type=float, nargs="+", default=[0.40, 0.52, 0.70], help="W/m.K")
ap.add_argument("--T_target", type=float, default=600.0)
a = ap.parse_args()
h, dt = (0.008, 60.0) if a.fast else (0.005, 30.0)

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT_DIR, exist_ok=True)
geo = Geometry()
lines, res = [], {}
t0 = time.time()
for D, k in itertools.product(a.diameters, a.k):
    M = Materials(); M.ceramic.k = k
    f = studies.furnace_section(M, Test(T_target=a.T_target), geo, h=h, dt=dt, diameter=D / 1000)
    p = studies.wall_section_iso834(M, geo, h=h, dt=dt, diameter=D / 1000, t_max_min=480)
    res[(D, k)] = (f, p)
    lines.append([D, k, f"{f['net_area_ratio']:.3f}", f"{f['T_min']:.1f}", f"{f['T_mean']:.1f}",
                   f"{f['s1_max_MPa']:.2f}", f"{f['s1_p99_MPa']:.2f}", f"{p['t_insulation_min']:.1f}",
                   f"{p.get('s1_p99_MPa', np.nan):.2f}"])
    print(f"D = {D:3.0f} mm, k = {k:.2f}: T_min {f['T_min']:5.0f} °C, tension p99 {f['s1_p99_MPa']:.2f} MPa, "
          f"insulation ISO 834 {p['t_insulation_min']:.0f} min  ({time.time() - t0:.0f} s)")

hdr = ["D_hole_mm", "k_W_mK", "net_area_ratio", "furnace_T_min_C", "furnace_T_mean_C",
       "furnace_tension_max_MPa", "furnace_tension_p99_MPa", "ISO834_insulation_min", "ISO834_tension_p99_30min_MPa"]
save_csv(os.path.join(OUT_DIR, "section_parametric.csv"), lines, hdr)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 2, figsize=(10, 7))
panels = [(axs[0, 0], lambda f, p: f["T_min"], f"Coldest point at loading (test {a.T_target:.0f} °C), °C"),
           (axs[0, 1], lambda f, p: f["s1_p99_MPa"], "Thermal tension, 99th percentile (test), MPa"),
           (axs[1, 0], lambda f, p: p["t_insulation_min"], "Insulation under ISO 834 (wall), min"),
           (axs[1, 1], lambda f, p: f["net_area_ratio"], "Net / gross area")]
colours = plt.cm.Blues(np.linspace(0.45, 0.95, len(a.k)))
for ax, fn, tit in panels:
    for c, k in zip(colours, a.k):
        ax.plot(a.diameters, [fn(*res[(D, k)]) for D in a.diameters], "o-", color=c, label=f"k = {k:.2f}")
    ax.set_title(tit, fontsize=9); ax.set_xlabel("hole diameter (mm)"); ax.grid(alpha=0.3)
axs[0, 0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "section_parametric.png"), dpi=130)
print("Results in", os.path.abspath(OUT_DIR), f"({time.time() - t0:.0f} s)")
