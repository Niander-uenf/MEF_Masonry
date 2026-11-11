"""Simulates the UFV tests on the 3D prism: heating at 5 °C/min, a 30 min
dwell and the maximum load of the test applied at the end of the dwell.

Usage: python examples/03_ufv_test_3d.py [--temperatures 600 800 1000] [--h 8] [--dz 12]
If files data/thermocouples_<T>.csv exist (columns: time_min, T_C),
they are overlaid on the chart for validation.
"""
import argparse
import csv
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, studies
from mef_masonry.parameters import Geometry, Test, LOADS_UFV_KGF
from mef_masonry.thermal import Materials
from mef_masonry.output import save_csv

ap = argparse.ArgumentParser()
ap.add_argument("--temperatures", type=float, nargs="+", default=[600, 800, 1000])
ap.add_argument("--h", type=float, default=8.0, help="element size in the section, mm")
ap.add_argument("--dz", type=float, default=12.0, help="vertical step, mm")
ap.add_argument("--dt", type=float, default=30.0, help="time step, s")
ap.add_argument("--hc", type=float, default=10.0, help="convection in the furnace, W/m².K")
ap.add_argument("--emissivity", type=float, default=0.8)
ap.add_argument("--no_vtk", action="store_true")
a = ap.parse_args()

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT_DIR = os.path.join(BASE, "results")
os.makedirs(OUT_DIR, exist_ok=True)
geo = Geometry()
M = Materials()
m = msh.prism_mesh(geo, h=a.h / 1000, dz=a.dz / 1000, t_joint=M.joint.thickness)
print("Mesh:", msh.summary(m))
print(f"Joint E = {M.joint.E / 1e6:.1f} MPa (recalibrate with example 02 if the mesh changes)\n")

all_, histories = [], {}
for T_target in a.temperatures:
    t0 = time.time()
    tst = Test(T_target=T_target, load_kgf=LOADS_UFV_KGF.get(int(T_target), 10000.0),
                 h_convection=a.hc, emissivity=a.emissivity)
    r = studies.simulate_test_3d(m, M, tst, geo, dt=a.dt,
                                  vtk_folder=None if a.no_vtk else os.path.join(OUT_DIR, "vtk"))
    ind = r["indicators"]
    all_.append(ind)
    histories[T_target] = r
    print(f"Test at {T_target:.0f} °C (load {tst.load_kgf:.1f} kgf), {time.time() - t0:.0f} s")
    print(f"  temperature at loading: surface {ind['T_surface']:.0f} °C, "
          f"centre of the web {ind['T_web']:.0f} °C, "
          f"mean of the middle third {ind['T_mean_ceramic_third']:.0f} °C")
    print(f"  thermal expansion between platens {ind['expansion_between_plates_mm']:.3f} mm; "
          f"shortening under the load {ind['load_shortening_mm']:.3f} mm")
    print(f"  maximum tension in the ceramic: thermal {ind['tension_max_thermal_MPa']:.2f} MPa "
          f"(p99 {ind['tension_p99_thermal_MPa']:.2f}); total {ind['tension_max_total_MPa']:.2f} MPa "
          f"(p99 {ind['tension_p99_total_MPa']:.2f}); central courses {ind['tension_max_total_central_courses_MPa']:.2f} MPa")
    print(f"  ceramic with tension > 1 MPa: {100 * ind['ceramic_fraction_tension_above_1MPa']:.0f}% of the volume\n")

    rt = r["thermal"]
    lines = [[f"{t / 60:.2f}", f"{g:.1f}"] + [f"{rt.probes[k][i]:.1f}" for k in r["probes"]]
              for i, (t, g) in enumerate(zip(rt.t, rt.T_gas))]
    save_csv(os.path.join(OUT_DIR, f"history_{int(T_target)}C.csv"), lines,
               ["time_min", "furnace_C"] + r["probes"])

keys = list(all_[0].keys())
save_csv(os.path.join(OUT_DIR, "ufv_test_3d.csv"), [[str(i[k]) for k in keys] for i in all_], keys)

# figure: temperatures during the test, with thermocouples if available
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, len(a.temperatures), figsize=(5 * len(a.temperatures), 3.8), squeeze=False)
for ax, T_target in zip(axs[0], a.temperatures):
    r = histories[T_target]; rt = r["thermal"]
    ax.plot(rt.t / 60, rt.T_gas, "k--", lw=1, label="furnace")
    for k, sty in zip(r["probes"], ("-", "-", ":")):
        ax.plot(rt.t / 60, rt.probes[k], sty, lw=1.5, label=k)
    fname = os.path.join(BASE, "data", f"thermocouples_{int(T_target)}.csv")
    if os.path.exists(fname):
        d = np.genfromtxt(fname, delimiter=",", names=True)
        ax.plot(d[d.dtype.names[0]], d[d.dtype.names[1]], "o", ms=3, label="thermocouple (test)")
    ax.axvline(r["t_end_ramp"] / 60, color="0.6", lw=0.8)
    ax.set_title(f"Test at {T_target:.0f} °C"); ax.set_xlabel("time (min)"); ax.set_ylabel("°C")
    ax.grid(alpha=0.3)
axs[0][0].legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "ufv_test_temperatures.png"), dpi=130)
print("Results in", os.path.abspath(OUT_DIR))
