"""Convergence studies of the FEM (Table 5 of the paper).

  mesh_nl  : prism at 30 °C (or --T) with the calibrated non-linear model, on
             several meshes; peak load, displacement at the peak and chord E.
             The calibrated parameters (gap, effective fc) stay fixed.
  step_nl  : the same analysis on the default mesh with different load increments.
  tol_nl   : the same analysis with different convergence tolerances.
  thermal  : heated test (thermal + thermoelastic), on several meshes and
             time steps; temperatures at the moment of loading and thermal tension.
             The joint stiffness is recalibrated on each mesh (E_prism = 1,100 MPa).

Usage:
  python examples/09_convergence_study.py mesh_nl --meshes 12x20 10x15 8x12
  python examples/09_convergence_study.py mesh_nl --element H8 --meshes 12x20   (or --element W15 --meshes 16x30)
  python examples/09_convergence_study.py step_nl --steps 20 40 80
  python examples/09_convergence_study.py tol_nl
  python examples/09_convergence_study.py thermal --meshes 12x20 8x12 --dts 15 30 60 --T 600
Each study appends lines to results/convergence_<study>.csv.
"""
import argparse
import csv
import os
import sys
import time
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, studies
from mef_masonry.parameters import Geometry, Test, CeramicDamage, FrictionJoint, LOADS_UFV_KGF, CALIBRATION_30C
from mef_masonry.thermal import Materials

ap = argparse.ArgumentParser()
ap.add_argument("study", choices=["mesh_nl", "step_nl", "tol_nl", "thermal"])
ap.add_argument("--meshes", nargs="+", default=["12x20"], help="h x dz in mm, for example 12x20 8x12")
ap.add_argument("--steps", type=int, nargs="+", default=[40])
ap.add_argument("--dts", type=float, nargs="+", default=[30.0], help="time steps of the thermal analysis, s")
ap.add_argument("--T", type=float, default=30.0, help="target temperature of the test, °C")
ap.add_argument("--delta_max", type=float, default=0.8, help="mm")
ap.add_argument("--element", default="W6", choices=["W6", "H8", "W15"], help="3D element")
a = ap.parse_args()

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT_DIR = os.path.join(BASE, "results")
os.makedirs(OUT_DIR, exist_ok=True)
fname = os.path.join(OUT_DIR, f"convergence_{a.study}.csv")
geo = Geometry()
P30 = LOADS_UFV_KGF[30] * 9.80665


def write_file(line: dict):
    new = not os.path.exists(fname)
    with open(fname, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(line))
        if new:
            w.writeheader()
        w.writerow(line)
    print("  ->", line, flush=True)


def mesh(txt):
    h, dz = (float(v) / 1000 for v in txt.lower().split("x"))
    return msh.prism_mesh(geo, h=h, dz=dz, t_joint=Materials().joint.thickness, element=a.element), h, dz


def nl_test(m, n_steps=40, tol_u=1e-4, tol_d=5e-3):
    M = Materials()
    Pd = CeramicDamage(fc=CALIBRATION_30C["fc"])
    Pj = FrictionJoint()
    fields = studies.heating_fields(m, M, Test(T_target=a.T)) if a.T > M.ceramic.T_ref + 1 else None
    t0 = time.time()
    mod, (d, P) = studies.nonlinear_test(m, M, Pd, Pj, fields=fields, delta_max=a.delta_max / 1000,
                                            n_steps=n_steps, geo=geo, verbose=True, tol_u=tol_u, tol_d=tol_d)
    i = int(np.argmax(P))
    hist = [p for p in mod.history if p.phase == "load"]
    return {"T_C": a.T, "nodes": m.n_nodes, "free_dofs_n": int(np.count_nonzero(~mod.base_fixed)),
            "ceramic_elements": len(mod.ic), "peak_kN": round(P[i] / 1e3, 2), "delta_peak_mm": round(d[i] * 1e3, 3),
            "E_chord_MPa": round(studies.chord_modulus(mod, P30, geo), 0),
            "contact_at_peak_pct": round(100 * hist[i - 1].contact, 1) if i >= 1 else "",
            "increments": len(hist), "not_converged": sum(1 for p in hist if not p.converged),
            "time_min": round((time.time() - t0) / 60, 1)}


if a.study == "mesh_nl":
    for txt in a.meshes:
        m, h, dz = mesh(txt)
        print(f"Mesh {txt}: {msh.summary(m)}", flush=True)
        write_file({"element": a.element, "mesh_mm": txt, **nl_test(m)})

elif a.study == "step_nl":
    m, h, dz = mesh(a.meshes[0])
    for n in a.steps:
        print(f"{n} increments up to {a.delta_max} mm", flush=True)
        write_file({"mesh_mm": a.meshes[0], "increment_mm": round(a.delta_max / n, 4), **nl_test(m, n_steps=n)})

elif a.study == "tol_nl":
    m, h, dz = mesh(a.meshes[0])
    for tu, td in ((1e-4, 5e-3), (2e-5, 1e-3)):
        print(f"tol_u = {tu}, tol_d = {td}", flush=True)
        write_file({"mesh_mm": a.meshes[0], "tol_u": tu, "tol_d": td, **nl_test(m, tol_u=tu, tol_d=td)})

elif a.study == "thermal":
    for txt in a.meshes:
        m, h, dz = mesh(txt)
        M = Materials()
        Ej, _ = studies.calibrate_joint(m, M, 1100.0, geo, verbose=False)
        for dt in a.dts:
            t0 = time.time()
            ind = studies.simulate_test_3d(m, M, Test(T_target=a.T), geo, dt=dt)["indicators"]
            write_file({"mesh_mm": txt, "dt_s": dt, "nodes": m.n_nodes, "E_joint_MPa": round(Ej / 1e6, 1),
                    "T_surface_C": round(ind["T_surface"], 1), "T_mean_third_C": round(ind["T_mean_ceramic_third"], 1),
                    "T_web_C": round(ind["T_web"], 1),
                    "tension_max_thermal_MPa": round(ind["tension_max_thermal_MPa"], 3),
                    "tension_p99_thermal_MPa": round(ind["tension_p99_thermal_MPa"], 3),
                    "time_min": round((time.time() - t0) / 60, 1)})
