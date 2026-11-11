"""Calibration of the non-linear model on the prism at 30 °C, in two independent
steps:

  1. Stiffness: with the contact as stiff as the ceramic, finds the
     standard deviation of the bed-joint gaps that reproduces the chord
     modulus of the prism between 5% and 30% of the peak load (E = 1,100 MPa, BPQ), as the
     mean of three seeds of the random field (the stiffness varies by about 20%
     between seeds). Uses only the initial part of the curve.
  2. Strength: with that gap, finds the effective compressive strength
     of the ceramic (fc) that reproduces the peak load of the test at 30 °C.

The calibrated values go into mef_masonry/parameters.py (FrictionJoint.gap_std
and CALIBRATION_30C) and serve, unchanged, to predict the heated tests
(07_ufv_test_nonlinear.py --calibrated).

Usage: python examples/08_calibrate_imperfection.py [--E_prism 1100] [--P30 98.7]
       [--gaps 0.024 0.032] [--seeds 1 2 3] [--h 12] [--dz 20]
Time: about 15 min with the default mesh (2 cores).
"""
import argparse
import os
import sys
import time
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, studies
from mef_masonry.parameters import Geometry, CeramicDamage, FrictionJoint, LOADS_UFV_KGF
from mef_masonry.thermal import Materials
from mef_masonry.output import save_csv

ap = argparse.ArgumentParser()
ap.add_argument("--E_prism", type=float, default=1100.0, help="chord modulus 5-30%% of the prism, MPa")
ap.add_argument("--P30", type=float, default=LOADS_UFV_KGF[30] * 9.80665 / 1e3, help="peak load at 30 °C, kN")
ap.add_argument("--gaps", type=float, nargs="+", default=[0.024, 0.032], help="standard deviations to test, mm")
ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3], help="seeds of the gap field")
ap.add_argument("--h", type=float, default=12.0)
ap.add_argument("--dz", type=float, default=20.0)
a = ap.parse_args()

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT_DIR, exist_ok=True)
geo, M, Pd = Geometry(), Materials(), CeramicDamage()
m = msh.prism_mesh(geo, h=a.h / 1000, dz=a.dz / 1000, t_joint=M.joint.thickness)
P30 = a.P30 * 1e3
Pj = FrictionJoint(seed=a.seeds[0])
print(f"Mesh: {msh.summary(m)}")
print(f"Targets: chord E (5-30%) = {a.E_prism:.0f} MPa; peak load = {a.P30:.1f} kN; "
      f"E_contact = {Pj.E_contact / 1e6:.0f} MPa\n")

t0 = time.time()
print("1. Joint gap from the prism modulus")
sg, pts_E = studies.calibrate_gap(m, M, Pd, Pj, P30, a.E_prism, gaps=[f * 1e-3 for f in a.gaps],
                                   seeds=a.seeds, geo=geo)
print(f"   gap_std = {sg * 1e3:.4f} mm ({time.time() - t0:.0f} s)\n")
sg = round(sg * 1e6) / 1e6          # 0,001 mm

t0 = time.time()
print("2. Effective compressive strength of the ceramic from the peak load at 30 °C")
Pj = replace(Pj, gap_std=sg)
fc, Pmax, pts_fc = studies.calibrate_fc(m, M, Pd, Pj, P30, geo=geo)
print(f"   effective fc = {fc / 1e6:.2f} MPa (literature: {Pd.fc / 1e6:.2f} MPa); "
      f"peak load {Pmax / 1e3:.1f} kN ({time.time() - t0:.0f} s)\n")

lines = [["stage 1", f"{s * 1e3:.4f}", "", f"{E:.1f}", ""] for s, E, _ in pts_E]
lines += [["stage 2", f"{sg * 1e3:.4f}", f"{f / 1e6:.3f}", "", f"{P / 1e3:.2f}"] for f, P in pts_fc]
lines.append(["calibrated", f"{sg * 1e3:.4f}", f"{fc / 1e6:.3f}", f"{a.E_prism:.0f}", f"{Pmax / 1e3:.2f}"])
save_csv(os.path.join(OUT_DIR, "calibration_30C.csv"), lines,
           ["stage", "gap_std_mm", "fc_MPa", "E_chord_MPa", "max_load_kN"])
print("Put in mef_masonry/parameters.py:")
print(f"  FrictionJoint.gap_std = {sg:.3e}")
print(f"  CALIBRATION_30C = {{'fc': {fc:.3e}}}")
