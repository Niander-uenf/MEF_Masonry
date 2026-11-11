"""Calibrates the stiffness of the dry joints so that the prism reproduces the measured modulus,
and compares experimental x numerical x formulas (series and NBR 16868).

Usage: python examples/02_calibrate_joints.py [--E_target 1100] [--h 8] [--dz 12]
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, studies
from mef_masonry.parameters import Geometry, E_DRY_PRISM_MPA, LOADS_UFV_KGF
from mef_masonry.thermal import Materials
from mef_masonry.output import save_csv

ap = argparse.ArgumentParser()
ap.add_argument("--E_target", type=float, default=1100.0, help="modulus of the dry prism, MPa (BPQ: 1,001 to 1,216)")
ap.add_argument("--h", type=float, default=8.0, help="element size in the section, mm")
ap.add_argument("--dz", type=float, default=12.0, help="vertical step, mm")
a = ap.parse_args()

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT_DIR, exist_ok=True)
geo = Geometry()
M = Materials()
t0 = time.time()
m = msh.prism_mesh(geo, h=a.h / 1000, dz=a.dz / 1000, t_joint=M.joint.thickness)
print("Mesh:", msh.summary(m))

print(f"\nCalibrating the joint E for E_prism = {a.E_target:.0f} MPa (middle third, gross area)")
Ej, hist = studies.calibrate_joint(m, M, a.E_target, geo)
r = studies.prism_modulus(m, M, geo=geo)
base = (r["z2"] - r["z1"]) * 1000
Ej_serie = studies.series_estimate(M, a.E_target, base, 2)

fpk = LOADS_UFV_KGF[30] * 9.80665 / (geo.length * geo.width) / 1e6
lines = [
    ["E of the unit (gross area)", "experimental (BPQ, n = 6)", f"{M.ceramic.E_unit_gross / 1e6:.0f}"],
    ["E of the ceramic material (net area)", "derived: E_unit / 0.717", f"{M.ceramic.E_material / 1e6:.0f}"],
]
for k, v in E_DRY_PRISM_MPA.items():
    lines.append([f"E of the dry prism {k}", "experimental (BPQ, n = 4)", f"{v:.0f}"])
lines += [
    ["E of the prism = 600 fpk", f"NBR 16868-1, fpk = {fpk:.2f} MPa (UFV, 30 °C)", f"{600 * fpk:.0f}"],
    ["E of the joint: series association", f"theoretical, base {base:.0f} mm with 2 joints", f"{Ej_serie / 1e6:.1f}"],
    ["E of the joint: calibrated in the FEM", "numerical 3D", f"{Ej / 1e6:.1f}"],
    ["E of the prism in the middle third", "numerical 3D, calibrated joint", f"{r['E_extensometer']:.0f}"],
    ["E of the prism between platens", "numerical 3D, calibrated joint", f"{r['E_plates']:.0f}"],
    ["Normal stiffness of the joint k_n = E_j/t", "numerical 3D (N/mm³)", f"{Ej / M.joint.thickness / 1e9:.1f}"],
]
print("\nSummary (MPa, unless stated otherwise):")
for l in lines:
    print(f"  {l[0]:<40s} {l[1]:<45s} {l[2]:>8s}")
save_csv(os.path.join(OUT_DIR, "joint_calibration.csv"), lines, ["quantity", "origin", "value"])
save_csv(os.path.join(OUT_DIR, "calibration_iterations.csv"), [[f"{x:.3f}", f"{y:.2f}"] for x, y in hist],
           ["E_joint_MPa", "E_prism_MPa"])
print(f"\nUse E_joint = {Ej / 1e6:.1f} MPa in the other examples (parameters.DryJoint.E). "
      f"Time: {time.time() - t0:.0f} s")
