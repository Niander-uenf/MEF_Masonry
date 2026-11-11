"""UFV test with the non-linear model: damage in the ceramic and dry joints with
friction. Heats the prism field by field (5 °C/min ramp and dwell) and then
compresses it with the top platen under displacement control, beyond the peak.

The bed joints have initial gaps (irregularity of the faces), with the
standard deviation calibrated from the prism modulus at 30 °C (script 08; default in
FrictionJoint). With --gap 0, contact is perfect and the joint stiffness is
recalibrated from the prism E.

Strengths: those of the literature (CeramicDamage) by default; with --calibrated, the
effective compressive strength that reproduces the peak load at 30 °C
(CALIBRATION_30C, script 08).

Usage: python examples/07_ufv_test_nonlinear.py [--temperatures 30 600 800 1000]
       [--calibrated] [--fc 5.76] [--ft 0.66] [--gap 0.028] [--seed 1]
       [--h 12] [--dz 20] [--delta_max 0.8] [--steps 40]
Time: about 10 min per temperature with the default mesh.
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, studies
from mef_masonry.parameters import Geometry, Test, CeramicDamage, FrictionJoint, LOADS_UFV_KGF, CALIBRATION_30C
from mef_masonry.thermal import Materials
from mef_masonry.nonlinear import NonlinearModel
from mef_masonry.output import save_csv, save_vtk

ap = argparse.ArgumentParser()
ap.add_argument("--temperatures", type=float, nargs="+", default=[30, 600, 800, 1000])
ap.add_argument("--h", type=float, default=12.0)
ap.add_argument("--dz", type=float, default=20.0)
ap.add_argument("--ft", type=float, default=CeramicDamage().ft / 1e6, help="tensile strength of the ceramic, MPa")
ap.add_argument("--fc", type=float, default=None, help="compressive strength of the ceramic, MPa")
ap.add_argument("--calibrated", action="store_true", help="uses the effective fc calibrated at 30 °C")
ap.add_argument("--E_prism", type=float, default=1100.0, help="target prism modulus at 30 °C, MPa")
ap.add_argument("--delta_max", type=float, default=0.8, help="maximum platen displacement, mm")
ap.add_argument("--steps", type=int, default=40)
ap.add_argument("--gap", type=float, default=FrictionJoint().gap_std * 1e3,
                help="standard deviation of the bed-joint gaps, mm (0 = perfect contact)")
ap.add_argument("--E_contact", type=float, default=(FrictionJoint().E_contact or 0) / 1e6,
                help="E of the joints where there is contact, MPa")
ap.add_argument("--seed", type=int, default=FrictionJoint().seed, help="seed of the gap field")
ap.add_argument("--label", default="")
a = ap.parse_args()

BASE = os.path.join(os.path.dirname(__file__), "..")
OUT_DIR = os.path.join(BASE, "results")
os.makedirs(os.path.join(OUT_DIR, "vtk"), exist_ok=True)
geo = Geometry()
M = Materials()
fc = a.fc * 1e6 if a.fc else (CALIBRATION_30C["fc"] if a.calibrated else CeramicDamage().fc)
Pd = CeramicDamage(ft=a.ft * 1e6, fc=fc)
Pj = FrictionJoint(gap_std=a.gap * 1e-3, seed=a.seed,
                 E_contact=(a.E_contact * 1e6 or None) if a.gap > 0 else None)
P30 = LOADS_UFV_KGF[30] * 9.80665
m = msh.prism_mesh(geo, h=a.h / 1000, dz=a.dz / 1000, t_joint=M.joint.thickness)
print("Mesh:", msh.summary(m))
print(f"Ceramic: ft = {Pd.ft / 1e6:.2f} MPa, fc = {Pd.fc / 1e6:.2f} MPa, G_f = {Pd.G_f:.1f} N/m, "
      f"G_c = {Pd.G_c:.0f} N/m; joints: mu = {Pj.mu}")

t0 = time.time()
if Pj.gap_std <= 0:
    mod0 = NonlinearModel(m, M, Pd, Pj, verbose=False)
    print("Perfect contact. Calibrating the joint stiffness for E_prism =", a.E_prism, "MPa")
    Ej = studies.calibrate_joint_nonlinear(mod0, a.E_prism, geo)
    M.joint.E = Ej; M.capping.E = Ej
    del mod0
    print(f"  E_joint = {Ej / 1e6:.1f} MPa ({time.time() - t0:.0f} s)\n")
else:
    print(f"Joints with gaps: std {Pj.gap_std * 1e3:.3f} mm, initial contact {100 * Pj.initial_contact:.0f}%, "
          f"correlation {Pj.corr_length * 1e3:.0f} mm, E_contact {Pj.E_contact / 1e6:.0f} MPa, seed {Pj.seed}\n")

summary, pd_curves = [], {}
for T_target in a.temperatures:
    t0 = time.time()
    print(f"Test at {T_target:.0f} °C")
    fields = studies.heating_fields(m, M, Test(T_target=T_target)) if T_target > M.ceramic.T_ref + 1 else None
    mod, (deltas, loads) = studies.nonlinear_test(m, M, Pd, Pj, fields=fields, delta_max=a.delta_max / 1000,
                                                      n_steps=a.steps, geo=geo, verbose=True)
    vc = mod.vol[mod.ic]
    frac_aq = float(vc[mod.d_heating > 0.5].sum() / vc.sum())
    d_heat = np.full(len(m.conn3), np.nan); d_heat[mod.id_glob[mod.ic]] = mod.d_heating
    i = int(np.argmax(loads))
    P_exp = LOADS_UFV_KGF.get(int(T_target), np.nan) * 9.80665
    summary.append([T_target, f"{loads[i] / 1e3:.1f}", f"{P_exp / 1e3:.1f}", f"{(loads[i] / P_exp - 1) * 100:+.0f}%",
                   f"{deltas[i] * 1e3:.3f}", f"{100 * frac_aq:.2f}", f"{studies.chord_modulus(mod, P30, geo):.0f}",
                   sum(1 for p in mod.history if p.phase == "load" and not p.converged), f"{time.time() - t0:.0f}"])
    pd_curves[T_target] = (deltas, loads)
    hist = [p for p in mod.history if p.phase == "load"]
    save_csv(os.path.join(OUT_DIR, f"nl_load_displacement_{int(T_target)}C{a.label}.csv"),
               [[f"{p.label * 1e3:.4f}", f"{p.load_N / 1e3:.3f}", f"{-p.eps_ext * 1e3:.4f}", f"{100 * p.contact:.1f}",
                 f"{100 * p.damaged_vol:.2f}"] for p in hist],
               ["delta_mm", "load_kN", "middle_third_strain_permil", "joint_contact_pct", "ceramica_d_maior_0.5_pct"])
    gap_el = np.full(len(m.conn3), np.nan); gap_el[mod.id_glob[mod.ij]] = mod.gap * 1e3
    save_vtk(os.path.join(OUT_DIR, "vtk", f"nl_prism_{int(T_target)}C{a.label}.vtk"), m,
               node_data={"T_C": mod.T, "u_mm": mod.u.reshape(-1, 3) * 1e3},
               elem_data={"material": m.mat.astype(float), "final_damage": mod.global_damage(),
                           "damage_after_heating": d_heat, "initial_gap_mm": gap_el})
    print(f"  peak {loads[i] / 1e3:.1f} kN at delta = {deltas[i] * 1e3:.3f} mm; test {P_exp / 1e3:.1f} kN; "
          f"ceramic with d > 0.5 after heating: {100 * frac_aq:.2f}%  ({time.time() - t0:.0f} s)\n")

hdr = ["T_target_C", "peak_model_kN", "peak_test_kN", "difference", "delta_at_peak_mm",
       "damaged_ceramic_after_heating_pct", "E_chord_5_30_MPa", "increments_not_converged", "time_s"]
save_csv(os.path.join(OUT_DIR, f"nl_summary{a.label}.csv"), summary, hdr)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.5, 4.3))
colours = plt.cm.Blues(np.linspace(0.35, 0.95, len(a.temperatures)))
for c, T in zip(colours, a.temperatures):
    d, P = pd_curves[T]
    ax.plot(d * 1e3, P / 1e3, "-", color=c, lw=1.8, label=f"model, {T:.0f} °C")
    Pe = LOADS_UFV_KGF.get(int(T))
    if Pe:
        ax.axhline(Pe * 9.80665 / 1e3, color=c, ls=":", lw=1)
ax.set_xlabel("platen displacement after heating (mm)"); ax.set_ylabel("load (kN)")
ax.set_title("BCEPQ prism: non-linear model (lines) and peak load of the test (dotted)", fontsize=9)
ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, f"nl_load_displacement{a.label}.png"), dpi=130)
print("Summary:")
for l in summary:
    print("  ", dict(zip(hdr, l)))
