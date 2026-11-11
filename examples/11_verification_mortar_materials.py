"""Verification of phase 1: temperature-dependent specific heat, mortar
joints and library.

V1  Variable heat capacity c(T) (peak between 100 and 200 °C): body of
    very high conductivity (Biot << 1) heated by convection; compared with the
    ODE  rho c(T) A dT/dt = h P (Tg - T), integrated with a tolerance of 1e-10.
V2  Column of bricks and mortar (linear 3D) with equal nu/E in both materials
    (same lateral strain, hence uniform uniaxial stress): shortening =
    sum of h_i / E_i (series association, exact).
V3  Same column with a larger nu/E in the mortar (0.30/1 GPa against 0.15/3 GPa): the
    mortar tries to expand more than the brick and pulls it laterally
    (Hilsdorf mechanism): lateral tensile stress in the brick next to the joint.
V4  Non-linear column with equal nu/E: without confinement, failure comes from the
    (weaker) mortar, with peak load = fa × area. With equal Poisson's ratios and
    different E, the mortar is confined and the load rises (recorded separately).
V5  Properties per material: conductivity, capacity and thermal strain
    of the mortar reach the solvers.
V6  Library: the items that run pass validation; those of phase 2 onwards
    are rejected.

Usage: python examples/11_verification_mortar_materials.py
Writes results/verification_mortar_materials.csv.
"""
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, thermal as th, mechanical as mech, project as pj, library as lib, curves
from mef_masonry.mesh import CERAMIC, MORTAR
from mef_masonry.parameters import Geometry, CeramicDamage
from mef_masonry.thermal import Materials, Boundary
from mef_masonry.nonlinear import NonlinearModel
from mef_masonry.output import save_csv

lines = []


def record(case, model, ref, tol_rel, unit_str=""):
    error = abs(model - ref) / max(abs(ref), 1e-30)
    ok = error <= tol_rel
    lines.append([case, round(model, 6), round(ref, 6), unit_str, f"{100 * error:.3f}%", f"{100 * tol_rel:.1f}%",
                   "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: model {model:.6g} {unit_str}, reference {ref:.6g} ({100 * error:.3f}%)")


# ---------------------------------------------------------------------------- V1
a = 0.05
m2 = msh.rectangle_mesh(a, a, a / 10)
M = Materials()
M.ceramic.k, M.ceramic.rho, M.ceramic.c = 2000.0, 1000.0, 1000.0
tab_c = [(20, 1000.0), (100, 1000.0), (115, 3000.0), (200, 1000.0), (1200, 1000.0)]
M.ceramic.table_c = tab_c
M.ceramic.T_ref = 20.0
h, Tg = 25.0, 500.0
comp = th.boundary_lengths_2d(m2)
ctr = [Boundary(comp, curves.constant(Tg), h, 0.0)]
inst = [600.0, 1200.0, 2400.0]
r = th.solve_2d(m2, M, ctr, 2400.0, dt=1.0, T_initial=20.0, field_instants=inst,
                   mat=np.full(len(m2.tri), CERAMIC, dtype=np.int8))
A, P = a * a, 4 * a
cT = lambda T: np.interp(T, [x[0] for x in tab_c], [x[1] for x in tab_c])  # noqa: E731
sol = solve_ivp(lambda t, T: h * P * (Tg - T) / (1000.0 * cT(T) * A), (0, 2400), [20.0], rtol=1e-10, atol=1e-10,
                dense_output=True, max_step=1.0)
for t in inst:
    record(f"V1 temperature with c(T) at {t / 60:.0f} min", float(r.fields[t].mean()), float(sol.sol(t)[0]), 0.005,
              "°C")

# ---------------------------------------------------------------------------- V2 and V3
geo = Geometry(unit_length=0.23, width=0.10, course_height=0.07, hole_diameter=0.0, holes_per_unit=0,
                units_per_course=1, n_courses=3, bond="stack_bond", joint_type="mortar")
mm = msh.prism_mesh(geo, h=0.01, dz=0.005, t_joint=0.01, symmetry=True, plates=False, capping=False)
Mm = Materials()
Mm.ceramic.E_unit_gross, Mm.ceramic.net_area_ratio, Mm.ceramic.nu = 3000e6, 1.0, 0.15
Mm.mortar.E_unit_gross, Mm.mortar.nu = 1000e6, 0.05
P_N = 50e3
R = mech.solve_3d(mm, Mm, load_N=P_N, include_thermal=False)["load"]
H = geo.prism_height
d_fem = -mech.shortening_between_plates(mm, R, H)
hs = {CERAMIC: 0.0, MORTAR: 0.0}
for c in mm.layers:
    hs[MORTAR if c.kind == "joint" else CERAMIC] += c.z1 - c.z0
Aarea = geo.length * geo.width
d_ref = P_N / Aarea * (hs[CERAMIC] / 3000e6 + hs[MORTAR] / 1000e6)
record("V2 shortening of the brick + mortar column (series)", d_fem * 1e6, d_ref * 1e6, 0.002, "µm")

Mm.mortar.nu, Mm.ceramic.nu = 0.30, 0.15
R3 = mech.solve_3d(mm, Mm, load_N=P_N, include_thermal=False)["load"]
active = np.nonzero(R3.active)[0]
mat_at = mm.mat[active]
zc = mm.xyz[mm.conn3[active]].mean(1)[:, 2]
zj = [0.5 * (c.z0 + c.z1) for c in mm.layers if c.kind == "joint"][0]
near = (mat_at == CERAMIC) & (np.abs(zc - zj) < 0.012)
sx_brick = float(np.median(R3.sig[near, 0]))
sx_mor = float(np.median(R3.sig[mat_at == MORTAR, 0]))
ok3 = sx_brick > 0 and sx_mor < 0
lines.append(["V3 brick in tension and mortar in lateral compression (Hilsdorf)",
               round(sx_brick / 1e6, 4), round(sx_mor / 1e6, 4), "MPa (σx brick / mortar)", "—", "signs",
               "ok" if ok3 else "FAILED"])
print(f"{'ok ' if ok3 else 'XX '} V3 σx in the brick next to the joint {sx_brick / 1e6:+.4f} MPa; in the mortar {sx_mor / 1e6:+.4f} MPa")

# ---------------------------------------------------------------------------- V4
Mn = Materials()
mn = msh.prism_mesh(geo, h=0.02, dz=0.01, t_joint=0.01, symmetry=True, plates=False, capping=False)
Pd = CeramicDamage(ft=5e6, fc=20e6)
Pa = CeramicDamage(ft=2e6, fc=4e6, du_t=0.065e-3, du_c=1.6e-3)
P4max = {}
for name, nu_b, nu_m in (("equal nu/E", 0.15, 0.05), ("equal Poisson", 0.18, 0.18)):
    Mn = Materials()
    Mn.ceramic.E_unit_gross, Mn.ceramic.net_area_ratio, Mn.ceramic.nu = 3000e6, 1.0, nu_b
    Mn.mortar.E_unit_gross, Mn.mortar.nu = 1000e6, nu_m
    mod = NonlinearModel(mn, Mn, Pd, None, verbose=False, mortar_damage=Pa)
    mod.heat({0.0: np.full(mn.n_nodes, Mn.ceramic.T_ref)}, pre_load_N=0.0)
    d4, P4 = mod.apply_load(1.0e-3, 50)
    P4max[name] = float(P4.max())
record("V4 peak load of the column = fa × area (equal nu/E)", P4max["equal nu/E"] / 1e3, 4e6 * Aarea / 1e3, 0.02,
          "kN")
ok4 = P4max["equal Poisson"] > P4max["equal nu/E"]
lines.append(["V4 confined mortar (equal Poisson's ratios, smaller E) resists more", round(P4max["equal Poisson"] / 1e3, 3),
               round(P4max["equal nu/E"] / 1e3, 3), "kN", "—", "largest", "ok" if ok4 else "FAILED"])
print(f"{'ok ' if ok4 else 'XX '} V4 with equal Poisson's ratios: {P4max['equal Poisson'] / 1e3:.1f} kN "
      f"(> {P4max['equal nu/E'] / 1e3:.1f} kN without confinement)")

# ---------------------------------------------------------------------------- V5
M5 = Materials()
M5.mortar.k, M5.mortar.rho, M5.mortar.c, M5.mortar.alpha = 1.15, 2000.0, 1000.0, 10e-6
M5.mortar.T_ref = 20.0
mat5 = np.array([MORTAR], dtype=np.int8)
record("V5 mortar conductivity in the thermal solver", float(M5.k(mat5, np.array([300.0]))[0]), 1.15, 1e-12,
          "W/m·K")
record("V5 volumetric capacity of the mortar", float(M5.rho_c(mat5)[0]) / 1e6, 2.0, 1e-12, "MJ/m³·K")
eps5 = mech._thermal_strain_nodes(mat5, np.array([[520.0] * 6]), M5)[0, 0]
record("V5 thermal strain of the mortar at 520 °C", eps5 * 1e3, 10e-6 * 500 * 1e3, 1e-12, "‰")

# ---------------------------------------------------------------------------- V6
valid = rejected = 0
for k, item in lib.UNITS.items():
    p = pj.complete(lib.apply(pj.default_project(), "unit", k))
    errors = pj.validate(p)
    if item.get("phase", 0) <= lib.CURRENT_PHASE:
        valid += not errors
        if errors:
            print("   ", k, errors)
    else:
        rejected += bool(errors)
n_ok = sum(1 for v in lib.UNITS.values() if v.get("phase", 0) <= lib.CURRENT_PHASE)
n_fut = len(lib.UNITS) - n_ok
record("V6 library items that run and pass validation", valid, n_ok, 0.0, "items")
record("V6 items of future phases rejected by validation", rejected, n_fut, 0.0, "items")

output = os.path.join(os.path.dirname(__file__), "..", "results", "verification_mortar_materials.csv")
save_csv(output, lines, ["case", "model", "reference", "unit", "error", "tolerance", "status"])
nf = sum(1 for l in lines if l[-1] != "ok")
print(f"\n{len(lines) - nf}/{len(lines)} cases within tolerance. Written to {os.path.abspath(output)}")
sys.exit(1 if nf else 0)
