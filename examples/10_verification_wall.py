"""Verification of the wall model (wall.py) against known solutions.

V1  Amplification of the sinusoidal imperfection (elastic pin-ended column):
    deflection = a0 (P/Pcr') / (1 - P/Pcr'), with Pcr' = pi² EI / H'² and H' shortened by P/EA
V2  Free thermal bowing (linear gradient): arc of radius t / (alpha dT)
V3  Eccentric load at both ends (secant formula): e [sec(pi/2 sqrt(P/Pcr)) - 1]
V4  Fibre section from the 2D mesh: net/gross area ratio and inertia against the geometry
V5  Short wall at 20 °C, centred load: capacity = prism strength (gross area)
V6  Short wall heated uniformly from 20 to 600 °C under 1% of the load: free expansion and
    capacity = prism strength at 600 °C
V7  Convergence with the number of elements (3 m wall): cold capacity and fire
    resistance under ISO 834 at 40% of the capacity

Usage: python examples/10_verification_wall.py [--fast]   (--fast skips V7)
Writes results/verification_wall.csv.
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, thermal as th, curves
from mef_masonry.parameters import Geometry
from mef_masonry.wall import FibreSection, LinearLaw, WallLaw, WallModel
from mef_masonry.wall_data import PRISM_CURVES_BCEPQ, H_REF, A_GROSS_REF
from mef_masonry.thermal import Materials, Boundary
from mef_masonry.output import save_csv

ap = argparse.ArgumentParser()
ap.add_argument("--fast", action="store_true")
a = ap.parse_args()

lines = []


def record(case, model, ref, tol_rel, unit_str=""):
    error = abs(model - ref) / abs(ref)
    ok = error <= tol_rel
    lines.append([case, round(model, 6), round(ref, 6), unit_str, f"{100 * error:.3f}%", f"{100 * tol_rel:.1f}%",
                   "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: model {model:.6g} {unit_str}, reference {ref:.6g} ({100 * error:.3f}%)")


zero = lambda T: 0.0 * np.asarray(T, float)  # noqa: E731
t, E, H = 0.14, 1e9, 3.0
rect = FibreSection.rectangle(t, 200)
EI, EA = E * rect.inertia, E * rect.area
Pcr = math.pi ** 2 * EI / H ** 2

# V1
a0 = H / 1000
for r in (0.5, 0.9):
    m = WallModel(H, rect, LinearLaw(E), zero, n_el=24, imperfection=a0)
    m.apply_load(r * Pcr, 0, 0, n_steps=20)
    Pc = Pcr / (1 - r * Pcr / EA) ** 2
    ref = a0 * (r * Pcr / Pc) / (1 - r * Pcr / Pc)
    record(f"V1 additional deflection with P = {r:.1f} Pcr", -m.mid_deflection * 1e3, ref * 1e3, 0.005, "mm")

# V2
alpha, dT = 6e-6, 800.0
m = WallModel(H, rect, LinearLaw(E), lambda T: alpha * (np.asarray(T) - 20.0), n_el=24, T_initial=20.0)
Tf = 20.0 + dT * (0.5 - rect.y / t)
d, info, _ = m.equilibrate(np.zeros(m.nd), Tf)
m.accept(d, Tf, info)
R = t / (alpha * dT)
record("V2 free thermal deflection", -m.mid_deflection * 1e3, (R - math.sqrt(R * R - H * H / 4)) * 1e3, 0.005, "mm")
Lm = H * (1 + alpha * dT / 2)                     # elongated mean axis, bent into an arc of radius R
chord = 2 * R * math.sin(Lm / (2 * R))
record("V2 top displacement (chord of the arc)", m.top_displacement * 1e3, (chord - H) * 1e3, 0.005, "mm")

# V3
e = 0.010
m = WallModel(H, rect, LinearLaw(E), zero, n_el=24)
m.apply_load(0.5 * Pcr, e, e, n_steps=20)
Pc = Pcr / (1 - 0.5 * Pcr / EA) ** 2
record("V3 deflection with an eccentricity of 10 mm, P = 0.5 Pcr", m.mid_deflection * 1e3,
          e * (1 / math.cos(math.pi / 2 * math.sqrt(0.5 * Pcr / Pc)) - 1) * 1e3, 0.005, "mm")

# V4
geo1 = Geometry(units_per_course=1.0)
m2 = msh.section_mesh(geo1, h=0.005, symmetry=False, bands=())
sec = FibreSection.from_mesh(m2)
Rf = geo1.hole_diameter / 2
nfu = len(geo1.hole_centres)
ratio = 1 - nfu * math.pi * Rf ** 2 / (geo1.length * t)
I_an = t ** 3 / 12 - nfu * math.pi * Rf ** 4 / 4 / geo1.length
record("V4 net/gross area ratio of the fibres", sec.area / t, ratio, 0.002)
record("V4 inertia per metre of the fibres", sec.inertia * 1e6, I_an * 1e6, 0.002, "×10⁻⁶ m⁴/m")

# V5 and V6
A_net_ref = A_GROSS_REF * (1 - 3 * math.pi * 0.042 ** 2 / A_GROSS_REF)
law = WallLaw.from_prism(PRISM_CURVES_BCEPQ, H_REF, A_net_ref)
M = Materials(); M.ceramic.T_ref = 20.0
f_gross = {T: max(c[1] for c in v) * 1e3 / A_GROSS_REF / 1e6 for T, v in PRISM_CURVES_BCEPQ.items()}
m = WallModel(0.3, sec, law, M.ceramic.thermal_strain, n_el=4, T_initial=20.0)
P = m.capacity(0, 0, tol=1e-4)
record("V5 capacity of the short wall at 20 °C (gross area)", P / t / 1e6, f_gross[30], 0.003, "MPa")
m = WallModel(0.3, sec, law, M.ceramic.thermal_strain, n_el=4, T_initial=20.0)
P0 = 0.01 * f_gross[30] * 1e6 * t                 # 1% of the capacity, so that the equilibrium is unique
m.apply_load(P0, 0, 0)
u0 = m.top_displacement
fields = [(60.0 * k, np.full(len(sec.y), 20.0 + 29.0 * k)) for k in range(21)]      # 20 -> 600 °C
tc, _ = m.fire(P0, 0, 0, fields)
record("V6 free expansion from 20 to 600 °C under 1% of the load", (m.top_displacement - u0) * 1e3,
          0.3 * float(M.ceramic.thermal_strain(600.0)) * 1e3, 0.01, "mm")
P = m.capacity(0, 0, tol=1e-4)
record("V6 capacity at 600 °C (gross area)", P / t / 1e6, f_gross[600], 0.003, "MPa")

# V7
if not a.fast:
    Ti = 20.0
    exp = th.boundary_lengths_2d(m2, faces=("y0",)); nexp = th.boundary_lengths_2d(m2, faces=("y1",))
    ctr = [Boundary(exp, curves.iso834(Ti), 25, 0.8), Boundary(nexp, curves.constant(Ti), 4.0, 0.8)]
    r = th.solve_2d(m2, M, ctr, 120 * 60, dt=30, T_initial=Ti, field_instants=list(np.arange(60.0, 7201, 60.0)))
    fields = [(0.0, np.full(len(sec.y), Ti))] + [(tt, sec.temperatures(m2, r.fields[tt])) for tt in sorted(r.fields)]
    res = {}
    for n in (16, 32, 64):
        mdl = WallModel(H, sec, law, M.ceramic.thermal_strain, n_el=n, imperfection=H / 450, T_initial=Ti)
        Pm = mdl.capacity(0, 0)
        mdl = WallModel(H, sec, law, M.ceramic.thermal_strain, n_el=n, imperfection=H / 450, T_initial=Ti)
        mdl.apply_load(0.4 * Pm, 0, 0)
        tc, _ = mdl.fire(0.4 * Pm, 0, 0, fields)
        res[n] = (Pm, tc / 60 if tc else np.inf)
        print(f"    n = {n}: capacity {Pm / 1e3:.2f} kN/m, fire resistance at 40%: {res[n][1]:.1f} min")
    for n in (16, 32):
        record(f"V7 cold capacity with {n} elements (ref. 64)", res[n][0] / 1e3, res[64][0] / 1e3, 0.01, "kN/m")
        record(f"V7 fire resistance with {n} elements (ref. 64)", res[n][1], res[64][1],
                  0.02 if n == 32 else 0.20, "min")

output = os.path.join(os.path.dirname(__file__), "..", "results", "verification_wall.csv")
save_csv(output, lines, ["case", "model", "reference", "unit", "error", "tolerance", "status"])
nf = sum(1 for l in lines if l[-1] != "ok")
print(f"\n{len(lines) - nf}/{len(lines)} cases within tolerance. Written to {os.path.abspath(output)}")
sys.exit(1 if nf else 0)
