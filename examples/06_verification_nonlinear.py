"""Verification of the non-linear model (damage in the ceramic and dry joints).

  N1. Uniaxial tension with a weakened layer: peak strength and
      dissipated energy = G_f x area, for two element sizes
      (mesh objectivity).
  N2. Homogeneous uniaxial compression: peak stress = fc.
  N3. Dry joint: transmits no tension; in compression, stiffness equal to that of the
      linear model; under shear with compression, the force saturates at mu x N.
  N4. Joint with initial gap g: almost no force until the gap closes; then
      force = series stiffness x (displacement - g).
  N5. Crack closure: a bar cracked in tension and then compressed
      recovers the intact stiffness in compression.

Usage: python examples/06_verification_nonlinear.py [--element H8|W15]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, mechanical as mech
from mef_masonry.mesh import Layer, CERAMIC
from mef_masonry.thermal import Materials
from mef_masonry.nonlinear import NonlinearModel
from mef_masonry.parameters import CeramicDamage, FrictionJoint
from mef_masonry.output import save_csv

ap = argparse.ArgumentParser()
ap.add_argument("--element", default="W6", choices=["W6", "H8", "W15"])
EL = ap.parse_args().element
_extrude = msh.extrude


def extrude(*a, **k):
    return _extrude(*a, kind=EL, **k)


msh.extrude = extrude
print(f"3D element: {EL}\n")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT_DIR, exist_ok=True)
lines = []


def record(case, quantity, fem, ref, error, tol, un):
    ok = abs(error) <= tol
    lines.append([case, quantity, f"{fem:.5g}", f"{ref:.5g}", f"{error:.3g}", f"{tol:.3g}", un, "OK" if ok else "FAILED"])
    print(f"  {quantity:<46s} FEM {fem:11.5g}  reference {ref:11.5g}  error {error:8.3g} {un}  {'OK' if ok else 'FAILED'}")


M = Materials()
Pd = CeramicDamage()
Pj0 = FrictionJoint(gap_std=0.0, E_contact=None)   # perfect contact, joint E from Materials
B = 0.02
A = B * B

print("N1. Uniaxial tension: peak and dissipated energy (mesh objectivity)")
for dz in (0.010, 0.005):
    m2 = msh.rectangle_mesh(B, B, 0.01)
    L = 0.10
    m3 = msh.extrude(m2, [Layer(0.0, L, "uniform", material=CERAMIC)], dz=dz)
    zc = m3.xyz[m3.conn3].mean(1)[:, 2]
    weak = np.abs(zc - (L / 2 - dz / 2)) < dz / 4      # one layer in the middle
    factor = np.where(weak, 0.95, 1.0)
    mod = NonlinearModel(m3, M, Pd, Pj0, support="platen", ft_factor_elem=factor, verbose=False)
    presc = 3 * mod.top + 2
    ds = np.linspace(0, 0.12e-3, 121)
    Ps = []
    for dlt in ds[1:]:
        P, it, conv = mod._step(presc, np.full(len(presc), dlt), np.zeros(mod.ndof))
        Ps.append(-P)
    Ps = np.r_[0.0, Ps]
    W = np.trapezoid(Ps, ds) if hasattr(np, "trapezoid") else np.trapz(Ps, ds)
    Gf = Pd.du_t * Pd.ft * 0.95
    record(f"N1 dz={dz * 1e3:.0f}mm", "peak stress (MPa)", Ps.max() / A / 1e6, 0.95 * Pd.ft / 1e6,
              (Ps.max() / A - 0.95 * Pd.ft) / (0.95 * Pd.ft) * 100, 2.0, "%")
    # 10% tolerance: in the band, the damaged layer has its lateral contraction
    # restrained by its neighbours, which raises the effective principal stress by up to ~9%
    record(f"N1 dz={dz * 1e3:.0f}mm", "dissipated energy / area (N/m)", W / A, Gf, (W / A - Gf) / Gf * 100, 10.0, "%")

print("\nN5. Crack closure: the bar of N1, cracked in tension, is compressed")
E0 = float(M.ceramic.E_T(np.array([M.ceramic.T_ref]))[0])
k_el = E0 * A / L
dlt_c = 2e-6
Pc5 = mod._step(presc, np.full(len(presc), -dlt_c), np.zeros(mod.ndof))[0]
record("N5", "compressive stiffness after cracking / intact", Pc5 / (k_el * dlt_c), 1.0,
          (Pc5 / (k_el * dlt_c) - 1) * 100, 2.0, "%")

print("\nN2. Homogeneous uniaxial compression: peak = fc")
m2 = msh.rectangle_mesh(B, B, 0.01)
m3 = msh.extrude(m2, [Layer(0.0, 0.04, "uniform", material=CERAMIC)], dz=0.01)
mod = NonlinearModel(m3, M, Pd, Pj0, support="platen", verbose=False)
presc = 3 * mod.top + 2
Ps = []
for dlt in np.linspace(0, 0.12e-3, 61)[1:]:
    Ps.append(mod._step(presc, np.full(len(presc), -dlt), np.zeros(mod.ndof))[0])
record("N2", "peak stress (MPa)", max(Ps) / A / 1e6, Pd.fc / 1e6, (max(Ps) / A - Pd.fc) / Pd.fc * 100, 2.0, "%")

print("\nN3. Dry joint between two units (wide, low body, as in a direct shear test)")
Bj, tb, tj = 0.10, 0.01, 0.002
lay = [Layer(0.0, tb, "uniform", material=CERAMIC), Layer(tb, tb + tj, "joint"),
       Layer(tb + tj, 2 * tb + tj, "uniform", material=CERAMIC)]
m3 = msh.extrude(msh.rectangle_mesh(Bj, Bj, 0.02), lay, dz=0.005)
Aj = Bj * Bj
mod = NonlinearModel(m3, M, Pd, Pj0, support="platen", verbose=False)
presc = 3 * mod.top + 2
d0 = 2e-7
Pc = mod._step(presc, np.full(len(presc), -d0), np.zeros(mod.ndof))[0]
mod = NonlinearModel(m3, M, Pd, Pj0, support="platen", verbose=False)
Pt = mod._step(presc, np.full(len(presc), +d0), np.zeros(mod.ndof))[0]
k_an = Aj / (2 * tb / M.ceramic.E_material + tj / M.joint.E)      # series: ceramic + thin layer
record("N3", "compressive stiffness (kN/mm)", Pc / d0 / 1e6, k_an / 1e6, (Pc / d0 - k_an) / k_an * 100, 1.0, "%")
record("N3", "force in tension / force in compression", abs(Pt) / Pc, 0.0, abs(Pt) / Pc, 1e-3, "-")

mod = NonlinearModel(m3, M, Pd, Pj0, support="fixed_end", verbose=False)
N = 0.5e6 * Aj                     # compression of 0.5 MPa
F = np.zeros(mod.ndof)
F[3 * mod.top + 2] = -N * mod.top_area[mod.top] / mod.top_area[mod.top].sum()
presc = 3 * mod.top
for s in (1e-6, 3e-6, 1e-5, 3e-5, 1e-4):
    mod._step(presc, np.full(len(presc), s), F)
K, Fth = mod._assemble(mod.d, mod.En, mod.Gs)
V = (K @ mod.u - Fth - F)[presc].sum()
mu = FrictionJoint().mu
record("N3", "final shear / (mu N)", V / (mu * N), 1.0, (V / (mu * N) - 1) * 100, 2.0, "%")

print("\nN4. Joint with uniform gap g: almost zero force until it closes; then, series stiffness")
g = 4e-6
mod = NonlinearModel(m3, M, Pd, Pj0, support="platen", verbose=False)
mod.gap[mod.bed] = g; mod.eg = mod.gap / mod.dzj
mod.closed = mod.eg <= 0; mod.En = np.where(mod.closed, mod.Ej, mod.Pj.open_factor * mod.Ej)
presc = 3 * mod.top + 2
P1 = mod._step(presc, np.full(len(presc), -0.5 * g), np.zeros(mod.ndof))[0]
d1 = 1e-6
P2 = mod._step(presc, np.full(len(presc), -(g + d1)), np.zeros(mod.ndof))[0]
record("N4", "force at 0.5 g / (k x g)", P1 / (k_an * g), 0.0, P1 / (k_an * g), 1e-3, "-")
record("N4", "force at g + d / (k x d)", P2 / (k_an * d1), 1.0, (P2 / (k_an * d1) - 1) * 100, 1.0, "%")

save_csv(os.path.join(OUT_DIR, "verification_nonlinear.csv" if EL == "W6" else f"verification_nonlinear_{EL}.csv"), lines,
           ["case", "quantity", "MEF", "reference", "error", "tolerance", "unit", "status"])
failures = sum(1 for l in lines if l[-1] != "OK")
print(f"\n{len(lines) - failures} of {len(lines)} checks OK.")
sys.exit(1 if failures else 0)
