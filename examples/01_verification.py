"""Verification of the code against analytical solutions.

Run before any study:  python examples/01_verification.py
Writes results/verification.csv. All cases must appear as OK.
Other 3D elements:  python examples/01_verification.py --element H8   (or W15)
-> results/verification_H8.csv
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, thermal as th, mechanical as mech
from mef_masonry.mesh import Layer, CERAMIC, JOINT
from mef_masonry.thermal import Materials, Boundary
from mef_masonry.parameters import Geometry
from mef_masonry import curves
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


def record(case, quantity, fem, analytic, error, tolerance, unit_str):
    ok = abs(error) <= tolerance
    lines.append([case, quantity, f"{fem:.6g}", f"{analytic:.6g}", f"{error:.3g}", f"{tolerance:.3g}",
                   unit_str, "OK" if ok else "FAILED"])
    print(f"  {quantity:<42s} FEM {fem:12.5g}  analytical {analytic:12.5g}  error {error:9.3g} {unit_str}"
          f"  {'OK' if ok else 'FAILED'}")


def plate_1d_centre(Fo, terms=200):
    n = np.arange(terms)
    return np.sum(4 * (-1) ** n / ((2 * n + 1) * np.pi) * np.exp(-((2 * n + 1) * np.pi / 2) ** 2 * Fo))


# ---------------------------------------------------------------------------
print("V1. Transient conduction in a plate with imposed surface temperature")
M = Materials()
k, rc = M.ceramic.k, M.ceramic.rho * M.ceramic.c
esp = 0.14; l = esp / 2; a = k / rc
Ti, Ts = 30.0, 600.0
Fos = [0.05, 0.1, 0.2, 0.5, 1.0]
times = [F * l ** 2 / a for F in Fos]
for name, kind in ((f"3D ({EL})", "3d"), ("2D (triangles)", "2d")):
    for dt in (60.0, 15.0):
        m2 = msh.rectangle_mesh(esp, 0.01, 0.0035, hy=0.005)
        if kind == "3d":
            m3 = msh.extrude(m2, [Layer(0.0, 0.01, "uniform", material=CERAMIC)], dz=0.01)
            area = th.side_areas_3d(m3, faces=("x0", "x1"))
            ctr = [Boundary(area, curves.constant(Ts), 1e8, 0.0)]
            centre = np.nonzero(np.isclose(m3.xyz[:, 0], l))[0]
            r = th.solve_3d(m3, M, ctr, max(times), dt=dt, T_initial=Ti, probes={"c": centre},
                               record_every=dt, field_instants=times)
            T_c = [r.fields[t][centre].mean() for t in times]
        else:
            area = th.boundary_lengths_2d(m2, faces=("x0", "x1"))
            ctr = [Boundary(area, curves.constant(Ts), 1e8, 0.0)]
            centre = np.nonzero(np.isclose(m2.xy[:, 0], l))[0]
            mat = np.full(len(m2.tri), CERAMIC, dtype=np.int8)
            r = th.solve_2d(m2, M, ctr, max(times), mat=mat, dt=dt, T_initial=Ti,
                               record_every=dt, field_instants=times)
            T_c = [r.fields[t][centre].mean() for t in times]
        print(f" {name}, dt = {dt:.0f} s")
        for F, Tc in zip(Fos, T_c):
            Tan = Ts + (Ti - Ts) * plate_1d_centre(F)
            tol = 6.0 if dt == 60.0 else 2.0
            record(f"V1 {name} dt={dt:.0f}s", f"T at the centre, Fo = {F}", Tc, Tan, Tc - Tan, tol, "°C")

# ---------------------------------------------------------------------------
print("\nV2. Uniaxial compression of a solid block (3D)")
M = Materials()
E = M.ceramic.E_material; nu = M.ceramic.nu
B, H, P = 0.10, 0.20, 10e3
m2 = msh.rectangle_mesh(B, B, 0.02)
m3 = msh.extrude(m2, [Layer(0.0, H, "uniform", material=CERAMIC)], dz=0.04)
sig = P / B ** 2
for mode_, label in (("rigid", "rigid platen"), ("steel", "pressure on the top")):
    R = mech.solve_3d(m3, M, load_N=P, include_thermal=False, plates=mode_)["load"]
    uz_top = mech.mean_uz_plane(m3, R.u, H)
    record(f"V2 {label}", "top displacement (mm)", uz_top * 1e3, -sig * H / E * 1e3,
              (uz_top + sig * H / E) / (sig * H / E) * 100, 0.1, "%")
    dux = R.u[np.isclose(m3.xyz[:, 0], B), 0].mean() - R.u[np.isclose(m3.xyz[:, 0], 0), 0].mean()
    record(f"V2 {label}", "lateral expansion (mm)", dux * 1e3, nu * sig / E * B * 1e3,
              (dux - nu * sig / E * B) / (nu * sig / E * B) * 100, 0.1, "%")
    record(f"V2 {label}", "mean szz stress (MPa)", np.nanmean(R.sig[:, 2]) / 1e6, -sig / 1e6,
              (np.nanmean(R.sig[:, 2]) + sig) / sig * 100, 0.1, "%")

# ---------------------------------------------------------------------------
print("\nV3. Uniform free thermal expansion (3D)")
dT = 500.0
T = np.full(m3.n_nodes, M.ceramic.T_ref + dT)
R = mech.solve_3d(m3, M, T_nodes=T)["thermal"]
uz_top = mech.mean_uz_plane(m3, R.u, H)
expected = M.ceramic.alpha * dT * H
record("V3", "free elongation (mm)", uz_top * 1e3, expected * 1e3,
          (uz_top - expected) / expected * 100, 0.1, "%")
record("V3", "largest |stress| (MPa)", np.nanmax(np.abs(R.sig)) / 1e6, 0.0,
          np.nanmax(np.abs(R.sig)) / 1e6, 1e-3, "MPa")

# ---------------------------------------------------------------------------
print("\nV4. Layered column: ceramic + thin joints, nu = 0 (series association)")
M = Materials()
M.ceramic.nu = 0.0; M.joint.nu = 0.0; M.joint.E = 100e6
lay, z = [], 0.0
for i in range(5):
    lay.append(Layer(z, z + 0.058, "uniform", material=CERAMIC)); z += 0.058
    if i < 4:
        lay.append(Layer(z, z + 0.002, "uniform", material=JOINT)); z += 0.002
m3 = msh.extrude(m2, lay, dz=0.02)
R = mech.solve_3d(m3, M, load_N=P, include_thermal=False)["load"]
flex = 5 * 0.058 / M.ceramic.E_material + 4 * 0.002 / M.joint.E
uz_top = mech.mean_uz_plane(m3, R.u, z)
record("V4", "total shortening (mm)", -uz_top * 1e3, sig * flex * 1e3,
          (-uz_top - sig * flex) / (sig * flex) * 100, 0.1, "%")
Eeq_fem = sig * z / -uz_top; Eeq_an = z / flex
record("V4", "equivalent modulus (MPa)", Eeq_fem / 1e6, Eeq_an / 1e6,
          (Eeq_fem - Eeq_an) / Eeq_an * 100, 0.1, "%")

# ---------------------------------------------------------------------------
print("\nV5. Free plate with a parabolic temperature profile through the thickness (3D)")
M = Materials()
E = M.ceramic.E_material; nu = M.ceramic.nu; alpha = M.ceramic.alpha
Lx, t, Hz, dT0 = 0.24, 0.04, 0.24, 300.0
m2 = msh.rectangle_mesh(Lx, t, 0.01, hy=0.0025)
# the quadratic wedge uses 2 cm layers (same number of nodes over the height; the
# direct factorisation of W15 with 1 cm layers runs out of memory without pypardiso)
m3 = msh.extrude(m2, [Layer(0.0, Hz, "uniform", material=CERAMIC)], dz=0.02 if EL == "W15" else 0.01)
y = m3.xyz[:, 1]
T = M.ceramic.T_ref + dT0 * (2 * y / t - 1) ** 2
t0 = time.time()
R = mech.solve_3d(m3, M, T_nodes=T)["thermal"]
cnt = m3.xyz[m3.conn3].mean(1)
near = (np.abs(cnt[:, 0] - Lx / 2) < 0.011) & (np.abs(cnt[:, 2] - Hz / 2) < 0.011)
yc = cnt[near, 1]
an = E * alpha * (dT0 / 3 - dT0 * (2 * yc / t - 1) ** 2) / (1 - nu)
ref = E * alpha * dT0
for comp, name in ((0, "sxx"), (2, "szz")):
    err = np.max(np.abs(R.sig[near, comp] - an)) / ref * 100
    record("V5", f"{name}: largest error / (E alpha dT0)", np.max(np.abs(R.sig[near, comp])) / 1e6,
              np.max(np.abs(an)) / 1e6, err, 3.0, "%")
print(f"  ({len(m3.conn3)} elements {EL}, {time.time() - t0:.1f} s)")

# ---------------------------------------------------------------------------
print("\nV6. The same parabolic profile in the 2D section (generalised plane strain)")
m2 = msh.rectangle_mesh(Lx, t, 0.0025)
T = M.ceramic.T_ref + dT0 * (2 * m2.xy[:, 1] / t - 1) ** 2
mat = np.full(len(m2.tri), CERAMIC, dtype=np.int8)
s = mech.solve_section_gps(m2, M, T, mat=mat)
cnt = m2.xy[m2.tri].mean(1)
near = np.abs(cnt[:, 0] - Lx / 2) < 0.006
an = E * alpha * (dT0 / 3 - dT0 * (2 * cnt[near, 1] / t - 1) ** 2) / (1 - nu)
for comp in ("sx", "sz"):
    err = np.max(np.abs(s[comp][near] - an)) / ref * 100
    record("V6", f"{comp}: largest error / (E alpha dT0)", np.max(np.abs(s[comp][near])) / 1e6,
              np.max(np.abs(an)) / 1e6, err, 3.0, "%")

# ---------------------------------------------------------------------------
print("\nV7. Section with holes: 2D thermal x extruded 3D (must coincide)")
M = Materials()
geo = Geometry()
f, t_end, _ = curves.ramp_dwell(T_target=600.0)
m2 = msh.section_mesh(geo, h=0.006, symmetry=True)
ctr2 = [Boundary(th.boundary_lengths_2d(m2), f, 10.0, 0.8)]
r2 = th.solve_2d(m2, M, ctr2, t_end, dt=30.0)
m3 = msh.extrude(m2, [Layer(0.0, 0.02, "unit", course=0)], dz=0.01, geo=geo)
m3.mat[m3.mat == JOINT] = CERAMIC
ctr3 = [Boundary(th.side_areas_3d(m3), f, 10.0, 0.8)]
r3 = th.solve_3d(m3, M, ctr3, t_end, dt=30.0)
T3 = r3.T_final[:m2.n_nodes]          # corner nodes of the plane z = 0 (they come first in all meshes)
for name, v2, v3 in (("section minimum", r2.T_final.min(), T3.min()),
                     ("nodal mean", r2.T_final.mean(), T3.mean())):
    # H8 and W15: the 3D section (quadrilaterals, or quadratic triangles) is not the T3 mesh
    # of the 2D model; both converge to the same solution, but do not coincide node by node
    record("V7", f"T {name} at the end (°C)", v3, v2, v3 - v2, 0.5 if EL == "W6" else 1.0, "°C")

save_csv(os.path.join(OUT_DIR, "verification.csv" if EL == "W6" else f"verification_{EL}.csv"), lines,
           ["case", "quantity", "MEF", "analytical/reference", "error", "tolerance", "unit", "status"])
failures = sum(1 for l in lines if l[-1] != "OK")
print(f"\n{len(lines) - failures} of {len(lines)} checks OK.")
sys.exit(1 if failures else 0)
