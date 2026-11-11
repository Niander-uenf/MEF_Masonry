"""Verification of phase 2: rectangular cells and grout.

V1  Net area of the rectangular cells: the area of the unit elements in the
    2D mesh, course by course (with running bond the cells move), equals
    the gross area minus the cells computed from the geometry. Includes a block with
    hollow shells (small cells in the face shells).
V2  Grout volume in the 3D prism (stack bond, mortar, all cells
    grouted): grouted area × height between the cappings.
V3  Stiffness in parallel: one course of units with alternate grouted cells,
    with the same Poisson's ratio in both materials (same lateral strain, no
    interaction): shortening = P H / (Eb Ab + Eg Ag), exact.
V4  The grout properties reach the solvers: conductivity, capacity,
    modulus and thermal strain, including the EN 1992-1-2 tables.
V5  2D section and wall: the grout enters as a material in the thermal analysis and as fibres
    in the wall; the fibre area is the load-bearing area (unit + grout).
V6  Grout increases the strength of the non-linear prism (concrete block B1 of
    Camacho et al., 2015) and an informative comparison with the tests: hollow prism
    5.63 MPa and grouted with G1 11.29 MPa (gross area). Only with --full
    (about 10 minutes).
V7  Library: the units already released pass validation; those not yet released
    (soil-cement brick, without measured holes) are rejected.

Usage: python examples/12_verification_cells_grout.py [--full]
Writes results/verification_cells_grout.csv.
"""
import os
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, thermal as th, mechanical as mech, project as pj, library as lib
from mef_masonry.mesh import CERAMIC, GROUT, AIR
from mef_masonry.parameters import Geometry, EN1992_E_FACTOR, EN1992_THERMAL_STRAIN
from mef_masonry.thermal import Materials
from mef_masonry.wall import FibreSection
from mef_masonry.output import save_csv

FULL = "--full" in sys.argv
lines = []


def record(case, model, ref, tol_rel, unit_str=""):
    error = abs(model - ref) / max(abs(ref), 1e-30)
    ok = error <= tol_rel
    lines.append([case, round(model, 6), round(ref, 6), unit_str, f"{100 * error:.3f}%", f"{100 * tol_rel:.1f}%",
                   "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: model {model:.6g} {unit_str}, reference {ref:.6g} ({100 * error:.3f}%)")


def compare(case, model, ref, unit_str=""):
    """Comparison with a test (informative: does not count as a failure)."""
    error = (model - ref) / ref
    lines.append([case, round(model, 4), round(ref, 4), unit_str, f"{100 * error:+.1f}%", "—", "comparison"])
    print(f"--  {case}: model {model:.4g} {unit_str}, test {ref:.4g} ({100 * error:+.1f}%)")


def condition(case, ok, a, b, unit_str, criterion):
    lines.append([case, round(a, 6), round(b, 6), unit_str, "—", criterion, "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: {a:.4g} and {b:.4g} {unit_str} ({criterion})")


# ---------------------------------------------------------------------------- V1
tj = 0.010
cases = {
    "concrete 14×19×39, one and a half units, running bond": Geometry(
        hole_shape="rectangular", unit_length=0.40, width=0.14, course_height=0.20, holes_per_unit=2,
        units_per_course=1.5, n_courses=2, bond="running_bond", joint_thickness=tj, joint_type="mortar",
        face_shell=0.032, web=0.030, end_shell=0.030),
    "clay block with hollow shells": Geometry(
        hole_shape="rectangular", unit_length=0.30, width=0.14, course_height=0.20, holes_per_unit=2,
        units_per_course=1, n_courses=1, bond="stack_bond", joint_thickness=tj, joint_type="mortar",
        face_shell=0.030, web=0.008, end_shell=0.008, wall_cells=6, min_shell=0.008),
}
for name, g in cases.items():
    assert not g.validate(), g.validate()
    m2 = msh.rect_section_mesh(g, h=0.008, symmetry=False, bands=g.joint_positions, t_joint=tj)
    a = np.abs(m2.areas())
    for i in range(g.n_courses):
        hole, _ = msh.course_states(m2, g, i)
        ref = g.length * g.width - g.hole_area(i, t_joint=tj)
        record(f"V1 unit area in course {i + 1}, {name}", a[~hole].sum() * 1e4, ref * 1e4, 1e-9, "cm²")

# ---------------------------------------------------------------------------- V2
gv = Geometry(hole_shape="rectangular", unit_length=0.30, width=0.14, course_height=0.20, holes_per_unit=2,
               units_per_course=1, n_courses=3, bond="stack_bond", joint_thickness=tj, joint_type="mortar",
               face_shell=0.028, web=0.028, end_shell=0.028, grouting="all")
t_cap = 0.003
mv = msh.prism_mesh(gv, h=0.01, dz=0.02, t_joint=tj, symmetry=False, plates=False, capping=True)
xyz = mv.xyz[mv.conn3]
p = xyz[:, :3, :2]
air = 0.5 * np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
vol = air * (xyz[:, 3, 2] - xyz[:, 0, 2])
A_g = gv.hole_area(0, t_joint=tj) - gv.hole_area(0, grout=True, t_joint=tj)
record("V2 grout volume in the 3-course prism", vol[mv.mat == GROUT].sum() * 1e6,
          A_g * (gv.prism_height - 2 * t_cap) * 1e6, 1e-9, "cm³")
v_air = float(vol[mv.mat == AIR].sum() * 1e6)
condition("V2 no void in the fully grouted prism with full bedding", v_air == 0.0, v_air, 0.0, "cm³", "zero")

# ---------------------------------------------------------------------------- V3
g3 = replace(gv, n_courses=1, grouting="alternate", holes_per_unit=4, units_per_course=1)
m3 = msh.prism_mesh(g3, h=0.01, dz=0.02, t_joint=tj, symmetry=True, plates=False, capping=False)
M3 = Materials()
M3.ceramic.E_unit_gross, M3.ceramic.net_area_ratio, M3.ceramic.nu = 9000e6, 1.0, 0.2
M3.grout.E_unit_gross, M3.grout.nu = 26000e6, 0.2
P_N = 100e3
R3 = mech.solve_3d(m3, M3, load_N=P_N, include_thermal=False)["load"]
H3 = g3.prism_height
d3 = -mech.shortening_between_plates(m3, R3, H3)
A_b = g3.length * g3.width - g3.hole_area(0, t_joint=tj)
A_g3 = g3.hole_area(0, t_joint=tj) - g3.hole_area(0, grout=True, t_joint=tj)
record("V3 shortening in parallel unit + grout (equal ν)", d3 * 1e6,
          P_N * H3 / (9000e6 * A_b + 26000e6 * A_g3) * 1e6, 1e-6, "µm")

# ---------------------------------------------------------------------------- V4
M4 = Materials()
g = M4.grout
g.k, g.rho, g.c, g.T_ref = 1.75, 2400.0, 1000.0, 20.0
g.E_unit_gross, g.net_area_ratio = 26.1e9, 1.0
mat4 = np.array([GROUT], dtype=np.int8)
record("V4 grout conductivity", float(M4.k(mat4, np.array([300.0]))[0]), 1.75, 1e-12, "W/m·K")
record("V4 volumetric capacity of the grout", float(M4.rho_c(mat4)[0]) / 1e6, 2.4, 1e-12, "MJ/m³·K")
g.table_E_factor, g.table_thermal_strain = list(EN1992_E_FACTOR), list(EN1992_THERMAL_STRAIN)
E4, _ = mech._properties(mat4, np.array([600.0]), M4)
record("V4 grout E at 600 °C (EN 1992-1-2: 0.45 × 0.0025/0.025)", E4[0] / 1e9, 26.1 * 0.045, 1e-12, "GPa")
eps4 = mech._thermal_strain_nodes(mat4, np.array([[400.0] * 6]), M4)[0, 0]
record("V4 grout thermal strain at 400 °C (siliceous)", eps4 * 1e3,
          (-1.8e-4 + 9e-6 * 400 + 2.3e-11 * 400 ** 3) * 1e3, 1e-3, "‰")

# ---------------------------------------------------------------------------- V5
p5 = pj.complete(lib.apply(pj.default_project(), "unit", "concrete_m15_29"))
geo5 = pj.geometry(p5)
geo1 = replace(geo5, units_per_course=1.0)
m25 = msh.section_mesh(geo1, h=0.005, symmetry=False, bands=())
mat5 = th.section_material(m25, geo1)
a5 = np.abs(m25.areas())
A_g5 = geo1.hole_area(0) - geo1.hole_area(0, grout=True)
record("V5 grout area in the 2D section (thermal material)", a5[mat5 == GROUT].sum() * 1e4, A_g5 * 1e4, 1e-9,
          "cm²")
sec = FibreSection.from_mesh(m25, geo1)
record("V5 wall fibre area = load-bearing area (unit + grout)", sec.area * sec.length * 1e4,
          geo1.length * geo1.width * pj.solid_area_ratio(geo1) * 1e4, 1e-9, "cm²")
record("V5 grout fibres in the wall", float((sec.A[sec.grout] * sec.length).sum() * 1e4), A_g5 * 1e4,
          1e-9, "cm²")

# ---------------------------------------------------------------------------- V6
if FULL:
    from mef_masonry.nonlinear import NonlinearModel
    res = {}
    for grt in ("none", "all"):
        q = pj.complete(lib.apply(pj.default_project(), "unit", "concrete_m15_29"))
        q["geometry"]["grouting"] = grt
        q["model"]["analysis"] = "nonlinear"
        q["test"].update({"kind": "ambient", "T_initial": 20})
        q["mesh"].update({"h": 12, "dz": 20})
        geo = pj.geometry(q)
        M = pj.materials(q, geo)
        m = pj.mesh_3d(q, geo)
        mod = NonlinearModel(m, M, pj.damage(q), pj.friction_joint(q), verbose=False, mortar_damage=pj.damage_mor(q),
                              grout_damage=pj.damage_grt(q))
        mod.heat({0.0: np.full(m.n_nodes, M.ceramic.T_ref)})
        d, P = mod.apply_load(0.006 * geo.prism_height, 60)
        res[grt] = float(P.max()) / (geo.length * geo.width) / 1e6
        print(f"    prism {grt}: {res[grt]:.2f} MPa (gross area)")
    condition("V6 grout increases the prism strength (MPa, gross area)", res["all"] > 1.5 * res["none"],
             res["all"], res["none"], "MPa", "grouted > 1.5 × hollow")
    compare("V6 hollow prism B1 against the test (Camacho et al., 2015)", res["none"], 5.63, "MPa")
    compare("V6 grouted prism B1-G1 against the test", res["all"], 11.29, "MPa")

# ---------------------------------------------------------------------------- V7
valid = rejected = 0
for k, item in lib.UNITS.items():
    q = pj.complete(lib.apply(pj.default_project(), "unit", k))
    errors = pj.validate(q)
    if item.get("phase", 0) <= lib.CURRENT_PHASE:
        valid += not errors
        if errors:
            print("   ", k, errors)
    else:
        rejected += bool(errors)
n_ok = sum(1 for v in lib.UNITS.values() if v.get("phase", 0) <= lib.CURRENT_PHASE)
record("V7 library units that run and pass validation", valid, n_ok, 0.0, "items")
record("V7 units of future phases rejected", rejected, len(lib.UNITS) - n_ok, 0.0, "items")

output = os.path.join(os.path.dirname(__file__), "..", "results", "verification_cells_grout.csv")
save_csv(output, lines, ["case", "model", "reference", "unit", "error", "tolerance", "status"])
nf = sum(1 for l in lines if l[-1] == "FAILED")
nv = sum(1 for l in lines if l[-1] != "comparison")
print(f"\n{nv - nf}/{nv} cases within tolerance"
      + (f" (plus {len(lines) - nv} comparisons with tests)" if len(lines) > nv else "")
      + f". Written to {os.path.abspath(output)}")
sys.exit(1 if nf else 0)
