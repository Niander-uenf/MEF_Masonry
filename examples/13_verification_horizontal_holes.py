"""Verification of phase 3: units with horizontal holes (infill).

V1  Hole volume in the 3D prism (one and a half units in running bond, 3 courses,
    mortar joints): holes × unit length of each course (the
    head joints close the holes), exact.
V2  Minimum horizontal section: in a row of holes, the solid area is the gross
    area minus the columns of holes (outside the head joints), exact.
V3  Vertical section (2D): area of the holes and of the half joints, exact.
V4  Thermal: the vertical 2D section reproduces the 3D model extruded along the holes
    (one course heated by ISO 834 on one face, 60 min; same properties
    in the unit and the joint so that both meshes solve the same problem).
V5  Wall: the fibres lie in the continuous walls (faces and vertical webs),
    with an area per metre equal to the width minus the columns of holes, and the
    centroid on the mid-plane.
V6  Radiation in the holes: the dimension used is the hole width in the direction of the
    heat flow.
V7  Non-linear prism of the 8-hole block at 20 °C: runs to the peak, and the peak
    is that of the walls that run down continuously (faces and vertical web), fc × their
    area, with a 5% tolerance (the horizontal webs carry almost nothing).
    Only with --full (about 1 min).
V8  Library: the units with horizontal holes pass validation; the soil-cement
    brick (hole diameter still without a source) is rejected.

Usage: python examples/13_verification_horizontal_holes.py [--full]
Writes results/verification_horizontal_holes.csv.
"""
import os
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, thermal as th, project as pj, library as lib, curves
from mef_masonry.mesh import AIR, MORTAR, CERAMIC
from mef_masonry.parameters import Geometry
from mef_masonry.thermal import Materials, Boundary
from mef_masonry.wall import FibreSection
from mef_masonry.output import save_csv

FULL = "--full" in sys.argv
lines = []


def record(case, model, ref, tol_rel, unit_str="", tol_abs=0.0):
    error = abs(model - ref) / max(abs(ref), 1e-30)
    ok = error <= tol_rel or abs(model - ref) <= tol_abs
    tol = f"{100 * tol_rel:.1f}%" if tol_abs <= 0 else f"{tol_abs:g} {unit_str}"
    lines.append([case, round(model, 6), round(ref, 6), unit_str, f"{100 * error:.3f}%", tol,
                   "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: model {model:.6g} {unit_str}, reference {ref:.6g} ({100 * error:.3f}%)")


def condition(case, ok, a, b, unit_str, criterion):
    lines.append([case, round(a, 6), round(b, 6), unit_str, "—", criterion, "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: {a:.4g} and {b:.4g} {unit_str} ({criterion})")


def volumes(m):
    xyz = m.xyz[m.conn3]
    p = xyz[:, :3, :2]
    a = 0.5 * np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                     - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
    return a * (xyz[:, 3, 2] - xyz[:, 0, 2])


tj = 0.010
g = Geometry(hole_shape="horizontal", width=0.090, course_height=0.200, unit_length=0.200,
              joint_thickness=tj, joint_type="mortar", face_shell=0.008, web=0.007, bed_shell=0.008,
              holes_width=2, holes_height=4, units_per_course=1.5, bond="running_bond", n_courses=3)
assert not g.validate(), g.validate()
wy, hz = g.hole_width_h, g.hole_height_h(tj)

# ---------------------------------------------------------------------------- V1 and V2
m = msh.prism_mesh(g, h=0.008, dz=0.010, t_joint=tj, symmetry=False, plates=False, capping=True)
vol = volumes(m)
ref = sum(g.holes_width * g.holes_height * wy * hz * (g.length - tj * len(j)) for j in g.joints_per_course)
record("V1 volume of the horizontal holes in the 3D prism", vol[m.mat == AIR].sum() * 1e6, ref * 1e6, 1e-9, "cm³")

c_row = [k for k, c in enumerate(m.layers) if c.kind == "unit" and c.hole_h == 1 and c.course == 1][0]
sel = m.layer == c_row
dz = (m.xyz[m.conn3[sel][:, 3], 2] - m.xyz[m.conn3[sel][:, 0], 2])
nsub = len(np.unique(np.round(m.xyz[m.conn3[sel][:, 0], 2], 9)))
a_sol = (vol[sel][m.mat[sel] != AIR] / dz[m.mat[sel] != AIR]).sum() / nsub
ref2 = g.length * g.width - g.holes_width * wy * (g.length - tj * len(g.joints_per_course[1]))
record("V2 solid area in a row of holes (course 2)", a_sol * 1e4, ref2 * 1e4, 1e-9, "cm²")

# ---------------------------------------------------------------------------- V3
m2 = msh.section_mesh(g, h=0.004, symmetry=False, bands=())
a2 = np.abs(m2.areas())
mat2 = th.section_material(m2, g)
record("V3 hole area in the vertical section", a2[mat2 == AIR].sum() * 1e4,
          g.holes_width * g.holes_height * wy * hz * 1e4, 1e-9, "cm²")
record("V3 area of the half joints in the vertical section", a2[mat2 == MORTAR].sum() * 1e4, g.width * tj * 1e4,
          1e-9, "cm²")

# ---------------------------------------------------------------------------- V4
g4 = replace(g, unit_length=0.040, units_per_course=1, bond="stack_bond", n_courses=1)
M4 = Materials()
M4.ceramic.k, M4.ceramic.rho, M4.ceramic.c = 0.90, 1600.0, 920.0
M4.mortar.k, M4.mortar.rho, M4.mortar.c = 0.90, 1600.0, 920.0
M4.hole_diameter = g4.hole_dimension
Ti, tf, h4 = 20.0, 3600.0, 0.004
f_iso = curves.iso834(Ti)
m2v = msh.section_mesh(g4, h=h4, symmetry=False, bands=())
ctr2 = [Boundary(th.boundary_lengths_2d(m2v, faces=("y0",)), f_iso, 25.0, 0.8),
        Boundary(th.boundary_lengths_2d(m2v, faces=("y1",)), curves.constant(Ti), 4.0, 0.8)]
face2 = np.nonzero(th.boundary_lengths_2d(m2v, faces=("y1",)) > 0)[0]
r2 = th.solve_2d(m2v, M4, ctr2, tf, dt=15.0, T_initial=Ti, mat=th.section_material(m2v, g4))
m3 = msh.prism_mesh(g4, h=h4, dz=h4, t_joint=tj, symmetry=False, plates=False, capping=False)
ctr3 = [Boundary(th.side_areas_3d(m3, faces=("y0",)), f_iso, 25.0, 0.8),
        Boundary(th.side_areas_3d(m3, faces=("y1",)), curves.constant(Ti), 4.0, 0.8)]
r3 = th.solve_3d(m3, M4, ctr3, tf, dt=15.0, T_initial=Ti)
face3 = np.nonzero(np.isclose(m3.xyz[:, 1], g4.width))[0]
record("V4 unexposed face at 60 min: vertical 2D section × 3D", float(r2.T_final[face2].mean()),
          float(r3.T_final[face3].mean()), 0.0, "°C", tol_abs=1.0)
mid2 = np.nonzero(np.isclose(m2v.xy[:, 1], g4.width / 2))[0]
mid3 = np.nonzero(np.isclose(m3.xyz[:, 1], g4.width / 2))[0]
record("V4 vertical web (mid-thickness) at 60 min: 2D × 3D", float(r2.T_final[mid2].mean()),
          float(r3.T_final[mid3].mean()), 0.0, "°C", tol_abs=1.0)

# ---------------------------------------------------------------------------- V5
sec = FibreSection.from_mesh(m2, g)
record("V5 fibre area per metre of wall = width - columns of holes", sec.area * 1e3,
          (g.width - g.holes_width * wy) * 1e3, 1e-9, "mm²/mm")
cy = float((sec.A * sec.y).sum() / sec.area)
condition("V5 centroid of the fibres on the mid-plane", abs(cy) < 1e-9, cy * 1e3, 0.0, "mm", "|y| < 1e-6 mm")
cnt = m2.xy[m2.tri[sec.tri]].mean(1)
inside = sum(((cnt[:, 1] > y0) & (cnt[:, 1] < y1)).sum() for y0, y1 in g.holes_h_y())
condition("V5 no fibre in the columns of holes", inside == 0, float(inside), 0.0, "fibres", "zero")

# ---------------------------------------------------------------------------- V6
p6 = pj.complete(lib.apply(pj.default_project(), "unit", "infill_8_holes"))
geo6 = pj.geometry(p6)
record("V6 hole dimension in the radiation = hole width", pj.materials(p6, geo6).hole_diameter * 1e3,
          geo6.hole_width_h * 1e3, 1e-12, "mm")

# ---------------------------------------------------------------------------- V7
if FULL:
    from mef_masonry.nonlinear import NonlinearModel
    q = pj.complete(lib.apply(pj.default_project(), "unit", "infill_8_holes"))
    q["model"]["analysis"] = "nonlinear"
    q["test"].update({"kind": "ambient", "T_initial": 20})
    q["mesh"].update({"h": 12, "dz": 20})
    geo = pj.geometry(q)
    M = pj.materials(q, geo)
    m7 = pj.mesh_3d(q, geo)
    mod = NonlinearModel(m7, M, pj.damage(q), pj.friction_joint(q), verbose=False, mortar_damage=pj.damage_mor(q),
                          grout_damage=pj.damage_grt(q))
    mod.heat({0.0: np.full(m7.n_nodes, M.ceramic.T_ref)})
    d, P = mod.apply_load(0.006 * geo.prism_height, 60)
    i = int(np.argmax(P))
    A_walls = geo.length * (geo.width - geo.holes_width * geo.hole_width_h)
    ceiling = q["model"]["fc"] * 1e6 * A_walls
    condition("V7 8-hole prism: peak reached (the load falls afterwards)", i < len(P) - 1, P[i] / 1e3, P[-1] / 1e3,
             "kN", "peak before the end")
    record("V7 peak of the 8-hole prism = fc × area of the continuous walls", P[i] / 1e3, ceiling / 1e3, 0.05, "kN")

# ---------------------------------------------------------------------------- V8
ok8 = []
for k in ("infill_8_holes", "infill_6_holes"):
    ok8.append(not pj.validate(pj.complete(lib.apply(pj.default_project(), "unit", k))))
reject = bool(pj.validate(pj.complete(lib.apply(pj.default_project(), "unit", "soil_cement_modular"))))
record("V8 units with horizontal holes that pass validation", sum(ok8), 2, 0.0, "items")
condition("V8 soil-cement brick (without hole diameter) rejected", reject, float(reject), 1.0, "", "rejected")

output = os.path.join(os.path.dirname(__file__), "..", "results", "verification_horizontal_holes.csv")
save_csv(output, lines, ["case", "model", "reference", "unit", "error", "tolerance", "status"])
nf = sum(1 for l in lines if l[-1] == "FAILED")
print(f"\n{len(lines) - nf}/{len(lines)} cases within tolerance. Written to {os.path.abspath(output)}")
sys.exit(1 if nf else 0)
