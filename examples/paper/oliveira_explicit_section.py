"""Section 7.3 of the paper: where does the extra heat of the explicit section of the concrete block
of Oliveira et al. (2021) come from? Varies the cavity model and the conduction through the webs and measures the
insulation time. Base: project "oliveira" (lower-limit k, 1.5%, emissivity 0.7).

The cavity is a material with equivalent conductivity k_air + 4 sigma T^3 D eps/(2 - eps),
the exchange between infinite parallel plates at a distance D. In a cell of finite dimensions
only part of the radiation from the hot shell reaches the cold shell directly: the view factor between the
opposite faces of the cell. Scaling the radiative term by this factor is equivalent to using an
effective emissivity eps' with eps'/(2 - eps') = F eps/(2 - eps)."""
import os, sys, json
from math import atan, log, sqrt, pi
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oliveira import project, HERE
from mef_masonry import project as pj, thermal as th, mesh as msh, curves
from mef_masonry.mesh import CERAMIC, MORTAR
from mef_masonry.thermal import Boundary


def F_rectangles(a, b, c):
    """View factor between two aligned parallel rectangles a x b at a distance c."""
    X, Y = a / c, b / c
    return 2 / (pi * X * Y) * (log(sqrt((1 + X * X) * (1 + Y * Y) / (1 + X * X + Y * Y)))
        + X * sqrt(1 + Y * Y) * atan(X / sqrt(1 + Y * Y)) + Y * sqrt(1 + X * X) * atan(Y / sqrt(1 + X * X))
        - X * atan(X) - Y * atan(Y))


def F_strips(w, c):
    """2D view factor between two parallel strips of width w at a distance c."""
    return sqrt(1 + (c / w) ** 2) - c / w


def effective_eps(F, eps=0.9):
    r = F * eps / (2 - eps)
    return 2 * r / (1 + r)


def run_case(name, eps_holes=0.9, web_k_factor=1.0):
    p = project("oliveira")
    p["test"]["hole_emissivity"] = eps_holes
    geo = pj.geometry(p)
    M = pj.materials(p, geo)
    ts, msh_ = p["test"], p["mesh"]
    m2 = msh.section_mesh(geo, h=msh_["h"] / 1e3, symmetry=False, bands=())
    mat = th.section_material(m2, geo).copy()
    if web_k_factor != 1.0:
        yc = m2.xy[m2.tri][:, :, 1].mean(1)
        W, tc = geo.width, p["geometry"]["face_shell"] / 1e3
        web = (mat == CERAMIC) & (yc > tc) & (yc < W - tc)
        mat[web] = MORTAR
        c = M.ceramic
        M.mortar.k, M.mortar.rho, M.mortar.c = c.k * web_k_factor, c.rho, c.c
        M.mortar.table_k = [(T, v * web_k_factor) for T, v in c.table_k]
        M.mortar.table_rho, M.mortar.table_c = c.table_rho, c.table_c
    f, tf, _ = pj.thermal_programme(p)
    nexp = th.boundary_lengths_2d(m2, faces=("y1",))
    ctr = [Boundary(th.boundary_lengths_2d(m2, faces=("y0",)), f, ts["h_convection"], ts["emissivity"], "exposed"),
           Boundary(nexp, curves.constant(ts["T_initial"]), 4.0, ts["emissivity"], "unexposed")]
    face = np.nonzero(nexp > 0)[0]; w = nexp[face] / nexp[face].sum(); Ti = ts["T_initial"]
    r = th.solve_2d(m2, M, ctr, tf, dt=msh_["dt"], T_initial=Ti, mat=mat, record_every=60.0,
                       probes={"face": face}, stop_if=lambda t, T: (T[face] @ w - Ti >= 140) or (T[face].max() - Ti >= 180))
    I = r.stopped_at / 60 if r.stopped_at else None
    print(f"{name:58s} eps_holes = {eps_holes:.3f}  k webs x {web_k_factor:.2f}  ->  insulation = {I} min", flush=True)
    return {"eps_holes": eps_holes, "web_k_factor": web_k_factor, "I_min": I}


if __name__ == "__main__":
    p = project("oliveira"); g = p["geometry"]
    D = g["width"] - 2 * g["face_shell"]                                   # 68 mm, thickness of the cell
    a = (g["unit_length"] - 2 * g["end_shell"] - (g["holes_per_unit"] - 1) * g["web"]) / g["holes_per_unit"]
    b = g["course_height"] - 2 * g["bed_shell"] if g.get("bed_shell") else g["course_height"]
    F2, F3 = F_strips(a, D), F_rectangles(a, b, D)
    print(f"cell: {a:.1f} x {b:.1f} mm, opposite faces at {D:.0f} mm; 2D view factor = {F2:.3f}, 3D = {F3:.3f}")
    out = {"cell_mm": [a, b, D], "view_factor_2D": F2, "view_factor_3D": F3,
           "effective_eps_2D": effective_eps(F2), "effective_eps_3D": effective_eps(F3)}
    out["reference"] = run_case("reference: infinite parallel plates, eps = 0.9")
    out["F_2D"] = run_case("direct radiation between the faces, 2D view factor", effective_eps(F2))
    out["F_3D"] = run_case("direct radiation between the faces, 3D view factor", effective_eps(F3))
    out["no_radiation"] = run_case("no radiation in the cavities (air conduction only)", 0.0)
    out["webs_half"] = run_case("webs with half the conductivity", 0.9, 0.5)
    out["F_3D_webs_half"] = run_case("3D view factor and webs with half the conductivity", effective_eps(F3), 0.5)
    out["no_radiation_webs_half"] = run_case("no radiation and webs with half the conductivity", 0.0, 0.5)
    json.dump(out, open(os.path.join(HERE, "oliveira_explicit_section.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
