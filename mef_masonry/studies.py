"""High-level functions: joint calibration, 3D simulation of the UFV test
and quick analyses of the 2D section (furnace, and wall under standard fire)."""
import copy
import os
from dataclasses import replace
from typing import Optional

import numpy as np

from . import mesh as msh, thermal as th, mechanical as mech, curves
from .mesh import AIR, CERAMIC, JOINT, CAPPING, STEEL, HOLE
from .parameters import Geometry, Test
from .thermal import Materials, Boundary


# ----------------------------------------------------------------------------
# Prism stiffness and calibration of the dry joints
# ----------------------------------------------------------------------------
def prism_modulus(m: msh.Mesh3D, M: Materials, load_N=1e5, geo: Geometry = Geometry(),
                  middle_third=True):
    """Secant modulus of the prism on the gross area (MPa): from the extensometer in the middle
    third of the outer face and from the shortening between plates."""
    R = mech.solve_3d(m, M, load_N=load_N, include_thermal=False)["load"]
    A_gross = geo.length * geo.width
    sig = load_N / A_gross
    H = geo.prism_height
    eps, z1, z2 = mech.extensometer_strain(m, R, H / 3, 2 * H / 3)
    d = mech.shortening_between_plates(m, R, H)
    return {"E_extensometer": sig / -eps / 1e6, "E_plates": sig / (-d / H) / 1e6,
            "z1": z1, "z2": z2, "result": R}


def series_estimate(M: Materials, E_target_MPa, base_mm, n_joints):
    """Joint E (material, Pa) from the series association on the gross area:
    base/E_target = (base - n t)/E_unit + n t/(E_j r), r = net/gross area."""
    t = M.joint.thickness
    r = M.ceramic.net_area_ratio
    base = base_mm / 1000
    flex_j = base / (E_target_MPa * 1e6) - (base - n_joints * t) / M.ceramic.E_unit_gross
    if flex_j <= 0:
        raise ValueError("Target stiffer than the unit itself; check the data.")
    return n_joints * t / flex_j / r


def calibrate_joint(m: msh.Mesh3D, M: Materials, E_target_MPa=1100.0, geo: Geometry = Geometry(),
                   tol=0.002, max_iter=8, verbose=True):
    """Adjusts M.joint.E (and the capping, if equal) so that the middle-third modulus
    equals E_target_MPa. Uses the secant method on 1/E_joint. Returns (E_joint, history)."""
    base = prism_modulus(m, M, geo=geo)
    n_j = int(sum(1 for c in m.layers if c.kind == "joint"
                  and base["z1"] <= 0.5 * (c.z0 + c.z1) <= base["z2"]))
    Ej = series_estimate(M, E_target_MPa, (base["z2"] - base["z1"]) * 1000, max(n_j, 1))
    hist = []
    x_prev = f_prev = None
    for it in range(max_iter):
        M.joint.E = Ej; M.capping.E = Ej
        Ep = prism_modulus(m, M, geo=geo)["E_extensometer"]
        hist.append((Ej / 1e6, Ep))
        if verbose:
            print(f"  iteration {it}: E_joint = {Ej / 1e6:8.2f} MPa -> E_prism = {Ep:8.1f} MPa")
        f = 1 / Ep - 1 / E_target_MPa
        if abs(Ep - E_target_MPa) / E_target_MPa < tol:
            break
        x = 1 / Ej
        if x_prev is None:
            x_new = x * (1 + 0.3 * np.sign(f))
        else:
            x_new = x - f * (x - x_prev) / (f - f_prev)
        x_prev, f_prev = x, f
        Ej = 1 / x_new
    return Ej, hist


# ----------------------------------------------------------------------------
# UFV test in 3D
# ----------------------------------------------------------------------------
def simulate_test_3d(m: msh.Mesh3D, M: Materials, test: Test, geo: Geometry = Geometry(),
                      dt=30.0, h_ends=0.0, vtk_folder: Optional[str] = None, verbose=False):
    """Heats the prism with the test protocol and applies the peak load at the end
    of the dwell. Returns a dictionary with histories and indicators."""
    f, t_end, t_ramp = curves.ramp_dwell(test.T_initial, test.rate_C_min,
                                             test.T_target, test.dwell_s)
    boundaries = [Boundary(th.side_areas_3d(m), f, test.h_convection, test.emissivity, "sides")]
    if h_ends > 0:
        boundaries.append(Boundary(th.top_base_areas_3d(m), curves.constant(test.T_initial),
                                  h_ends, 0.0, "top and base"))
    H = geo.prism_height
    zm = m.zs[np.argmin(np.abs(m.zs - H / 2))]
    L, Wm = m.m2d.length, m.m2d.model_width
    probes = {
        "surface (mid-height, x = 140 mm)": msh.nearest_node(m.xyz, (0.140, 0.0, zm)),
        "web centre (x = 280 mm, mid-height)": msh.nearest_node(m.xyz, (0.280, Wm, zm)),
    }
    cer_mid = np.unique(m.conn3[(m.mat == CERAMIC) & (np.abs(m.xyz[m.conn3[:, 0], 2] - H / 2) < H / 6)])
    probes["mean of the ceramic in the middle third"] = cer_mid
    rt = th.solve_3d(m, M, boundaries, t_end, dt=dt, T_initial=test.T_initial, probes=probes,
                        record_every=60.0, field_instants=[t_ramp], verbose=verbose)
    T = rt.T_final
    res = mech.solve_3d(m, M, T_nodes=T, load_N=test.load_N)
    Rt, Rc = res["thermal"], res["load"]
    Rtot = mech.combine(Rt, Rc)
    A_gross = geo.length * geo.width
    eps_c, _, _ = mech.extensometer_strain(m, Rc, H / 3, 2 * H / 3)
    ind = {
        "T_target": test.T_target, "time_to_load_min": t_end / 60,
        "T_surface": rt.probes[list(probes)[0]][-1],
        "T_web": rt.probes[list(probes)[1]][-1],
        "T_mean_ceramic_third": rt.probes[list(probes)[2]][-1],
        "T_min_ceramic": T[np.unique(m.conn3[m.mat == CERAMIC])].min(),
        "expansion_between_plates_mm": mech.shortening_between_plates(m, Rt, H) * 1e3,
        "load_shortening_mm": -mech.shortening_between_plates(m, Rc, H) * 1e3,
        "E_middle_third_MPa": test.load_N / A_gross / -eps_c / 1e6,
    }
    for name, R in (("thermal", Rt), ("load", Rc), ("total", Rtot)):
        s, p99, pos = mech.max_tension(m, R)
        sm, p99m, posm = mech.max_tension(m, R, zmin=geo.course_height, zmax=H - geo.course_height)
        ind[f"tension_max_{name}_MPa"] = s / 1e6
        ind[f"tension_p99_{name}_MPa"] = p99 / 1e6
        ind[f"tension_max_{name}_central_courses_MPa"] = sm / 1e6
        ind[f"position_tension_max_{name}_mm"] = tuple(np.round(pos * 1000, 1))
    s1 = Rtot.princ[:, 2]
    vol = m.elem.volumes(m.xyz, m.conn3)
    cer = m.mat == CERAMIC
    ind["ceramic_fraction_tension_above_1MPa"] = float(vol[cer & (s1 > 1e6)].sum() / vol[cer].sum())
    if vtk_folder:
        from .output import save_vtk
        os.makedirs(vtk_folder, exist_ok=True)
        save_vtk(os.path.join(vtk_folder, f"prism_{int(test.T_target)}C.vtk"), m,
                   node_data={"T_C": T, "u_total_mm": Rtot.u * 1e3, "u_thermal_mm": Rt.u * 1e3},
                   elem_data={"material": m.mat.astype(float),
                               "s1_total_MPa": Rtot.princ[:, 2] / 1e6,
                               "s1_thermal_MPa": Rt.princ[:, 2] / 1e6,
                               "s3_total_MPa": Rtot.princ[:, 0] / 1e6,
                               "szz_total_MPa": Rtot.sig[:, 2] / 1e6})
    return {"thermal": rt, "mechanical": res, "total": Rtot, "indicators": ind, "probes": list(probes),
            "t_end_ramp": t_ramp}


# ----------------------------------------------------------------------------
# 2D section: furnace test and wall under standard fire
# ----------------------------------------------------------------------------
def _geo_hole(geo: Geometry, diameter):
    g = copy.copy(geo)
    g.hole_diameter = diameter
    return g


def furnace_section(M: Materials, test: Test, geo: Geometry = Geometry(), h=0.005, dt=30.0,
                diameter=None, y_holes=None):
    """Prism section in the furnace (all faces exposed) to the end of the dwell;
    thermal stresses in generalised plane strain."""
    g = _geo_hole(geo, geo.hole_diameter if diameter is None else diameter)
    holes = None if y_holes is None else [(x, y_holes) for x in g.hole_centres]
    m2 = msh.section_mesh(g, h=h, symmetry=False, bands=(), holes=holes)
    M2 = copy.deepcopy(M); M2.hole_diameter = g.hole_diameter
    f, t_end, _ = curves.ramp_dwell(test.T_initial, test.rate_C_min, test.T_target, test.dwell_s)
    ctr = [Boundary(th.boundary_lengths_2d(m2), f, test.h_convection, test.emissivity)]
    r = th.solve_2d(m2, M2, ctr, t_end, dt=dt, T_initial=test.T_initial)
    s = mech.solve_section_gps(m2, M2, r.T_final)
    sol = m2.mark != HOLE
    solid_nodes = np.unique(m2.tri[sol])
    a = s["area"]
    return {"T_min": r.T_final[solid_nodes].min(), "T_mean": float((r.T_final[m2.tri].mean(1)[sol] * a[sol]).sum() / a[sol].sum()),
            "s1_max_MPa": np.nanmax(s["s1"]) / 1e6, "s1_p99_MPa": np.nanpercentile(s["s1"][sol], 99) / 1e6,
            "fraction_s1_above_1MPa": float(a[sol & (s["s1"] > 1e6)].sum() / a[sol].sum()),
            "net_area_ratio": float(a[sol].sum() / (g.length * g.width)),
            "mesh": m2, "T": r.T_final, "stresses": s}


def wall_section_iso834(M: Materials, geo: Geometry = Geometry(), h=0.005, dt=30.0, t_max_min=240.0,
                        diameter=None, y_holes=None, h_exposed=25.0, eps_exposed=0.8,
                        h_unexposed=4.0, eps_unexposed=0.8, T_amb=20.0, t_min_stress=30.0):
    """BCEPQ wall (140 mm thick) with the ISO 834 standard fire on face
    y = 0 and ambient conditions on face y = W; faces x = 0 and x = L adiabatic
    (repetition along the wall). Insulation: mean rise of 140 °C or
    maximum of 180 °C on the unexposed face (NBR 5628 / ISO 834-1).
    Holes treated as closed cavities (no vertical convection)."""
    g = _geo_hole(geo, geo.hole_diameter if diameter is None else diameter)
    holes = None if y_holes is None else [(x, y_holes) for x in g.hole_centres]
    m2 = msh.section_mesh(g, h=h, symmetry=False, bands=(), holes=holes)
    M2 = copy.deepcopy(M); M2.hole_diameter = g.hole_diameter
    M2.ceramic.T_ref = T_amb
    exp = th.boundary_lengths_2d(m2, faces=("y0",))
    nexp = th.boundary_lengths_2d(m2, faces=("y1",))
    ctr = [Boundary(exp, curves.iso834(T_amb), h_exposed, eps_exposed, "exposed"),
           Boundary(nexp, curves.constant(T_amb), h_unexposed, eps_unexposed, "unexposed")]
    face = np.nonzero(nexp > 0)[0]
    w = nexp[face] / nexp[face].sum()

    def criterion(t, T):
        return (T[face] @ w - T_amb >= 140.0) or (T[face].max() - T_amb >= 180.0)
    r = th.solve_2d(m2, M2, ctr, t_max_min * 60, dt=dt, T_initial=T_amb, stop_if=criterion,
                       probes={"unexposed face (mean)": face}, record_every=60.0,
                       field_instants=[t_min_stress * 60])
    out = {"t_insulation_min": (r.stopped_at / 60 if r.stopped_at else np.inf),
           "net_area_ratio": float(np.abs(m2.areas())[m2.mark != HOLE].sum() / (g.length * g.width)),
           "mesh": m2, "result": r}
    if t_min_stress * 60 in r.fields:
        s = mech.solve_section_gps(m2, M2, r.fields[t_min_stress * 60])
        out["s1_max_MPa"] = np.nanmax(s["s1"]) / 1e6
        out["s1_p99_MPa"] = np.nanpercentile(s["s1"][m2.mark != HOLE], 99) / 1e6
    return out


# ----------------------------------------------------------------------------
# Non-linear analysis of the prism
# ----------------------------------------------------------------------------
def nonlinear_modulus(mod, geo: Geometry = Geometry(), delta=2e-6):
    """Prism modulus (MPa, gross area, middle third of the outer face) in the
    elastic range of the non-linear model, with the joints closed."""
    from .nonlinear import NonlinearModel  # noqa
    m = mod.m
    state = mod.state()
    presc = 3 * mod.top + 2
    P, _, _ = mod._step(presc, mod.u[presc] - delta, np.zeros(mod.ndof))
    u = mod.u.reshape(-1, 3)
    H = geo.prism_height
    z1 = m.zs[np.argmin(np.abs(m.zs - H / 3))]; z2 = m.zs[np.argmin(np.abs(m.zs - 2 * H / 3))]
    face = np.nonzero(np.isclose(m.xyz[:, 1], 0.0))[0]
    uz = lambda z: u[face[np.isclose(m.xyz[face, 2], z)], 2].mean()
    eps = (uz(z2) - uz(z1)) / (z2 - z1)
    mod.restore(state)
    return P / (geo.length * geo.width) / -eps / 1e6


def calibrate_joint_nonlinear(mod, E_target_MPa=1100.0, geo: Geometry = Geometry(), tol=0.002, verbose=True):
    """Adjusts the E of the thin-layer joints (non-linear model) so that the prism
    reproduces E_target in the middle third. Secant on 1/E_joint."""
    def apply(Ej):
        mod.Ej[:] = Ej
        mod.En = mod.Ej.copy(); mod.Gs = mod.Ej / (2 * (1 + mod.nuj)); mod.closed[:] = True
        mod._pc = None
    Ej = float(mod.Ej[0])
    x_prev = f_prev = None
    for it in range(8):
        apply(Ej)
        Ep = nonlinear_modulus(mod, geo)
        if verbose:
            print(f"  iteration {it}: E_joint = {Ej / 1e6:7.2f} MPa -> E_prism = {Ep:7.1f} MPa")
        f = 1 / Ep - 1 / E_target_MPa
        if abs(Ep - E_target_MPa) / E_target_MPa < tol:
            break
        x = 1 / Ej
        x_new = x * (1 + 0.3 * np.sign(f)) if x_prev is None else x - f * (x - x_prev) / (f - f_prev)
        x_prev, f_prev = x, f
        Ej = 1 / x_new
    apply(Ej)
    return Ej


def heating_fields(m, M: Materials, test: Test, n_ramp=6, n_dwell=2, dt=30.0):
    """Runs the 3D thermal analysis and returns {time: field} at instants of the ramp and the dwell."""
    f, t_end, t_ramp = curves.ramp_dwell(test.T_initial, test.rate_C_min, test.T_target, test.dwell_s)
    times = list(np.linspace(0, t_ramp, n_ramp + 1)[1:]) + list(np.linspace(t_ramp, t_end, n_dwell + 1)[1:])
    times = [round(t / dt) * dt for t in times]
    ctr = [Boundary(th.side_areas_3d(m), f, test.h_convection, test.emissivity)]
    r = th.solve_3d(m, M, ctr, t_end, dt=dt, T_initial=test.T_initial, field_instants=times)
    return r.fields


def chord_modulus(mod, P_ref_N, geo: Geometry = Geometry(), f1=0.05, f2=0.30):
    """Chord modulus of the prism (MPa, gross area) between f1 and f2 of P_ref, with the
    strain of the virtual extensometer (mod.set_extensometer)."""
    h = [p for p in mod.history if p.phase == "load"]
    P = np.r_[0.0, [p.load_N for p in h]]; e = np.r_[0.0, [-p.eps_ext for p in h]]
    if len(P) < 3 or P.max() < f2 * P_ref_N:
        return float("nan")
    i = int(np.argmax(P >= f2 * P_ref_N)) + 1
    e1, e2 = np.interp([f1 * P_ref_N, f2 * P_ref_N], P[:i], e[:i])
    return (f2 - f1) * P_ref_N / (geo.length * geo.width) / (e2 - e1) / 1e6


def nonlinear_test(m, M: Materials, Pd, Pj, fields=None, delta_max=0.9e-3, n_steps=30,
                      max_load=None, geo: Geometry = Geometry(), verbose=True, tol_u=1e-4, tol_d=5e-3):
    """One prism test with the non-linear model: heating (fields
    {time: T_nodes}, or uniform temperature T_ref) under the seating load
    and then displacement-controlled loading beyond the peak. Returns the
    model (with the history and the middle-third extensometer) and (deltas, loads)."""
    from .nonlinear import NonlinearModel
    mod = NonlinearModel(m, M, Pd, Pj, verbose=verbose)
    mod.tol_u, mod.tol_d = tol_u, tol_d
    H = geo.prism_height
    mod.set_extensometer(H / 3, 2 * H / 3)
    if fields is None:
        fields = {0.0: np.full(m.n_nodes, M.ceramic.T_ref)}
    mod.heat(fields)
    deltas, loads = mod.apply_load(delta_max, n_steps, max_load=max_load)
    return mod, (deltas, loads)


def calibrate_gap(m, M: Materials, Pd, Pj, P_ref_N, E_target_MPa=1100.0, gaps=(0.024e-3, 0.032e-3),
                   seeds=(1, 2, 3), geo: Geometry = Geometry(), verbose=True):
    """Standard deviation of the bed-joint gaps that reproduces, as the mean over the
    seeds of the random field, the 5-30% chord modulus of P_ref (E_target).
    Each case loads the prism to 0.32 P_ref. The stiffness varies by about 20%
    between seeds, hence the mean. Returns (std, [(std, mean E, [E per seed])])."""
    from dataclasses import replace
    pts = []
    for sg in gaps:
        Es = []
        for sd in seeds:
            mod, _ = nonlinear_test(m, M, Pd, replace(Pj, gap_std=sg, seed=sd), delta_max=0.45e-3,
                                       n_steps=30, max_load=0.32 * P_ref_N, geo=geo, verbose=False)
            Es.append(chord_modulus(mod, P_ref_N, geo))
        pts.append((sg, float(np.mean(Es)), Es))
        if verbose:
            print(f"    gap {sg * 1e3:.3f} mm -> E_chord = {np.mean(Es):7.1f} MPa (seeds: "
                  + ", ".join(f"{E:.0f}" for E in Es) + ")", flush=True)
    a = np.array([(p[0], p[1]) for p in pts])
    if len(a) == 1:
        return float(a[0, 0]), pts
    c = np.polyfit(a[:, 1], a[:, 0], 1)                 # std linear in E over the tested range
    if not (a[:, 1].min() <= E_target_MPa <= a[:, 1].max()) and verbose:
        print("    warning: E_target outside the tested range; value extrapolated", flush=True)
    return float(np.polyval(c, E_target_MPa)), pts


def calibrate_fc(m, M: Materials, Pd, Pj, P_target_N, geo: Geometry = Geometry(), tol=0.03, max_it=3, verbose=True):
    """Effective compressive strength of the ceramic (fc) that reproduces the peak
    load P_target at 30 °C. The peak load of the model is almost proportional to fc;
    the iteration uses this proportion. Returns (fc, peak, [(fc, peak)])."""
    from dataclasses import replace
    fc = Pd.fc
    pts = []
    for it in range(max_it):
        mod, (d, P) = nonlinear_test(m, M, replace(Pd, fc=fc), Pj, geo=geo, verbose=False)
        pts.append((fc, P.max()))
        if verbose:
            print(f"    fc = {fc / 1e6:5.2f} MPa -> peak load {P.max() / 1e3:6.1f} kN", flush=True)
        if abs(P.max() / P_target_N - 1) < tol:
            break
        fc = fc * P_target_N / P.max()
    return pts[-1][0], pts[-1][1], pts
