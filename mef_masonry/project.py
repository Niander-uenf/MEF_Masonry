"""Project of the program: JSON inputs -> model -> analysis -> results.

A project is a dictionary with six groups (geometry, ceramic, model, test,
mesh, sensitivity). Inputs are in laboratory units (mm, MPa, °C, min); the
conversion to SI is done here. SCHEMA describes each field (label, unit,
help, options) and serves both the interface and the validation.

    from mef_masonry.project import default_project, run
    p = default_project()
    p["test"]["kind"] = "hot"; p["test"]["T_target"] = 600
    res = run(p, "results/app/test")
"""
import copy
import json
import math
import os
import time
from dataclasses import replace
from typing import Callable, Dict, List

import numpy as np

from . import mesh as msh, thermal as th, mechanical as mech, curves, studies
from .mesh import AIR, CERAMIC, JOINT, CAPPING, STEEL, MORTAR, GROUT
from .parameters import Geometry, CeramicDamage, FrictionJoint, CALIBRATION_30C
from .thermal import Materials, Boundary
from .output import save_vtk, save_csv
from . import visual
from . import library as _lib
from . import legacy, __version__

VERSION = "2.0"   # format of the project file (2.0: English names; 1.0: Portuguese, converted by legacy.py)

# ----------------------------------------------------------------------------
# Input schema
# ----------------------------------------------------------------------------
def _n(id, label, default, unit_str="", help_text="", minimum=None, maximum=None, step=None, **kw):
    d = {"id": id, "label": label, "kind": "number", "default": default, "unit": unit_str, "help": help_text}
    if minimum is not None: d["min"] = minimum
    if maximum is not None: d["max"] = maximum
    if step is not None: d["step"] = step
    d.update(kw)
    return d


def _o(id, label, default, options, help_text="", **kw):
    return {"id": id, "label": label, "kind": "option", "default": default, "options": options, "help": help_text, **kw}


def _b(id, label, default, help_text="", **kw):
    return {"id": id, "label": label, "kind": "boolean", "default": default, "help": help_text, **kw}


def _t(id, label, columns, help_text="", **kw):
    return {"id": id, "label": label, "kind": "table", "default": [], "columns": columns, "help": help_text, **kw}


NL = {"model.analysis": ["nonlinear"]}
LIN = {"model.analysis": ["linear"]}
AQ = {"test.kind": ["hot", "residual", "transient"]}
NLP = {"model.analysis": ["nonlinear", "wall", "sensitivity"]}
CRP = {**NLP, "model.creep": [True]}
WAL = {"model.analysis": ["wall", "sensitivity"]}
WALLS = ("wall", "sensitivity")
P3D = {"model.analysis": ["thermal", "linear", "nonlinear"]}
MOR = {"geometry.joint_type": ["mortar"]}
DRY = {"geometry.joint_type": ["dry"]}
CIRC = {"geometry.hole_shape": ["circular"]}
RECT = {"geometry.hole_shape": ["rectangular"]}
HOR = {"geometry.hole_shape": ["horizontal"]}
VERT = {"geometry.hole_shape": ["circular", "rectangular"]}
RECT_HOR = {"geometry.hole_shape": ["rectangular", "horizontal"]}
GRT = {"geometry.grouting": ["all", "alternate"]}

_HELP_SENS = {
    "k": "Measurement uncertainty and variation between batches; also applies to the k(T) table, if any",
    "c": "Includes the uncertainty of the moisture peak; also applies to the c(T) table",
    "eps_holes": "Fired clay: about 0.85 to 0.95; the lower limit covers smoother surfaces",
    "h": "EN 1991-1-2 uses 25 W/m²·K on the face exposed to the standard fire",
    "eps": "EN 1991-1-2 recommends 0.8 for the surface of the member",
    "alpha": "Fired clay: about 4 to 8 ×10⁻⁶/°C; replaces α and scales the thermal strain table",
    "ktr": "From no creep (0) to the upper limit for concrete (2.35; Anderberg and Thelandersson, 1976)",
    "strength": "Scatter of the prism strength; multiplies the stresses of the wall law",
    "stiffness": "Scatter of the prism modulus; multiplies the modulus of the wall law",
    "imp": "Initial deflection H/n at mid-height; the amplitude varies linearly between the two limits",
    "exc": "Applied equally at top and base; positive towards the exposed face",
}


def _sensitivity_fields():
    from .sensitivity import PARAMETERS
    cs = [
        _n("trajectories", "Trajectories r", 10, "", "10 to 20 is usual; each trajectory costs (parameters + 1) "
           "wall runs", minimum=2, maximum=100, step=1),
        _n("levels", "Grid levels p", 4, "", "Even; with 4 levels, each step covers 2/3 of the range",
           minimum=4, maximum=10, step=1),
        _n("seed", "Seed for drawing the trajectories", 1, "", "The same seed repeats the same plan",
           minimum=0, step=1),
        _n("processes", "Parallel processes", 0, "", "0 = automatic (computer cores minus one)",
           minimum=0, maximum=64, step=1),
        _b("fast_mesh", "Fast mesh in the screening (10 mm section and 16 beam elements)", True,
           "About 6 times faster than the recommended mesh (5 mm and 32 elements). For the BCEPQ wall, the "
           "times changed by less than 7 min and the order of the parameters was the same with both meshes"),
    ]
    for d in PARAMETERS:
        ch = d["key"]
        sset = {"sensitivity.var_" + ch: [True]}
        cs.append(_b("var_" + ch, f"Vary: {d['name'][0].lower() + d['name'][1:]}", True, _HELP_SENS.get(ch, ""),
                     line_break=True))
        if d["kind"] == "rel":
            cs.append(_n("range_" + ch, f"{d['name']}: range", d["default"][0], "± %", minimum=0, maximum=90,
                         show_if=sset))
        else:
            lim = ({"minimum": 0, "maximum": 1} if ch in ("eps_holes", "eps") else
                   {"minimum": 0} if ch in ("h", "alpha", "ktr", "imp") else {})
            cs.append(_n("min_" + ch, f"{d['name']}: from", d["default"][0], d["unit"], show_if=sset, **lim))
            cs.append(_n("max_" + ch, f"{d['name']}: to", d["default"][1], d["unit"], show_if=sset, **lim))
    return cs


SCHEMA = [
    {"id": "geometry", "title": "Geometry", "description":
        "Specimen of solid units, with vertical round holes (aligned at mid-width), with "
        "vertical rectangular cells (concrete and structural clay blocks, with or without grout) or with horizontal "
        "holes (infill blocks): single unit (1 unit, 1 course), prism (default) or small wall (more "
        "units per course and more courses), with dry or mortar joints.", "fields": [
        _o("unit_library", "Library unit", "bcepq", _lib.options(_lib.UNITS),
           "Fills in the geometry, the joint and the materials; then change what you have measured"),
        _n("unit_length", "Unit length", 280, "mm", "Centre to centre: with mortar, unit + joint",
           minimum=40),
        _n("width", "Unit width", 140, "mm", minimum=40),
        _n("course_height", "Course height", 60, "mm", "Centre to centre: with mortar, unit + joint", minimum=20),
        _o("hole_shape", "Perforation", "circular",
           [["circular", "Vertical round holes (or solid unit)"],
            ["rectangular", "Vertical rectangular cells (concrete, structural clay)"],
            ["horizontal", "Horizontal holes (infill blocks)"]],
           "The horizontal holes run through the unit along its length"),
        _n("hole_diameter", "Hole diameter", 84, "mm", "0 = solid unit", minimum=0, show_if=CIRC),
        _n("holes_per_unit", "Holes or cells per unit", 2, "", "Even when there is a half unit or running bond; "
           "0 = solid unit", minimum=0, step=1, show_if=VERT),
        _n("holes_width", "Holes across the unit thickness", 2, "", "Columns of horizontal holes", minimum=1, step=1,
           show_if=HOR),
        _n("holes_height", "Holes over the unit height", 4, "", "Rows of horizontal holes (8 holes = 2 × 4)",
           minimum=1, step=1, show_if=HOR),
        _n("face_shell", "Thickness of the face shells", 25, "mm",
           "Longitudinal walls of the unit", minimum=3, show_if=RECT_HOR),
        _n("web", "Web thickness", 25, "mm", "Internal walls between the cells or between the holes",
           minimum=3, show_if=RECT_HOR),
        _n("bed_shell", "Thickness of the bed shells", 8, "mm",
           "Top and bottom walls, on which the mortar is laid", minimum=3, show_if=HOR),
        _n("end_shell", "Thickness of the end shells", 25, "mm", "Transverse walls at the ends of the unit",
           minimum=3, show_if=RECT),
        _n("wall_cells", "Small cells per unit in each face shell", 0, "",
           "Clay blocks with hollow shells; 0 = solid shells", minimum=0, step=1, show_if=RECT),
        _n("min_shell", "Wall around the small cells", 8, "mm", minimum=3,
           show_if={"geometry.hole_shape": ["rectangular"]}),
        _o("grouting", "Grout", "none",
           [["none", "No grout"], ["all", "All holes or cells"],
            ["alternate", "Alternate holes or cells (one yes, one no)"]],
           "The grout fills the holes or the main cells over the full height", show_if=VERT),
        _n("units_per_course", "Units per course", 1.5, "", "Multiple of 0.5 (1.5 = one and a half units)", minimum=0.5, step=0.5),
        _n("n_courses", "Number of courses", 5, "", minimum=1, step=1),
        _o("bond", "Bond", "running_bond",
           [["running_bond", "Running bond (staggered head joints)"], ["stack_bond", "Stack bond"]]),
        _o("joint_type", "Joint type", "dry",
           [["dry", "Dry (thin layer with contact, gaps and friction)"],
            ["mortar", "Mortar (continuous layer of real thickness, bonded to the unit)"]]),
        _n("joint_thickness", "Joint thickness", 2, "mm",
           "Dry: thin layer (2 mm). Mortar: real thickness (about 10 mm)", minimum=0.5),
        _o("bedding", "Mortar bedding", "total",
           [["total", "Over the whole face of the unit"], ["lateral", "On the face shells only (holes and webs without mortar)"]],
           show_if=MOR),
        _b("filled_head_joints", "Head joints filled with mortar", True, show_if=MOR),
        _b("capping", "Capping at top and base", True),
        _n("plate_thickness", "Thickness of the steel plates", 25, "mm", "Used only in the thermal analysis", minimum=0),
    ]},
    {"id": "ceramic", "title": "Materials", "description":
        "Properties of the unit and the mortar. The temperature-dependent tables are optional.", "fields": [
        _o("material_library", "Unit material (library)", "bcepq", _lib.options(_lib.MATERIALS),
           "Fills in the properties with their source; the others stay as they are"),
        _n("k", "Thermal conductivity", 0.52, "W/m·K", minimum=0.01),
        _n("rho", "Density", 1210, "kg/m³", minimum=100),
        _n("c", "Specific heat", 1285, "J/kg·K", minimum=100),
        _n("E_unit_gross", "Unit modulus (gross area)", 2801, "MPa",
           "The material modulus is obtained by dividing by the net-to-gross area ratio of the geometry", minimum=1),
        _n("nu", "Poisson's ratio", 0.18, "", minimum=0, maximum=0.49),
        _n("alpha", "Thermal expansion coefficient", 6.0, "×10⁻⁶/°C", "Provisional for the BCEPQ", minimum=0),
        _t("table_k", "k as a function of T", ["T (°C)", "k (W/m·K)"]),
        _t("table_E_factor", "Factor on E as a function of T", ["T (°C)", "factor"]),
        _t("table_thermal_strain", "Free thermal strain (replaces α)", ["T (°C)", "strain (‰)"]),
        _t("table_c", "Specific heat as a function of T", ["T (°C)", "c (J/kg·K)"],
           "Use for the moisture peak of concrete and soil-cement"),
        _t("table_rho", "Density as a function of T", ["T (°C)", "ρ (kg/m³)"]),
        _o("mortar_library", "Mortar (library)", "1:1:6", _lib.options(_lib.MORTARS),
           "Mohamad, Lourenço and Roman (2006); thermal properties from NBR 15220-2", show_if=MOR),
        _n("mor_fc", "Compressive strength of the mortar", 4.2, "MPa", minimum=0.1, show_if=MOR),
        _n("mor_ft", "Tensile strength of the mortar", 0.78, "MPa", "Estimate: 0.3 fc^(2/3)", minimum=0.01,
           show_if=MOR),
        _n("mor_E", "Mortar modulus", 4033, "MPa", "Initial tangent", minimum=1, show_if=MOR),
        _n("mor_nu", "Mortar Poisson's ratio", 0.20, "", minimum=0, maximum=0.49, show_if=MOR),
        _n("mor_k", "Thermal conductivity of the mortar", 1.15, "W/m·K", minimum=0.01, show_if=MOR),
        _n("mor_rho", "Density of the mortar", 2000, "kg/m³", minimum=100, show_if=MOR),
        _n("mor_c", "Specific heat of the mortar", 1000, "J/kg·K", minimum=100, show_if=MOR),
        _n("mor_alpha", "Thermal expansion of the mortar", 10.0, "×10⁻⁶/°C", "Provisional", minimum=0, show_if=MOR),
        _o("grout_library", "Grout (library)", "G1", _lib.options(_lib.GROUTS),
           "Camacho et al. (2015); E from EN 1992-1-1; thermal properties from NBR 15220-2", show_if=GRT),
        _n("grt_fc", "Compressive strength of the grout", 17.7, "MPa", minimum=0.5, show_if=GRT),
        _n("grt_ft", "Tensile strength of the grout", 2.04, "MPa", "Estimate: 0.3 fc^(2/3)", minimum=0.05,
           show_if=GRT),
        _n("grt_E", "Grout modulus", 26100, "MPa", "Estimate: 22 (fc/10)^0.3 GPa (EN 1992-1-1)", minimum=1,
           show_if=GRT),
        _n("grt_nu", "Grout Poisson's ratio", 0.20, "", minimum=0, maximum=0.49, show_if=GRT),
        _n("grt_k", "Thermal conductivity of the grout", 1.75, "W/m·K", minimum=0.01, show_if=GRT),
        _n("grt_rho", "Density of the grout", 2400, "kg/m³", minimum=100, show_if=GRT),
        _n("grt_c", "Specific heat of the grout", 1000, "J/kg·K", minimum=100, show_if=GRT),
        _b("grt_en1992", "Hot grout as the siliceous concrete of EN 1992-1-2", True,
           "Thermal strain and reduction of E, fc and ft with temperature (Table 3.1)", show_if=GRT),
        _n("grt_alpha", "Thermal expansion of the grout", 10.0, "×10⁻⁶/°C", "Provisional", minimum=0,
           show_if={**GRT, "ceramic.grt_en1992": [False]}),
    ]},
    {"id": "model", "title": "Models and elements", "description":
        "Analysis type, finite element and constitutive models.", "fields": [
        _o("analysis", "Analysis", "nonlinear", [
            ["thermal", "3D thermal (temperature fields)"],
            ["linear", "3D linear thermo-elastic (stresses and strains)"],
            ["nonlinear", "3D non-linear to failure (damage and friction joints)"],
            ["section_2d", "2D section (thermal and generalised plane-strain stresses)"],
            ["wall", "Wall: vertical strip with second order (2D section + fibre beams)"],
            ["sensitivity", "Wall sensitivity under ISO 834 (Morris: which inputs weigh most)"]],
           "Wall: the height comes from the number of courses and the length from the units per course; the section "
           "of one unit (with its holes) is heated in 2D and each ceramic triangle becomes a fibre. "
           "Sensitivity: runs the wall many times, varying the uncertain inputs over the ranges of the "
           "Sensitivity page, and ranks those that change the fire resistance and insulation times most"),
        _o("element", "Finite element", "W6", [
            ["W6", "6-node linear wedge (W6)"],
            ["H8", "8-node linear hexahedron (H8)"],
            ["W15", "15-node quadratic wedge (W15)"],
            ["T3", "3-node linear triangle (2D section)"]],
           "W6: default of the paper. H8: hexahedra on the quadrilaterals of the O-grid mesh, half the "
           "elements. W15: quadratic wedge, about 4.7 times more nodes than W6 on the same mesh; "
           "use coarser meshes (for example, 16 × 30 mm). In the 2D section the element is always T3",
           show_if={"model.analysis": ["thermal", "linear", "nonlinear", "section_2d"]}),
        _o("wall_law", "Wall material", "bcepq", [
            ["bcepq", "BCEPQ prism curves at 30, 600, 800 and 1000 °C (calibrated non-linear model)"],
            ["prism", "Calibrate now with the non-linear prism model at 30 °C (a few more minutes)"]],
           "Homogenised stress–strain law without tension (the joints open). With the prism calibration, "
           "the 30 °C curve comes from the current inputs and the other temperatures follow the BCEPQ ratios",
           show_if=WAL),
        _n("wall_strength_factor", "Factor on the masonry strength (wall law)", 1.0, "",
           "Multiplies the stresses of the law at all temperatures; 1 = law without adjustment", minimum=0.1, show_if=WAL),
        _n("wall_stiffness_factor", "Factor on the masonry stiffness (wall law)", 1.0, "",
           "Multiplies the modulus of the law (the peak strains change in inverse ratio); 1 = law without adjustment",
           minimum=0.1, show_if=WAL),
        _n("E_joint", "E of the joint layer (linear model)", 84.4, "MPa",
           "Calibrated for E_prism = 1,100 MPa with the 8 × 12 mm mesh", minimum=0.1, show_if={"model.analysis": ["linear"], "geometry.joint_type": ["dry"]}),
        _b("calibrate_joint", "Calibrate the joint E from the prism modulus", False, show_if={"model.analysis": ["linear"], "geometry.joint_type": ["dry"]}),
        _n("E_prism_target", "Prism modulus for the calibration", 1100, "MPa", minimum=10, show_if={"model.analysis": ["linear"], "geometry.joint_type": ["dry"]}),
        _o("ceramic_model", "Damage model of the unit, mortar and grout", "mazars_damage", [
            ["mazars_damage", "Tensile and compressive damage with crack closure (Mazars)"],
            ["isotropic_damage", "Damage without crack closure"]], show_if=NL),
        _n("ft", "Tensile strength of the unit ft", 0.66, "MPa", "Gooch et al. (2023)", minimum=0.01, show_if=NL),
        _n("fc", "Compressive strength of the unit fc", 2.48, "MPa",
           "2.48 = effective, calibrated at 30 °C; 5.76 = derived from the unit (BPQ)", minimum=0.1, show_if=NLP),
        _n("du_t", "Tensile ductility index G_f/ft", 0.018, "mm", minimum=0.001, show_if=NL),
        _n("du_c", "Compressive ductility index G_c/fc", 0.33, "mm", minimum=0.01, show_if=NL),
        _t("table_ft_factor", "Factor on ft as a function of T", ["T (°C)", "factor"], show_if=NL),
        _t("table_fc_factor", "Factor on fc as a function of T", ["T (°C)", "factor"], show_if=NL),
        _n("mu", "Friction of the dry joint μ", 0.60, "", "Zero cohesion", minimum=0, show_if={"model.analysis": ["nonlinear"], "geometry.joint_type": ["dry"]}),
        _n("gap_std", "Standard deviation of the joint gaps", 0.028, "mm", "0 = perfect contact", minimum=0,
           show_if={"model.analysis": ["nonlinear"], "geometry.joint_type": ["dry"]}),
        _n("initial_contact", "Area in contact without load", 20, "%", minimum=1, maximum=100, show_if={"model.analysis": ["nonlinear"], "geometry.joint_type": ["dry"]}),
        _n("corr_length", "Correlation length of the gaps", 40, "mm", minimum=5, show_if={"model.analysis": ["nonlinear"], "geometry.joint_type": ["dry"]}),
        _n("E_contact", "E of the dry joint and of the capping where there is contact", 3907, "MPa", minimum=1, show_if=NL),
        _n("seed", "Seed of the gap field", 1, "", minimum=0, step=1, show_if={"model.analysis": ["nonlinear"], "geometry.joint_type": ["dry"]}),
        _b("creep", "Transient creep (load-induced thermal strain)", True,
           "On first heating under compression, ε_tr = φ(T) σ / fc20, with φ(T) = k_tr ε_th(T) "
           "(Anderberg and Thelandersson, 1976) or a table measured in the transient test", show_if=NLP),
        _n("k_tr_unit", "Coefficient k_tr of the unit", 0.0, "",
           "Fired clay: no published data (0). Concrete with the EN 1992-1-2 curves: 0, because "
           "creep is already included in them. Concrete without those curves: 1.8 to 2.35", minimum=0, show_if=CRP),
        _t("table_creep_unit", "φ(T) of the unit, from the transient test (replaces k_tr)", ["T (°C)", "φ (‰)"],
           "φ = (transient test strain − free − elastic) × fc20 / σ", show_if=CRP),
        _n("k_tr_mortar", "Coefficient k_tr of the mortar", 2.0, "",
           "Prakash et al. (2020), cement mortar", minimum=0, show_if={**CRP, **MOR}),
        _n("k_tr_grout", "Coefficient k_tr of the grout", 0.0, "",
           "0 with the EN 1992-1-2 curves (creep already implicit)", minimum=0, show_if={**CRP, **GRT}),
    ]},
    {"id": "test", "title": "Test and load", "description":
        "Thermal programme of the furnace and loading.", "fields": [
        _o("kind", "Test type", "hot", [
            ["ambient", "Compression at room temperature"],
            ["hot", "Hot compression (ramp, dwell and loading in the furnace)"],
            ["residual", "Residual compression (heat, cool, then load)"],
            ["iso834", "ISO 834 standard fire on one face, under constant load"],
            ["transient", "Transient test: constant load and heating in the furnace to failure"]]),
        _n("T_initial", "Initial temperature", 30, "°C"),
        _n("rate", "Heating rate", 5, "°C/min", minimum=0.1, show_if=AQ),
        _n("T_target", "Target furnace temperature", 600, "°C", show_if=AQ),
        _n("dwell", "Dwell at the target temperature", 30, "min", minimum=0, show_if=AQ),
        _n("cooling_rate", "Cooling rate", 2, "°C/min", minimum=0.1, show_if={"test.kind": ["residual"]}),
        _n("final_dwell", "Wait after cooling", 60, "min", minimum=0, show_if={"test.kind": ["residual"]}),
        _n("iso_duration", "Duration of the ISO 834 exposure", 120, "min", minimum=5, show_if={"test.kind": ["iso834"]}),
        _n("h_convection", "Furnace convection coefficient", 10, "W/m²·K", "ISO 834: 25 on the exposed face", minimum=0),
        _n("emissivity", "Surface emissivity", 0.8, "", minimum=0, maximum=1),
        _n("hole_emissivity", "Emissivity of the hole walls", 0.9, "", minimum=0, maximum=1),
        _n("load", "Load", 97.9, "kN",
           "Linear: load applied at the end of the thermal programme. ISO 834 and transient test: constant load. "
           "Non-linear (other tests): only the reference for the chord modulus; failure comes from the analysis. "
           "Wall: total load over the wall length (ISO 834); at 30 °C, only the reference", minimum=0),
        _n("seating_load", "Seating load during heating", 2, "kN",
           "In the transient test and in the 3D ISO 834, the constant load applies", minimum=0, show_if=P3D),
        _o("base_support", "Base support", "pinned",
           [["pinned", "Pinned (free rotation)"], ["fixed", "Fixed (rotation restrained)"]],
           show_if=WAL),
        _o("top_support", "Top support", "pinned",
           [["pinned", "Pinned (free rotation)"], ["fixed", "Fixed (rotation restrained)"]],
           "The top always moves vertically under the load", show_if=WAL),
        _n("top_eccentricity", "Load eccentricity at the top", 0, "mm",
           "Positive towards the exposed face", show_if=WAL),
        _n("base_eccentricity", "Eccentricity of the reaction at the base", 0, "mm",
           "Positive towards the exposed face", show_if=WAL),
        _n("imperfection", "Initial imperfection: height divided by", 450, "",
           "Sinusoidal initial deflection at mid-height, towards the exposed face (H/450 is a usual value); "
           "0 = perfect wall", minimum=0, show_if=WAL),
        _n("load_levels", "Load levels for the time × load curve", 5, "",
           "Repeats the ISO 834 run at 10% to 90% of the cold capacity; 0 = the given load only", minimum=0, maximum=12,
           step=1, show_if={"model.analysis": ["wall"], "test.kind": ["iso834"]}),
    ]},
    {"id": "mesh", "title": "Mesh and solution", "description":
        "Discretisation and control of the numerical solution.", "fields": [
        _n("h", "Element size in plan", 12, "mm", "8 mm in the thermal and linear analyses of the paper; 12 mm in the non-linear one; "
           "5 mm in the 2D section and the wall",
           minimum=2),
        _n("dz", "Maximum height of the element layers", 20, "mm", minimum=2, show_if=P3D),
        _n("n_elements", "Beam elements over the height", 32, "",
           "Even number; with 32, the fire resistance is within 1% of that obtained with 64", minimum=4, step=1,
           show_if=WAL),
        _n("dt", "Thermal time step", 30, "s", minimum=1),
        _n("delta_max", "Maximum platen displacement", 0.8, "mm", minimum=0.05, show_if=NL),
        _n("increment", "Displacement increment", 0.02, "mm", minimum=0.001, show_if=NL),
        _n("tol_u", "Displacement tolerance", 1e-4, "", minimum=1e-8, show_if=NL),
        _n("tol_d", "Damage tolerance", 5e-3, "", minimum=1e-6, show_if=NL),
        _n("n_fields_ramp", "Thermal fields during the ramp", 6, "", minimum=1, step=1, show_if=NL),
        _n("n_fields_dwell", "Thermal fields during the dwell", 2, "", minimum=1, step=1, show_if=NL),
        _n("iso_interval", "Interval between fields under ISO 834", 5, "min", minimum=1, show_if=NL),
    ]},
    {"id": "sensitivity", "title": "Sensitivity", "description":
        "Global sensitivity analysis of the wall under ISO 834 by the Morris method (elementary effects). "
        "Each selected parameter varies over the given range, around the values of the other pages; the program "
        "runs the wall r × (parameters + 1) times and ranks the parameters by the mean effect μ* on the fire "
        "resistance and insulation times. The default ranges are starting points: replace them with the uncertainty "
        "you know for each input.", "fields": _sensitivity_fields()},
]


def default_project() -> Dict:
    p = {"name": "BCEPQ prism", "version": VERSION}
    for g in SCHEMA:
        p[g["id"]] = {c["id"]: copy.deepcopy(c["default"]) for c in g["fields"]}
    return p


def _unit_preset(key, extra):
    """Ready-made project from a library unit and the extra changes."""
    q = copy.deepcopy(_lib.UNITS[key]["changes"])
    q.setdefault("geometry", {})["unit_library"] = key
    for g, d in extra.items():
        q.setdefault(g, {}).update(d)
    return q


PRESETS = {
    "BCEPQ — 30 °C, non-linear": {"test": {"kind": "ambient", "load": 98.7}},
    "BCEPQ — 600 °C hot, non-linear": {"test": {"kind": "hot", "T_target": 600, "load": 97.9}},
    "BCEPQ — 800 °C hot, non-linear": {"test": {"kind": "hot", "T_target": 800, "load": 98.9}},
    "BCEPQ — 1000 °C hot, non-linear": {"test": {"kind": "hot", "T_target": 1000, "load": 114.4}},
    "BCEPQ — 600 °C, thermal and stresses (linear, mesh of the paper)": {
        "model": {"analysis": "linear"}, "test": {"kind": "hot", "T_target": 600, "load": 97.9},
        "mesh": {"h": 8, "dz": 12}},
    "BCEPQ — residual after 600 °C, non-linear": {"test": {"kind": "residual", "T_target": 600}},
    "BCEPQ wall — ISO 834, 2D section (insulation)": {
        "model": {"analysis": "section_2d", "element": "T3"},
        "test": {"kind": "iso834", "T_initial": 20, "h_convection": 25, "iso_duration": 240}, "mesh": {"h": 5}},
    "Solid brick with 1:1:6 mortar — prism at 30 °C, non-linear": _unit_preset("solid_brick_10x6x22", {
        "model": {"analysis": "nonlinear"}, "test": {"kind": "ambient", "load": 50},
        "mesh": {"h": 8, "dz": 10, "delta_max": 1.2, "increment": 0.02}}),
    "Solid brick with 1:1:6 mortar — 600 °C, thermal and stresses (linear)": _unit_preset("solid_brick_10x6x22", {
        "model": {"analysis": "linear"}, "test": {"kind": "hot", "T_target": 600, "load": 50},
        "mesh": {"h": 8, "dz": 10}}),
    "Concrete block 14 × 19 × 29 (Camacho et al., 2015) — hollow prism at 20 °C, non-linear": _unit_preset(
        "concrete_m15_29", {"geometry": {"grouting": "none"}, "model": {"analysis": "nonlinear"},
                            "test": {"kind": "ambient", "T_initial": 20, "load": 180},
                            "mesh": {"h": 12, "dz": 20, "delta_max": 1.5, "increment": 0.03}}),
    "Concrete block 14 × 19 × 29 (Camacho et al., 2015) — grouted prism (G1) at 20 °C, non-linear": _unit_preset(
        "concrete_m15_29", {"geometry": {"grouting": "all"}, "model": {"analysis": "nonlinear"},
                            "test": {"kind": "ambient", "T_initial": 20, "load": 470},
                            "mesh": {"h": 12, "dz": 20, "delta_max": 2.0, "increment": 0.03}}),
    "Grouted concrete block 14 × 19 × 39 — ISO 834, 2D section (insulation)": _unit_preset(
        "concrete_m15_39", {"geometry": {"grouting": "alternate", "units_per_course": 1},
                            "ceramic": {**_lib.GROUTS["G1"]["changes"]["ceramic"], "grout_library": "G1"},
                            "model": {"analysis": "section_2d", "element": "T3"},
                            "test": {"kind": "iso834", "T_initial": 20, "h_convection": 25, "iso_duration": 240},
                            "mesh": {"h": 5}}),
    "Infill block with 8 holes — ISO 834, vertical 2D section (insulation)": _unit_preset(
        "infill_8_holes", {"model": {"analysis": "section_2d", "element": "T3"},
                            "test": {"kind": "iso834", "T_initial": 20, "h_convection": 25, "iso_duration": 240},
                            "mesh": {"h": 4}}),
    "Infill block with 8 holes — prism at 20 °C, non-linear": _unit_preset(
        "infill_8_holes", {"model": {"analysis": "nonlinear"},
                            "test": {"kind": "ambient", "T_initial": 20, "load": 30},
                            "mesh": {"h": 12, "dz": 20, "delta_max": 1.2, "increment": 0.03}}),
    "BCEPQ — transient test under 49 kN (50% of the 30 °C load), non-linear": {
        "test": {"kind": "transient", "load": 49, "T_target": 800, "rate": 5, "dwell": 30}},
    "Solid brick with 1:1:6 mortar — transient test under 46 kN (40%), non-linear": _unit_preset(
        "solid_brick_10x6x22", {"model": {"analysis": "nonlinear"},
                           "test": {"kind": "transient", "load": 46, "T_target": 800, "rate": 5, "dwell": 30,
                                      "T_initial": 20}}),
    "BCEPQ wall 2.94 × 3.00 m — ISO 834 under load": {
        "geometry": {"units_per_course": 10.5, "n_courses": 50},
        "model": {"analysis": "wall"},
        "test": {"kind": "iso834", "T_initial": 20, "h_convection": 25, "iso_duration": 240, "load": 150},
        "mesh": {"h": 5}},
    "BCEPQ wall 2.94 × 3.00 m — capacity at 20 °C": {
        "geometry": {"units_per_course": 10.5, "n_courses": 50},
        "model": {"analysis": "wall"},
        "test": {"kind": "ambient", "T_initial": 20, "load": 150}, "mesh": {"h": 5}},
    "BCEPQ wall 2.94 × 3.00 m — sensitivity under ISO 834 (Morris, 11 parameters)": {
        "geometry": {"units_per_course": 10.5, "n_courses": 50},
        "model": {"analysis": "sensitivity"},
        "test": {"kind": "iso834", "T_initial": 20, "h_convection": 25, "iso_duration": 300, "load": 150},
        "mesh": {"h": 5}},
}


def apply(p: Dict, changes: Dict) -> Dict:
    """Copy of project p with the changes applied (a dict {group: {field: value}})."""
    q = copy.deepcopy(p)
    for g, d in legacy.translate_project(changes).items():
        if isinstance(d, dict):
            q.setdefault(g, {}).update(d)
        else:
            q[g] = d
    return q


def complete(p: Dict) -> Dict:
    """Fills missing fields with the defaults (projects of earlier versions, including those of
    version 1.x, with Portuguese names)."""
    return apply(default_project(), {k: v for k, v in legacy.translate_project(p).items()})


# ----------------------------------------------------------------------------
# Conversion to the model objects
# ----------------------------------------------------------------------------
def _tab(t, fy=1.0):
    if not t:
        return None
    a = [(float(r[0]), float(r[1]) * fy) for r in t if r is not None and len(r) >= 2
         and r[0] not in ("", None) and r[1] not in ("", None)]
    a.sort()
    return a or None


def geometry(p) -> Geometry:
    g = p["geometry"]
    return Geometry(unit_length=g["unit_length"] / 1e3, width=g["width"] / 1e3,
                     course_height=g["course_height"] / 1e3, hole_diameter=g["hole_diameter"] / 1e3,
                     n_courses=int(g["n_courses"]), plate_thickness=g["plate_thickness"] / 1e3,
                     holes_per_unit=int(g["holes_per_unit"]), units_per_course=float(g["units_per_course"]),
                     bond=g["bond"], joint_type=g.get("joint_type", "dry"),
                     bedding=g.get("bedding", "total"),
                     filled_head_joints=bool(g.get("filled_head_joints", True)),
                     joint_thickness=g.get("joint_thickness", 2) / 1e3,
                     hole_shape=g.get("hole_shape", "circular"), face_shell=g.get("face_shell", 25) / 1e3,
                     web=g.get("web", 25) / 1e3, end_shell=g.get("end_shell", 25) / 1e3,
                     wall_cells=int(g.get("wall_cells", 0)), min_shell=g.get("min_shell", 8) / 1e3,
                     grouting=("none" if g.get("hole_shape") == "horizontal"
                                   else g.get("grouting", "none")),
                     holes_width=int(g.get("holes_width", 2)), holes_height=int(g.get("holes_height", 4)),
                     bed_shell=g.get("bed_shell", 8) / 1e3)


def net_area_ratio(geo: Geometry) -> float:
    """Unit area (without holes) / gross area, in course 0: converts E_unit_gross into the material E."""
    return 1.0 - geo.hole_area(0) / (geo.length * geo.width)


def solid_area_ratio(geo: Geometry) -> float:
    """Load-bearing area (unit + grout) / gross area, in course 0."""
    return 1.0 - geo.hole_area(0, grout=True) / (geo.length * geo.width)


def materials(p, geo: Geometry) -> Materials:
    c, mdl, ts = p["ceramic"], p["model"], p["test"]
    M = Materials()
    M.ceramic.k, M.ceramic.rho, M.ceramic.c = c["k"], c["rho"], c["c"]
    M.ceramic.E_unit_gross = c["E_unit_gross"] * 1e6
    M.ceramic.net_area_ratio = net_area_ratio(geo)
    M.ceramic.nu, M.ceramic.alpha = c["nu"], c["alpha"] * 1e-6
    M.ceramic.T_ref = ts["T_initial"]
    M.ceramic.table_k = _tab(c.get("table_k"))
    M.ceramic.table_E_factor = _tab(c.get("table_E_factor"))
    eps = _tab(c.get("table_thermal_strain"), 1e-3)
    M.ceramic.table_thermal_strain = eps
    M.ceramic.table_c = _tab(c.get("table_c"))
    M.ceramic.table_rho = _tab(c.get("table_rho"))
    a = M.mortar
    a.k, a.rho, a.c = c["mor_k"], c["mor_rho"], c["mor_c"]
    a.E_unit_gross, a.net_area_ratio = c["mor_E"] * 1e6, 1.0
    a.nu, a.alpha, a.T_ref = c["mor_nu"], c["mor_alpha"] * 1e-6, ts["T_initial"]
    grt = M.grout
    grt.k, grt.rho, grt.c = c.get("grt_k", 1.75), c.get("grt_rho", 2400), c.get("grt_c", 1000)
    grt.E_unit_gross, grt.net_area_ratio = c.get("grt_E", 26100) * 1e6, 1.0
    grt.nu, grt.alpha, grt.T_ref = c.get("grt_nu", 0.2), c.get("grt_alpha", 10) * 1e-6, ts["T_initial"]
    if c.get("grt_en1992", True):
        from .parameters import EN1992_E_FACTOR, EN1992_THERMAL_STRAIN
        grt.table_E_factor, grt.table_thermal_strain = list(EN1992_E_FACTOR), list(EN1992_THERMAL_STRAIN)
    t = p["geometry"]["joint_thickness"] / 1e3
    for j in (M.joint, M.capping):
        j.thickness, j.E, j.nu = t, mdl["E_joint"] * 1e6, c["nu"]
    M.air.wall_emissivity = ts["hole_emissivity"]
    M.hole_diameter = max(geo.hole_dimension, 1e-6)
    return M


def damage_mor(p) -> CeramicDamage:
    """Mortar damage: given fc and ft; ductility indices of Dauda et al. (2020)."""
    from .parameters import mortar_damage
    c, mdl = p["ceramic"], p["model"]
    d = mortar_damage(fc=c["mor_fc"] * 1e6, ft=c["mor_ft"] * 1e6)
    d.crack_closure = mdl["ceramic_model"] == "mazars_damage"
    return d


def damage_grt(p) -> CeramicDamage:
    """Grout damage: given fc and ft; same ductility indices as the mortar."""
    from .parameters import grout_damage
    c, mdl = p["ceramic"], p["model"]
    d = grout_damage(fc=c.get("grt_fc", 17.7) * 1e6, ft=c.get("grt_ft", 2.04) * 1e6)
    d.crack_closure = mdl["ceramic_model"] == "mazars_damage"
    if c.get("grt_en1992", True):
        from .parameters import EN1992_FC_FACTOR, EN1992_FT_FACTOR
        d.table_fc_factor, d.table_ft_factor = list(EN1992_FC_FACTOR), list(EN1992_FT_FACTOR)
    return d


def damage(p) -> CeramicDamage:
    mdl = p["model"]
    return CeramicDamage(ft=mdl["ft"] * 1e6, fc=mdl["fc"] * 1e6, du_t=mdl["du_t"] * 1e-3, du_c=mdl["du_c"] * 1e-3,
                        crack_closure=(mdl["ceramic_model"] == "mazars_damage"),
                        table_ft_factor=_tab(mdl.get("table_ft_factor")), table_fc_factor=_tab(mdl.get("table_fc_factor")))


def friction_joint(p) -> FrictionJoint:
    mdl, ts = p["model"], p["test"]
    return FrictionJoint(mu=mdl["mu"], gap_std=mdl["gap_std"] * 1e-3, initial_contact=mdl["initial_contact"] / 100,
                       corr_length=mdl["corr_length"] * 1e-3, E_contact=mdl["E_contact"] * 1e6,
                       seed=int(mdl["seed"]), pre_load_N=ts["seating_load"] * 1e3)


def validate(p) -> List[str]:
    errors = []
    try:
        geo = geometry(p)
        errors += geo.validate()
    except Exception as e:  # noqa
        return [f"Invalid geometry: {e}"]
    mdl, ts, msh_ = p["model"], p["test"], p["mesh"]
    if mdl["analysis"] == "section_2d":
        if ts["kind"] == "ambient":
            errors.append("The 2D section computes thermal effects only; choose a test with heating.")
    elif mdl["analysis"] in WALLS:
        if mdl["analysis"] == "sensitivity":
            errors += _validate_sensitivity(p)
        elif ts["kind"] not in ("ambient", "iso834"):
            errors.append("For the wall, choose compression at room temperature (capacity) or ISO 834 under load.")
        if ts["kind"] == "iso834" and ts["load"] <= 0:
            errors.append("For the wall under ISO 834, give the constant load (total over the wall length).")
        t2 = p["geometry"]["width"] / 2
        if max(abs(ts["top_eccentricity"]), abs(ts["base_eccentricity"])) >= t2:
            errors.append("The eccentricity must be less than half the wall thickness.")
        if msh_["h"] > 10:
            errors.append("For the wall, use elements of at most 10 mm in the section (recommended: 5 mm).")
        if msh_["n_elements"] < 4:
            errors.append("Use at least 4 beam elements over the height.")
    else:
        if mdl["element"] == "T3":
            errors.append("The T3 triangle is for the 2D section only; in 3D use W6, H8 or W15.")
    if ts["kind"] in ("hot", "residual", "transient") and ts["T_target"] <= ts["T_initial"]:
        errors.append("The target temperature must be higher than the initial one.")
    if ts["kind"] == "transient":
        if mdl["analysis"] not in ("nonlinear", "thermal"):
            errors.append("The transient test belongs to the 3D non-linear model (or to the thermal one only).")
        elif mdl["analysis"] == "nonlinear" and ts["load"] <= 0:
            errors.append("For the transient test, give the constant load.")
    tab = mdl.get("table_creep_unit") or []
    if any(len(r) != 2 or r[1] < 0 for r in tab):
        errors.append("The φ(T) table of the unit needs (T, φ) pairs with φ ≥ 0.")
    if mdl["analysis"] == "linear" and ts["kind"] != "iso834" and ts["load"] <= 0:
        errors.append("The linear analysis needs a positive load.")
    if mdl["analysis"] == "nonlinear" and ts["kind"] == "iso834" and ts["load"] <= 0:
        errors.append("For ISO 834 with the non-linear analysis, give the constant load.")
    g = p["geometry"]
    if g.get("joint_type") == "mortar" and g["joint_thickness"] > 0.3 * g["course_height"]:
        errors.append("The mortar joint must be at most 30% of the course height (centre to centre).")
    if g.get("unit_library") in _lib.UNITS and _lib.UNITS[g["unit_library"]].get("phase", 0) > _lib.CURRENT_PHASE:
        errors.append("The unit chosen in the library does not run in this version yet; choose another or 'Custom'.")
    if mdl["analysis"] not in WALLS and msh_["h"] > p["geometry"]["width"] / 3:
        errors.append("Element size too large for the unit width.")
    if mdl["analysis"] == "nonlinear" and msh_["increment"] >= msh_["delta_max"]:
        errors.append("The increment must be smaller than the maximum displacement.")
    return errors


# ----------------------------------------------------------------------------
# Thermal programme
# ----------------------------------------------------------------------------
def thermal_programme(p):
    """Returns (f(t), t_final in s, milestones {name: t}) of the furnace gas."""
    ts = p["test"]
    Ti = ts["T_initial"]
    if ts["kind"] in ("hot", "transient"):
        f, tf, tr = curves.ramp_dwell(Ti, ts["rate"], ts["T_target"], ts["dwell"] * 60)
        return f, tf, {"end of ramp": tr, "load": tf}
    if ts["kind"] == "residual":
        f, tf, (t1, t2, t3) = curves.ramp_dwell_cooling(Ti, ts["rate"], ts["T_target"], ts["dwell"] * 60,
                                                               ts["cooling_rate"], ts["final_dwell"] * 60)
        return f, tf, {"end of ramp": t1, "end of dwell": t2, "end of cooling": t3, "load": tf}
    if ts["kind"] == "iso834":
        return curves.iso834(Ti), ts["iso_duration"] * 60, {"end": ts["iso_duration"] * 60}
    return curves.constant(Ti), 0.0, {"load": 0.0}


def _boundaries_3d(p, m):
    ts = p["test"]
    f, tf, _ = thermal_programme(p)
    if ts["kind"] == "iso834":
        return [Boundary(th.side_areas_3d(m, faces=("y0",)), f, ts["h_convection"], ts["emissivity"], "exposed"),
                Boundary(th.side_areas_3d(m, faces=("y1",)), curves.constant(ts["T_initial"]), 4.0,
                         ts["emissivity"], "unexposed")]
    return [Boundary(th.side_areas_3d(m), f, ts["h_convection"], ts["emissivity"], "sides")]


def _probes_3d(m, geo):
    H = geo.prism_height
    zm = m.zs[np.argmin(np.abs(m.zs - H / 2))]
    W, Wm = m.m2d.width, m.m2d.model_width
    cf = geo.hole_centres
    x_web = 0.5 * (cf[0] + cf[1]) if len(cf) >= 2 else (cf[0] if cf else geo.length / 2)
    if geo.rectangular:
        cell = sorted(c for c in geo.cells(0, m.m2d.band_width) if c[4])
        if len(cell) >= 2:
            x_web = 0.5 * (cell[0][1] + cell[1][0])
    y_mid = min(W / 2, Wm)
    if len(cf) == 1:            # one hole: core of the shell between the hole and the face
        y_mid = 0.5 * (W / 2 - geo.hole_diameter / 2)
    s = {"surface (y = 0, mid-height)": msh.nearest_node(m.xyz, (geo.length / 3, 0.0, zm)),
         "web core (mid-height)": msh.nearest_node(m.xyz, (x_web, y_mid, zm))}
    if not m.m2d.symmetry:
        s["opposite face (y = W, mid-height)"] = msh.nearest_node(m.xyz, (geo.length / 3, W, zm))
    cer = np.unique(m.conn3[(m.mat == CERAMIC) & (np.abs(m.xyz[m.conn3[:, 0], 2] - H / 2) < H / 6)])
    s["mean of the ceramic in the middle third"] = cer
    return s


# ----------------------------------------------------------------------------
# Execution
# ----------------------------------------------------------------------------
class Cancelled(Exception):
    pass


class _Output:
    """Results accumulated during a run."""

    def __init__(self, folder):
        self.folder = folder
        self.summary, self.tables, self.figures, self.files, self.notices = [], [], [], [], []
        self.visual = None

    def pos(self, m, states, description="", material=None):
        visual.export(os.path.join(self.folder, "visual.json"), m, states, description, material)
        self.visual = "visual.json"

    def item(self, label, value, unit_str=""):
        if isinstance(value, (float, np.floating)):
            value = float(f"{value:.4g}")
        self.summary.append([label, value, unit_str])

    def table(self, title, columns, lines, file=None):
        self.tables.append({"title": title, "columns": columns, "lines": lines})
        if file:
            save_csv(os.path.join(self.folder, file), lines, columns)
            self.files.append({"name": file, "description": title})

    def figure_(self, fig, file, title):
        fig.savefig(os.path.join(self.folder, file), dpi=130, facecolor=fig.get_facecolor())
        self.figures.append({"file": file, "title": title})
        import matplotlib.pyplot as plt
        plt.close(fig)

    def vtk(self, file, m, **kw):
        save_vtk(os.path.join(self.folder, file), m, **kw)
        self.files.append({"name": file, "description": "3D fields for ParaView"})

    def dict(self):
        return {"summary": self.summary, "tables": self.tables, "figures": self.figures,
                "files": self.files, "notices": self.notices, "visual": self.visual}


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e4e3df", "axes.edgecolor": "#52514e"})
    return plt


COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6"]


ELEMENTS_3D = {"W6": "W6 wedges", "H8": "H8 hexahedra", "W15": "W15 wedges"}


def mesh_3d(p, geo=None):
    """3D mesh of the prism with the chosen element."""
    geo = geometry(p) if geo is None else geo
    msh_ = p["mesh"]
    el = p["model"]["element"] if p["model"]["element"] in ELEMENTS_3D else "W6"
    return msh.prism_mesh(geo, h=msh_["h"] / 1e3, dz=msh_["dz"] / 1e3, t_joint=p["geometry"]["joint_thickness"] / 1e3,
                           symmetry=(p["test"]["kind"] != "iso834"), capping=p["geometry"]["capping"],
                           element=el)


def mesh_figure(p, folder=None, file="mesh.png"):
    """Figure of the section and of the vertical cut of the mesh; returns (path or bytes, summary)."""
    _plt()
    from matplotlib.figure import Figure          # no pyplot: may run alongside an analysis
    from matplotlib.collections import PolyCollection
    geo = geometry(p)
    msh_, ts, mdl = p["mesh"], p["test"], p["model"]
    if mdl["analysis"] in WALLS:
        return _wall_figure(p, geo, folder, file)
    if mdl["analysis"] == "section_2d":
        m2 = msh.section_mesh(geo, h=msh_["h"] / 1e3, symmetry=False, bands=())
        info = msh.summary(m2)
    else:
        m = mesh_3d(p, geo)
        m2 = m.m2d
        info = msh.summary(m)
    if m2.vertical:                          # vertical section (horizontal holes), thickness × height
        fig = Figure(figsize=(10, 4.2))
        a0 = fig.subplots(1, 1)
        mat2 = th.section_material(m2, geo)
        colours2 = {AIR: "#cde2fb", CERAMIC: "#d9c7b0", MORTAR: "#7f8fa6", JOINT: "#eb6834"}
        pc = PolyCollection(m2.xy[m2.tri][:, :, ::-1] * 1e3, facecolors=[colours2[int(k)] for k in mat2],
                            edgecolors="#52514e", linewidths=0.25)
        a0.add_collection(pc); a0.autoscale(); a0.set_aspect("equal"); a0.grid(False)
        a0.set_title("Vertical section of one course (valid over the whole length)\n"
                     "holes in blue, half joints in blue-grey", fontsize=9)
        a0.set_xlabel("y, thickness (mm; y = 0 is the exposed face)"); a0.set_ylabel("height in the course (mm)")
        fig.tight_layout()
        if folder:
            path = os.path.join(folder, file)
            fig.savefig(path, dpi=130)
            return path, info
        import io
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=110)
        return buf.getvalue(), info
    fig = Figure(figsize=(10, 3.6))
    ax = fig.subplots(1, 1 if mdl["analysis"] == "section_2d" else 2,
                      gridspec_kw=None if mdl["analysis"] == "section_2d" else {"width_ratios": [2.2, 1]})
    ax = np.atleast_1d(ax)
    pol, mark = m2.tri, m2.mark
    quad = mdl["analysis"] != "section_2d" and m.kind == "H8"
    if quad:
        pol, mark = m2.quad, m2.mark_quad
    _, gr0 = msh.course_states(m2, geo, 0, quad=quad)
    colour = np.where(gr0, "#5b5f66", np.where(mark == 0, "#cde2fb", np.where(mark > 0, "#eb6834", "#d9c7b0")))
    pc = PolyCollection(m2.xy[pol] * 1e3, facecolors=colour, edgecolors="#52514e", linewidths=0.25)
    ax[0].add_collection(pc); ax[0].autoscale(); ax[0].set_aspect("equal")
    ax[0].set_title("Plan (head-joint bands in orange; the horizontal holes appear in the cut)"
                    if geo.horizontal else
                    "Section at course 1 (holes in blue, grout in grey, head-joint bands in orange)"
                    if gr0.any() else "Section at course 1 (holes in blue, head-joint bands in orange)",
                    fontsize=9, loc="left")
    ax[0].set_xlabel("x (mm)"); ax[0].set_ylabel("y (mm)"); ax[0].grid(False)
    if mdl["analysis"] != "section_2d":
        colours = {AIR: "#cde2fb", CERAMIC: "#d9c7b0", JOINT: "#eb6834", CAPPING: "#eda100", STEEL: "#8a8984",
                 MORTAR: "#7f8fa6", GROUT: "#5b5f66"}
        # vertical cut: one column of elements at y ~ 0
        quads = []
        cs = []
        mid = geo.rectangular or geo.grouting != "none"     # cut through the holes
        # horizontal holes: face x = 0 (thickness × height), where the rows of holes appear
        axis, coords = (0, [1, 2]) if geo.horizontal else (1, [0, 2])
        v_cut = 0.0 if geo.horizontal else (geo.width / 2 if mid else 0.0)
        for f in m.elem.faces_vis:
            if len(f) != 4:
                continue
            fc = m.conn3[:, f]
            sel = np.isclose(m.xyz[fc, axis], v_cut).all(1)
            if geo.horizontal:                        # first course only, so that the holes appear
                sel &= (m.xyz[fc, 2] <= geo.course_height + 1e-9).all(1)
            for e in np.nonzero(sel)[0]:
                quads.append(m.xyz[fc[e]][:, coords])
                cs.append(colours[int(m.mat[e])])
        pc2 = PolyCollection(np.array(quads) * 1e3, facecolors=cs, edgecolors="#52514e", linewidths=0.2)
        ax[1].add_collection(pc2); ax[1].autoscale(); ax[1].set_aspect("equal")
        ax[1].set_title("Face x = 0, course 1" if geo.horizontal else "Cut at y = W/2" if mid
                        else "Face y = 0 (joints and plates)", fontsize=9, loc="left")
        ax[1].set_xlabel("y (mm)" if geo.horizontal else "x (mm)"); ax[1].set_ylabel("z (mm)"); ax[1].grid(False)
    fig.tight_layout()
    if folder:
        path = os.path.join(folder, file)
        fig.savefig(path, dpi=130)
        return path, info
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    return buf.getvalue(), info


def run(p: Dict, folder: str, log: Callable[[str], None] = print) -> Dict:
    """Runs the project and writes everything to `folder`. Returns the results dictionary."""
    p = complete(p)
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "project.json"), "w", encoding="utf-8") as f:
        json.dump(p, f, ensure_ascii=False, indent=2)
    errors = validate(p)
    if errors:
        raise ValueError(" ".join(errors))
    t0 = time.time()
    S = _Output(folder)
    geo = geometry(p)
    M = materials(p, geo)
    an = p["model"]["analysis"]
    log(f"mef_masonry {__version__}: analysis '{an}', test '{p['test']['kind']}'")
    S.item("Net-to-gross area ratio", net_area_ratio(geo))
    if geo.grouting != "none":
        S.item("Load-bearing (unit + grout) to gross area ratio", solid_area_ratio(geo))
    if an not in WALLS:
        S.item("E of the material (net area)", M.ceramic.E_material / 1e6, "MPa")
    if an == "section_2d":
        _run_2d(p, geo, M, S, log)
    elif an == "wall":
        mesh_figure(p, folder)
        S.figures.append({"file": "mesh.png", "title": "Fibre section and wall layout"})
        _run_wall(p, geo, M, S, log)
    elif an == "sensitivity":
        mesh_figure(p, folder)
        S.figures.append({"file": "mesh.png", "title": "Fibre section and wall layout (nominal values)"})
        _run_sensitivity(p, geo, S, log)
    else:
        msh_ = p["mesh"]
        m = mesh_3d(p, geo)
        log("Mesh: " + msh.summary(m))
        S.item("Mesh nodes", m.n_nodes)
        S.item(f"Elements ({ELEMENTS_3D[m.kind]})", len(m.conn3))
        if m.kind == "W15" and m.n_nodes > (40000 if an == "nonlinear" else 80000):
            try:
                import pypardiso  # noqa: F401
                solver = ""
            except ImportError:
                solver = " and install pypardiso"
            S.notices.append(f"Large mesh for W15 ({m.n_nodes} nodes): the analysis may take hours and a lot of memory. "
                            f"Consider the 16 × 30 mm mesh{solver}.")
            log("Warning: " + S.notices[-1])
        mesh_figure(p, folder)
        S.figures.append({"file": "mesh.png", "title": "Mesh"})
        if an == "thermal":
            _run_thermal(p, geo, M, m, S, log)
        elif an == "linear":
            _run_linear(p, geo, M, m, S, log)
        else:
            _run_nonlinear(p, geo, M, m, S, log)
    S.item("Computation time", (time.time() - t0) / 60, "min")
    res = S.dict()
    res["project"] = p
    with open(os.path.join(folder, "results.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2, default=float)
    log("Done.")
    return res


def _thermal_3d(p, geo, M, m, S, log, instants=()):
    f, tf, milestones = thermal_programme(p)
    if tf <= 0:
        return None, {}, milestones
    probes = _probes_3d(m, geo)
    log(f"Thermal: {tf / 60:.0f} min of test, Δt = {p['mesh']['dt']:.0f} s")
    rt = th.solve_3d(m, M, _boundaries_3d(p, m), tf, dt=p["mesh"]["dt"], T_initial=p["test"]["T_initial"],
                        probes=probes, record_every=60.0, field_instants=sorted(set(instants)))
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 3.8))
    t = rt.t / 60
    ax.plot(t, rt.T_gas, color="#9a9893", ls=(0, (4, 2)), lw=1.5, label="furnace gas")
    for k, (name, v) in enumerate(rt.probes.items()):
        ax.plot(t, v, color=COLOURS[k % len(COLOURS)], lw=1.8, label=name)
    ax.set_xlabel("time (min)"); ax.set_ylabel("temperature (°C)"); ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    S.figure_(fig, "temperatures.png", "Temperatures during the test")
    names = list(rt.probes)
    S.table("Temperatures during the test", ["time (min)", "gas (°C)"] + [f"{n} (°C)" for n in names],
             [[round(a, 2), round(b, 1)] + [round(rt.probes[n][i], 1) for n in names]
              for i, (a, b) in enumerate(zip(t, rt.T_gas))][:: max(1, len(t) // 60)],
             file="temperatures.csv")
    for n in names:
        S.item(f"T final — {n}", rt.probes[n][-1], "°C")
    return rt, probes, milestones


def _run_thermal(p, geo, M, m, S, log):
    _, tf, _ = thermal_programme(p)
    dt = p["mesh"]["dt"]
    inst = [round(tf * k / 12 / dt) * dt for k in range(1, 12)] if tf > 0 else []
    rt, _, _ = _thermal_3d(p, geo, M, m, S, log, instants=inst)
    if rt is None:
        S.notices.append("Test at room temperature: there is no thermal field to compute.")
        return
    S.vtk("final_temperature.vtk", m, node_data={"T_C": rt.T_final}, elem_data={"material": m.mat.astype(float)})
    sty = [{"label": f"t = {t / 60:.0f} min", "fields": {"T": rt.fields[t]}} for t in sorted(rt.fields)]
    sty.append({"label": f"t = {tf / 60:.0f} min (end)", "fields": {"T": rt.T_final}})
    S.pos(m, sty, "3D thermal")


def _run_linear(p, geo, M, m, S, log):
    ts, mdl = p["test"], p["model"]
    if mdl["calibrate_joint"]:
        log(f"Calibrating the joint E for E_prism = {mdl['E_prism_target']:.0f} MPa")
        Ej, _ = studies.calibrate_joint(m, M, mdl["E_prism_target"], geo, verbose=False)
        S.item("Calibrated joint E", Ej / 1e6, "MPa")
    if ts["kind"] == "iso834":
        n = max(2, int(ts["iso_duration"] // 30))
        inst = [ts["iso_duration"] * 60 * (k + 1) / n for k in range(n)]
        rt, _, _ = _thermal_3d(p, geo, M, m, S, log, instants=inst)
        lines, iso_states = [], []
        for t in sorted(rt.fields):
            R = mech.solve_3d(m, M, T_nodes=rt.fields[t], load_N=ts["load"] * 1e3)
            Rt = mech.combine(R["thermal"], R["load"]) if ts["load"] > 0 else R["thermal"]
            s, p99, _ = mech.max_tension(m, Rt)
            d = mech.shortening_between_plates(m, Rt, geo.prism_height) * 1e3
            lines.append([round(t / 60, 1), round(s / 1e6, 3), round(p99 / 1e6, 3), round(d, 3)])
            iso_states.append({"label": f"t = {t / 60:.0f} min", "u": Rt.u,
                                "fields": {"T": rt.fields[t], **visual.u_fields(Rt.u), "s1": Rt.princ[:, 2] / 1e6,
                                           "s3": Rt.princ[:, 0] / 1e6, "szz": Rt.sig[:, 2] / 1e6}})
            log(f"  t = {t / 60:.0f} min: max tension {s / 1e6:.2f} MPa")
        S.table("Stresses and strain during ISO 834 (thermo-elastic)",
                 ["time (min)", "maximum tension (MPa)", "p99 tension (MPa)", "change in height (mm)"], lines,
                 file="iso834_linear.csv")
        S.vtk("final_field.vtk", m, node_data={"T_C": rt.T_final, "u_mm": Rt.u * 1e3},
              elem_data={"material": m.mat.astype(float), "s1_MPa": Rt.princ[:, 2] / 1e6})
        S.pos(m, iso_states, "Thermo-elastic, ISO 834")
        return
    rt, probes, milestones = _thermal_3d(p, geo, M, m, S, log)
    T = rt.T_final if rt is not None else None
    log("Thermo-elastic: thermal and load cases")
    R = mech.solve_3d(m, M, T_nodes=T, load_N=ts["load"] * 1e3, include_thermal=T is not None)
    Rc = R["load"]
    Rt = R.get("thermal")
    Rtot = mech.combine(Rt, Rc) if Rt is not None else Rc
    H = geo.prism_height
    A = geo.length * geo.width
    eps, _, _ = mech.extensometer_strain(m, Rc, H / 3, 2 * H / 3)
    S.item("Specimen modulus, middle third", ts["load"] * 1e3 / A / -eps / 1e6, "MPa")
    S.item("Shortening under load between plates", -mech.shortening_between_plates(m, Rc, H) * 1e3, "mm")
    lines = []
    for name, Rr in (("thermal", Rt), ("load", Rc), ("total", Rtot)):
        if Rr is None:
            continue
        s, p99, pos = mech.max_tension(m, Rr)
        x, y, z = (float(v) for v in np.round(pos * 1e3, 1))
        lines.append([name, round(s / 1e6, 3), round(p99 / 1e6, 3), x, y, z])
    S.table("Maximum principal tension in the ceramic",
             ["case", "maximum (MPa)", "99th percentile (MPa)", "x (mm)", "y (mm)", "z (mm)"],
             lines, file="max_tension.csv")
    if Rt is not None:
        S.item("Thermal expansion between plates", mech.shortening_between_plates(m, Rt, H) * 1e3, "mm")
    dn = {"u_total_mm": Rtot.u * 1e3}
    if T is not None:
        dn["T_C"] = T
    S.vtk("stresses.vtk", m, node_data=dn, elem_data={"material": m.mat.astype(float),
                                                      "s1_total_MPa": Rtot.princ[:, 2] / 1e6,
                                                      "s3_total_MPa": Rtot.princ[:, 0] / 1e6})
    sty = []
    for name, Rr in (("thermal only", Rt), ("load only", Rc), ("thermal + load", Rtot)):
        if Rr is None:
            continue
        c = {**visual.u_fields(Rr.u), "s1": Rr.princ[:, 2] / 1e6, "s3": Rr.princ[:, 0] / 1e6, "szz": Rr.sig[:, 2] / 1e6}
        if T is not None:
            c["T"] = T
        sty.append({"label": name, "fields": c, "u": Rr.u})
    S.pos(m, sty, "3D thermo-elastic")


def _nl_fields(p, geo, M, m, S, log):
    """Temperature fields for heating the non-linear model."""
    ts, msh_ = p["test"], p["mesh"]
    f, tf, milestones = thermal_programme(p)
    if tf <= 0:
        return None
    if ts["kind"] in ("hot", "transient"):
        tr = milestones["end of ramp"]
        inst = list(np.linspace(0, tr, int(msh_["n_fields_ramp"]) + 1)[1:]) + \
            list(np.linspace(tr, tf, int(msh_["n_fields_dwell"]) + 1)[1:])
    else:  # residual
        t1, t2, t3 = milestones["end of ramp"], milestones["end of dwell"], milestones["end of cooling"]
        nr, npat = int(msh_["n_fields_ramp"]), int(msh_["n_fields_dwell"])
        inst = list(np.linspace(0, t1, nr + 1)[1:]) + list(np.linspace(t1, t2, npat + 1)[1:]) + \
            list(np.linspace(t2, t3, nr + 1)[1:]) + list(np.linspace(t3, tf, npat + 1)[1:])
    dt = msh_["dt"]
    inst = sorted({round(t / dt) * dt for t in inst})
    rt, _, _ = _thermal_3d(p, geo, M, m, S, log, instants=inst)
    return rt.fields


def _creep_coef(p):
    """(enabled, k_unit, k_mortar, k_grout, table_unit [(T, φ)]) of the transient creep."""
    mdl, g = p["model"], p["geometry"]
    tab = [(float(T), float(v) * 1e-3) for T, v in (mdl.get("table_creep_unit") or [])]
    kb = float(mdl.get("k_tr_unit", 0.0))
    km = float(mdl.get("k_tr_mortar", 2.0)) if g.get("joint_type") == "mortar" else 0.0
    kg = float(mdl.get("k_tr_grout", 0.0)) if g.get("grouting", "none") != "none" else 0.0
    enabled = bool(mdl.get("creep", True)) and (kb > 0 or km > 0 or kg > 0 or bool(tab))
    return enabled, kb, km, kg, tab


def _enable_creep(mod, p, S=None):
    enabled, kb, km, kg, tab = _creep_coef(p)
    if enabled:
        mod.set_creep(kb, km, kg, tab or None)
        if S is not None:
            S.item("Transient creep: k_tr of the unit, the mortar and the grout",
                   f"{'table' if tab else f'{kb:g}'}; {km:g}; {kg:g}")
    elif S is not None:
        S.item("Transient creep", "off" if not p["model"].get("creep", True) else "zero (k_tr = 0)")
    return enabled


def _run_nonlinear(p, geo, M, m, S, log):
    from .nonlinear import NonlinearModel
    ts, msh_ = p["test"], p["mesh"]
    Pd, Pj = damage(p), friction_joint(p)
    H = geo.prism_height
    if ts["kind"] == "iso834":
        return _nonlinear_iso(p, geo, M, m, S, log, Pd, Pj)
    if ts["kind"] == "transient":
        return _nonlinear_transient(p, geo, M, m, S, log, Pd, Pj)
    fields = _nl_fields(p, geo, M, m, S, log)
    n_steps = max(2, int(round(msh_["delta_max"] / msh_["increment"])))
    mod = NonlinearModel(m, M, Pd, Pj, verbose=True, mortar_damage=damage_mor(p), grout_damage=damage_grt(p))
    mod.tol_u, mod.tol_d = msh_["tol_u"], msh_["tol_d"]
    mod.instantaneous = []
    mod.set_extensometer(H / 3, 2 * H / 3)
    _enable_creep(mod, p, S)
    if fields is None:
        fields = {0.0: np.full(m.n_nodes, M.ceramic.T_ref)}
    log(f"Heating: {len(fields)} field(s) under {ts['seating_load']:.1f} kN")
    mod.heat(fields)
    vc = mod.vol[mod.ic]
    frac = float(vc[mod.d_heating > 0.5].sum() / vc.sum())
    log(f"Loading to {msh_['delta_max']} mm in increments of {msh_['increment']} mm")
    d, P = mod.apply_load(msh_["delta_max"] / 1e3, n_steps)
    i = int(np.argmax(P))
    hist = [q for q in mod.history if q.phase == "load"]
    S.item("Peak load", P[i] / 1e3, "kN")
    S.item("Compressive strength (gross area)", P[i] / (geo.length * geo.width) / 1e6, "MPa")
    S.item("Platen displacement at peak", d[i] * 1e3, "mm")
    S.item("Ceramic with damage > 0.5 after heating", 100 * frac, "%")
    if ts["load"] > 0:
        S.item("Chord modulus 5–30% of the reference load", studies.chord_modulus(mod, ts["load"] * 1e3, geo), "MPa")
    nc = sum(1 for q in hist if not q.converged)
    S.item("Increments accepted without full convergence", nc)
    if nc:
        S.notices.append(f"{nc} increment(s) accepted at the minimum step without full convergence (see log).")
    if P[-1] > 0.95 * P[i] and i == len(P) - 1:
        S.notices.append("The load was still rising at the end: increase the maximum displacement to capture the peak.")
    S.table("Load × displacement curve",
             ["displacement (mm)", "load (kN)", "strain in the middle third (‰)", "contact in the joints (%)",
              "ceramic with damage > 0.5 (%)"],
             [[round(q.label * 1e3, 4), round(q.load_N / 1e3, 3), round(-q.eps_ext * 1e3, 4),
               round(100 * q.contact, 1), round(100 * q.damaged_vol, 2)] for q in hist],
             file="load_displacement.csv")
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(d * 1e3, P / 1e3, color=COLOURS[0], lw=2)
    ax.plot(d[i] * 1e3, P[i] / 1e3, "o", color=COLOURS[0], ms=7, mec="white", mew=1.5)
    ax.annotate(f"{P[i] / 1e3:.1f} kN", (d[i] * 1e3, P[i] / 1e3), xytext=(6, 6), textcoords="offset points", fontsize=8)
    ax.set_xlabel("platen displacement (mm)"); ax.set_ylabel("load (kN)"); ax.set_xlim(0, None); ax.set_ylim(0, None)
    fig.tight_layout()
    S.figure_(fig, "load_displacement.png", "Load × displacement curve")
    gap = np.full(len(m.conn3), np.nan); gap[mod.id_glob[mod.ij]] = mod.gap * 1e3
    d_heat = np.full(len(m.conn3), np.nan); d_heat[mod.id_glob[mod.ic]] = mod.d_heating
    S.vtk("damage.vtk", m, node_data={"T_C": mod.T, "u_mm": mod.u.reshape(-1, 3) * 1e3},
          elem_data={"material": m.mat.astype(float), "effective_damage": mod.global_damage("d"),
                      "tensile_damage": mod.global_damage("dt"), "compressive_damage": mod.global_damage("dc"),
                      "damage_after_heating": d_heat, "initial_gap_mm": gap})
    sty = visual.nonlinear_states(mod)
    if sty:
        sty[0]["fields"]["gap"] = gap
    S.pos(m, sty, "3D non-linear")


def _nonlinear_iso(p, geo, M, m, S, log, Pd, Pj):
    """ISO 834 on one face under constant load; failure by the axial
    deformation limits (EN 1363-1: C = h/100 and dC/dt = 3h/1000 per minute,
    h in mm) or by loss of equilibrium of the model."""
    from .nonlinear import NonlinearModel
    ts, msh_ = p["test"], p["mesh"]
    H = geo.prism_height
    dtc = msh_["iso_interval"] * 60
    inst = list(np.arange(dtc, ts["iso_duration"] * 60 + 1e-6, dtc))
    rt, _, _ = _thermal_3d(p, geo, M, m, S, log, instants=inst)
    mod = NonlinearModel(m, M, Pd, Pj, verbose=True, mortar_damage=damage_mor(p), grout_damage=damage_grt(p))
    mod.tol_u, mod.tol_d = msh_["tol_u"], msh_["tol_d"]
    mod.instantaneous = []
    _enable_creep(mod, p, S)
    P = ts["load"] * 1e3
    log(f"Applying the constant load of {ts['load']:.1f} kN")
    mod.heat({0.0: np.full(m.n_nodes, M.ceramic.T_ref)}, pre_load_N=P)
    u0 = float(mod.u[3 * mod.top + 2].mean())
    hmm = H * 1e3
    C_lim, v_lim = hmm / 100, 3 * hmm / 1000
    lines, t_fail, prev = [[0.0, 0.0, 0.0, 0.0]], None, (0.0, 0.0)
    for t in sorted(rt.fields):
        mod.heat({t: rt.fields[t]}, pre_load_N=P)
        C = -(float(mod.u[3 * mod.top + 2].mean()) - u0) * 1e3          # contraction positive, mm
        v = (C - prev[1]) / ((t - prev[0]) / 60)
        frac = float(mod.vol[mod.ic][np.maximum(mod.dt, mod.dc) > 0.5].sum() / mod.vol[mod.ic].sum())
        lines.append([round(t / 60, 1), round(C, 4), round(v, 4), round(100 * frac, 2)])
        log(f"  t = {t / 60:.0f} min: contraction {C:.3f} mm, rate {v:.3f} mm/min")
        if t_fail is None and C >= C_lim and v >= v_lim:
            t_fail = t / 60
            break
        if not mod.history[-1].converged and abs(C) > 10 * C_lim:
            t_fail = t / 60
            S.notices.append("The model lost equilibrium (unbounded displacements).")
            break
        prev = (t, C)
    S.item("Axial contraction limit C = h/100", C_lim, "mm")
    S.item("Contraction rate limit 3h/1000", v_lim, "mm/min")
    S.item("Time to loss of load-bearing capacity", t_fail if t_fail else f"> {ts['iso_duration']:.0f}", "min")
    if mod.creep is not None:
        S.item("Mean transient strain at the end (vertical)", -mod.mean_transient_strain * 1e3, "‰")
    S.table("Axial contraction during ISO 834 (positive = shortening)",
             ["time (min)", "contraction (mm)", "rate (mm/min)", "ceramic with damage > 0.5 (%)"], lines,
             file="iso834_contraction.csv")
    plt = _plt()
    a = np.array(lines, float)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(a[:, 0], a[:, 1], color=COLOURS[0], lw=2, label="model contraction")
    ax.axhline(C_lim, color="#9a9893", ls=(0, (4, 2)), lw=1, label="limit h/100")
    ax.set_xlabel("exposure time (min)"); ax.set_ylabel("axial contraction (mm)"); ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    S.figure_(fig, "iso834_contraction.png", "Axial contraction under constant load")
    S.notices.append("Failure criterion by the deformation limits of EN 1363-1 for vertical members; "
                    "without tables of E(T), ft(T) and fc(T), the model tends not to fail.")
    S.vtk("damage.vtk", m, node_data={"T_C": mod.T, "u_mm": mod.u.reshape(-1, 3) * 1e3},
          elem_data={"material": m.mat.astype(float), "tensile_damage": mod.global_damage("dt"),
                      "compressive_damage": mod.global_damage("dc")})
    S.pos(m, visual.nonlinear_states(mod), "3D non-linear, ISO 834")


def _nonlinear_transient(p, geo, M, m, S, log, Pd, Pj):
    """Transient test: constant load applied cold and heating in the furnace
    (ramp and dwell) to failure. Failure is the loss of equilibrium of the model
    or the axial deformation limits of EN 1363-1 (C = h/100 and dC/dt = 3h/1000
    per minute). The strain × temperature curve is the classical output of the test:
    its difference from the free expansion is the load-induced strain."""
    from .nonlinear import NonlinearModel
    ts, msh_ = p["test"], p["mesh"]
    H = geo.prism_height
    A = geo.length * geo.width
    fields = _nl_fields(p, geo, M, m, S, log)
    f_gas, tf, _ = thermal_programme(p)
    mod = NonlinearModel(m, M, Pd, Pj, verbose=True, mortar_damage=damage_mor(p), grout_damage=damage_grt(p))
    mod.tol_u, mod.tol_d = msh_["tol_u"], msh_["tol_d"]
    mod.instantaneous = []
    mod.set_extensometer(H / 3, 2 * H / 3)
    _enable_creep(mod, p, S)
    P = ts["load"] * 1e3
    S.item("Constant load", ts["load"], "kN")
    S.item("Stress on the gross area", P / A / 1e6, "MPa")
    log(f"Transient test: {ts['load']:.1f} kN cold and heating at {ts['rate']:g} °C/min up to "
        f"{ts['T_target']:.0f} °C")
    mod.heat({0.0: np.full(m.n_nodes, M.ceramic.T_ref)}, pre_load_N=P)
    u0 = float(mod.u[3 * mod.top + 2].mean())
    e0 = mod.extensometer_strain()
    vc = mod.vol[mod.ic]
    bl = ~(mod.mor | mod.grt_)

    def T_unit():
        Tc = mod.T[mod.conn3[mod.ic]].mean(1)
        return float((Tc * vc)[bl].sum() / vc[bl].sum())

    hmm = H * 1e3
    C_lim, v_lim = hmm / 100, 3 * hmm / 1000
    Ti = ts["T_initial"]
    lines = [[0.0, round(Ti, 1), round(T_unit(), 1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    fail, prev = None, (0.0, 0.0)
    for t in sorted(fields):
        mod.heat({t: fields[t]}, pre_load_N=P)
        conv = mod.history[-1].converged
        Tb = T_unit()
        C = -(float(mod.u[3 * mod.top + 2].mean()) - u0) * 1e3            # contraction positive, mm
        v = (C - prev[1]) / ((t - prev[0]) / 60)
        eps = -(mod.extensometer_strain() - e0) * 1e3                   # shortening positive, ‰
        free = -float(M.ceramic.thermal_strain(np.array([Tb]))[0]) * 1e3
        tr = -mod.mean_transient_strain * 1e3
        frac = float(vc[np.maximum(mod.dt, mod.dc) > 0.5].sum() / vc.sum())
        lines.append([round(t / 60, 2), round(float(f_gas(t)), 1), round(Tb, 1), round(eps, 4), round(free, 4),
                       round(eps - free, 4), round(tr, 4), round(C, 4), round(100 * frac, 2)])
        log(f"  t = {t / 60:.0f} min: unit {Tb:.0f} °C, strain {eps:+.3f} ‰ (free {free:+.3f} ‰), "
            f"contraction {C:.3f} mm")
        if (C >= C_lim and v >= v_lim) or (not conv and abs(C) > C_lim):
            fail = (t, Tb, float(f_gas(t)))
            break
        prev = (t, C)
    if fail:
        S.item("Time to failure", fail[0] / 60, "min")
        S.item("Furnace temperature at failure", fail[2], "°C")
        S.item("Mean unit temperature at failure", fail[1], "°C")
    else:
        S.item("Time to failure", f"> {tf / 60:.0f} (no failure)", "min")
    if mod.creep is not None:
        S.item("Mean transient strain at the end (vertical)", lines[-1][6], "‰")
    S.table("Transient test (shortening strains positive, in the middle third)",
             ["time (min)", "furnace (°C)", "unit, mean (°C)", "total strain (‰)", "free (‰)",
              "load-induced (‰)", "transient creep (‰)", "top contraction (mm)",
              "unit with damage > 0.5 (%)"], lines, file="transient_test.csv")
    plt = _plt()
    a = np.array(lines, float)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(a[:, 2], a[:, 3], color=COLOURS[0], lw=2, label=f"total, under {ts['load']:g} kN")
    ax.plot(a[:, 2], a[:, 4], color="#9a9893", lw=1.5, ls=(0, (4, 2)), label="free expansion")
    if mod.creep is not None:
        ax.plot(a[:, 2], a[:, 6], color=COLOURS[1], lw=1.5, label="transient creep")
    ax.axhline(0, color="#d6d4cf", lw=0.8)
    ax.set_xlabel("mean unit temperature (°C)"); ax.set_ylabel("shortening (‰)")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    S.figure_(fig, "transient_test.png", "Strain × temperature in the transient test")
    S.notices.append("Strain in the middle third of the face y = 0; shortening positive and expansion negative.")
    S.vtk("damage.vtk", m, node_data={"T_C": mod.T, "u_mm": mod.u.reshape(-1, 3) * 1e3},
          elem_data={"material": m.mat.astype(float), "tensile_damage": mod.global_damage("dt"),
                      "compressive_damage": mod.global_damage("dc")})
    S.pos(m, visual.nonlinear_states(mod), "3D non-linear, transient test")


def _run_2d(p, geo, M, S, log):
    from .mesh import HOLE
    ts, msh_ = p["test"], p["mesh"]
    m2 = msh.section_mesh(geo, h=msh_["h"] / 1e3, symmetry=False, bands=())
    mat2 = th.section_material(m2, geo)
    log("2D mesh: " + msh.summary(m2))
    S.item("Mesh nodes", m2.n_nodes)
    S.item("Elements (T3 triangles)", len(m2.tri))
    f, tf, _ = thermal_programme(p)
    face = None
    if ts["kind"] == "iso834":
        exp = th.boundary_lengths_2d(m2, faces=("y0",))
        nexp = th.boundary_lengths_2d(m2, faces=("y1",))
        ctr = [Boundary(exp, f, ts["h_convection"], ts["emissivity"], "exposed"),
               Boundary(nexp, curves.constant(ts["T_initial"]), 4.0, ts["emissivity"], "unexposed")]
        face = np.nonzero(nexp > 0)[0]
        w = nexp[face] / nexp[face].sum()
        Ti = ts["T_initial"]
        probes = {"unexposed face (mean)": face}

        def stop(t, T):
            return (T[face] @ w - Ti >= 140.0) or (T[face].max() - Ti >= 180.0)
    else:
        ctr = [Boundary(th.boundary_lengths_2d(m2), f, ts["h_convection"], ts["emissivity"])]
        cnt = m2.xy[m2.tri].mean(1)
        sol = np.unique(m2.tri[mat2 == CERAMIC])
        probes = {"ceramic minimum": sol[np.argsort(np.linalg.norm(m2.xy[sol] - [geo.length / 2, geo.width / 2], axis=1))[:1]]}
        stop = None
    log(f"2D thermal: {tf / 60:.0f} min")
    r = th.solve_2d(m2, M, ctr, tf, dt=msh_["dt"], T_initial=ts["T_initial"], probes=probes, record_every=60.0, mat=mat2,
                       stop_if=stop)
    if ts["kind"] == "iso834":
        S.item("Insulation time (140/180 °C on the unexposed face)",
               r.stopped_at / 60 if r.stopped_at else f"> {tf / 60:.0f}", "min")
    s = mech.solve_section_gps(m2, M, r.T_final, mat=mat2)
    sol = mat2 != AIR
    S.item("Maximum tension (generalised plane strain)", np.nanmax(s["s1"]) / 1e6, "MPa")
    S.item("p99 tension", np.nanpercentile(s["s1"][sol], 99) / 1e6, "MPa")
    from .output import section_figure
    Te = r.T_final[m2.tri].mean(1)
    section_figure(m2, Te, "Temperature at the end (°C)", os.path.join(S.folder, "section_T.png"), "°C")
    S.figures.append({"file": "section_T.png", "title": "Temperature in the section at the end"})
    section_figure(m2, s["s1"] / 1e6, "Maximum principal stress (MPa)", os.path.join(S.folder, "section_s1.png"), "MPa",
                 mask=sol, cmap="RdBu_r")
    S.figures.append({"file": "section_s1.png", "title": "Maximum principal stress in the section"})
    names = list(r.probes)
    S.table("Temperatures during the test", ["time (min)", "gas (°C)"] + [f"{n} (°C)" for n in names],
             [[round(a / 60, 2), round(b, 1)] + [round(float(np.mean(r.probes[n][i])), 1) for n in names]
              for i, (a, b) in enumerate(zip(r.t, r.T_gas))], file="temperatures_2d.csv")
    S.vtk("section.vtk", m2, node_data={"T_C": r.T_final}, elem_data={"s1_MPa": s["s1"] / 1e6})
    S.pos(m2, [{"label": f"end, t = {r.t[-1] / 60:.0f} min", "fields": {"T": r.T_final, "s1": s["s1"] / 1e6}}],
          "2D section", material=mat2)


# ----------------------------------------------------------------------------
# Wall (vertical strip with second order)
# ----------------------------------------------------------------------------
def _wall_section(p, geo):
    """2D mesh of one unit length (period of the section) and the fibres."""
    from .wall import FibreSection
    geo1 = replace(geo, units_per_course=1.0)
    m2 = msh.section_mesh(geo1, h=p["mesh"]["h"] / 1e3, symmetry=False, bands=())
    return geo1, m2, FibreSection.from_mesh(m2, geo1)


def _n_elements(p):
    n = max(4, int(p["mesh"]["n_elements"]))
    return n + (n % 2)


def _wall_figure(p, geo, folder=None, file="mesh.png"):
    _plt()
    from matplotlib.figure import Figure
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Rectangle, FancyArrowPatch
    ts = p["test"]
    geo1, m2, sec = _wall_section(p, geo)
    H, t, L = geo.prism_height, geo.width, geo.length
    n = _n_elements(p)
    info = (f"Wall {L * 1e3:.0f} × {H * 1e3:.0f} × {t * 1e3:.0f} mm (H/t = {H / t:.1f}); "
            + ("vertical section of one course" if m2.vertical else "section of one unit")
            + f" with {len(sec.y)} fibres; {n} beam elements")
    fig = Figure(figsize=(10, 3.8))
    ax = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1.3, 1] if m2.vertical else [2.4, 1]})
    if m2.vertical:
        # thickness horizontal (exposed face on the left), course height vertical
        fib = np.zeros(len(m2.tri), bool); fib[sec.tri] = True
        mat2 = th.section_material(m2, geo1)
        colour = np.where(mat2 == AIR, "#cde2fb", np.where(~fib, "#efe6da",
                       np.where(mat2 == CERAMIC, "#d9c7b0", "#7f8fa6")))
        pc = PolyCollection(m2.xy[m2.tri][:, :, ::-1] * 1e3, facecolors=colour, edgecolors="#52514e", linewidths=0.25)
        ax[0].add_collection(pc); ax[0].autoscale(); ax[0].set_aspect("equal"); ax[0].grid(False)
        ax[0].axvline(0, color="#eb6834", lw=3)
        ax[0].set_title("Vertical section of one course (fibres: continuous shells)", fontsize=9, loc="left")
        ax[0].set_xlabel("y, thickness (mm; y = 0 is the face " + ("exposed)" if ts["kind"] == "iso834" else "y = 0)"))
        ax[0].set_ylabel("height in the course (mm)")
    else:
        _, gr0 = msh.course_states(m2, geo1, 0)
        colour = np.where(gr0, "#5b5f66", np.where(m2.mark == 0, "#cde2fb", "#d9c7b0"))
        pc = PolyCollection(m2.xy[m2.tri] * 1e3, facecolors=colour, edgecolors="#52514e", linewidths=0.25)
        ax[0].add_collection(pc); ax[0].autoscale(); ax[0].set_aspect("equal"); ax[0].grid(False)
        ax[0].axhline(0, color="#eb6834", lw=3)
        ax[0].text(m2.length * 1e3 / 2, -6, "exposed face (ISO 834)" if ts["kind"] == "iso834" else "face y = 0",
                   color="#eb6834", ha="center", va="top", fontsize=8)
        ax[0].set_ylim(-18, t * 1e3 + 4)
        ax[0].set_title("Fibre section: one unit length (unit in beige, holes in blue"
                        + (", grout in grey)" if gr0.any() else ")"), fontsize=9, loc="left")
        ax[0].set_xlabel("x (mm)"); ax[0].set_ylabel("y (mm)")
    # schematic elevation: y horizontal (exposed on the left), height vertical
    a = ax[1]; a.grid(False)
    imp = H / ts["imperfection"] if ts["imperfection"] > 0 else 0.0
    xs = np.linspace(0, H, n + 1)
    scl = 0.25 * t / max(imp, 1e-9) if imp > 0 else 0
    ys = -imp * np.sin(np.pi * xs / H) * scl
    a.add_patch(Rectangle((-t / 2 * 1e3, 0), t * 1e3, H * 1e3, facecolor="#efe6da", edgecolor="#9a9893", lw=0.8))
    a.plot(ys * 1e3, xs * 1e3, "-o", color="#2a78d6", ms=2.5, lw=1.2)
    if ts["kind"] == "iso834":
        a.axvline(-t / 2 * 1e3, color="#eb6834", lw=3)
    for z, sup, ex in ((0.0, ts["base_support"], ts["base_eccentricity"]), (H, ts["top_support"], ts["top_eccentricity"])):
        yy = -ex
        if sup == "fixed":
            a.plot([-t / 2 * 1e3 * 1.6, t / 2 * 1e3 * 1.6], [z * 1e3] * 2, color="#52514e", lw=3)
        else:
            a.plot(yy, z * 1e3, marker="^" if z == 0 else "v", color="#52514e", ms=9)
    yt = -ts["top_eccentricity"]
    a.add_patch(FancyArrowPatch((yt, H * 1e3 * 1.16), (yt, H * 1e3 * 1.02), arrowstyle="-|>", mutation_scale=12,
                                color="#eb6834" if ts["kind"] == "iso834" else "#52514e", lw=1.5))
    a.text(yt, H * 1e3 * 1.18, "P", ha="center", va="bottom", fontsize=9)
    a.set_xlim(-t * 1e3 * 1.4, t * 1e3 * 1.4); a.set_ylim(-0.06 * H * 1e3, 1.26 * H * 1e3)
    a.set_title("Vertical strip (axis and imperfection magnified)", fontsize=9, loc="left")
    a.set_xlabel("y (mm)"); a.set_ylabel("height (mm)")
    fig.tight_layout()
    if folder:
        path = os.path.join(folder, file)
        fig.savefig(path, dpi=130)
        return path, info
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    return buf.getvalue(), info


def _prism_curve_30(p, log):
    """Runs the non-linear prism model (5 courses of one and a half units) at 30 °C with
    the current inputs; returns (delta mm, load kN), the prism H and the net area."""
    from .nonlinear import NonlinearModel
    q = complete(copy.deepcopy(p))
    q["geometry"].update({"units_per_course": 1.5, "n_courses": 5, "bond": "running_bond"})
    q["model"].update({"analysis": "nonlinear", "element": "W6"})
    q["test"].update({"kind": "ambient"})
    q["mesh"].update({"h": 12, "dz": 20, "delta_max": 0.8, "increment": 0.02})
    geo = geometry(q)
    M = materials(q, geo)
    m = mesh_3d(q, geo)
    log("Calibration from the prism at 30 °C: " + msh.summary(m))
    mod = NonlinearModel(m, M, damage(q), friction_joint(q), verbose=False, mortar_damage=damage_mor(q), grout_damage=damage_grt(q))
    mod.heat({0.0: np.full(m.n_nodes, M.ceramic.T_ref)})
    if geo.joint_type == "mortar":       # mortared prism: later peak
        d, P = mod.apply_load(0.006 * geo.prism_height, 60)
    else:
        d, P = mod.apply_load(0.8e-3, 40)
    i = int(np.argmax(P))
    log(f"  prism: peak {P[i] / 1e3:.1f} kN at {d[i] * 1e3:.3f} mm")
    A_net = geo.length * geo.width * solid_area_ratio(geo)
    return np.c_[d * 1e3, P / 1e3], geo.prism_height, A_net, P[i] / (geo.length * geo.width)


def _wall_law(p, geo, log, S):
    from .wall import WallLaw
    from .wall_data import PRISM_CURVES_BCEPQ, H_REF, A_GROSS_REF
    ref = WallLaw.from_prism(PRISM_CURVES_BCEPQ, H_REF, A_GROSS_REF * net_area_ratio(Geometry()))
    if p["model"]["wall_law"] != "prism":
        return ref
    curve, H_p, A_net, f_gross = _prism_curve_30(p, log)
    new_law = WallLaw.from_prism({30.0: curve}, H_p, A_net)
    k_s = new_law.s_peak[0] / ref.s_peak[0]
    k_e = new_law.e_peak[0] / ref.e_peak[0]
    S.item("Strength of the calibrated prism at 30 °C (gross area)", f_gross / 1e6, "MPa")
    log(f"  ratios to the BCEPQ: stress {k_s:.3f}, strain {k_e:.3f}")
    return ref.scaled(k_s, k_e).with_curve(30.0, new_law.e[0], new_law.s[0])


def _wall_creep(p, geo, M, sec, S=None):
    """c(T) per fibre (1/Pa) of the transient creep in the wall, or None.

    Unit and bed joint are in series over the course height: c = (1 - f) φ_b/fc_b
    + f φ_a/fc_a, with f = joint thickness / course height (mortar joints only).
    The grout is continuous: c = φ_g/fc_g."""
    enabled, kb, km, kg, tab = _creep_coef(p)
    if S is not None:
        S.item("Transient creep: k_tr of the unit, the mortar and the grout",
               f"{'table' if tab else f'{kb:g}'}; {km:g}; {kg:g}" if enabled else
               ("off" if not p["model"].get("creep", True) else "zero (k_tr = 0)"))
    if not enabled:
        return None
    c, mdl = p["ceramic"], p["model"]
    fc_b = mdl["fc"] * 1e6
    fc_a = c.get("mor_fc", 4.2) * 1e6
    fc_g = c.get("grt_fc", 17.7) * 1e6
    f_a = geo.joint_thickness / geo.course_height if geo.joint_type == "mortar" else 0.0
    grt = np.asarray(sec.grout, bool)

    def phi_unit(T):
        if tab:
            Ts, vs = zip(*sorted(tab))
            return np.interp(T, Ts, vs)
        return kb * np.maximum(M.ceramic.thermal_strain(T), 0.0)

    def coef(T):
        T = np.asarray(T, float)
        cb = (1 - f_a) * phi_unit(T) / fc_b + f_a * km * np.maximum(M.mortar.thermal_strain(T), 0.0) / fc_a
        cg = kg * np.maximum(M.grout.thermal_strain(T), 0.0) / fc_g
        return np.where(grt, cg, cb)
    return coef


def _run_wall(p, geo, M, S, log, law=None, fast=False):
    """Wall under ISO 834 or cold. `law`: ready-made wall law (without the project
    factors); `fast`: no figures and no load levels, returns R, I and Pmax."""
    from .wall import WallModel
    from .output import section_figure
    ts, msh_ = p["test"], p["mesh"]
    geo1, m2, sec = _wall_section(p, geo)
    H, t, L = geo.prism_height, geo.width, geo.length
    Ti = ts["T_initial"]
    n = _n_elements(p)
    log(f"Wall {L * 1e3:.0f} × {H * 1e3:.0f} × {t * 1e3:.0f} mm; section with {len(sec.y)} fibres; "
        f"{n} beam elements")
    S.item("Wall height", H * 1e3, "mm")
    S.item("Esbeltez H/t", H / t)
    S.item("Fibres in the section", len(sec.y))
    if law is None:
        law = _wall_law(p, geo, log, S)
    fr = p["model"].get("wall_strength_factor", 1.0)
    fk = p["model"].get("wall_stiffness_factor", 1.0)
    if fr != 1.0 or fk != 1.0:
        law = law.scaled(fr, fr / fk)
        S.item("Factors on the strength and stiffness of the law", f"{fr:.3g}; {fk:.3g}")
    if not fast and p["model"]["wall_law"] != "prism" and (geo.joint_type != "dry"
                                                  or p["ceramic"].get("material_library", "bcepq") != "bcepq"):
        S.notices.append("The wall law used is that of the BCEPQ prism; for other masonry, choose 'Calibrate now "
                        "with the non-linear prism model' under Models and elements.")
    ratio = sec.area / t
    S.item(f"Masonry strength at {Ti:.0f} °C (gross area)", law.strength(Ti) * ratio / 1e6, "MPa")
    imp = H / ts["imperfection"] if ts["imperfection"] > 0 else 0.0
    e_t, e_b = ts["top_eccentricity"] / 1e3, ts["base_eccentricity"] / 1e3
    grt_fib, mor_fib = sec.grout, sec.mortar

    def eth(T):
        """Free thermal strain per fibre (unit, grout or joint)."""
        return np.where(grt_fib, M.grout.thermal_strain(T),
                        np.where(mor_fib, M.mortar.thermal_strain(T), M.ceramic.thermal_strain(T)))

    flu = _wall_creep(p, geo, M, sec, S)

    def new():
        return WallModel(H, sec, law, eth, n_el=n, base_support=ts["base_support"], top_support=ts["top_support"],
                            imperfection=imp, T_initial=Ti, creep=flu)

    if imp == 0 and e_t == 0 and e_b == 0:
        S.notices.append("Perfect wall and centred load: without imperfection, buckling is not detected when cold.")
    # ---------------------------------------------------------------- cold capacity
    log(f"Capacity at {Ti:.0f} °C (increasing load until equilibrium is lost)")
    m0 = new()
    curve = [(0.0, 0.0, 0.0)]
    Pmax = m0.capacity(e_t, e_b, record=lambda P: curve.append((P, m0.mid_deflection, m0.top_displacement)))
    N_sec = law.strength(Ti) * sec.area
    S.item(f"Wall capacity at {Ti:.0f} °C", Pmax * L / 1e3, "kN")
    S.item(f"Capacity per metre at {Ti:.0f} °C", Pmax / 1e3, "kN/m")
    S.item("Reduction factor (capacity / section strength)", Pmax / N_sec)
    log(f"  capacity {Pmax * L / 1e3:.1f} kN ({Pmax / 1e3:.1f} kN/m); reduction factor {Pmax / N_sec:.3f}")
    plt = _plt()

    def stress_field(mod):
        a = np.full(len(m2.tri), np.nan)
        sg = mod.mid_height_fibres()
        if sg is not None:
            a[sec.tri] = sg / 1e6
        return a

    if ts["kind"] != "iso834":
        c = np.array(curve)
        S.table("Load × deflection, cold", ["load (kN)", "load (kN/m)", "mid-height deflection (mm)", "top displacement (mm)"],
                 [[round(a * L / 1e3, 2), round(a / 1e3, 3), round(b * 1e3, 3), round(d * 1e3, 4)] for a, b, d in c],
                 file="wall_capacity.csv")
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.8))
        ax[0].plot(-c[:, 1] * 1e3, c[:, 0] * L / 1e3, color=COLOURS[0], lw=2)
        ax[0].set_xlabel("mid-height deflection, towards y < 0 (mm)"); ax[0].set_ylabel("load (kN)")
        ax[0].axhline(Pmax * L / 1e3, color="#9a9893", ls=(0, (4, 2)), lw=1)
        ax[0].annotate(f"{Pmax * L / 1e3:.0f} kN", (0, Pmax * L / 1e3), xytext=(4, 4), textcoords="offset points", fontsize=8)
        ax[1].plot(m0.deflection * 1e3, m0.x0 * 1e3, color=COLOURS[0], lw=2)
        ax[1].axvline(0, color="#9a9893", lw=0.8)
        ax[1].set_xlabel("horizontal displacement at peak (mm)"); ax[1].set_ylabel("height (mm)")
        fig.tight_layout()
        S.figure_(fig, "wall_capacity.png", "Load × deflection and deformed shape at peak")
        S.pos(m2, [{"label": "peak load (mid-height)", "fields": {"szz": stress_field(m0)}}],
              "Wall: vertical stress in the fibres at mid-height", material=th.section_material(m2, geo1))
        S.notices.append("Vertical strip model: the vertical edges of the wall do not restrain bending. "
                        "Material homogenised from the prism, without tension (the joints open).")
        return
    # ---------------------------------------------------------------- 2D thermal
    f, tf, _ = thermal_programme(p)
    exp = th.boundary_lengths_2d(m2, faces=("y0",))
    nexp = th.boundary_lengths_2d(m2, faces=("y1",))
    ctr = [Boundary(exp, f, ts["h_convection"], ts["emissivity"], "exposed"),
           Boundary(nexp, curves.constant(Ti), 4.0, ts["emissivity"], "unexposed")]
    face = np.nonzero(nexp > 0)[0]
    w = nexp[face] / nexp[face].sum()
    face_exp = np.nonzero(exp > 0)[0]
    we = exp[face_exp] / exp[face_exp].sum()
    t_iso = []

    def insulation(tt, T):
        if not t_iso and (T[face] @ w - Ti >= 140.0 or T[face].max() - Ti >= 180.0):
            t_iso.append(tt)
        return False

    inst = list(np.arange(60.0, tf + 1e-6, 60.0))
    log(f"2D thermal analysis of the section: {tf / 60:.0f} min of ISO 834, Δt = {msh_['dt']:.0f} s")
    r = th.solve_2d(m2, M, ctr, tf, dt=msh_["dt"], T_initial=Ti, record_every=60.0, field_instants=inst,
                       mat=th.section_material(m2, geo1),
                       stop_if=insulation)
    times = [0.0] + sorted(r.fields)
    Tnodes = {0.0: np.full(m2.n_nodes, Ti), **r.fields}
    fields = [(tt, sec.temperatures(m2, Tnodes[tt])) for tt in times]
    S.item("Insulation time (140/180 °C on the unexposed face)", t_iso[0] / 60 if t_iso else f"> {tf / 60:.0f}", "min")
    # ---------------------------------------------------------------- fire under load
    P = ts["load"] * 1e3 / L
    S.item("Applied load", ts["load"], "kN")
    S.item("Applied load per metre", P / 1e3, "kN/m")
    S.item(f"Load / capacity at {Ti:.0f} °C", 100 * P / Pmax, "%")
    snapshot_step = max(15.0, round(tf / 60 / 8 / 15) * 15) * 60
    snapshots = {}
    hist = []
    mod = new()

    def record(tt):
        hist.append((tt, mod.mid_deflection, mod.max_deflection, mod.top_displacement))
        k = round(tt / 60)
        if abs(tt - 60 * k) < 1e-6 and (tt % snapshot_step < 1e-6) and tt > 0:
            snapshots[tt] = (mod.deflection.copy(), stress_field(mod))
        if abs(tt - 60 * k) < 1e-6 and k % 10 == 0 and k > 0:
            log(f"  t = {k} min: deflection {mod.max_deflection * 1e3:.1f} mm, top {mod.top_displacement * 1e3:+.2f} mm")

    log(f"ISO 834 under {ts['load']:.1f} kN ({P / 1e3:.1f} kN/m, {100 * P / Pmax:.0f}% of the cold capacity)")
    if P >= Pmax or not mod.apply_load(P, e_t, e_b):
        t_col, t_last = 0.0, 0.0
        S.notices.append("The given load is larger than the cold capacity of the wall.")
    else:
        hist.append((0.0, mod.mid_deflection, mod.max_deflection, mod.top_displacement))
        snapshots[0.0] = (mod.deflection.copy(), stress_field(mod))
        t_col, t_last = mod.fire(P, e_t, e_b, fields, record=record)
        snapshots[t_last] = (mod.deflection.copy(), stress_field(mod))
    R = t_col / 60 if t_col is not None else None
    S.item("Fire resistance time (loss of equilibrium)", R if R is not None else f"> {tf / 60:.0f}", "min")
    # axial deformation criterion of EN 1363-1 (h in mm)
    hmm = H * 1e3
    C_lim, v_lim = hmm / 100, 3 * hmm / 1000
    t_C = None
    if hist:
        hm = [(a, d) for a, b, c, d in hist if abs(a / 60 - round(a / 60)) < 1e-6]
        u0 = hm[0][1]
        for (ta, ua), (tb, ub) in zip(hm[:-1], hm[1:]):
            C = -(ub - u0) * 1e3
            v = -(ub - ua) * 1e3 / ((tb - ta) / 60)
            if C >= C_lim and v >= v_lim:
                t_C = tb / 60
                break
    S.item("Axial contraction criterion (EN 1363-1)", t_C if t_C else "not reached", "min" if t_C else "")
    if fast:
        return {"R": R if R is not None else tf / 60, "R_censored": R is None,
                "I": t_iso[0] / 60 if t_iso else tf / 60, "I_censored": not t_iso,
                "Pmax": Pmax * L / 1e3, "deflection": hist[-1][2] * 1e3 if hist else float("nan"),
                "capacity_load": 100 * P / Pmax, "tf": tf / 60}
    if hist:
        S.item("Maximum deflection at the end", hist[-1][2] * 1e3, "mm")
        xc, rc = mod.critical_compression()
        S.item("Height of the most compressed fibre at the end", xc * 1e3, "mm")
    if R is not None:
        log(f"  loss of equilibrium at {R:.1f} min")
    a = np.array(hist) if hist else np.zeros((1, 4))
    S.table("History of ISO 834 under load",
             ["time (min)", "mid-height deflection (mm)", "maximum deflection (mm)", "top displacement (mm)"],
             [[round(x / 60, 3), round(b * 1e3, 3), round(c * 1e3, 3), round(d * 1e3, 4)] for x, b, c, d in a
              if abs(x / 60 - round(x / 60)) < 1e-6 or x == a[-1, 0]],
             file="wall_history.csv")
    # ---------------------------------------------------------------- load levels
    levels = []
    k = int(ts["load_levels"])
    if k > 0:
        log(f"Fire resistance × load curve: {k} levels")
        for frac in np.linspace(0.1, 0.9, k):
            mk = new()
            Pk = frac * Pmax
            if not mk.apply_load(Pk, e_t, e_b):
                levels.append((frac, Pk, 0.0)); continue
            tc, _ = mk.fire(Pk, e_t, e_b, fields)
            Rk = tc / 60 if tc is not None else np.inf
            levels.append((frac, Pk, Rk))
            log(f"  {100 * frac:.0f}% ({Pk * L / 1e3:.0f} kN): " + (f"{Rk:.1f} min" if np.isfinite(Rk) else f"> {tf / 60:.0f} min"))
        S.table("Fire resistance × load", ["load / cold capacity (%)", "load (kN)", "load (kN/m)",
                                                  "fire resistance (min)"],
                 [[round(100 * fr, 1), round(Pk * L / 1e3, 1), round(Pk / 1e3, 2),
                   round(Rk, 1) if np.isfinite(Rk) else f"> {tf / 60:.0f}"] for fr, Pk, Rk in levels],
                 file="wall_time_load.csv")
    # ---------------------------------------------------------------- figures
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    ax[0].plot(a[:, 0] / 60, -a[:, 2] * 1e3, color=COLOURS[0], lw=2)
    ax[0].set_xlabel("exposure time (min)"); ax[0].set_ylabel("maximum deflection towards the fire (mm)")
    ax[1].plot(a[:, 0] / 60, a[:, 3] * 1e3, color=COLOURS[1], lw=2)
    ax[1].axhline(0, color="#9a9893", lw=0.8)
    ax[1].set_xlabel("exposure time (min)"); ax[1].set_ylabel("vertical displacement of the top (mm)")
    for x in ax:
        if R is not None:
            x.axvline(R, color="#9a9893", ls=(0, (4, 2)), lw=1)
    if R is not None:
        ax[0].annotate(f"loss of equilibrium: {R:.0f} min", (R, 0), xytext=(-4, 6), textcoords="offset points",
                       ha="right", fontsize=8)
    fig.tight_layout()
    S.figure_(fig, "wall_deflection.png", "Deflection and top displacement during ISO 834")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.2))
    tf_ord = sorted(snapshots)
    cmap = plt.get_cmap("inferno")
    for i, tt in enumerate(tf_ord):
        colour = cmap(0.1 + 0.75 * i / max(1, len(tf_ord) - 1))
        ax[0].plot(snapshots[tt][0] * 1e3, mod.x0 * 1e3, color=colour, lw=1.6, label=f"{tt / 60:.0f} min")
    ax[0].axvline(0, color="#9a9893", lw=0.8)
    ax[0].set_xlabel("horizontal displacement (mm; < 0 = towards the fire)"); ax[0].set_ylabel("height (mm)")
    if ax[0].get_legend_handles_labels()[0]:
        ax[0].legend(fontsize=7, frameon=False)
    # mean temperature through the thickness
    yb = np.linspace(-t / 2, t / 2, 29)
    ib = np.clip(np.digitize(sec.y, yb) - 1, 0, len(yb) - 2)
    yc = 0.5 * (yb[1:] + yb[:-1])
    Ab = np.bincount(ib, sec.A, len(yc))
    for i, tt in enumerate([x for x in tf_ord if x > 0 and x in Tnodes] or [times[-1]]):
        Tf = sec.temperatures(m2, Tnodes[tt])
        Tm = np.bincount(ib, sec.A * Tf, len(yc)) / np.maximum(Ab, 1e-30)
        ok = Ab > 0
        colour = cmap(0.1 + 0.75 * i / max(1, len(tf_ord) - 1))
        ax[1].plot((yc[ok] + t / 2) * 1e3, Tm[ok], color=colour, lw=1.6, label=f"{tt / 60:.0f} min")
    ax[1].set_xlabel("distance from the exposed face (mm)"); ax[1].set_ylabel("mean ceramic temperature (°C)")
    ax[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    S.figure_(fig, "wall_deformed.png", "Deformed shape and temperatures through the thickness")
    if levels:
        fr = np.array([x[0] for x in levels]); Rs = np.array([x[2] for x in levels])
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        ok = np.isfinite(Rs)
        ax.plot(Rs[ok], 100 * fr[ok], "o-", color=COLOURS[0], lw=2, mec="white", mew=1.2)
        if (~ok).any():
            ax.plot(np.full((~ok).sum(), tf / 60), 100 * fr[~ok], ">", color=COLOURS[0], ms=8)
        if R is not None:
            ax.plot([R], [100 * P / Pmax], "s", color=COLOURS[1], ms=7, mec="white", label="given load")
            ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel("fire resistance (min)"); ax.set_ylabel("load / cold capacity (%)")
        ax.set_xlim(0, tf / 60 * 1.05); ax.set_ylim(0, 100)
        fig.tight_layout()
        S.figure_(fig, "wall_time_load.png", "Fire resistance × load level")
    Te = Tnodes[times[-1]][m2.tri].mean(1)
    section_figure(m2, Te, f"Temperature at {times[-1] / 60:.0f} min (°C)", os.path.join(S.folder, "section_T.png"), "°C")
    S.figures.append({"file": "section_T.png", "title": "Temperature in the section at the end of ISO 834"})
    sty = []
    for tt in tf_ord:
        tn = min(Tnodes, key=lambda x: abs(x - tt))
        sty.append({"label": f"t = {tt / 60:.0f} min", "fields": {"T": Tnodes[tn], "szz": snapshots[tt][1]}})
    S.pos(m2, sty, "Wall: temperature and vertical stress in the fibres at mid-height", material=th.section_material(m2, geo1))
    S.notices.append("Vertical strip model: the vertical edges of the wall do not restrain bending. Material "
                    "homogenised from the curves of the prism heated at 5 °C/min (steady-state test); "
                    + ("with transient creep from the given k_tr; " if flu is not None else
                       "without transient creep (k_tr = 0 or switched off); ")
                    + "spalling and integrity (E) are not included.")
    if M.ceramic.table_thermal_strain is None:
        S.notices.append(f"Thermal strain from α = {M.ceramic.alpha * 1e6:.1f}×10⁻⁶/°C (provisional): the thermal "
                        "deflection and the fire resistance depend directly on this value.")
    if "fixed" in (ts["base_support"], ts["top_support"]):
        S.notices.append("With fixed supports, the expansion of the exposed face is restrained and the ceramic of that face "
                        "crushes next to the support"
                        + ("; without transient creep, the result tends to be conservative." if flu is None else
                           ". Transient creep relieves part of that compression."))


# ----------------------------------------------------------------------------
# Wall sensitivity under ISO 834 (Morris)
# ----------------------------------------------------------------------------
class _SilentOutput(_Output):
    """Output that keeps only the summary (runs of the sensitivity study)."""

    def __init__(self):
        super().__init__(None)

    def pos(self, *a, **k):
        pass

    def table(self, *a, **k):
        pass

    def figure_(self, fig, *a, **k):
        import matplotlib.pyplot as plt
        plt.close(fig)

    def vtk(self, *a, **k):
        pass


def evaluate_wall(q: Dict, law=None) -> Dict:
    """One wall run under ISO 834 without figures: R and I (min) and Pmax (kN)."""
    q = complete(q)
    geo = geometry(q)
    M = materials(q, geo)
    return _run_wall(q, geo, M, _SilentOutput(), lambda s: None, law=law, fast=True)


def _validate_sensitivity(p) -> List[str]:
    from . import sensitivity as sens
    ts, sset = p["test"], p["sensitivity"]
    errors = []
    if ts["kind"] != "iso834":
        errors.append("The sensitivity study is for the wall under ISO 834: choose that test under Test and load.")
    if ts["load"] <= 0:
        errors.append("Give the constant load of the wall (total over the length).")
    pars = sens.active(sset)
    if len(pars) < 2:
        errors.append("Select at least two parameters on the Sensitivity page.")
    t2 = p["geometry"]["width"] / 2
    for d in pars:
        lo, hi, ch = d["lo"], d["hi"], d["key"]
        if d["kind"] == "rel":
            if not 0 < hi < 90:
                errors.append(f"{d['name']}: the range must be between 0 and 90%.")
            continue
        if lo == hi:
            errors.append(f"{d['name']}: the two limits of the range are equal.")
        if ch == "imp" and min(lo, hi) <= 0:
            errors.append("Imperfection: use H divided by positive values (for example, 1000 and 250).")
        if ch == "exc" and max(abs(lo), abs(hi)) >= t2:
            errors.append("Eccentricity: the range must lie within half the wall thickness.")
        if ch in ("eps", "eps_holes") and not (0 <= min(lo, hi) and max(lo, hi) <= 1):
            errors.append(f"{d['name']}: use values between 0 and 1.")
        if ch in ("h", "alpha", "ktr") and min(lo, hi) < 0:
            errors.append(f"{d['name']}: use positive values.")
    if int(sset["trajectories"]) < 2:
        errors.append("Use at least 2 trajectories.")
    if int(sset["levels"]) % 2:
        errors.append("The number of grid levels must be even (4, 6 or 8).")
    return errors


def _sensitivity_base(p):
    """Nominal project of the runs: wall under ISO 834, no load levels, fast mesh."""
    b = complete(copy.deepcopy(p))
    b["model"]["analysis"] = "wall"
    b["test"].update({"kind": "iso834", "load_levels": 0})
    if b["sensitivity"].get("fast_mesh", True):
        b["mesh"]["h"] = max(b["mesh"]["h"], 10)
        b["mesh"]["n_elements"] = min(int(b["mesh"]["n_elements"]), 16)
    return b


def _range_text(d):
    if d["kind"] == "rel":
        return f"±{d['hi']:g}%"
    if d["key"] == "imp":
        return f"H/{d['lo']:g} to H/{d['hi']:g}"
    u = f" {d['unit']}" if d["unit"] else ""
    return f"{d['lo']:g} to {d['hi']:g}{u}"


def _behaviour(ms, sg, mmax):
    """Reading of the pair (μ*, σ): Garcia Sanchez et al. (2014)."""
    if not np.isfinite(ms) or ms < 0.05 * mmax:
        return "small effect"
    r = sg / ms if ms > 0 else np.inf
    return ("almost linear" if r < 0.1 else "monotonic" if r < 0.5 else "almost monotonic" if r < 1.0
            else "non-linear or with interactions")


def _run_sensitivity(p, geo, S, log):
    from . import sensitivity as sens
    sset = p["sensitivity"]
    base = _sensitivity_base(p)
    pars = sens.active(sset)
    k, r, nv = len(pars), int(sset["trajectories"]), int(sset["levels"])
    X = sens.trajectories(k, r, nv, int(sset["seed"]))
    projects = [base] + [sens.apply(base, pars, x) for x in X.reshape(-1, k)]
    procs = int(sset["processes"]) or sens.default_processes()
    msh_ = base["mesh"]
    log(f"Sensitivity (Morris): {k} parameters, {r} trajectories, grid of {nv} levels; "
        f"{len(projects)} wall runs (h = {msh_['h']:g} mm, {_n_elements(base)} beam elements), "
        f"{procs} process(es)")
    S.item("Varied parameters", k)
    S.item("Wall runs", len(projects))
    geo_b = geometry(base)
    law = _wall_law(base, geo_b, log, S)
    t0 = time.time()
    res = sens.run_many(projects, law, procs, log)
    S.item("Mean time per run (with the parallel processes)", (time.time() - t0) / len(projects), "s")
    failures = [(i, x["error"]) for i, x in enumerate(res) if "error" in x]
    nom = res[0]
    if "error" in nom:
        raise RuntimeError("The nominal wall run failed: " + nom["error"])
    body = res[1:]

    def col(ch):
        return np.array([x.get(ch, np.nan) for x in body], float).reshape(r, k + 1)

    Y = {"R": col("R"), "I": col("I"), "Pmax": col("Pmax")}
    cens = {"R": int(sum(x.get("R_censored", False) for x in body)),
            "I": int(sum(x.get("I_censored", False) for x in body))}
    ef = {s: sens.effects(X, Y[s]) for s in Y}
    tf = nom["tf"]
    # ------------------------------------------------------------ summary
    S.item("Nominal fire resistance (screening mesh)", nom["R"] if not nom["R_censored"] else f"> {tf:.0f}",
           "min")
    S.item("Nominal insulation time", nom["I"] if not nom["I_censored"] else f"> {tf:.0f}", "min")
    S.item("Nominal cold capacity", nom["Pmax"], "kN")
    for s, lbl in (("R", "fire resistance"), ("I", "insulation")):
        e = ef[s]
        if np.all(~np.isfinite(e["mu_star"])):
            continue
        i = int(np.nanargmax(e["mu_star"]))
        S.item(f"Largest effect on {lbl}", f"{pars[i]['name']} (μ* = {e['mu_star'][i]:.0f} min)")
    vR = Y["R"][np.isfinite(Y["R"])]
    if vR.size:
        S.item("Fire resistance over the runs (minimum to maximum)", f"{vR.min():.0f} to {vR.max():.0f}", "min")
    # ------------------------------------------------------------ tables and files
    names = [d["name"] for d in pars]
    for s, tit, u in sens.OUTPUTS:
        e = ef[s]
        mmax = np.nanmax(e["mu_star"]) if np.any(np.isfinite(e["mu_star"])) else 0.0
        order = np.argsort(-np.nan_to_num(e["mu_star"], nan=-1))
        lines = [[names[i], _range_text(pars[i]), round(float(e["mu_star"][i]), 2),
                   f"{e['ic_lo'][i]:.1f} to {e['ic_hi'][i]:.1f}", round(float(e["mu"][i]), 2),
                   round(float(e["sigma"][i]), 2),
                   round(float(e["sigma"][i] / e["mu_star"][i]), 2) if e["mu_star"][i] > 1e-9 else "—",
                   _behaviour(e["mu_star"][i], e["sigma"][i], mmax)] for i in order]
        S.table(f"Morris effects on the {tit} ({u})",
                 ["parameter", "band", f"μ* ({u})", "95% CI of μ*", f"μ ({u})", f"σ ({u})", "σ/μ*", "behaviour"],
                 lines, file=f"sensitivity_{s}.csv")
    hdr = ["run", "trajectory", "step"] + [f"{d['short']} ({d['unit']})" if d["unit"] else d["short"]
                                                for d in pars] + ["R (min)", "I (min)", "Pmax (kN)",
                                                                  "deflection at the end (mm)", "remark"]
    lines = []
    for j, (x, rr) in enumerate(zip([None] + list(X.reshape(-1, k)), res)):
        tj, ps = ("nominal", "") if j == 0 else ((j - 1) // (k + 1) + 1, (j - 1) % (k + 1))
        vals = (["nominal"] * k if x is None else
                [round(sens.value(d, xi), 4) if d["kind"] == "abs" else f"{sens.value(d, xi):+.1f}%"
                 for d, xi in zip(pars, x)])
        obs = rr.get("error", "") or "; ".join(n for n, c in (("no failure until the end", rr.get("R_censored")),
                                                               ("no loss of insulation", rr.get("I_censored"))) if c)
        lines.append([j, tj, ps] + vals + [round(rr.get("R", np.nan), 2), round(rr.get("I", np.nan), 2),
                                             round(rr.get("Pmax", np.nan), 1), round(rr.get("deflection", np.nan), 1),
                                             obs])
    S.table("Runs of the study (Morris plan)", hdr, lines, file="sensitivity_runs.csv")
    # ------------------------------------------------------------ figures
    _sensitivity_figures(S, pars, ef)
    # ------------------------------------------------------------ notices
    S.notices.append("Morris is a screening method: it ranks the inputs by the mean effect μ*, the mean change of the "
                    "output when the parameter spans its whole range (μ with sign: + when the output grows from the "
                    "first to the second limit). A high σ indicates a non-linear effect or interaction with other "
                    "inputs. The results hold for the given ranges.")
    if base["sensitivity"].get("fast_mesh", True):
        S.notices.append("Screening with the fast mesh: the times differ by a few minutes from those of the recommended mesh "
                        "(5 mm and 32 elements). For the BCEPQ wall, the order of the parameters was the same with both "
                        "meshes; for another wall, confirm by repeating the study without the fast mesh (about 6 "
                        "times longer).")
    if cens["R"]:
        S.notices.append(f"In {cens['R']} run(s) the wall did not lose equilibrium within {tf:.0f} min; the fire "
                        f"resistance entered as {tf:.0f} min (a lower bound), which reduces the effects on R. "
                        "Increase the ISO 834 duration to measure them in full.")
    if cens["I"]:
        S.notices.append(f"In {cens['I']} run(s) the unexposed face did not reach 140/180 °C within {tf:.0f} min; the "
                        f"insulation time entered as {tf:.0f} min.")
    if failures:
        S.notices.append(f"{len(failures)} run(s) failed and were left out of the effects (first: run "
                        f"{failures[0][0]}, {failures[0][1]}).")
    if any(d["key"] == "ktr" for d in pars) and not p["model"].get("creep", True):
        S.notices.append("Transient creep was switched off in the project; it was switched on in the runs to "
                        "vary the k_tr of the unit.")


def _sensitivity_figures(S, pars, ef):
    plt = _plt()
    from . import sensitivity as sens
    ink, ink2, grid = "#2b2a28", "#6b6a66", "#e4e3df"
    short = [d["short"] for d in pars]
    names = [d["name"] for d in pars]
    oR = np.argsort(np.nan_to_num(ef["R"]["mu_star"], nan=-1))       # order of R, from smallest to largest
    # ---------------------------------------------- μ* with CI (small multiples, same order)
    fig, ax = plt.subplots(1, 3, figsize=(11, 0.34 * len(pars) + 1.6), sharey=True,
                           gridspec_kw={"width_ratios": [1.15, 1.15, 0.9]})
    y = np.arange(len(pars))
    for a, (s, tit, u) in zip(ax, sens.OUTPUTS):
        e = ef[s]
        ms = np.nan_to_num(e["mu_star"][oR])
        a.barh(y, ms, height=0.56, color=COLOURS[0], zorder=2)
        lo, hi = np.nan_to_num(e["ic_lo"][oR]), np.nan_to_num(e["ic_hi"][oR])
        a.hlines(y, lo, hi, color=ink, lw=1.0, zorder=3)
        top = np.argsort(-ms)[:3]
        xmax = max(hi.max(), ms.max(), 1e-9)
        for i in top:
            if ms[i] > 0.02 * xmax:
                a.text(hi[i] + 0.02 * xmax, y[i], f"{ms[i]:.0f}" if xmax >= 10 else f"{ms[i]:.1f}",
                       va="center", fontsize=8, color=ink)
        a.set_xlim(0, xmax * 1.18)
        a.set_title(f"{tit} ({u})", fontsize=9, loc="left", color=ink)
        a.grid(axis="y", visible=False); a.grid(axis="x", color=grid)
        a.tick_params(axis="y", length=0)
        a.set_xlabel(f"μ* ({u})", color=ink2)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([names[i] for i in oR], fontsize=8)
    fig.tight_layout()
    S.figure_(fig, "sensitivity_mu.png", "Mean effect μ* of each parameter, with 95% interval (bootstrap); "
             "same order in the three panels, by fire resistance")
    # ---------------------------------------------- μ*–σ plane
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for a, s, tit in ((ax[0], "R", "fire resistance"), (ax[1], "I", "insulation time")):
        e = ef[s]
        ms, sg = np.nan_to_num(e["mu_star"]), np.nan_to_num(e["sigma"])
        lim = max(ms.max(), sg.max(), 1e-9) * 1.15
        xs = np.array([0, lim])
        for rz, lbl in ((1.0, "σ/μ* = 1"), (0.5, "0.5"), (0.1, "0.1")):
            a.plot(xs, rz * xs, color="#b9b7b1", lw=0.9, ls=(0, (4, 3)), zorder=1)
            xe = lim if rz <= 1 else lim / rz
            a.text(xe * 0.98, rz * xe * 0.98, lbl, fontsize=7, color=ink2, ha="right", va="bottom")
        a.scatter(ms, sg, s=46, color=COLOURS[0], edgecolor="white", linewidth=1.5, zorder=3)
        for i in range(len(pars)):
            if ms[i] >= 0.08 * ms.max():
                a.annotate(short[i], (ms[i], sg[i]), xytext=(5, 3), textcoords="offset points", fontsize=8,
                           color=ink)
        a.set_xlim(0, lim); a.set_ylim(0, lim)
        a.set_xlabel("μ* (min)", color=ink2); a.set_ylabel("σ (min)", color=ink2)
        a.set_title(tit, fontsize=9, loc="left", color=ink)
    fig.tight_layout()
    S.figure_(fig, "sensitivity_mu_sigma.png", "μ*–σ plane: to the right, larger effect; above the line σ/μ* = 1, "
             "non-linear effect or with interactions (small parameters without label)")
