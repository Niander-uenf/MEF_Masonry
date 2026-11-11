"""Library of units, materials and mortars.

Each item carries the changes it applies to the project (group -> fields, in
the units of the interface), its source and a note on what remains to be checked. Items with
"phase" > CURRENT_PHASE do not run yet (geometry under development) and appear
disabled in the interface. Sources and discussion: document "Library of
blocks and bricks of MEF".

Dimension convention: unit_length and course_height are centre to centre
(unit + joint); with a dry joint, the joint is a thin layer inside the unit.
"""
from typing import Dict

from . import legacy
from .parameters import EN1992_FC_FACTOR, EN1992_FT_FACTOR, EN1992_E_FACTOR, EN1992_THERMAL_STRAIN

CURRENT_PHASE = 3


def _ft_mortar(fc):
    """ft = 0.3 fc^(2/3) in MPa (Eurocode 2; used by Dauda, Iuorio and Lourenço, 2020)."""
    return round(0.3 * fc ** (2 / 3), 2)


# ---------------------------------------------------------------------------
# Mortars (Mohamad, Lourenço and Roman, 2006: fa and initial tangent Ea, block M1)
# ---------------------------------------------------------------------------
_THERM_MOR = {"mor_k": 1.15, "mor_rho": 2000, "mor_c": 1000}     # NBR 15220-2, via LabEEE/UFSC

MORTARS: Dict[str, Dict] = {
    "1:0.25:3": {"name": "1 : 0.25 : 3 (cement : lime : sand), fa = 19.9 MPa",
                 "changes": {"ceramic": {"mor_fc": 19.9, "mor_E": 11230, **_THERM_MOR}}},
    "1:0.5:4.5": {"name": "1 : 0.5 : 4.5, fa = 8.6 MPa",
                  "changes": {"ceramic": {"mor_fc": 8.63, "mor_E": 6409, **_THERM_MOR}}},
    "1:1:6": {"name": "1 : 1 : 6, fa = 4.2 MPa",
              "changes": {"ceramic": {"mor_fc": 4.20, "mor_E": 4033, **_THERM_MOR}}},
    "1:2:9": {"name": "1 : 2 : 9, fa = 2.3 MPa",
              "changes": {"ceramic": {"mor_fc": 2.29, "mor_E": 2042, **_THERM_MOR}}},
    "custom": {"name": "Custom (values given below)", "changes": {}},
}
for _k, _v in MORTARS.items():
    c = _v["changes"].get("ceramic")
    if c:
        c["mor_ft"] = _ft_mortar(c["mor_fc"])
        c["mor_nu"] = 0.20
        _v["source"] = "Mohamad, Lourenço and Roman (2006); thermal properties from NBR 15220-2; ft = 0.3 fc^(2/3)"
        _v["note"] = ("fa and Ea measured on 5 × 10 cm cylinders (block M1). ν = 0.20 is nominal and the mortar α "
                      "is provisional.")
    else:
        _v["source"], _v["note"] = "", "Give the measured properties of the mortar."


# ---------------------------------------------------------------------------
# Grouts (Camacho et al., 2015: mixes by mass, slump 200 ± 30 mm)
# ---------------------------------------------------------------------------
def _E_concrete(fc):
    """Ecm = 22 (fcm/10)^0.3 GPa (EN 1992-1-1), in MPa."""
    return round(22000 * (fc / 10) ** 0.3, -2)


_THERM_GRT = {"grt_k": 1.75, "grt_rho": 2400, "grt_c": 1000}       # concrete, NBR 15220-2
GROUTS: Dict[str, Dict] = {
    "G1": {"name": "G1: 1 : 1.83 : 2.17 (cement : sand : gravel, by mass), w/c 0.67, fgm = 17.7 MPa",
           "changes": {"ceramic": {"grt_fc": 17.73, **_THERM_GRT}}},
    "G2": {"name": "G2: 1 : 1.20 : 1.80, w/c 0.55, fgm = 26.7 MPa",
           "changes": {"ceramic": {"grt_fc": 26.66, **_THERM_GRT}}},
    "custom": {"name": "Custom (values given below)", "changes": {}},
}
for _k, _v in GROUTS.items():
    c = _v["changes"].get("ceramic")
    if c:
        c["grt_ft"] = _ft_mortar(c["grt_fc"])
        c["grt_E"] = _E_concrete(c["grt_fc"])
        c["grt_nu"] = 0.20
        _v["source"] = "Camacho et al. (2015); E from EN 1992-1-1; ft = 0.3 fc^(2/3); thermal properties from NBR 15220-2"
        _v["note"] = ("Only fgm was measured; E, ft and ν are estimates. Hot, the grout follows the siliceous concrete "
                      "of EN 1992-1-2 (option below).")
    else:
        _v["source"], _v["note"] = "", "Give the measured properties of the grout."


# ---------------------------------------------------------------------------
# Unit materials
# ---------------------------------------------------------------------------
def _concrete_k_upper():
    """Conductivity of concrete, upper limit of EN 1992-1-2 (W/m.K)."""
    return [[T, round(2 - 0.2451 * (T / 100) + 0.0107 * (T / 100) ** 2, 3)] for T in range(20, 1201, 100)]


def _concrete_c_dry():
    """Specific heat of dry concrete, EN 1992-1-2 (J/kg.K), without the moisture peak."""
    return [[20, 900], [100, 900], [200, 1000], [400, 1100], [1200, 1100]]


def _tab(t, f=1.0):
    """Table [(T, v)] -> [[T, v*f]] (interface format)."""
    return [[T, round(v * f, 6)] for T, v in t]


MATERIALS: Dict[str, Dict] = {
    "bcepq": {"name": "BCEPQ ceramic (measured)",
              "changes": {"ceramic": {"k": 0.52, "rho": 1210, "c": 1285, "E_unit_gross": 2801, "nu": 0.18,
                                        "alpha": 6.0, "table_k": [], "table_c": [], "table_rho": [],
                                        "table_E_factor": [], "table_thermal_strain": []},
                           "model": {"table_ft_factor": [], "table_fc_factor": [], "k_tr_unit": 0.0,
                                      "table_creep_unit": []}},
              "source": "BCEPQ manuscript and Cerâmica 2026",
              "note": "α still provisional (no dilatometry)."},
    "clay_nbr15220": {"name": "Fired clay, NBR 15220-2 (1600 kg/m³)",
                          "changes": {"ceramic": {"k": 0.90, "rho": 1600, "c": 920, "table_k": [],
                                                    "table_c": [], "table_rho": [], "table_E_factor": [],
                                                    "table_thermal_strain": []},
                                       "model": {"table_ft_factor": [], "table_fc_factor": [], "k_tr_unit": 0.0,
                                      "table_creep_unit": []}},
                          "source": "NBR 15220-2, via the LabEEE/UFSC catalogue",
                          "note": "Only the thermal properties come from the standard: give E, ν and α of the unit or calibrate."},
    "concrete_nbr15220": {"name": "Concrete (block), NBR 15220-2 + EN 1992-1-2",
                          "changes": {"ceramic": {"k": 1.75, "rho": 2400, "c": 1000, "nu": 0.20, "alpha": 10.0,
                                                    "table_k": _concrete_k_upper(),
                                                    "table_c": _concrete_c_dry(), "table_rho": [],
                                                    "table_E_factor": _tab(EN1992_E_FACTOR),
                                                    "table_thermal_strain": _tab(EN1992_THERMAL_STRAIN, 1e3)},
                                       "model": {"table_fc_factor": _tab(EN1992_FC_FACTOR),
                                                  "table_ft_factor": _tab(EN1992_FT_FACTOR),
                                                  "k_tr_unit": 0.0, "table_creep_unit": []}},
                          "source": ("NBR 15220-2 (20 °C); EN 1992-1-2, siliceous aggregate: k (upper limit), "
                                    "dry c, fc(θ) and εc1(θ) from Table 3.1, ft(θ) and thermal strain"),
                          "note": ("Curves of structural concrete, applied to the unit as an approximation. The specific heat "
                                   "peak of moisture (100 to 200 °C) is not included. The E factor is fc(θ)/εc1(θ) "
                                   "relative to 20 °C. The thermal strain was not checked against an open copy of the "
                                   "standard. As the EN curves already include transient creep, k_tr of the unit is "
                                   "0. Give E, fc and ft of the unit.")},
    "custom": {"name": "Custom (values given below)", "changes": {},
                      "source": "", "note": "Give the measured properties of the unit."},
}


# ---------------------------------------------------------------------------
# Units (blocks and bricks)
# ---------------------------------------------------------------------------
def _geo(L, H, C, joint, **kw):
    """Centre-to-centre geometry from the actual dimensions (mm) and the joint."""
    dry = joint == 0
    g = {"width": L, "course_height": H + joint, "unit_length": C + joint,
         "joint_type": "dry" if dry else "mortar", "joint_thickness": 2 if dry else joint}
    g.update(kw)
    return g


_CONCRETE = {**MATERIALS["concrete_nbr15220"]["changes"]["ceramic"]}
_STRUCTURAL_CLAY = {**MATERIALS["clay_nbr15220"]["changes"]["ceramic"], "nu": 0.15, "alpha": 6.0}
# ductility indices: concrete from Model Code 2010 (G_f = 73 fcm^0.18 N/m; G_f/ft ≈ 0.1 mm) and G_c/fc
# as for the mortar; fired clay with the values of clay units (Dauda et al., 2020)
_CONCRETE_MODEL = {"du_t": 0.10, "du_c": 1.6, **MATERIALS["concrete_nbr15220"]["changes"]["model"]}
_CLAY_MODEL = {"du_t": 0.018, "du_c": 0.33, "table_ft_factor": [], "table_fc_factor": [],
                    "k_tr_unit": 0.0, "table_creep_unit": []}

UNITS: Dict[str, Dict] = {
    "bcepq": {
        "name": "BCEPQ (interlocking unit, UFV): 140 × 60 × 280 mm, 2 holes of 84 mm, dry joint",
        "phase": 0,
        "changes": {"geometry": {"unit_length": 280, "width": 140, "course_height": 60, "hole_diameter": 84,
                                   "holes_per_unit": 2, "units_per_course": 1.5, "n_courses": 5,
                                   "bond": "running_bond", "joint_type": "dry", "joint_thickness": 2,
                                   "capping": True, "hole_shape": "circular", "grouting": "none"},
                     "ceramic": {**MATERIALS["bcepq"]["changes"]["ceramic"], "material_library": "bcepq"},
                     "model": {**_CLAY_MODEL, "ft": 0.66, "fc": 2.48}},
        "source": "UFV tests and BCEPQ manuscript",
        "note": "fc = 2.48 MPa is the effective value calibrated on the prism at 30 °C."},
    "solid_brick_10x6x22": {
        "name": "Solid clay brick 100 × 60 × 220 mm, 1:1:6 mortar of 10 mm",
        "phase": 1,
        "changes": {"geometry": _geo(100, 60, 220, 10, hole_diameter=0, holes_per_unit=0, units_per_course=1,
                                       n_courses=5, bond="stack_bond", bedding="total", hole_shape="circular",
                                       grouting="none",
                                       filled_head_joints=True, capping=True),
                     "ceramic": {**MATERIALS["clay_nbr15220"]["changes"]["ceramic"],
                                  **MORTARS["1:1:6"]["changes"]["ceramic"],
                                  "material_library": "clay_nbr15220", "mortar_library": "1:1:6"},
                     "model": {**_CLAY_MODEL, "ft": 0.66, "fc": 5.76}},
        "source": "Dimensions and thermal properties: LabEEE/UFSC catalogue (NBR 15220-2); mortar: Mohamad et al. (2006)",
        "note": ("E, ν, α, fc and ft of the brick kept the values of the BCEPQ ceramic (material fc = 5.76 MPa, "
                 "ft = 0.66 MPa from Gooch et al., 2023): replace them with those of your batch.")},
    "soil_cement_modular": {
        "name": "Modular soil-cement brick 125 × 62.5 × 250 mm, 2 holes — under development (measure the holes)",
        "phase": 9, "changes": {}, "source": "Floresta/Sahara; Garcia Jr. (2025)",
        "note": "The hole diameter and the properties of the soil-cement have no open source."},
    "concrete_m15_39": {
        "name": "Concrete block 140 × 190 × 390 mm, 2 cells, 1:0.5:4.5 mortar (Mohamad et al., 2006)",
        "phase": 2,
        "changes": {"geometry": _geo(140, 190, 390, 10, hole_shape="rectangular", holes_per_unit=2,
                                       face_shell=32, web=30, end_shell=30, wall_cells=0,
                                       hole_diameter=0, units_per_course=1, n_courses=3, bond="stack_bond",
                                       bedding="total", filled_head_joints=True, capping=True,
                                       grouting="none"),
                     "ceramic": {**_CONCRETE, "E_unit_gross": 9160,
                                  **MORTARS["1:0.5:4.5"]["changes"]["ceramic"],
                                  "material_library": "concrete_nbr15220", "mortar_library": "1:0.5:4.5"},
                     "model": {**_CONCRETE_MODEL, "ft": 1.14, "fc": 18.5}},
        "source": "Mohamad, Lourenço and Roman (2006); NBR 6136; thermal properties from NBR 15220-2 and EN 1992-1-2",
        "note": ("fb = 18.5 MPa on the net area (ratio 0.58). The shells are 25 mm at the top and tapered: the model uses "
                 "equivalent prismatic shells (32 mm on the faces, 30 mm on the webs and ends), which give the same net "
                 "area. ft = 0.1137 fb^0.789 (same source). Unit E not measured: the net E of 15.8 GPa comes from the "
                 "prism modulus (15.15 GPa) through the series association. Check the number of courses and the "
                 "bedding of the prism in the source.")},
    "concrete_m15_29": {
        "name": "Concrete block 140 × 190 × 290 mm (B1), 3-block prism, grout G1 (Camacho et al., 2015)",
        "phase": 2,
        "changes": {"geometry": _geo(140, 190, 290, 10, hole_shape="rectangular", holes_per_unit=2,
                                       face_shell=28, web=28, end_shell=28, wall_cells=0,
                                       hole_diameter=0, units_per_course=1, n_courses=3, bond="stack_bond",
                                       bedding="total", filled_head_joints=True, capping=True,
                                       grouting="all"),
                     "ceramic": {**_CONCRETE, "E_unit_gross": 5300,
                                  **MORTARS["1:0.5:4.5"]["changes"]["ceramic"],
                                  "mor_fc": 6.10, "mor_ft": _ft_mortar(6.10), "mor_E": 5050,
                                  **GROUTS["G1"]["changes"]["ceramic"],
                                  "material_library": "concrete_nbr15220", "mortar_library": "custom",
                                  "grout_library": "G1"},
                     "model": {**_CONCRETE_MODEL, "ft": 0.98, "fc": 15.39}},
        "source": "Camacho, Logullo, Parsekian and Soudais (2015); NBR 6136",
        "note": ("Block B1: fb = 8.64 MPa gross and 15.39 MPa net (ratio 0.57; equivalent shells of 28 mm). "
                 "Measured prisms (gross area): hollow 5.63 MPa; with grout G1 11.29 MPa; with G2 10.50 MPa. Mortar "
                 "1:0.5:4.5 with 6.10 MPa (E = 5050 MPa, interpolated between the mortars of Mohamad et al., 2006). Unit E estimated from the "
                 "failure strain (3.3‰): net E ≈ 2 fb/εu = 9.3 GPa. ft = 0.1137 fb^0.789.")},
    "concrete_m20_39": {
        "name": "Concrete block M-20 190 × 190 × 390 mm, class A (NBR 6136), 2 cells",
        "phase": 2,
        "changes": {"geometry": _geo(190, 190, 390, 10, hole_shape="rectangular", holes_per_unit=2,
                                       face_shell=32, web=25, end_shell=25, wall_cells=0,
                                       hole_diameter=0, units_per_course=1, n_courses=3, bond="stack_bond",
                                       bedding="total", filled_head_joints=True, capping=True,
                                       grouting="none"),
                     "ceramic": {**_CONCRETE, "E_unit_gross": 7340,
                                  **MORTARS["1:0.5:4.5"]["changes"]["ceramic"],
                                  "material_library": "concrete_nbr15220", "mortar_library": "1:0.5:4.5"},
                     "model": {**_CONCRETE_MODEL, "ft": 1.07, "fc": 17.2}},
        "source": "NBR 6136 (minimum shells of class A); other values estimated",
        "note": ("Shells at the minimum thicknesses (32 mm on the faces, 25 mm on the webs and ends; net area 0.46). "
                 "fc = minimum fbk of class A (8 MPa gross) on the net area; ft = 0.1137 fb^0.789; net E of "
                 "15.8 GPa as for the M-15 of Mohamad et al. Replace with the values of your batch.")},
    "concrete_c_9_39": {
        "name": "Non-load-bearing concrete block class C 90 × 190 × 390 mm, 2 cells (NBR 6136)",
        "phase": 2,
        "changes": {"geometry": _geo(90, 190, 390, 10, hole_shape="rectangular", holes_per_unit=2,
                                       face_shell=18, web=18, end_shell=18, wall_cells=0,
                                       hole_diameter=0, units_per_course=1, n_courses=3, bond="stack_bond",
                                       bedding="total", filled_head_joints=True, capping=True,
                                       grouting="none"),
                     "ceramic": {**_CONCRETE, "E_unit_gross": 3000,
                                  **MORTARS["1:1:6"]["changes"]["ceramic"],
                                  "material_library": "concrete_nbr15220", "mortar_library": "1:1:6"},
                     "model": {**_CONCRETE_MODEL, "ft": 0.48, "fc": 6.2}},
        "source": "NBR 6136 (18 mm shells, class C); other values estimated",
        "note": ("fc = minimum fbk of class C (3 MPa gross) on the net area (ratio 0.48); ft = 0.1137 fb^0.789; "
                 "net E ≈ 1000 fb (order of magnitude). Replace with the values of your batch.")},
    "struct_clay_14x19x29_solid": {
        "name": "Structural clay block 140 × 190 × 290 mm, solid shells, 2 cells (NBR 15270)",
        "phase": 2,
        "changes": {"geometry": _geo(140, 190, 290, 10, hole_shape="rectangular", holes_per_unit=2,
                                       face_shell=20, web=8, end_shell=20, wall_cells=0,
                                       hole_diameter=0, units_per_course=1, n_courses=3, bond="stack_bond",
                                       bedding="total", filled_head_joints=True, capping=True,
                                       grouting="none"),
                     "ceramic": {**_STRUCTURAL_CLAY, "E_unit_gross": 6000,
                                  **MORTARS["1:1:6"]["changes"]["ceramic"],
                                  "material_library": "clay_nbr15220", "mortar_library": "1:1:6"},
                     "model": {**_CLAY_MODEL, "ft": 1.49, "fc": 14.9}},
        "source": "NBR 15270-2 (face shells ≥ 20 mm, webs ≥ 8 mm); thermal properties from NBR 15220-2",
        "note": ("Geometry at the minimum thicknesses (net area 0.40); real blocks have thicker webs: measure "
                 "yours. fc of a 6 MPa gross block on the net area; ft = 0.10 fc and net E ≈ 1000 fb are "
                 "order-of-magnitude estimates.")},
    "struct_clay_14x19x29_hollow": {
        "name": "Structural clay block 140 × 190 × 290 mm, hollow shells (NBR 15270)",
        "phase": 2,
        "changes": {"geometry": _geo(140, 190, 290, 10, hole_shape="rectangular", holes_per_unit=2,
                                       face_shell=30, web=8, end_shell=8, wall_cells=6,
                                       min_shell=8, hole_diameter=0, units_per_course=1, n_courses=3,
                                       bond="stack_bond", bedding="total", filled_head_joints=True,
                                       capping=True, grouting="none"),
                     "ceramic": {**_STRUCTURAL_CLAY, "E_unit_gross": 4700,
                                  **MORTARS["1:1:6"]["changes"]["ceramic"],
                                  "material_library": "clay_nbr15220", "mortar_library": "1:1:6"},
                     "model": {**_CLAY_MODEL, "ft": 1.91, "fc": 19.1}},
        "source": "NBR 15270-2 (shells ≥ 8 mm, webs ≥ 7 mm); thermal properties from NBR 15220-2",
        "note": ("Illustrative: face shells of 30 mm with a row of 6 cells of 14 mm, internal walls of "
                 "8 mm (net area 0.31). fc of a 6 MPa gross block on the net area; ft = 0.10 fc and net E "
                 "≈ 1000 fb are estimates. Measure your block.")},
    "struct_clay_14x19x39": {
        "name": "Structural clay block 140 × 190 × 390 mm, solid shells, 3 cells (typical layout)",
        "phase": 2,
        "changes": {"geometry": _geo(140, 190, 390, 10, hole_shape="rectangular", holes_per_unit=3,
                                       face_shell=20, web=8, end_shell=20, wall_cells=0,
                                       hole_diameter=0, units_per_course=1, n_courses=3, bond="stack_bond",
                                       bedding="total", filled_head_joints=True, capping=True,
                                       grouting="none"),
                     "ceramic": {**_STRUCTURAL_CLAY, "E_unit_gross": 6000,
                                  **MORTARS["1:1:6"]["changes"]["ceramic"],
                                  "material_library": "clay_nbr15220", "mortar_library": "1:1:6"},
                     "model": {**_CLAY_MODEL, "ft": 1.49, "fc": 14.9}},
        "source": "NBR 15270-2; 3-cell layout to be checked in the manufacturer's catalogue",
        "note": ("Minimum thicknesses of the standard; with 3 cells the running-bond joints do not align the cells "
                 "(the model accepts it, without grout). Properties estimated as for the 29 cm block.")},
    "infill_8_holes": {
        "name": "Clay infill block 90 × 190 × 190 mm, 8 horizontal holes (2 × 4), 1:2:9 mortar",
        "phase": 3,
        "changes": {"geometry": _geo(90, 190, 190, 10, hole_shape="horizontal", holes_width=2, holes_height=4,
                                       face_shell=8, web=7, bed_shell=8, hole_diameter=0,
                                       holes_per_unit=0, units_per_course=1, n_courses=3, bond="stack_bond",
                                       bedding="total", filled_head_joints=True, capping=True,
                                       grouting="none"),
                     "ceramic": {**_STRUCTURAL_CLAY, "E_unit_gross": 1500,
                                  **MORTARS["1:2:9"]["changes"]["ceramic"],
                                  "material_library": "clay_nbr15220", "mortar_library": "1:2:9"},
                     "model": {**_CLAY_MODEL, "ft": 0.59, "fc": 5.87}},
        "source": "SINAPI 7271 (dimensions and holes); NBR 15270-1 (minimum shells and strength); NBR 15220-2",
        "note": ("Face and bed shells of 8 mm and webs of 7 mm (above the minima of 7 and 6 mm): holes of "
                 "33.5 × 38.3 mm. The walls that run down continuously (faces and vertical web) add up to 26% of the width; "
                 "fc comes from the minimum strength of 1.5 MPa on the gross area, brought to that area. ft = 0.10 fc and net E "
                 "≈ 1000 fb are estimates. Measure your block.")},
    "infill_6_holes": {
        "name": "Clay infill block 90 × 140 × 240 mm, 6 horizontal holes (2 × 3), 1:2:9 mortar (Guedes et al., 2019)",
        "phase": 3,
        "changes": {"geometry": _geo(90, 140, 240, 10, hole_shape="horizontal", holes_width=2, holes_height=3,
                                       face_shell=7.8, web=7.1, bed_shell=7.8, hole_diameter=0,
                                       holes_per_unit=0, units_per_course=1, n_courses=3, bond="stack_bond",
                                       bedding="total", filled_head_joints=True, capping=True,
                                       grouting="none"),
                     "ceramic": {**_STRUCTURAL_CLAY, "E_unit_gross": 1500,
                                  **MORTARS["1:2:9"]["changes"]["ceramic"],
                                  "material_library": "clay_nbr15220", "mortar_library": "1:2:9"},
                     "model": {**_CLAY_MODEL, "ft": 0.60, "fc": 5.95}},
        "source": "Guedes et al. (2019): shells and webs measured in three brickworks; NBR 15270-1; NBR 15220-2",
        "note": ("Face shells of 7.8 mm and webs of 7.1 mm (means of the measured ranges, 7.4 to 8.3 and 6.3 to "
                 "8.0 mm); bed shells equal to the face shells. fc comes from the minimum strength of 1.5 MPa gross, "
                 "brought to the area of the continuous walls (25%); ft = 0.10 fc and net E ≈ 1000 fb are estimates.")},
    "custom": {"name": "Custom (geometry given below)", "phase": 0, "changes": {},
                      "source": "", "note": "Give the geometry of the unit."},
}


def options(table: Dict[str, Dict]):
    return [[k, v["name"]] for k, v in table.items()]


def for_interface():
    """Changes and notes of each item, by selector field (sent to the interface)."""
    def pack(tab):
        return {k: {"changes": v["changes"], "note": v.get("note", ""), "source": v.get("source", "")}
                for k, v in tab.items()}
    return {"geometry.unit_library": pack(UNITS),
            "ceramic.material_library": pack(MATERIALS),
            "ceramic.mortar_library": pack(MORTARS),
            "ceramic.grout_library": pack(GROUTS)}


def apply(p: Dict, selector: str, key: str) -> Dict:
    """Applies an item to project p (copy), as the interface does."""
    import copy
    selector, key = legacy.library_item(selector, key)
    tab = {"unit": UNITS, "material": MATERIALS, "mortar": MORTARS, "grout": GROUTS}[selector]
    q = copy.deepcopy(p)
    for g, d in tab[key]["changes"].items():
        q.setdefault(g, {}).update(copy.deepcopy(d))
    field = {"unit": ("geometry", "unit_library"), "material": ("ceramic", "material_library"),
             "mortar": ("ceramic", "mortar_library"), "grout": ("ceramic", "grout_library")}[selector]
    q[field[0]][field[1]] = key
    return q
