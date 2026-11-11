"""Model parameters, with the origin of each value.

Internal units: SI (m, s, K or °C, Pa, W). Geometry inputs are also in
metres.
"""
from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np

Table = Optional[Sequence[Tuple[float, float]]]  # [(T in °C, value), ...]


def interpolate_(T, table: Table, default):
    """Interpolates a (T, value) table linearly; outside it, uses the end values.
    Without a table, returns the default (constant) value."""
    if table is None:
        return np.full(np.shape(T), default, dtype=float) if np.ndim(T) else float(default)
    t = np.asarray(table, dtype=float)
    return np.interp(T, t[:, 0], t[:, 1])


@dataclass
class Ceramic:
    """BCEPQ ceramic (extruded, pressed and fired unit)."""
    k: float = 0.52            # W/m.K   Cerâmica 2026 (mean of 20 measurements: 0.5188)
    rho: float = 1210.0        # kg/m³   BPQ manuscript (bulk density)
    c: float = 1285.0          # J/kg.K  derived: effusivity 898 W.s^0.5/m².K and k (Cerâmica 2026)
    E_unit_gross: float = 2801e6  # Pa  unit modulus (BPQ, n = 6), referred to the gross area
    net_area_ratio: float = 0.717  # net / gross area of the section with 84 mm holes
    nu: float = 0.18           # BPQ (0,182)
    alpha: float = 6e-6         # 1/K     PROVISIONAL, not measured (the paper uses 6e-6)
    T_ref: float = 30.0        # °C      temperature without thermal strain
    # Temperature dependence (optional). E.g.: table_k = [(20, 0.52), (1000, 0.60)]
    table_k: Table = None
    table_E_factor: Table = None   # multiplying factor on E: [(20, 1.0), (800, 0.6)]
    table_thermal_strain: Table = None  # total free thermal strain: [(30, 0), (573, 0.0035), ...]
    table_c: Table = None         # specific heat (J/kg.K) with temperature (moisture peak, for example)
    table_rho: Table = None       # density (kg/m³) with temperature (water loss)

    @property
    def E_material(self):
        """E of the ceramic material (net area). If the unit E was measured with
        stress on the gross area, the material is stiffer: E_mat = E_gross / (A_net/A_gross)."""
        return self.E_unit_gross / self.net_area_ratio

    def k_T(self, T):
        return interpolate_(T, self.table_k, self.k)

    def c_T(self, T):
        return interpolate_(T, self.table_c, self.c)

    def rho_T(self, T):
        return interpolate_(T, self.table_rho, self.rho)

    @property
    def variable_capacity(self):
        return self.table_c is not None or self.table_rho is not None

    def E_T(self, T):
        return self.E_material * interpolate_(T, self.table_E_factor, 1.0)

    def thermal_strain(self, T):
        if self.table_thermal_strain is not None:
            return interpolate_(T, self.table_thermal_strain, 0.0)
        return self.alpha * (np.asarray(T, float) - self.T_ref)


@dataclass
class DryJoint:
    """Dry joint between courses (equivalent thin layer).
    E_joint must be calibrated so that the prism reproduces the measured modulus
    (example 02). The default value comes from the calibration with E_prism = 1,100 MPa."""
    thickness: float = 0.002   # m
    E: float = 84.4e6          # Pa: calibrated (example 02) for E_prism = 1,100 MPa, mesh h = 8 mm, dz = 12 mm
    nu: float = 0.18
    k_factor: float = 1.0       # joint conductivity = k_factor * k of the ceramic


@dataclass
class Mortar(Ceramic):
    """Bedding mortar (continuous layer of real thickness, bonded to the
    unit). Thermal: NBR 15220-2 (via the LabEEE/UFSC catalogue). Mechanical, for
    the 1:1:6 mix (cement : lime : sand, by volume): Mohamad, Lourenço and Roman
    (2006), fa = 4.2 MPa and Ea = 4033 MPa (initial tangent modulus). E_material
    = E_unit_gross (net_area_ratio = 1)."""
    k: float = 1.15
    rho: float = 2000.0
    c: float = 1000.0
    E_unit_gross: float = 4033e6
    net_area_ratio: float = 1.0
    nu: float = 0.20            # nominal
    alpha: float = 10e-6         # PROVISIONAL: order of magnitude for cement mortars


# ---------------------------------------------------------------------------
# Hot siliceous-aggregate concrete (EN 1992-1-2). Table 3.1: fc(θ)/fck and
# εc1(θ); 3.2.2.2: kc,t(θ) = 1 up to 100 °C, falling linearly to 0 at 600 °C;
# 3.3.1: thermal strain -1.8e-4 + 9e-6 θ + 2.3e-11 θ³ up to 700 °C and 14e-3
# above. The E factor is the ratio of the secant moduli at the peak, fc(θ)/εc1(θ),
# with a floor of 0.001 (residual stiffness).
# ---------------------------------------------------------------------------
EN1992_FC_FACTOR = [(20, 1.0), (100, 1.0), (200, 0.95), (300, 0.85), (400, 0.75), (500, 0.60), (600, 0.45),
                   (700, 0.30), (800, 0.15), (900, 0.08), (1000, 0.04), (1100, 0.01), (1200, 0.0)]
_EN1992_EC1 = [0.0025, 0.0040, 0.0055, 0.0070, 0.0100, 0.0150, 0.0250, 0.0250, 0.0250, 0.0250, 0.0250, 0.0250, 0.0250]
EN1992_FT_FACTOR = [(20, 1.0), (100, 1.0), (600, 0.0), (1200, 0.0)]
EN1992_E_FACTOR = [(T, round(max(f * 0.0025 / e, 0.001), 4)) for (T, f), e in zip(EN1992_FC_FACTOR, _EN1992_EC1)]


def _eps_siliceous(T):
    return -1.8e-4 + 9e-6 * T + 2.3e-11 * T ** 3 if T <= 700 else 14e-3


EN1992_THERMAL_STRAIN = [(T, round(_eps_siliceous(T), 7)) for T in list(range(20, 701, 20)) + [1200]]
EN1992_THERMAL_STRAIN = sorted(dict(EN1992_THERMAL_STRAIN).items())


@dataclass
class Grout(Ceramic):
    """Grout of the holes or cells (fine fluid concrete). Default: grout G1 of
    Camacho et al. (2015), fgk = 17.7 MPa; E from EN 1992-1-2/EN 1992-1-1,
    Ecm = 22 (fcm/10)^0.3 GPa (estimate, not measured); thermal properties of normal
    concrete from NBR 15220-2. E_material = E_unit_gross (net_area_ratio = 1)."""
    k: float = 1.75
    rho: float = 2400.0
    c: float = 1000.0
    E_unit_gross: float = 26.1e9
    net_area_ratio: float = 1.0
    nu: float = 0.20
    alpha: float = 10e-6         # PROVISIONAL: siliceous concrete, order of magnitude


@dataclass
class Steel:
    k: float = 45.0
    rho: float = 7850.0
    c: float = 500.0
    E: float = 200e9
    nu: float = 0.30
    alpha: float = 12e-6


@dataclass
class HoleAir:
    """Air in the holes: conduction + radiation between the walls of the hole,
    represented by an equivalent conductivity k_rad = 4 σ T³ D ε/(2-ε)."""
    wall_emissivity: float = 0.9
    rho_c: float = 1.0e3       # J/m³.K (negligible)


@dataclass
class Geometry:
    """Parametric geometry of the prism of units with vertical round holes,
    aligned at mid-width. The default is the BCEPQ prism of UFV: five
    courses of one and a half units (420 x 140 mm), two 84 mm holes per unit,
    staggered head joints (even course: unit 0-280 and half unit
    280-420 mm; odd: half unit 0-140 and unit 140-420)."""
    unit_length: float = 0.280
    width: float = 0.140
    course_height: float = 0.060
    hole_diameter: float = 0.084       # 0 = solid unit
    n_courses: int = 5
    plate_thickness: float = 0.025
    holes_per_unit: int = 2
    units_per_course: float = 1.5      # prism length in units (multiple of 0.5)
    bond: str = "running_bond"     # "running_bond" (staggered head joints) or "stack_bond"
    # x positions of the head joints in each course; None = computed from the bond.
    # Accepts one number per course (old format) or one tuple per course.
    head_joints: Optional[Tuple] = None
    # Joints: "dry" (thin layer with contact, gaps and friction) or "mortar"
    # (continuous layer of real thickness). With mortar, unit_length and
    # course_height are centre-to-centre dimensions (unit + joint).
    joint_type: str = "dry"
    bedding: str = "total"         # mortar over the whole face ("total") or on the face shells only ("lateral")
    filled_head_joints: bool = True
    joint_thickness: float = 0.002      # m (defines the real extent of the unit in the rectangular cells)
    # Holes: "circular" (round vertical holes in line, hole_diameter),
    # "rectangular" (vertical cells: holes_per_unit main cells along
    # the length, face shells, webs and end shells; with
    # wall_cells > 0, each face shell has a row of small
    # cells, as in clay blocks with hollow shells) or "horizontal"
    # (grid of holes_width × holes_height rectangular holes running through
    # the unit along its length, as in infill blocks: face
    # shells, webs between the holes and bed shells at top and bottom).
    hole_shape: str = "circular"
    face_shell: float = 0.025       # m, longitudinal walls (faces of the unit)
    web: float = 0.025                # m, internal cross webs
    end_shell: float = 0.025          # m, transverse walls at the ends of the unit
    wall_cells: int = 0             # small cells per unit in each face shell
    min_shell: float = 0.008        # m, walls around the small cells
    grouting: str = "none"        # "none", "all" or "alternate" (main holes or cells)
    holes_width: int = 2              # horizontal holes across the unit thickness
    holes_height: int = 4               # horizontal holes over the unit height
    bed_shell: float = 0.008         # m, top and bottom walls of units with horizontal holes

    @property
    def rectangular(self):
        return self.hole_shape == "rectangular" and self.holes_per_unit > 0

    @property
    def horizontal(self):
        return self.hole_shape == "horizontal" and self.holes_width > 0 and self.holes_height > 0

    @property
    def has_holes(self):
        """Vertical holes (round or cells): those that receive grout."""
        return self.rectangular or (self.hole_shape == "circular" and self.hole_diameter > 0
                                   and self.holes_per_unit > 0)

    @property
    def shell_thickness(self):
        """Thickness of the face shells (for face-shell bedding)."""
        if self.rectangular or self.horizontal:
            return self.face_shell
        if self.hole_diameter <= 0:
            return self.width / 2
        return (self.width - self.hole_diameter) / 2

    @property
    def hole_dimension(self):
        """Hole dimension across the thickness (radiation in the air of the holes)."""
        if self.rectangular:
            return self.width - 2 * self.face_shell
        if self.horizontal:
            return self.hole_width_h
        return self.hole_diameter

    # ---------------------------------------------------------------- horizontal holes
    @property
    def hole_width_h(self):
        """Width (across the unit thickness) of each horizontal hole."""
        n = self.holes_width
        return (self.width - 2 * self.face_shell - (n - 1) * self.web) / n

    def hole_height_h(self, t_joint=None):
        """Height of each horizontal hole; the unit is course_height - joint high."""
        tj = self.joint_thickness if t_joint is None else t_joint
        n = self.holes_height
        return (self.course_height - tj - 2 * self.bed_shell - (n - 1) * self.web) / n

    def holes_h_y(self):
        """Bands (y0, y1) of the columns of horizontal holes."""
        if not self.horizontal:
            return []
        w, te, tw = self.hole_width_h, self.face_shell, self.web
        return [(te + k * (w + tw), te + k * (w + tw) + w) for k in range(self.holes_width)]

    def holes_h_z(self, zc, t_joint=None):
        """Bands (z0, z1) of the rows of horizontal holes of a unit centred at zc."""
        if not self.horizontal:
            return []
        tj = self.joint_thickness if t_joint is None else t_joint
        hz, tb, tw = self.hole_height_h(tj), self.bed_shell, self.web
        z0 = zc - (self.course_height - tj) / 2 + tb
        return [(z0 + k * (hz + tw), z0 + k * (hz + tw) + hz) for k in range(self.holes_height)]

    def course_units(self, i):
        """Intervals (a, b) of the units and half units of course i, centre to centre."""
        xs = [0.0] + list(self.joints_per_course[i]) + [self.length]
        return list(zip(xs[:-1], xs[1:]))

    def cells(self, i, t_joint=None):
        """Rectangular cells of course i: [(x0, x1, y0, y1, main)], in m."""
        if not self.rectangular:
            return []
        tj = self.joint_thickness if t_joint is None else t_joint
        W, Lt = self.width, self.length
        te, tw, tt = self.face_shell, self.web, self.end_shell
        res = []
        for a, b in self.course_units(i):
            x0 = a + (tj / 2 if a > 1e-9 else 0.0)
            x1 = b - (tj / 2 if b < Lt - 1e-9 else 0.0)
            n = max(1, int(round(self.holes_per_unit * (b - a) / self.unit_length)))
            Lc = (x1 - x0 - 2 * tt - (n - 1) * tw) / n
            for k in range(n):
                cx0 = x0 + tt + k * (Lc + tw)
                res.append((cx0, cx0 + Lc, te, W - te, True))
            if self.wall_cells > 0:
                ms = max(1, int(round(self.wall_cells * (b - a) / self.unit_length)))
                th = self.min_shell
                Ls = (x1 - x0 - (ms + 1) * th) / ms
                wy = te - 2 * th
                for k in range(ms):
                    sx0 = x0 + th + k * (Ls + th)
                    res.append((sx0, sx0 + Ls, th, th + wy, False))
                    res.append((sx0, sx0 + Ls, W - th - wy, W - th, False))
        return res

    def grouted_columns(self, t_joint=None):
        """x bands (x0, x1) of the grouted holes or main cells, in course 0."""
        if self.grouting == "none" or not self.has_holes:
            return []
        if self.rectangular:
            cols = sorted((c[0], c[1]) for c in self.cells(0, t_joint) if c[4])
        else:
            R = self.hole_diameter / 2
            cols = [(x - R, x + R) for x in self.hole_centres]
        if self.grouting == "alternate":
            cols = cols[::2]
        return cols

    def hole_area(self, i=0, grout=False, t_joint=None):
        """Hole area (m²) in course i: grout=False counts all holes (unit
        area = gross - holes); grout=True counts only the ungrouted holes."""
        if self.horizontal:          # minimum horizontal section: through the rows of holes
            return self.holes_width * self.hole_width_h * self.length
        if not self.has_holes:
            return 0.0
        cols = self.grouted_columns(t_joint) if grout else []

        def grouted(x0, x1):
            return any(min(x1, b) - max(x0, a) > 0.5 * (x1 - x0) for a, b in cols)
        if self.rectangular:
            return sum((c[1] - c[0]) * (c[3] - c[2]) for c in self.cells(i, t_joint)
                       if not (c[4] and grouted(c[0], c[1])))
        R = self.hole_diameter / 2
        return sum(np.pi * R * R for x in self.hole_centres if not grouted(x - R, x + R))

    @property
    def length(self):
        return self.units_per_course * self.unit_length

    @property
    def hole_spacing(self):
        return self.unit_length / self.holes_per_unit

    @property
    def hole_centres(self):
        if self.hole_shape != "circular" or self.hole_diameter <= 0 or self.holes_per_unit <= 0:
            return ()
        s = self.hole_spacing
        n = int(round(self.length / s))
        return tuple(s / 2 + k * s for k in range(n))

    @property
    def joints_per_course(self) -> Tuple[Tuple[float, ...], ...]:
        if self.head_joints is not None:
            return tuple(tuple(j) if isinstance(j, (tuple, list)) else (float(j),) for j in self.head_joints)
        Lb, Lt = self.unit_length, self.length
        res = []
        for i in range(self.n_courses):
            off = Lb / 2 if (self.bond == "running_bond" and i % 2 == 1) else 0.0
            xs = [round(off + k * Lb, 9) for k in range(int(Lt / Lb) + 2)]
            res.append(tuple(x for x in xs if 1e-6 < x < Lt - 1e-6))
        return tuple(res)

    @property
    def joint_positions(self) -> Tuple[float, ...]:
        return tuple(sorted({x for f in self.joints_per_course for x in f}))

    @property
    def prism_height(self):
        return self.n_courses * self.course_height

    def validate(self):
        """List of problems of the geometry (empty if everything is right)."""
        p = []
        if self.unit_length <= 0 or self.width <= 0 or self.course_height <= 0 or self.n_courses < 1:
            p.append("Dimensions and number of courses must be positive.")
        if abs(self.units_per_course * 2 - round(self.units_per_course * 2)) > 1e-9 or self.units_per_course < 0.5:
            p.append("Units per course must be a multiple of 0.5.")
        if self.bond not in ("running_bond", "stack_bond"):
            p.append("Bond must be 'running_bond' or 'stack_bond'.")
        if self.joint_type not in ("dry", "mortar"):
            p.append("Joint type must be 'dry' or 'mortar'.")
        if self.bedding not in ("total", "lateral"):
            p.append("Bedding must be 'total' or 'lateral'.")
        if self.grouting not in ("none", "all", "alternate"):
            p.append("Grouting must be 'none', 'all' or 'alternate'.")
        if self.grouting != "none" and not self.has_holes:
            p.append("Only vertical holes receive grout (solid unit or unit with horizontal holes).")
        if self.hole_shape not in ("circular", "rectangular", "horizontal"):
            p.append("Hole shape must be 'circular', 'rectangular' or 'horizontal'.")
        if self.hole_shape == "horizontal":
            if self.holes_width < 1 or self.holes_height < 1:
                p.append("With horizontal holes, use at least 1 hole across the thickness and 1 over the height.")
            elif self.hole_width_h < 0.003:
                p.append("Horizontal holes too narrow (less than 3 mm): check face shells and webs.")
            elif self.hole_height_h() < 0.003:
                p.append("Horizontal holes too low (less than 3 mm): check bed shells, webs and joint.")
        elif self.rectangular:
            n = self.holes_per_unit
            mid = abs(self.units_per_course % 1) > 1e-9
            if (self.bond == "running_bond" or mid) and n % 2 and self.grouting != "none":
                p.append("With grout and a half unit or running bond, the number of cells per unit must be even "
                         "(the cells must line up from one course to the next).")
            if 2 * self.face_shell >= 0.9 * self.width:
                p.append("The face shells take up almost the whole width of the unit.")
            if self.wall_cells > 0 and self.face_shell <= 2 * self.min_shell + 1e-4:
                p.append("With hollow shells, the face shell must be more than twice the minimum wall.")
            for i in range(self.n_courses):
                if any(c[1] - c[0] < 0.004 for c in self.cells(i)):
                    p.append("Cells too short (less than 4 mm): check webs, end shells and cells.")
                    break
        elif self.hole_diameter > 0 and self.hole_shape == "circular":
            n = self.holes_per_unit
            mid = abs(self.units_per_course % 1) > 1e-9
            if n < 1:
                p.append("With holes, use at least 1 hole per unit.")
            elif (self.bond == "running_bond" or mid) and n % 2:
                p.append("With a half unit or running bond, the number of holes per unit must be even "
                         "(the holes must line up between courses).")
            R, W = self.hole_diameter / 2, self.width
            if R >= 0.9 * W / 2:
                p.append("The hole is too wide for the unit width (wall less than 10% of the half width).")
            elif n >= 1 and R >= 0.9 * self.hole_spacing / 2:
                p.append("Holes too close to each other or to the joints; reduce the diameter or the number of holes.")
        return p


@dataclass
class Test:
    """UFV protocol (Jonathan's dissertation): furnace at 5 °C/min to the
    target temperature, 30 min dwell and compression still inside the furnace."""
    T_initial: float = 30.0
    rate_C_min: float = 5.0
    T_target: float = 600.0
    dwell_s: float = 1800.0
    h_convection: float = 10.0      # W/m².K (electric furnace, still air) — test 10 to 25
    emissivity: float = 0.8      # ceramic surface — test 0.8 to 0.9
    load_kgf: float = 9980.9      # peak load of the test

    @property
    def load_N(self):
        return self.load_kgf * 9.80665


@dataclass
class CeramicDamage:
    """Isotropic damage model of the ceramic (non-linear analysis).

    Rankine criterion in tension and von Mises in compression (effective
    stresses), exponential softening and regularisation by the crack
    band (Bažant and Oh, 1983). The fracture energies come from the ductility
    indices d_u = G/f used for clay units in micro-models."""
    ft: float = 0.66e6      # Pa: direct tension of pressed clay brick, net area (Gooch et al., 2023)
    fc: float = 5.76e6      # Pa: mean fb 4.13 MPa on gross area (BPQ) / 0.717 -> net area
    du_t: float = 0.018e-3  # m: G_f/ft of the clay unit (van der Pluijm 1992; Lourenço 1996; via Dauda et al., 2020)
    du_c: float = 0.33e-3   # m: G_c/fc of the clay unit (Angelillo et al. 2014; via Dauda et al., 2020)
    d_max: float = 0.9999   # maximum damage (minimum residual stiffness, for numerical stability)
    crack_closure: bool = True  # a compressed tensile crack transmits compression again (Mazars, 1984)
    table_ft_factor: Table = None  # factors on ft with temperature: [(T, factor), ...]
    table_fc_factor: Table = None

    @property
    def G_f(self):
        return self.du_t * self.ft   # N/m

    @property
    def G_c(self):
        return self.du_c * self.fc

    def ft_T(self, T):
        return self.ft * interpolate_(T, self.table_ft_factor, 1.0)

    def fc_T(self, T):
        return self.fc * interpolate_(T, self.table_fc_factor, 1.0)


def mortar_damage(fc: float = 4.2e6, ft: Optional[float] = None) -> "CeramicDamage":
    """Damage parameters of the mortar. ft = 0.3 fc^(2/3) (MPa; Eurocode 2
    formula used by Dauda, Iuorio and Lourenço, 2020). Ductility indices
    of the mortar: G_f/ft = 0.065 mm and G_c/fc = 1.6 mm (same source,
    after Lourenço, 1996)."""
    if ft is None:
        ft = 0.3 * (fc / 1e6) ** (2 / 3) * 1e6
    return CeramicDamage(ft=ft, fc=fc, du_t=0.065e-3, du_c=1.6e-3)


def grout_damage(fc: float = 17.7e6, ft: Optional[float] = None) -> "CeramicDamage":
    """Damage parameters of the grout: ft = 0.3 fc^(2/3) (MPa) and the same ductility
    indices as the mortar (G_f/ft = 0.065 mm, G_c/fc = 1.6 mm). For
    fc = 17.7 MPa, G_f = 0.13 N/mm, close to G_f = 73 fcm^0.18 N/m of the Model
    Code 2010 (0.12 N/mm)."""
    if ft is None:
        ft = 0.3 * (fc / 1e6) ** (2 / 3) * 1e6
    return CeramicDamage(ft=ft, fc=fc, du_t=0.065e-3, du_c=1.6e-3)


@dataclass
class FrictionJoint:
    """Dry joint as a thin non-linear layer: transmits no tension and
    transmits shear up to the Coulomb friction limit without cohesion."""
    mu: float = 0.60             # friction in sawn dry joints: 0.61-0.64 (Lourenço, Ramos and Vasconcelos, 2004)
    open_factor: float = 1e-4   # residual normal stiffness of the open joint (fraction of E_joint)
    plane_factor: float = 1e-3    # in-plane stiffness of the layer (for stability only)
    G_min_factor: float = 1e-3    # minimum shear stiffness when sliding
    eps_tol: float = 1e-6        # hysteresis band of the contact (normal strain) to avoid oscillation
    pre_load_N: float = 2e3     # seating load during heating (whole prism)
    # Irregularity of the bearing faces (bed joints). Without load, only
    # part of the nominal area touches: 23% in untreated units (Ngapeya and
    # Waldmann, 2020) and ~15% rising to ~95% near failure (Zahra et al.,
    # 2016). The rest has gaps that close with the load and concentrate stresses
    # (Ngapeya, Waldmann and Scholzen, 2018). The gaps are a Gaussian field
    # with standard deviation gap_std and correlation length corr_length.
    gap_std: float = 0.028e-3  # m; calibrated from the prism modulus (chord E 5-30% = 1,100 MPa as the
                                    # mean of 3 seeds; script 08 gives 0.029 mm; mesh h = 12, dz = 20 mm)
                                    # 0 = perfect contact
    E_contact: Optional[float] = 3907e6  # Pa: E of the joints and capping where there is contact; equal to that of the
                                         # ceramic (the joint flexibility comes from the gaps). None = Materials.joint.E
    corr_length: float = 0.04  # m (assumption: about half the hole diameter)
    initial_contact: float = 0.20   # fraction of the area in contact without load
    seed: int = 1                # seed of the random field


# Effective compressive strength of the ceramic that reproduces the peak load
# of the prism at 30 °C with the calibrated gaps (script 08, mesh h = 12 mm,
# dz = 20 mm). It is 0.43 times the value derived from the unit test. Likely
# explanations: the test of the low, wide unit (60 x 140 mm), capped and
# confined by platen friction, overestimates the material strength; the
# asperities smaller than the mesh are not in the model; the 30 °C load still
# needs confirmation. Use only as a calibration (option --calibrated of script 07).
CALIBRATION_30C = {"fc": 2.48e6}

# Peak loads of the UFV prisms (kgf); the 30 °C one needs confirmation.
LOADS_UFV_KGF = {30: 10066.6, 600: 9980.9, 800: 10087.3, 1000: 11662.1}

# Secant modulus of dry 5-course prisms (BPQ, n = 4 per group), MPa
E_DRY_PRISM_MPA = {"5ID (whole units)": 1216.0, "5HD (alternating half units)": 1001.0}

SIGMA = 5.670374419e-8  # Stefan-Boltzmann
