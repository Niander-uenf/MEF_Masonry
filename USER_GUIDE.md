# MEF_Masonry — user guide

*A short overview is in [README.md](README.md); the paper and [REPRODUCE.md](REPRODUCE.md) describe the studies.*

Thermal and thermo-mechanical finite elements for block masonry at high
temperature, without ANSYS. It only needs Python with NumPy, SciPy and
Matplotlib; the 3D results are written as VTK for viewing in ParaView (free),
and the program has a browser interface with its own 3D post-processor.

The package was born for the UFV test of the BCEPQ unit (a prism of 5 courses
of one and a half units, dry joints, furnace at 5 °C/min, a 30 min dwell and
loading at the end of the dwell). It now handles masonry of several units and,
in version 2.0, covers:

| | What it takes |
|---|---|
| Units | Solid; vertical round holes (BCEPQ); vertical rectangular cells (concrete NBR 6136, structural clay NBR 15270, including hollow shells), with or without grout; horizontal holes (infill) |
| Joints | Dry (contact, gaps and friction) or mortar (bonded layer with its own damage) |
| Specimens | Single unit, prism, small wall and full-scale wall |
| Analyses | 3D thermal; 3D thermo-elastic; 3D non-linear to failure; 2D section (horizontal, or vertical for horizontal holes); wall with second-order effects; transient creep in the 3D non-linear model and in the wall; global sensitivity of the wall (Morris) |
| Tests | Ambient, hot, residual, transient (constant load and heating to failure) and the ISO 834 standard fire |
| Library | 11 units, 3 materials, 4 mortars and 2 grouts, each with its source |
| Verification | 8 scripts, 240 cases in the fast run (01 and 06 run with the three elements), with one command: `examples/00_verify_all.py` |

The sections below describe each stage, from the most recent (version 2.0, in
English) to the oldest (Wall mode, 0.5).

## Version 2.0.0: the program in English

Up to version 1.0.1, the identifiers of the code, the interface, the messages,
the fields of the project files and this guide were in Portuguese. From 2.0.0
everything is in English. The solver did not change: the 240 verification
cases and the studies of the paper give the same results as in 1.0.1.

What changed name:

| Before (1.x) | Now (2.0) |
|---|---|
| Folders `exemplos/`, `exemplos/artigo/`, `resultados/`, `dados/` | `examples/`, `examples/paper/`, `results/`, `data/` |
| Modules `parametros`, `malha`, `elementos`, `elementos3d`, `termico`, `mecanico`, `naolinear`, `parede`, `dados_parede`, `sensibilidade`, `biblioteca`, `projeto`, `saida`, `curvas`, `estudos` | `parameters`, `mesh`, `elements`, `elements3d`, `thermal`, `mechanical`, `nonlinear`, `wall`, `wall_data`, `sensitivity`, `library`, `project`, `output`, `curves`, `studies` |
| `projeto_padrao()`, `aplicar()`, `executar()`, `ESQUEMA`, `PREDEFINIDOS`, `biblioteca.BLOCOS` | `default_project()`, `apply()`, `run()`, `SCHEMA`, `PRESETS`, `library.UNITS` |
| Groups of the project `geometria`, `ceramica`, `modelo`, `ensaio`, `malha`, `sensibilidade` | `geometry`, `ceramic`, `model`, `test`, `mesh`, `sensitivity` |
| Scripts `00_verificar_tudo.py`, `01_verificacao.py`, ... | `00_verify_all.py`, `01_verification.py`, ... (list under "Order of use") |
| Options `--completo`, `--rapido`, `--calibrado`, `--temperaturas`, `--elemento`, ... | `--full`, `--fast`, `--calibrated`, `--temperatures`, `--element`, ... |
| `LEIA-ME.md`, `Iniciar_MEF_Masonry.bat`, `Empacotar_Windows.bat`, `empacotar.py`, `iniciar.sh` | `USER_GUIDE.md`, `Start_MEF_Masonry.bat`, `Build_Windows.bat`, `build_exe.py`, `start.sh` |
| `projeto.json`, `resultados.json` in each run folder | `project.json`, `results.json` |

Numbers are written with a decimal point. In the interface, a decimal comma
is still accepted when typing a value.

What keeps working (`mef_masonry/legacy.py`):

- **Saved projects** (`.json`) of versions 1.x open as before: the group and
  field names and the option values are converted when the file is opened in
  the interface, and by `project.complete`, `project.apply` and `project.run`.
- **Your own scripts** written for 1.x: `from mef_masonry import projeto as pj`
  returns the module `project`, and `pj.projeto_padrao()`, `pj.aplicar()`,
  `pj.executar()` and the other old names resolve to the new functions. The
  library accepts the old selectors and keys
  (`library.apply(p, "bloco", "concreto_m15_39")`). The results come back with
  the names of version 2 (`res["summary"]`, not `res["resumo"]`).
- **The last project open in the browser** comes back the first time 2.0 opens.
- **History of analyses:** version 1.x stored the runs in `resultados_app/`;
  2.0 stores them in `results_app/`. Rename the old folder to `results_app`
  (or start the program with `python app.py --folder resultados_app`) and the
  old runs appear in Previous analyses, where their inputs can be reused; to see
  their results in the interface, run them again.

## Version 1.0.1: studies of the review of the paper

The solver did not change. Three scripts were added to `examples/paper/`, with
their results in `results/paper/`, and Figures 2 and 3 now state that they
are illustrative:

- `oliveira_explicit_section.py`: the infinite parallel-plate model of the
  cavities overestimates the radiative exchange in finite cells. With the view
  factor between the shells (effective emissivity ε′, with
  ε′/(2 − ε′) = F ε/(2 − ε)), the insulation of the wall of Oliveira et al.
  (2021) goes from 42.5 min to 58.5 min (2D factor) and 68.5 min (3D factor),
  against 65 to 82 min in the tests. To use the correction in a project, enter
  ε′ as the emissivity of the hole walls.
- `morris_checks.py`: the screening of the BCEPQ wall without the effects that
  involve runs without failure, and the nominal wall with k_tr = 2.35
  (224 min, against 211 min with k_tr = 0).
- `morris_narrow_alpha.py`: the same screening with α from 5 to 7 ×10⁻⁶/°C.

## Version 1.0: the program is renamed MEF_Masonry

The program was born for the BCEPQ, but it now computes solid units,
structural and infill clay units and concrete blocks with and without grout,
with dry or mortar joints. In version 1.0 it was therefore renamed
MEF_Masonry. The BCEPQ remains one of the units of the library and the default
project.

What changed name in 1.0:

| Before (up to 0.10) | From 1.0 |
|---|---|
| Folder `mef_bcepq` | Folder `MEF_Masonry` |
| Python package `mef_bcepq` | Package `mef_masonry` |
| `Iniciar_MEF_BCEPQ.bat` | `Iniciar_MEF_Masonry.bat` (`Start_MEF_Masonry.bat` in 2.0) |
| Executable `dist/MEF_BCEPQ/MEF_BCEPQ.exe` | `dist/MEF_Masonry/MEF_Masonry.exe` |

Scripts that do `import mef_bcepq` keep running: the old name is an alias of
`mef_masonry`, which warns about the change and will be removed in a future
version. Replace `mef_bcepq` with `mef_masonry` in your scripts.

## Installation

```
python -m pip install -r requirements.txt
```

Python 3.9 or newer. Optional: `pip install pypardiso` (Windows/Linux) speeds
up the 3D mechanical solver and reduces the memory.

## Program with interface

```
python app.py          (or double-click Start_MEF_Masonry.bat on Windows)
```

On Windows, `Start_MEF_Masonry.bat` looks for Python (`py` or `python`),
installs NumPy, SciPy and Matplotlib the first time and opens the program; if
something is missing, the window stays open with the explanation. On Linux or
macOS, `./start.sh` does the same.

It opens a local interface in the browser (nothing goes to the internet) with:

- **Parametric geometry:** dimensions of the unit, diameter and number of
  holes per unit, units per course (multiples of 0.5), number of courses, bond
  (running or stack), joint and capping. The mesh can be previewed.
- **Properties** of the ceramic, with optional tables against temperature
  (k, factor on E, thermal strain, ft and fc).
- **Models and elements:**
  - 3D thermal, 3D linear thermo-elastic, 3D non-linear to failure (Mazars
    damage with crack closure; joints with gaps and friction), 2D section or
    wall (vertical strip with second-order effects);
  - 3D elements W6 (linear wedge, the default of the paper), H8 (linear
    hexahedron on the quadrilaterals of the O-grid mesh) and W15 (quadratic
    wedge); T3 in the 2D section. See "Elements" below.
- **Tests:** compression at 30 °C, hot (ramp, dwell and loading), residual
  (heating, cooling and loading) and the ISO 834 standard fire on one face,
  under constant load. In the non-linear ISO 834 case, failure follows the
  axial deformation limits of EN 1363-1 (C = h/100 and dC/dt = 3h/1000 per
  minute); this case is experimental, because the load is applied as a
  pressure on the top, without a rigid platen.
- **Run** with a live log and cancellation.
- **Results:** indicators, charts, tables, CSV and VTK.
- **3D post-processor** (WebGL, works offline):
  - colour maps of temperature, displacements, principal and vertical
    stresses, tensile, compressive and effective damage, and joint gaps;
  - states of the heating and of the loading;
  - deformed shape with a scale factor;
  - cuts in x, y and z;
  - mirroring of the symmetry;
  - materials that can be hidden, value at a point (double click) and PNG image.
- **Projects** saved and opened as .json, and ready-made projects (the UFV
  tests).

The analyses are stored in `results_app/`, one folder per run, with
`project.json`, `results.json`, figures, CSV, VTK and `visual.json`.

**Executable (.exe):** on the target computer, run `python build_exe.py` (or
double-click `Build_Windows.bat`). PyInstaller builds
`dist/MEF_Masonry/MEF_Masonry.exe`, which can be copied to other machines
without Python. The build was not tested on Windows from this environment; if
a module is missing, add a `--hidden-import` to `build_exe.py`.

**From code:** `mef_masonry.project.run(p, folder)` runs the same project
without the interface (see `default_project()` and `SCHEMA`):

```python
from mef_masonry import project as pj, library as lib

p = pj.default_project()                       # all inputs at their defaults
p = lib.apply(p, "unit", "concrete_m15_39")    # unit from the library
p = pj.apply(p, {"model": {"analysis": "section_2d", "element": "T3"},
                 "test": {"kind": "iso834", "iso_duration": 120}})
res = pj.run(p, "results/example")             # CSV, PNG, VTK and JSON files
for item in res["summary"]:
    print(item)
```

## Global sensitivity of the wall (version 0.10)

Before calibrating the model and planning tests, it is worth knowing which
inputs change the response most. Version 0.10 brings a global sensitivity
analysis of the wall under ISO 834 by the **Morris method** (elementary
effects): the program runs the wall many times, varying the uncertain inputs
within ranges, and ranks which weigh most on the fire resistance (R) and
insulation (I) times and on the cold capacity.

**How to use it:** in Models and elements, choose the analysis "Wall
sensitivity under ISO 834 (Morris)". Page 6, Sensitivity, appears with the
parameters and the ranges. The rest of the project (geometry, materials, load,
supports) gives the nominal values around which everything varies. The
ready-made project "BCEPQ wall 2.94 × 3.00 m — sensitivity under ISO 834" is
already set up.

**The method:**

- Each parameter has a range, mapped to [0, 1] and divided into a grid of p
  levels (4 by default).
- A trajectory starts from a random point of the grid and changes one
  parameter at a time, in random order, by Δ = p / [2 (p − 1)] (2/3 of the
  range with 4 levels). With k parameters, each trajectory costs k + 1 runs.
- The elementary effect of parameter i is EE = [y(x + Δ eᵢ) − y(x)] / Δ, on
  the scale of the whole range. With r trajectories (10 by default), each
  parameter has r effects.
- μ* = mean of |EE| (Campolongo et al., 2007) is how much the output changes,
  on average, when the parameter spans its range: it is the measure that ranks
  the parameters. μ, the mean with sign, says which way the output moves. σ,
  the standard deviation, measures non-linearity and interactions.
- The 95% interval of μ* comes from a bootstrap over the trajectories.
- The ratio σ/μ* classifies the effect (Garcia Sanchez et al., 2014): almost
  linear (< 0.1), monotonic (< 0.5), almost monotonic (< 1) or non-linear and
  with interactions (≥ 1).

| Parameter | Default range | Why |
|---|---|---|
| Thermal conductivity of the unit | ±20% | Uncertainty of the measurement and variation between batches; also applies to the k(T) table, if any |
| Specific heat of the unit | ±15% | Includes the uncertainty of the moisture peak; also applies to the c(T) table |
| Emissivity of the hole walls | 0.7 to 0.95 | Ceramic: about 0.85 to 0.95; the low limit covers smoother surfaces |
| Convection on the exposed face | 15 to 35 W/m²·K | EN 1991-1-2 uses 25 W/m²·K on the face exposed to the standard fire |
| Surface emissivity | 0.6 to 0.95 | EN 1991-1-2 recommends 0.8 for the surface of the member |
| Thermal expansion of the unit | 4 to 8 ×10⁻⁶/°C | Red ceramic: about 4 to 8 ×10⁻⁶/°C; replaces α and scales the thermal strain table |
| Transient creep of the unit (k_tr) | 0 to 2.35 | From no creep (0) to the upper limit of concrete (2.35; Anderberg and Thelandersson, 1976) |
| Masonry strength (factor on the law) | ±15% | Scatter of the prism strength; multiplies the stresses of the wall law |
| Masonry stiffness (factor on the law) | ±25% | Scatter of the prism modulus; multiplies the modulus of the wall law |
| Initial imperfection (H divided by) | H/1000 to H/250 | Initial deflection H/n at mid-height; the amplitude varies linearly between the two limits |
| Load eccentricity (top and base) | −10 to 10 mm | Applied equally at the top and the base; positive towards the face exposed to fire |

**Cost:** r (k + 1) + 1 runs, with the nominal run. With the 11 parameters and
r = 10, that is 121 runs. With the fast mesh of the screening (10 mm section
and 16 beam elements), each run takes 10 to 30 s on one core (transient creep
is the most expensive); the recommended mesh is about 6 times slower. The runs
are executed in parallel (cores minus one, by default).

**Output:** in the summary, the nominal values and the parameter with the
largest effect on each time; a table of effects per output (μ*, interval, μ,
σ, σ/μ* and the reading); the figures of μ* and of the μ*–σ plane; and the
files `sensitivity_R.csv`, `sensitivity_I.csv`, `sensitivity_Pmax.csv` and
`sensitivity_runs.csv` (the complete plan, with the physical values and the
outputs of each run).

**BCEPQ wall** (ready-made project: 2.94 × 3.00 m, 150 kN, ISO 834 for 300 min,
11 parameters, r = 10, 121 runs with the fast mesh, 11.5 min on 2 cores; the
outputs are in `results/sensitivity_wall_bcepq/`). Nominal values: R = 211 min,
I = 112 min and a cold capacity of 554 kN.

| Fire resistance time R | μ* (min) | 95% CI | μ (min) | σ/μ* |
|---|---|---|---|---|
| Thermal expansion α (4 to 8 ×10⁻⁶/°C) | 253 | 181 to 328 | −253 | 0.50 |
| Load eccentricity (−10 to +10 mm) | 157 | 74 to 238 | +157 | 0.91 |
| Initial imperfection (H/1000 to H/250) | 55 | 21 to 102 | −55 | 1.28 |
| k_tr of the unit (0 to 2.35) | 24 | 4 to 51 | +24 | 1.72 |
| Masonry stiffness (±25%) | 20 | 4 to 40 | +20 | 1.59 |
| Others (k, ε of the face, c, ε of the holes, strength, h) | 0.3 to 10 | | | |

| Insulation time I | μ* (min) | 95% CI | μ (min) | σ/μ* |
|---|---|---|---|---|
| Conductivity k (±20%) | 42 | 39 to 47 | −42 | 0.15 |
| Specific heat c (±15%) | 27 | 26 to 29 | +27 | 0.12 |
| Emissivity of the faces (0.6 to 0.95) | 8.0 | 6.0 to 10.3 | +8.0 | 0.46 |
| Emissivity of the holes (0.7 to 0.95) | 7.6 | 6.8 to 8.4 | −7.6 | 0.18 |
| Convection on the exposed face (15 to 35 W/m²·K) | 0.8 | 0.6 to 1.1 | −0.8 | 0.52 |

On the cold capacity, the masonry strength and the eccentricity weigh most
(148 kN each), followed by the imperfection (87 kN) and the stiffness
(25 kN); the thermal inputs give zero effect, as they should, and the
mechanical inputs give zero on the insulation.

What this says for the test programme:

- **The α of the unit is the first measurement to make** (dilatometry up to
  1000 °C). The α = 6 ×10⁻⁶/°C of the BCEPQ is provisional and sits right on a
  threshold: no run with α ≥ 6.7 went beyond 300 min (medians of 90 and
  37 min), and about 70% of the runs with α ≤ 5.3 did. The mechanism is the
  thermal bowing towards the fire added to the second-order effects.
- **The eccentricity and the imperfection of the test must be measured.** A
  load shifted towards the exposed face offsets the bowing (median R above
  300 min with +10 mm, against 64 min with −10 mm). The 12 runs that failed
  before 30 min had α ≥ 6.7 and the load shifted towards the cold face.
- **k_tr and the stiffness come second**, with a non-linear effect and
  interactions (σ/μ* > 1): the transient test is worthwhile, but after α.
- **The masonry strength hardly changes R** (2 min), although it governs the
  cold capacity: the wall fails by instability, not by crushing.
- **On the insulation, k(T) and c(T) govern**; the emissivities weigh little
  and the convection nothing. In the calibration with the thermocouples, the
  focus is k and c (with the moisture peak); the boundary conditions weigh
  less.

In 59 of the 121 runs the wall did not fail within 300 min, and R entered as
300 min (a lower bound): the effects on R, above all that of α, are larger
than those of the table. The numbers hold for the ranges above, for the pinned
supports of the project and for the law of the BCEPQ prism.

**Verification:** `examples/15_verification_sensitivity.py`, 31 cases (32 with
`--full`, plus 4 informative results), all within tolerance:

- Morris plan (one parameter per step, step Δ, points on the grid);
- linear function (μ = aᵢ, σ = 0, exact), additive quadratic function (μ*
  equal to the separate calculation over the level pairs; with r = 2000,
  within 0.4% of the expected value), interaction x₁x₂ (σ > 0) and Sobol'
  g-function (known order of the parameters);
- a failed run is left out only in its two neighbouring effects; the same
  seed repeats the plan; the 95% interval contains μ*;
- physical ranges, k(T) and ε_th(T) tables scaled together, original project
  untouched, factors of the law (peak strength and modulus, exact);
- wall: the run of the study gives the same R, I and Pmax as the Wall
  analysis; in parallel and sequentially, identical results (`--full`).

With the recommended mesh (5 mm and 32 elements), the nominal R of the BCEPQ
drops from 210.9 to 205.1 min and I goes from 112.0 to 112.5 min; each run
takes about 6.6 times longer. In the first 3 trajectories of the study
(36 runs), also run with the recommended mesh, R changed by 0.9 min on average
(6.6 min at most) and I by 0.5 min at most; the first five parameters came out
in the same order with both meshes, with μ* within 5%.

## Transient creep and transient test (version 0.9)

In a hot test (heating and then loading), the material arrives hot and
unloaded. In a fire, the wall is loaded while it heats, and the material
shortens beyond its free expansion and elastic strain: this is **transient
creep**, or load-induced thermal strain. It is irreversible and only appears
on first heating under compression.

**The model** (page Models and elements, Transient creep):

- Anderberg and Thelandersson (1976): ε_tr = φ(T) σ / fc20, with
  φ(T) = k_tr ε_th(T). σ is the compressive stress and fc20 the compressive
  strength of the material at 20 °C. Instead of k_tr, the unit accepts a
  tabulated φ(T), taken from the transient test.
- Incremental: Δε_tr = Δφ σ / fc20 only when φ exceeds the largest value
  reached so far. It does not recover on cooling nor add up on reheating below
  the maximum temperature.
- In 3D, Δε_tr = Δφ / fc20 · [(1 + ν) σ⁻ − ν tr(σ⁻) I], with σ⁻ the
  compressive part of the stress: under uniaxial stress, creep shortens in the
  direction of the load and expands transversely with the material's own
  Poisson ratio. Each interval between thermal fields is crossed in substeps
  with a predictor and a corrector (trapezoidal rule), and the increment of
  ε_tr per substep stays below 2 × 10⁻⁴.
- In the wall, unit and bed joint are in series over the height of the course:
  c(T) = (1 − f) φ_b / fc_b + f φ_m / fc_m, with f the share of the joint in the
  height of the course. The grout, continuous, has its own coefficient.

**Default values** (creep is switched on):

| Material | k_tr | Why |
|---|---|---|
| Cement mortar | 2.0 | Prakash, Azenha, Pereira and Lourenço (2020), within the range of 1.8 to 2.5 of concrete |
| Clay unit (BCEPQ and others) | 0 | Fired clay has no cement paste and there is no published value; Prakash et al. (2020) also neglect creep in the clay unit |
| Concrete block and grout with the EN 1992-1-2 curves | 0 | The curves of the standard already include creep (Gernay and Franssen, 2012); adding more would count it twice |

With these defaults, the dry-jointed BCEPQ does not change: the results of the
earlier versions still hold. Creep appears in mortared masonry and wherever
you give k_tr or φ(T).

**Transient test** (Test = transient, only in the 3D non-linear model): the
constant load is applied cold and the furnace heats along the ramp and the
dwell until failure. Failure is the loss of equilibrium or the axial
deformation limits of EN 1363-1. The output is the strain × temperature curve,
with the free expansion and the creep part, and the temperatures at failure.
It is the test that measures φ(T): with the total strain ε, the free strain
ε_th and the elastic strain σ/E(T),

    φ(T) = (ε − ε_th − σ/E(T)) fc20 / σ

(shortening positive). The results of two new ready-made projects:

| Ready-made project | Without creep | With creep |
|---|---|---|
| Solid brick with 1:1:6 mortar, 46 kN (40%) up to 800 °C | 69% of the brick with damage > 0.5 at the end | 17% (creep of the mortar relieves the thermal stresses between brick and joint) |
| BCEPQ, 49 kN (50%) up to 800 °C, hypothetical k_tr of 2 in the unit | 54% of the unit with damage > 0.5 at the end | 44%, with 4.0‰ of mean creep; no failure in either (the BCEPQ does not lose strength when hot) |

**Verification:** `examples/14_verification_creep.py`, 19 cases (20 with
`--full`), all within tolerance:

- strain under constant stress and uniform heating, axial and lateral, with
  the three elements (exact);
- restrained heating against the analytical solution of the relaxation,
  σ = fc/k [exp(−E k ε_th / fc) − 1] (0.6% with the default step and 0.05%
  with 10⁻⁵);
- irreversibility on cooling and reheating, absence of creep under tension
  and tabulated φ(T) (exact);
- wall: shortening under constant load (exact) and series coefficient of unit
  and joint;
- BCEPQ wall: with k_tr = 0, switching creep on changes nothing.

**Sensitivity of the BCEPQ wall** (ready-made ISO 834 project under 150 kN, 32
elements, 240 min), with hypothetical k_tr in the unit:

| k_tr of the unit | 0 | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|
| Fire resistance time (min) | 205.1 | 179.4 | 167.6 | 176.2 | 200.9 |

The effect is not monotonic. Creep shortens the compressed fibres that heat
up: it takes stiffness from the section and, at the same time, reduces the
thermal deflection towards the fire. With k_tr near 1, the first effect
dominates and the time drops by 18%; with k_tr = 2, the two almost cancel out.
The deflection at failure stays between 43 and 46 mm in all cases, and the
failure mode is the same. As fired clay has no published k_tr, the value of the
BCEPQ must come from a transient test (loaded and heated prism), and the
program already reproduces that test for comparison.

## Horizontal holes and consolidation (version 0.8)

**Units with horizontal holes** (infill), with Perforation = Horizontal holes:

- **Geometry:** holes across the thickness and the height of the unit (the
  8-hole unit is 2 × 4), thickness of the face shells, of the webs (between the
  holes, in both directions) and of the bed shells (top and bottom). The holes
  run through the unit along its length; the mortar head joints close them.
- **3D:** the plan follows the edges of the columns of holes, and each course
  is divided into sublayers at the edges of the rows of holes; the holes open in
  the rows, outside the head joints. Radiation in the holes uses the width of
  the hole in the direction of the heat flow.
- **Vertical 2D section:** for these units, the 2D section is the vertical
  section of one course (thickness × height, with half a joint at the top and
  bottom). As the holes run through the unit, this section holds along the
  whole length of the wall. The top and bottom edges are periodic cuts (no heat
  exchange), and the generalised plane strain lies along the holes.
- **Wall:** the thermal analysis comes from the same vertical section, and the
  fibres are the walls that run down continuously through the course (faces
  and vertical webs, with the joints on top of them). The horizontal webs enter
  only the thermal analysis, because they carry no vertical compression.
- **Library:** 8-hole units of 90 × 190 × 190 mm (SINAPI 7271) and 6-hole units
  of 90 × 140 × 240 mm (walls measured by Guedes et al., 2019), with 1:2:9
  mortar. E, fc and ft are estimates from the minimum strength of 1.5 MPa of
  NBR 15270 (see the note of each item).
- **Ready-made projects:** ISO 834 insulation of the vertical section of the
  8-hole unit (24 min without rendering, with the properties of the library)
  and a non-linear prism at 20 °C.
- **Verification:** `examples/13_verification_horizontal_holes.py`, 12 cases
  (14 with `--full`), all within tolerance:
  - volume of the holes in the 3D prism, area of the minimum horizontal section
    and areas of the vertical section (exact);
  - thermal: the vertical 2D section reproduces the 3D model extruded along the
    holes (difference of 0.3 °C at 60 min of ISO 834);
  - fibres of the wall: area per metre equal to the width minus the holes, and
    centroid on the mid-plane;
  - non-linear 8-hole prism: the peak is that of the continuous walls (fc times
    their area, within 5%).

**Consolidation:**

- `examples/00_verify_all.py` runs all the verifications (01, 06 and 10 to 15,
  with the three elements where applicable) and writes
  `results/verification_summary.csv`. Without options it takes about 7 min;
  with `--full`, it includes the long cases (convergence of the wall, prisms of
  Camacho et al. and others) and takes about 35 min.
- The 2D section, the figures, the VTK and the post-processor show the
  vertical section upright (thickness horizontal, height vertical).

The results of the BCEPQ did not change: the prism at 30 °C still gives
98.68 kN and 1212 MPa.

## Rectangular cells and grout (version 0.7)

Besides round holes, the unit can have **rectangular cells**, like concrete
blocks (NBR 6136) and structural clay blocks (NBR 15270):

- **Geometry:** number of cells per unit, thickness of the face shells, of the
  webs and of the end shells. With "small cells per face shell" > 0, each face
  gets a row of small cells, as in clay blocks with hollow shells.
- **Mesh:** the structured grid follows all the edges of the cells of all the
  courses. With running bond, the cells move from one course to the next; each
  element is classified course by course (unit, hole or grout). The bed joints
  have mortar where there is unit on both sides and grout where the grouted
  column passes.
- **Grout:** in all the holes or cells, or alternate ones (one yes, one no),
  over the full height. It has its own damage, like the mortar, and its own
  properties. When hot, it follows by default the siliceous concrete of
  EN 1992-1-2: fc(θ) and εc1(θ) of Table 3.1 (the factor on E is
  fc(θ)/εc1(θ)), ft(θ) and the thermal strain. It also applies to round holes.
- **Where it enters:** in 3D (thermal, linear and non-linear), in the 2D section
  (the grout conducts heat in place of the air of the hole) and in Wall mode
  (the grout fibres add to the section, with their own thermal expansion).
- **Library:**
  - concrete blocks 14 × 19 × 39 (Mohamad et al., 2006), 14 × 19 × 29
    (Camacho et al., 2015), M-20 19 × 19 × 39 and infill class C 9 × 19 × 39
    (NBR 6136);
  - structural clay blocks 14 × 19 × 29 with solid shells and with hollow
    shells and 14 × 19 × 39 with 3 cells (NBR 15270);
  - grouts G1 (17.7 MPa) and G2 (26.7 MPa) of Camacho et al. (2015);
  - the material "concrete" now carries the EN 1992-1-2 tables (E, fc, ft and
    thermal strain against temperature).

  In the items without tests, E, fc and ft are estimates, and the note of each
  item says where they came from.
- **Verification:** `examples/12_verification_cells_grout.py`, 15 cases
  (16 with `--full`), all within tolerance:
  - area of the unit course by course (with running bond and hollow shells,
    exact);
  - volume of grout in the prism (exact);
  - stiffness in parallel unit + grout (exact);
  - grout properties in the solvers, with the EN 1992-1-2 tables;
  - grout in the 2D thermal analysis and in the fibres of the wall;
  - consistency of the library.

  With `--full`, it also runs the prism of Camacho et al. (block B1, 3 courses,
  about 10 minutes) and compares it with the tests, without calibrating
  anything:

  | Prism (gross area) | Model | Test | Difference |
  |---|---|---|---|
  | Hollow | 4.27 MPa | 5.63 MPa | −24% |
  | Grouted with G1 | 13.86 MPa | 11.29 MPa | +23% |

  In the hollow prism, the failure of the model comes from crushing of the
  mortar (6.1 MPa). The damage model has no dilatancy, and so it confines the
  mortar little. In the grouted prism, E, ft and ν of the grout are estimates,
  and the model includes neither the shrinkage of the grout nor its
  incompatibility with the unit. These two things explain much of the
  differences. For a study, calibrate the E of the joint and of the grout on a
  prism of your batch.

  The results of the BCEPQ did not change: the prism at 30 °C still gives
  98.68 kN and 1212 MPa, and the verifications 01, 06, 10 and 11 pass.

## Mortar joints and library (version 0.6)

From version 0.6, the joint can be **dry** (as in the BCEPQ: thin layer with
contact, gaps and friction) or **mortar**.

**Mortar joint:**

- It is a continuous layer of real thickness, with its own damage (the same
  Mazars model as the unit, with the parameters of the mortar) and perfectly
  bonded to the unit.
- The bedding can be over the whole face or on the face shells only (holes and
  webs without mortar). The head joints can be filled or empty.
- Length and height of the unit become centre-to-centre dimensions (unit +
  joint).
- When the mortar expands laterally more than the unit (larger ν/E), it pulls
  the unit in tension next to the joint (Hilsdorf mechanism), and the prism
  fails between the strength of the mortar and that of the unit.

**Other new features:**

- **Specific heat and density against temperature** (tables), by the apparent
  heat capacity method. They serve, for example, for the moisture peak of
  concrete and soil-cement.
- **Library** (`mef_masonry/library.py`), with the source of each item:
  - units: BCEPQ and solid brick 100 × 60 × 220 mm (those with rectangular
    cells came in version 0.7; those with horizontal holes appeared disabled
    until phase 3);
  - unit materials: BCEPQ, ceramic and concrete after NBR 15220-2, and concrete
    with k(T) and c(T) of EN 1992-1-2;
  - mortars: the mixes 1:0.25:3, 1:0.5:4.5, 1:1:6 and 1:2:9 (Mohamad,
    Lourenço and Roman, 2006).

  In the interface, choosing the item fills in the fields and shows the source
  and what remains to be checked. When a value is changed by hand, the item
  becomes "Custom".
- **Verification:** `examples/11_verification_mortar_materials.py`, 12 cases,
  all within tolerance:
  - c(T) against the exact solution of a heated body (error of 0.06%);
  - brick + mortar column in series (exact);
  - sign of the lateral Hilsdorf stresses;
  - failure of the mortar without confinement (0.4%) and gain with confinement;
  - properties of the mortar in the solvers;
  - consistency of the library.

  The results of the BCEPQ did not change: the prism at 30 °C still gives
  98.68 kN, and the verifications 01, 06 and 10 pass.

## Wall mode (version 0.5)

For full-scale walls, the solid 3D model becomes too heavy (a wall of
1.26 × 2.58 m has about 300,000 nodes) and lacks the effects that govern a
wall fire test: thermal bowing, second order and eccentricity. Wall mode
(`mef_masonry/wall.py`) solves this with a model of a vertical strip 1 m long:

- **Corotational 2D beam elements** (Crisfield 1991; Battini 2002), with large
  displacements and rotations: second order (P-delta), buckling and thermal
  bowing. The default is 32 elements over the height.
- **Fibre section:** the ceramic triangles of the 2D mesh of one unit length,
  with the holes. Each fibre has its own temperature, taken from the 2D thermal
  analysis with ISO 834 on one face (the same calculation as the 2D section).
- **Homogenised material:** stress (net area) × strain curve of the prism,
  computed by the calibrated 3D non-linear model at 30, 600, 800 and 1000 °C
  (`mef_masonry/wall_data.py`) and interpolated in temperature. Dry masonry
  carries no tension, because the joints open. Unloading is secant, and the
  post-peak branch is regularised by the element length. Optionally, the 30 °C
  curve comes from the prism model with the current inputs.
- **Free thermal strain** of the ceramic (α or table).
- **Inputs:** total load, pinned or fixed supports at the top and the base,
  eccentricities (positive towards the exposed side) and a sinusoidal initial
  imperfection (H/450 by default).
- **Outputs:**
  - cold capacity and reduction factor for slenderness;
  - fire resistance time (loss of equilibrium) and axial contraction criterion
    of EN 1363-1;
  - insulation time (140/180 °C);
  - deflection and top displacement over time;
  - deformed shapes;
  - temperature profiles;
  - fire resistance time × load level curve;
  - in the post-processor, the section with the temperature and the vertical
    stress in the fibres at mid-height.
- **Time:** the capacity takes about 5 s; ISO 834 for 4 h with 5 load levels,
  about 5 min.

Verification (`examples/10_verification_wall.py`, 14 cases, all within
tolerance):

- amplification of the Euler imperfection (0.1–0.3%);
- free thermal bowing (0.15%);
- secant formula (0.03%);
- area and inertia of the fibres (0.04%);
- capacity of the short wall equal to the strength of the prism at 20 and
  600 °C (0.04%);
- free expansion (0.15%);
- convergence with 16, 32 and 64 elements: 32 elements stay within 0.8% of the
  fire resistance time with 64.

Limitations of Wall mode:

- **Vertical edges:** they are free, without restraint to bending.
- **Spalling:** not included. Transient creep has been included since version
  0.9, but only with k_tr > 0 or a given φ(T) (fired clay has no published
  value).
- **Integrity (E):** not checked.
- **Temperature over the height:** the same over the whole wall.
- **α:** it is provisional, and the thermal deflection and the fire resistance
  time depend directly on it.
- **Fixed supports:** the restrained expansion makes the exposed face crush
  next to the support, and the model tends to lose equilibrium early (a
  conservative result).

## Order of use

| Step | Script | What it does | Time and memory* |
| --- | --- | --- | --- |
| 1 | `examples/01_verification.py` | Makes 36 comparisons with analytical and reference solutions; all must be OK | 11 s |
| 2 | `examples/02_calibrate_joints.py` | Fits the stiffness of the dry joints to the measured modulus of the prism and compares it with the series formula and NBR 16868 | 45 s, 1 GB |
| 3 | `examples/03_ufv_test_3d.py` | Simulates the 600, 800 and 1000 °C tests on the 3D prism | ~1 min per temperature, 1.6 GB |
| 4 | `examples/04_section_parametric.py` | Hole diameter × conductivity: temperature at loading, thermal tension and ISO 834 insulation | 1 min |
| 5 | `examples/05_section_optimisation.py` | Optimisation model (differential evolution) and Pareto front mass × insulation | 2 min |
| 6 | `examples/06_verification_nonlinear.py` | Verifies the non-linear model in 11 cases with a known solution | 2 min |
| 7 | `examples/08_calibrate_imperfection.py` | Calibrates the gap of the joints (from the E of the prism) and the effective fc (from the load at 30 °C) | 15 min |
| 8 | `examples/07_ufv_test_nonlinear.py` | Takes the prism to failure at 30, 600, 800 and 1000 °C (non-linear model) | ~5 min at 30 °C; 10–15 min heated |
| 9 | `examples/09_convergence_study.py` | Studies of mesh, element, increment, tolerances and Δt (Table 5 of the paper) | 5–60 min per case |
| 10 | `examples/10_verification_wall.py` | Verifies Wall mode in 14 cases (`--fast` skips the convergence) | 20 s; 4 min full |
| 11 | `examples/11_verification_mortar_materials.py` | Verifies c(T), the mortar joint and the library in 12 cases | 1 min |
| 12 | `examples/12_verification_cells_grout.py` | Verifies rectangular cells and grout in 15 cases (`--full` also runs the prisms of Camacho et al.) | 5 s; 10 min full |
| 13 | `examples/13_verification_horizontal_holes.py` | Verifies the horizontal holes in 12 cases (`--full` also runs the 8-hole prism) | 10 s; 1 min full |
| 14 | `examples/14_verification_creep.py` | Verifies transient creep in 19 cases (`--full` also runs the sensitivity of the BCEPQ wall and the transient test of the solid brick) | 2 min; 10 min full |
| 15 | `examples/15_verification_sensitivity.py` | Verifies the Morris method and the wall runs in 31 cases (`--full` also runs a small parallel study and the recommended mesh) | 1 min; 5 min full |
| — | `examples/00_verify_all.py` | Runs all the verifications and writes the summary | 7 min; 35 min full |

\* Measured on a computer with 2 cores. Default mesh of the linear 3D model:
h = 8 mm, dz = 12 mm (with h = 6 mm and dz = 10 mm, it needs about 5 GB of
memory). The non-linear model uses h = 12 mm and dz = 20 mm, about 3 times
faster; the calibration holds for that mesh.

The results go to `results/` (CSV, PNG and `vtk/`). Run the scripts from the
folder of the package, for example
`python examples/03_ufv_test_3d.py --temperatures 600`. Each script accepts
`--help`.

## The model

**Geometry.** The section of one and a half units (420 × 140 mm, three holes of
84 mm) has a structured mesh outside the holes and an O-grid mesh around each
hole, which follows the circumference. The 3D prism is the extrusion of the
section: steel plate, capping, five courses of 60 mm with dry joints of 2 mm,
capping and plate. The head joints between unit and half unit alternate
between x = 280 and 140 mm. The plane of symmetry y = 70 mm halves the model.

**Materials** (`mef_masonry/parameters.py`, with the sources):

| Property | Value | Source |
| --- | --- | --- |
| k of the ceramic | 0.52 W/m·K | Cerâmica 2026, n = 20 |
| ρ, c | 1210 kg/m³; 1285 J/kg·K | BPQ; c derived from the effusivity (898) and k |
| E of the unit (gross area) | 2801 MPa | BPQ, n = 6 |
| E of the material (net area) | 3907 MPa | 2801 / 0.717 |
| ν | 0.18 | BPQ |
| α | 6×10⁻⁶ /°C | **provisional**, not measured |
| E of the dry joint (2 mm) | 84.4 MPa | calibrated (step 2) for E_prism = 1100 MPa |

Temperature-dependent properties enter through tables (`table_k`,
`table_E_factor`, `table_thermal_strain`), for example the dilatometry curve
with the quartz jump near 573 °C.

**Thermal.** Transient conduction with implicit Euler. The exposed faces
exchange heat by convection (h = 10 W/m²·K) and radiation (ε = 0.8) with the
air of the furnace. The air of the holes conducts through an equivalent
radiative conductivity, k_rad = 4σT³Dε/(2−ε).

**Mechanics.** Linear thermo-elasticity with thermal strain. The platens of the
press are rigid frictionless supports: base with u_z = 0 and, under load, top
with a uniform displacement adjusted so that the reaction equals the load of
the test. The thermal case (free top) and the load case are added. Outputs:
principal stresses, expansion and shortening between plates, strain in the
middle third of the outer face (like an extensometer or DIC) and maximum
tension in the ceramic.

**Elements.** All isoparametric, with the full Jacobian
(`mef_masonry/elements3d.py`):

| Element | Nodes | Integration | On the same mesh | Use |
| --- | --- | --- | --- | --- |
| W6, linear wedge | 6 | 3 × 2 points | reference | default of the paper |
| H8, linear hexahedron | 8 | 2 × 2 × 2 | same nodes, half the elements | checking the influence of the element |
| W15, quadratic wedge | 15 | 6 × 3 points | about 4.7 times more nodes | bending and strong gradients, with a coarse mesh |

Damage and the state of the joints hold for the whole element and are
evaluated at its centre; the stiffness is integrated with the full rule. In the
W15, the width of the crack band uses the linear wedge of the corners, and the
thermal analysis uses consistent capacity and heat-exchange matrices (lumping
on the nodes of the quadratic faces distorts the incoming flux). Pressure loads
use the consistent vector. The W15 is heavy: on a 12 × 20 mm prism, 71,674 nodes
against 15,080 for the W6, and the linear model at 600 °C takes about 7 min and
4.6 GB of memory without pypardiso. For the non-linear model, use the coarser
mesh (16 × 30 mm, 35,869 nodes) and, if you can, install pypardiso; the
program warns when the W15 mesh exceeds 40,000 nodes in the non-linear model.

Choice of the element, non-linear test at 30 °C (Table 5 of the paper; 2 cores;
`examples/09_convergence_study.py mesh_nl --element ...`):

| Element | Mesh (mm) | Nodes | Peak (kN) | Chord E (MPa) | Time |
| --- | --- | --- | --- | --- | --- |
| W6 | 12 × 20 | 15,080 | 98.94 | 1212 | 3.4 min |
| H8 | 12 × 20 | 15,080 | 98.24 | 1211 | 3.1 min |
| W6 | 16 × 30 | 7,638 | 97.37 | 1338 | 0.7 min |
| W15 | 16 × 30 | 35,869 | 95.15 | 1331 | 29 min (1.7 GB) |

The element changes the peak by up to 2.3%; the mesh weighs more on the chord
modulus.

**2D section.** The same thermal analysis, with stresses in generalised plane
strain. It serves for quick studies: furnace with all faces exposed, or a wall
with ISO 834 on one face and ambient on the other. The insulation uses the
140/180 °C criterion on the unexposed face.

## Non-linear model

To reach the failure load, the ceramic and the joints stop being elastic
(`mef_masonry/nonlinear.py`).

**Ceramic: damage with crack closure** (Mazars, 1984; Mazars and
Pijaudier-Cabot, 1989). σ = (1 − d) D (ε − ε_t), with
1 − d = (1 − r·d_t)(1 − d_c):

- d_t (tension) starts with the Rankine criterion, when the effective principal
  stress reaches ft;
- d_c (compression) starts with the von Mises criterion of the compressive
  stresses, when it reaches fc;
- r is the weight of tension in the current state: 1 in pure tension, 0 in pure
  compression. A thermal crack compressed by the press closes and transmits
  compression again.

Both damage variables soften exponentially. The dissipated energy is
regularised by the crack band (Bažant and Oh, 1983), with the width of the
element measured normal to the crack (Oliver, 1989), so that the result does
not depend on the mesh.

**Dry joints and capping: contact with gaps and friction.** They are thin
layers (Desai et al., 1984) that transmit no tension and transmit shear up to
Coulomb friction without cohesion, μ|σ_n|. The bearing faces of the units are
not flat. Without load, only part of the nominal area is in contact: 23% in
untreated blocks (Ngapeya and Waldmann, 2020) and about 15% rising to 95% near
failure (Zahra et al., 2016). The gaps concentrate the load and govern the
cracking and the capacity of the dry prism (Ngapeya, Waldmann and Scholzen,
2018). In the model, each bed joint has a random field of initial gaps
(Gaussian, correlation of 40 mm, 20% of the area in contact without load).
Where there is contact, the joint is as stiff as the ceramic, so that the
flexibility of the joint comes only from the gaps that close. The gap enters as
an initial strain of the closed layer: σ_n = E (ε_n + g/t).

**Loading:** the prism is heated field by field, with the actual thermal fields
of the ramp and the dwell, under the seating load of 2 kN. Then the top platen
moves down under displacement control beyond the peak.

| Parameter | Value | Reference |
| --- | --- | --- |
| Tensile strength ft | 0.66 MPa | Gooch et al. (2023): direct tension of pressed clay brick, net area |
| Compressive strength fc (literature) | 5.76 MPa | BPQ: mean fb of 4.13 MPa on the gross area, divided by 0.717 |
| Tensile fracture energy G_f | 0.018 mm × ft ≈ 12 N/m | ductility index of clay units (van der Pluijm, 1992; Lourenço, 1996; via Dauda et al., 2020). It is consistent with the 0.05–0.085 N/mm measured for ft ≈ 3 MPa (Almeida et al., 2002) |
| Compressive fracture energy G_c | 0.33 mm × fc | compressive ductility index of clay units (Angelillo et al., 2014; via Dauda et al., 2020) |
| Friction in dry joints μ | 0.60, zero cohesion | 0.61–0.64 in sawn dry joints (Lourenço et al., 2004) |
| Area in contact without load | 20% | 23% (Ngapeya and Waldmann, 2020); ~15% (Zahra et al., 2016) |
| Standard deviation of the gaps | **0.028 mm** (calibrated) | chord E of the prism (5–30% of the peak load) = 1100 MPa as the mean of 3 seeds (BPQ: 1001–1216 MPa). Script 08 gives 0.029 mm; with 0.028 mm, the mean E is about 1130 MPa |
| Effective fc of the ceramic | **2.48 MPa** (calibrated) | peak load at 30 °C (98.7 kN); only with `--calibrated` |

**Calibration** (`08_calibrate_imperfection.py`), in two independent steps:
(1) the gap comes only from the initial part of the curve, from the modulus of
the prism, as the mean of three seeds of the random field (the stiffness varies
by about 20% between seeds); (2) with that gap, the effective fc comes from the
peak load at 30 °C. The script gave 0.029 mm and 2.48 MPa (default values:
0.028 mm and 2.48 MPa).

With the strengths of the literature, the model gives 230 kN at 30 °C, more
than twice the test. To reproduce the test, the effective fc must be 0.43 times
that of the unit. Likely explanations:

- the test of the low, wide unit (60 × 140 mm), capped and confined by the
  friction of the platens, overestimates the strength of the material;
- the asperities smaller than the mesh are not in the model;
- the load at 30 °C still has to be confirmed.

The heated tests are predictions, without a new fit.

**Results** (`results/nl_summary_calibrated.csv` and
`nl_load_displacement_calibrated.png`; gap of 0.028 mm, effective fc of
2.48 MPa; mesh h = 12 mm, dz = 20 mm):

| Test | Peak load of the model | Test (UFV) | Difference | Ceramic with d > 0.5 after heating | Chord E 5–30% |
| --- | --- | --- | --- | --- | --- |
| 30 °C | 98.9 kN | 98.7 kN | calibrated | 0.06% | 1212 MPa |
| 600 °C | 96.9 kN | 97.9 kN | −1% | 13% | 877 MPa |
| 800 °C | 97.4 kN | 98.9 kN | −2% | 17% | 960 MPa |
| 1000 °C | 99.6 kN | 114.4 kN | −13% | 18% | 1072 MPa |

With the fc of the literature (5.76 MPa), the same model gives 230 kN at 30 °C.

During heating, 13–18% of the ceramic cracks in tension. Even so, the peak load
hardly changes up to 800 °C, as in the tests, because the thermal cracks close
under the compression of the press. Without closure (isotropic damage, earlier
version, with fc of 2.58 MPa), the prism at 600 °C failed at 80 kN.

At 1000 °C the test gains strength, and the model does not follow because the
properties are constant with temperature. A likely explanation is the
continued sintering near the firing temperature.

The stiffness falls after heating, which is a prediction to be checked with the
DIC or the extensometers of the heated tests.

With other seeds of the gap field (2 and 3), the peak load at 30 °C was 99.5
and 98.8 kN (variation < 1%), and the chord E 964 and 1162 MPa. The stiffness
depends more on the distribution of the gaps than the strength does; the
scatter is similar to that of the BPQ prisms.

Time: about 3 min at 30 °C and 10–11 min in the heated tests (2 cores).
Increments accepted without full convergence (minimum step): 0, 5, 1 and 5,
respectively; the script reports this count.

**Temperature dependence.** The tables `table_ft_factor`, `table_fc_factor`
(CeramicDamage), `table_E_factor` and `table_thermal_strain` (Ceramic) take
the values. The default is constant, and α = 6×10⁻⁶ /°C is provisional. For
hollow clay bricks, the literature indicates that the modulus rises to about
1.75 times and falls above 800 °C, and that the expansion has a peak near
600 °C (Nguyen and Meftah, 2014, as reported by Prieler et al., 2025). The
numbers must come from the original or from Annex D of EN 1996-1-2 (paid
access); they were not invented here.

**Verification** (`examples/06_verification_nonlinear.py`, 11 cases):

- tensile and compressive peaks, with an error < 0.4%;
- fracture energy with two element sizes: error of 5–6% and difference between
  the meshes < 2%;
- a tensile crack under compression recovers the intact stiffness;
- dry joint: no tension, exact compressive stiffness, shear saturating at μN,
  and a gap that only transmits force after closing.

**Limitations:**

- scalar damage: the crack has no orientation of its own, and closure depends
  only on the sign of the stresses;
- friction treated as a secant stiffness, without sliding history;
- local damage with crack band (moderate mesh sensitivity);
- random gaps: the chord E varies by about 20% between seeds; the peak load,
  by less than 1%;
- transient creep by the model of Anderberg and Thelandersson, with k_tr (or
  φ(T)) per material and a transient Poisson ratio equal to the elastic one;
- effective fc calibrated, not measured.

### References of the non-linear model

- Al-Fakih, A.; Mohammed, B.S.; Liew, M.S. (2018). Behavior of the dry bed joint in the mortarless interlocking masonry system: an overview. *Civil Engineering Research Journal*, 4(3), 555639.
- Almeida, J.C.; Lourenço, P.B.; Barros, J.A. (2002). Characterization of brick and brick–mortar interface under uniaxial tension. Seminar, Belo Horizonte.
- Bažant, Z.P.; Oh, B.H. (1983). Crack band theory for fracture of concrete. *Materials and Structures*, 16, 155–177.
- Dauda, J.A.; Iuorio, O.; Lourenço, P.B. (2020). Numerical analysis and experimental characterisation of brick masonry. *Int. J. Masonry Research and Innovation*, 5(3), 321–347.
- Desai, C.S.; Zaman, M.M.; Lightner, J.G.; Siriwardane, H.J. (1984). Thin-layer element for interfaces and joints. *Int. J. Numer. Anal. Methods Geomech.*, 8, 19–43.
- Gooch, L.J.; Masia, M.J.; Stewart, M.G.; Lam, C.Y. (2023). Statistical assessment of tensile and shear properties of unreinforced clay brick masonry. *Construction and Building Materials*, 386, 131578.
- Jirásek, M.; Bažant, Z.P. (2002). *Inelastic Analysis of Structures*. Wiley.
- Lourenço, P.B.; Ramos, L.F.; Vasconcelos, G. (2004). On the cyclic behaviour of stone dry masonry joints. 13th Int. Brick and Block Masonry Conference, Amsterdam.
- Mazars, J. (1984). *Application de la mécanique de l'endommagement au comportement non linéaire et à la rupture du béton de structure*. Thesis, Université Paris VI.
- Mazars, J.; Pijaudier-Cabot, G. (1989). Continuum damage theory — application to concrete. *J. Eng. Mech.*, 115(2), 345–365.
- Ngapeya, G.G.C.; Waldmann, D.; Scholzen, F. (2018). Impact of the height imperfections of masonry blocks on the load bearing capacity of dry-stack masonry walls. *Construction and Building Materials*, 165, 898–913.
- Ngapeya, G.G.C.; Waldmann, D. (2020). Overcome of bed-joint imperfections and improvement of actual contact in dry-stacked masonry. *Construction and Building Materials*, 233, 117173.
- Nguyen, T.-D.; Meftah, F. (2014). Behavior of hollow clay brick masonry walls during fire. Part 2: 3D finite element modeling and spalling assessment. *Fire Safety Journal*, 66, 35–45.
- Oliver, J. (1989). A consistent characteristic length for smeared cracking models. *Int. J. Numer. Methods Eng.*, 28, 461–474.
- Prieler, R. et al. (2025). Thermo-mechanical analysis of masonry brick walls with embedded test specimen under fire exposure. *Fire Technology*, 61, 875–926.
- Zahra, T. et al. (2016), as reported by Al-Fakih et al. (2018): contact area of dry joints from ~15% to ~95% with the load.

## Verification

To check everything at once (for example, after installing or updating):

```
python examples/00_verify_all.py           (about 7 min)
python examples/00_verify_all.py --full    (about 35 min)
```

| Script | What it verifies | Cases |
|---|---|---|
| 01 (W6, H8, W15) | Thermal and linear thermo-elastic analyses against exact solutions | 36 per element |
| 06 (W6, H8, W15) | Damage and joints with friction | 11 per element |
| 10 | Wall mode | 10 (14 full) |
| 11 | c(T), mortar joint and library | 12 |
| 12 | Rectangular cells and grout | 15 (16 full, plus 2 comparisons with tests) |
| 13 | Horizontal holes | 12 (14 full) |
| 14 | Transient creep | 19 (20 full, plus 3 informative results) |
| 15 | Sensitivity (Morris) and wall runs | 31 (32 full, plus 4 informative results) |

`01_verification.py` writes `results/verification.csv`; with `--element H8` or
`--element W15`, it runs the same cases with that element
(`verification_H8.csv`, `verification_W15.csv`). All three pass the 36 cases.
`06_verification_nonlinear.py` accepts the same option, and all three pass the
11 cases of the non-linear model. Cases of `01`:

- transient conduction in a plate (Fourier series), in 2D and 3D, with two time steps;
- uniaxial compression, with a rigid platen and with a pressure;
- uniform free expansion;
- layered column, ceramic + joints (series association);
- plate with a parabolic temperature profile, in 3D and 2D;
- section with holes, comparing the 2D thermal analysis with the extruded 3D model.

Errors: below 0.5% in the stresses and 1.1 °C or less in the conduction with
dt = 15 s. In the last case, H8 and W15 have a tolerance of 1 °C (not 0.5 °C),
because the 3D section is not the same T3 mesh as the 2D model.

## Experimental data

Put files `thermocouples_600.csv` (and 800, 1000) with the columns
`time_min,T_C` in `data/`. Step 3 overlays these measurements on the curves of
the model.

## Limitations

- Linear model (steps 1 to 5): elastic material, without cracking or damage,
  and dry joints as equivalent thin layers. For the failure load, use the
  non-linear model (steps 6 to 8).
- Holes treated as closed cavities, without vertical convection.
- α, E(T) and the tensile strength still need to be measured.
- The calibrated stiffness of the joints depends on the mesh; recalibrate if h
  or dz change.
- Mortar joint perfectly bonded to the unit: no debonding or sliding at the
  interface. This is enough for the compression of prisms, but not for shear or
  out-of-plane bending governed by the bond. The α of the mortar is
  provisional.
- In the 2D section and in Wall mode, the section is that of the unit only (the
  mortar head joints do not enter the 2D thermal analysis). With rectangular
  cells and running bond, the 2D section is that of the first course.
- Grout perfectly bonded to the unit, without shrinkage. In the concrete and
  structural clay blocks of the library, E, fc and ft of the unit and of the
  grout are estimates (see the note of each item). The walls of the cells are
  prismatic (without the taper of real blocks); the library uses equivalent
  thicknesses that give the measured net area.
- Compressive damage has no dilatancy: prisms whose failure is governed by the
  crushing of confined mortar tend to fall below the test.
- Transient creep: the transient Poisson ratio is the elastic one (there are no
  data for another value); the k_tr of fired clay is unknown (default 0); basic
  (long-term) creep is not included.
- Sensitivity (Morris): it is a screening, which ranks the inputs but does not
  apportion the variance of the output (for that, Sobol' indices, with many
  more runs). The effects depend on the chosen ranges; a fire resistance
  without failure until the end of the exposure enters as the duration (a lower
  bound). Only the wall under ISO 834 has this study.
- Horizontal holes: holes open at the ends of the prism, without convection
  along them. In the wall, the vertical section is that of one course (the same
  over the whole height), and the fibres add up the continuous walls; the
  material law is still that of the prism (calibrate it on the prism of your
  unit). In the infill units of the library, E, fc and ft are estimates.

## For an experimental × numerical × theoretical paper

| Comparison | Experimental | Numerical | Theoretical |
| --- | --- | --- | --- |
| Thermal field | thermocouples | step 3 | conduction in a plate (verification) |
| Stiffness at 30 °C | E of the prism (BPQ, DIC) | step 2 | unit + joint in series; NBR 16868 (600 fpk) |
| Thermal stresses | observed cracking, ultrasound | steps 3 and 4 | plate with a parabolic profile |
| Geometry | — | steps 4 and 5 | — |
| Failure load × temperature | peak load of the tests | step 8 | — |
