# MEF_Masonry

*Português: veja [LEIA-ME.md](LEIA-ME.md), the full user documentation.*

MEF_Masonry is an open-source Python program for the thermal and
thermo-mechanical finite element analysis of block masonry in fire, from the
unit section to the full-scale wall. It was developed for the dry-stacked
pressed and fired clay block (BCEPQ) tested at the Federal University of
Viçosa and now covers solid units, units with round holes, vertically
perforated clay and concrete blocks with or without grout, and horizontally
perforated clay bricks, with dry or mortar joints.

| | |
|---|---|
| Version | 1.0.0 |
| Paper | [TO FILL: citation of the paper in Advances in Engineering Software, once accepted] |
| Archive | Zenodo, https://doi.org/10.5281/zenodo.23991397 (version 1.0.0) |
| Licence | BSD-3-Clause, see [LICENSE](LICENSE) |
| Contact | Niander Aguiar Cerqueira, Universidade Estadual do Norte Fluminense Darcy Ribeiro (UENF), niander@uenf.br |

## What it does

- **2D and 3D transient heat transfer** with radiation in the holes and cells
  (equivalent conductivity), convection and radiation on the faces, and
  tabulated k(T), c(T) and ρ(T).
- **3D thermo-mechanical prism model:** scalar damage with separate tensile and
  compressive variables and crack closure; dry joints as thin layers with a
  Gaussian random field of initial gaps, no tension and Coulomb friction;
  mortar joints and grout; transient creep (Anderberg and Thelandersson).
  Elements W6, H8 and W15.
- **Wall model:** corotational beam with a fibre section taken from the 2D
  thermal mesh of the unit, second-order effects, initial imperfection and end
  eccentricities, under ISO 834 or a furnace ramp.
- **Global sensitivity of the wall** by the Morris method, with μ*, μ, σ and
  bootstrap intervals, run in parallel.
- **Tests modelled:** ambient compression, hot compression, residual
  compression, transient test and one-sided ISO 834 exposure under load.
- **Browser interface** (local server, nothing leaves the computer) with a
  WebGL 3D post-processor, a library of units and materials with sources, and
  project files in JSON. Results are written as CSV, PNG, VTK (ParaView) and JSON.

The interface and the identifiers of the code are in Portuguese; this file
and the paper describe the program in English.

## Installation

Python 3.9 or newer (tested with 3.13).

```
python -m pip install -r requirements.txt
```

Optional: `pip install pypardiso` speeds up the 3D mechanical solver.

## Quick start

```
python app.py                 # opens the interface in the browser
```

On Windows, double-click `Iniciar_MEF_Masonry.bat`, which finds Python,
installs the dependencies on first use and opens the program. On Linux or
macOS, `./iniciar.sh` does the same.

Without the interface, the same project runs from code:

```python
from mef_masonry.projeto import projeto_padrao, executar
projeto = projeto_padrao()            # the BCEPQ prism of the UFV tests
executar(projeto, "resultados_app/exemplo")
```

A Windows executable that needs no Python is attached to each
[release](../../releases) and can be rebuilt with `python empacotar.py`.

## Verification

```
python exemplos/00_verificar_tudo.py              # 240 cases, about 7 min on two cores
python exemplos/00_verificar_tudo.py --completo   # with the long cases, about 45 min
```

Eight scripts compare the program with closed-form solutions, independent
calculations and consistency conditions (heat conduction and
thermo-elasticity with the three elements, damage and dry joints, wall model,
temperature-dependent properties and mortar joints, rectangular cells and
grout, horizontal perforations, transient creep, Morris screening). Each script
writes a CSV to `resultados/`, the summary goes to
`resultados/resumo_verificacoes.csv`, and the command exits with code 1 if any
case fails. The same suite runs on every push through GitHub Actions
(`.github/workflows/verify.yml`).

## Reproducing the paper

[REPRODUCE.md](REPRODUCE.md) lists the script or project behind every table
and figure of the paper, with its run time.

## Repository layout

| Path | Content |
|---|---|
| `mef_masonry/` | Solver package (meshes, elements, thermal, mechanical, non-linear, wall, sensitivity, library, project schema) |
| `app.py`, `mef_masonry/app/` | Local server and browser interface |
| `exemplos/` | Verification and study scripts; `exemplos/artigo/` holds the studies of the paper |
| `resultados/` | Verification records and the results of the paper studies |
| `dados/` | Place for measured thermocouple records; ANSYS input files of the model compared in Section 6 of the paper |
| `mef_bcepq/` | Alias of the old package name (versions up to 0.10), kept so that old scripts still run |
| `LEIA-ME.md` | Full documentation in Portuguese, version by version |

## Citing

Please cite the paper and the archived version of the software; see
[CITATION.cff](CITATION.cff).

## Limitations

The main limitations are listed in Section 9 of the paper and in `LEIA-ME.md`:
moisture transport is represented only through the specific heat, the
cavities through an equivalent conductivity, the wall law through the BCEPQ
prism curves scaled to other units, and the wall fire resistance has not yet
been compared with a loaded wall test.
