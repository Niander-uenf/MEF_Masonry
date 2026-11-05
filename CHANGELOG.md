# Changelog

Details of each version, with examples and verification cases, are in
[LEIA-ME.md](LEIA-ME.md) (Portuguese).

## 1.0.0
- The program is renamed MEF_Masonry (package `mef_masonry`); `mef_bcepq` stays as an alias for old scripts, and old projects and browser data still open.
- Scripts and results of the paper in `exemplos/artigo/` and `resultados/artigo/`, with `REPRODUCE.md` and the figure script.
- English README, citation metadata (CITATION.cff, .zenodo.json) and continuous verification on GitHub Actions.
- `mecanico.resolver_3d`: options `base_engastada` (base fixed in all directions) and `regiao_carga` (loaded part of the top face); comparison with the ANSYS model of the BCEPQ prism in `exemplos/artigo/benchmark_ansys.py`, with its input files in `dados/ansys_prisma_bcepq/`.

## 0.10
- Global sensitivity of the wall under ISO 834 by the Morris method, with bootstrap intervals and parallel runs.

## 0.9
- Transient creep (Anderberg and Thelandersson) in the 3D non-linear model and the wall; transient test (constant load, heating to failure).

## 0.8
- Horizontally perforated units (vertical section) and consolidation of the interface.

## 0.7
- Vertically perforated units with rectangular cells (concrete and structural clay blocks), with or without grout.

## 0.6
- Mortar joints with their own damage, temperature-dependent properties and the library of units and materials.

## 0.5
- Wall model: corotational beam with fibre sections fed by the 2D thermal field, second-order effects, ISO 834 under load.

## Earlier versions
- 2D and 3D heat transfer and thermo-elasticity, the non-linear 3D prism model with damage and dry joints, the calibration on the UFV tests and the browser interface.
