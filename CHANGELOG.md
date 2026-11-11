# Changelog

Details of each version, with examples and verification cases, are in
[USER_GUIDE.md](USER_GUIDE.md).

## Unreleased
- `examples/paper/figures.py`: the labels of Fig. 2 are rounded once from the stored μ* (253 and 87 instead of 254 and 88).

## 2.0.0
- Archived at Zenodo: https://doi.org/10.5281/zenodo.24147378.
- The whole program is in English: modules, functions, classes and variables, comments and docstrings, the browser
  interface and its messages, the fields and option values of the project files, the names of the scripts, results
  and data files, and the documentation (`USER_GUIDE.md`, formerly `LEIA-ME.md`). Main renamings: `exemplos/` ->
  `examples/` (`exemplos/artigo/` -> `examples/paper/`), `resultados/` -> `results/`, `dados/` -> `data/`;
  modules `parametros` -> `parameters`, `malha` -> `mesh`, `elementos` -> `elements`, `termico` -> `thermal`,
  `mecanico` -> `mechanical`, `naolinear` -> `nonlinear`, `parede` -> `wall`, `dados_parede` -> `wall_data`,
  `sensibilidade` -> `sensitivity`, `biblioteca` -> `library`, `projeto` -> `project`, `saida` -> `output`,
  `curvas` -> `curves`, `estudos` -> `studies`; `projeto_padrao` -> `default_project`, `aplicar` -> `apply`,
  `executar` -> `run`; script options in English (`--completo` -> `--full`, `--calibrado` -> `--calibrated`,
  `--temperaturas` -> `--temperatures`, ...). Numbers are written with a decimal point, also in the mortar mixes
  of the library (`1:0,25:3` -> `1:0.25:3`, `1:0,5:4,5` -> `1:0.5:4.5`).
- Compatibility (`mef_masonry/legacy.py`): project files of versions 1.x (Portuguese group, field and option names)
  are converted when opened, in the interface and by `project.complete`, `project.apply` and `project.run`; the old
  module names (`import mef_masonry.projeto`, `from mef_masonry import biblioteca`) and the old names of their
  functions, classes and constants still work; the library accepts the old selectors and keys
  (`library.apply(p, "bloco", "concreto_m15_39")`). Results are returned with the names of version 2.
  The interface lists the runs of versions 1.x and reuses their inputs.
- The solver is unchanged: the 240 verification cases and the studies of the paper were rerun and give the same
  results as the code of version 1.0.1. The prism results of versions 1.0.x (`convergencia_malha_nl.csv`,
  `comparacao_elementos.csv`, `nl_resumo_calibrado.csv`) had been produced by an earlier code; regenerated, they
  differ from them by up to 1.2% in the peak load (98.9 instead of 98.6 kN at 30 °C), and Tables 5 and 9 of the
  paper were updated.
- `examples/09_convergence_study.py`: option `--element` (W6, H8 or W15), to repeat the comparison of elements of
  Table 5 of the paper.

## 1.0.1
- Solver unchanged. Studies added during the review of the paper, in `exemplos/artigo/` with results in `resultados/artigo/` (in 2.0.0, `examples/paper/` and `results/paper/`):
  `oliveira_umidade_migrada.py` (moisture moved to the unexposed shell), `oliveira_secao_explicita.py` (view factors of the cells and conduction of the webs),
  `morris_verificacoes.py` (Morris effects without the runs that did not fail; nominal wall with k_tr = 2.35) and `morris_alfa_estreito.py` (screening with α from 5 to 7 × 10⁻⁶/°C).
- Figures 2 and 3 state that they are illustrative.

## 1.0.0
- The program is renamed MEF_Masonry (package `mef_masonry`); `mef_bcepq` stays as an alias for old scripts, and old projects and browser data still open.
- Scripts and results of the paper in `exemplos/artigo/` and `resultados/artigo/`, with `REPRODUCE.md` and the figure script.
- English README, citation metadata (CITATION.cff, .zenodo.json) and continuous verification on GitHub Actions.
- `mecanico.resolver_3d` (`mechanical.solve_3d` in 2.0.0): options `base_engastada` / `fixed_base` (base fixed in all directions) and `regiao_carga` / `load_region` (loaded part of the top face); comparison with the ANSYS model of the BCEPQ prism in `exemplos/artigo/benchmark_ansys.py`, with its input files in `dados/ansys_prisma_bcepq/`.

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
