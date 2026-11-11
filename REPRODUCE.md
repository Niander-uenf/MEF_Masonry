# Reproducing the paper

Every table and figure of the paper "MEF_Masonry: an open finite element
program for thermal and thermo-mechanical analysis of block masonry in fire"
comes from a script or a ready-made project of this repository. Run the
commands from the repository root; results go to `results/` (scripts) or
`results_app/` (projects run from the interface). Times are for two cores of
an Intel Xeon at 2.8 GHz.

| Paper | What | How to repeat it | Output | Time |
|---|---|---|---|---|
| Table 3, Section 4 | Verification suite, 240 cases | `python examples/00_verify_all.py` | `results/verification_summary.csv` and one CSV per script | 7 min |
| Section 4 | Long verification cases | `python examples/00_verify_all.py --full` | same | 45 min |
| Table 4, Section 5.1 | Convergence of the BCEPQ wall | `python examples/paper/wall_convergence.py` | `results/paper/wall_convergence.json` | 15 min |
| Table 5, Section 5.2 | Convergence of the dry-stacked prism | `python examples/09_convergence_study.py mesh_nl --meshes 12x20 10x15 8x12`; the other elements with `--element H8 --meshes 12x20` and `--element W15 --meshes 16x30`; increments and tolerances with the studies `step_nl` and `tol_nl` | `results/convergence_mesh_nl.csv`, `results/convergence_step_nl.csv`, `results/convergence_tol_nl.csv` | 1 h |
| Table 6, Section 5.3 | Run times and parallel speed-up | Run times: the "Computation time" item in the summary of each analysis; speed-up: `python examples/paper/speedup.py` | `results.json` of each analysis; speed-up printed | 4 min |
| Tables 7 and 8, Section 6 | Comparison with the ANSYS model of the BCEPQ prism (Sá, 2025) | `python examples/paper/benchmark_ansys.py` (ANSYS input files in `data/ansys_prism_bcepq/`) | `results/paper/benchmark_ansys.json` | 6 min |
| Table 9, Section 7.1 | Hot compression of BCEPQ prisms | `python examples/07_ufv_test_nonlinear.py --calibrated --temperatures 30 600 800 1000 --label _calibrated` | `results/nl_summary_calibrated.csv` | 40 min |
| Table 10, Section 7.3 | Insulation of the walls of Oliveira et al. (2021a) | `python examples/paper/oliveira.py`, `python examples/paper/oliveira_sensitivity.py`, `python examples/paper/oliveira_moved_moisture.py` (moisture moved to the unexposed shell) and `python examples/paper/oliveira_explicit_section.py` (view factors of the cells, conduction of the webs) | `results/paper/oliveira*.json` | 3 min |
| Fig. 2, Table 11, Section 8.1 | Morris screening of the BCEPQ wall | Interface: ready-made project "BCEPQ wall 2.94 × 3.00 m — sensitivity under ISO 834 (Morris, 11 parameters)", run into `results/sensitivity_wall_bcepq/` | `results/sensitivity_wall_bcepq/` | 12 min |
| Sections 3.5 and 8.1 | Checks of the screening: runs without failure excluded, transient creep coefficient 2.35 in the nominal wall, range of α narrowed to 5 to 7 × 10⁻⁶/°C | `python examples/paper/morris_checks.py` and `python examples/paper/morris_narrow_alpha.py` | `results/paper/morris_checks.json`, `results/paper/morris_alpha_5_7/` | 2 min and 40 min |
| Fig. 3, Section 8.2 | Illustrative fire-resistance charts | `python examples/paper/charts.py` | `results/paper/charts/` | 10 min |

Figures 1 to 3 are drawn from these results by `python examples/paper/figures.py`,
which writes SVG to `results/paper/figures/` (and PNG and PDF if `rsvg-convert`
is installed).

The wall analyses of Sections 8.1 and 8.2 use the fast mesh of the sensitivity
studies (10 mm section, 16 beam elements); Section 5.1 shows how far it is from
the recommended mesh.

Version 2.0.0 only translated the program into English: the solver is that of
version 1.0.1, and the results in `results/` were regenerated with 2.0.0.
