Studies of the paper on the program (Advances in Engineering Software)

Run them from the root of the repository; the results go to results/paper/.

  python examples/paper/wall_convergence.py           section 5.1: convergence of the BCEPQ wall (mesh, elements, step); ~15 min on 2 cores
  python examples/paper/benchmark_ansys.py            section 6: comparison with the ANSYS model of the BCEPQ prism (Sá, 2025); ~6 min
                                                      (--fast: only the 10 mm meshes, ~15 s)
  python examples/paper/speedup.py                    section 5.3: parallel speed-up (8 wall runs, 1 and 2 processes); ~4 min
  python examples/paper/oliveira.py                   section 7.3: insulation of the walls of Oliveira et al. (2021), three hypotheses; ~1 min
  python examples/paper/oliveira_sensitivity.py       section 7.3: moisture, emissivity, mesh, solid and homogeneous section; ~2 min
  python examples/paper/oliveira_moved_moisture.py    section 7.3: water moved to the cold shell (test of moisture transport); ~5 s
  python examples/paper/oliveira_explicit_section.py  section 7.3: view factor of the cells and conduction of the webs; ~5 s
  python examples/paper/morris_checks.py              sections 3.5 and 8.1: Morris without the truncated runs and nominal wall with k_tr = 2.35; ~2 min
  python examples/paper/morris_narrow_alpha.py        section 8.1: Morris screening with α from 5 to 7 ×10⁻⁶/°C; ~12 min on 2 cores
  python examples/paper/charts.py                     section 8.2: design charts (5 heights × 3 values of α × 9 load levels); ~10 min
  python examples/paper/figures.py                    figures 1 to 3, from the results above (SVG; PNG and PDF with rsvg-convert)

The sensitivity study of section 8.1 is the ready-made project "BCEPQ wall 2.94 × 3.00 m — sensitivity under ISO 834"
(results in results/sensitivity_wall_bcepq/).
