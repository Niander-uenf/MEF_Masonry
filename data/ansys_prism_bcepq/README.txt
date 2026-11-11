ANSYS model of the BCEPQ prism (Sá, 2025) — input files used in Section 6 of the paper

Origin: the ANSYS project "prejeto mini parede e prisma" (original name, in Portuguese: "mini wall and prism
project") of Bernard Silva Corrêa de Sá (co-author of the paper; undergraduate research, LECIV/UENF),
ANSYS Mechanical 2024 R1, solved on 06/01/2025. Report: Sá, B.S.C. (2025), "Análise e simulação numérica:
aplicação e resultados" (Analysis and numerical simulation: application and results).

  ds_thermal_600C.dat         transient thermal analysis (B5): temperature imposed on the side faces,
                              30 -> 600 °C in 720 s and dwell up to 2520 s; no specific heat, so the
                              result is the steady state (whole prism at the imposed temperature)
  ds_structural_600C.dat      static analysis (C5): SOLID187, 17 bonded contacts (MPC), fixed base,
                              imported temperature (600.0 °C at all 23,070 nodes, T_ref = 22 °C),
                              pressure of 2.6467 MPa on the ledge of the top course (97.8 kN)
  solve_structural_600C.out   solver output of the structural case (contact warnings included)

The files keep their original content; only their names were translated (ds_termico_600C.dat,
ds_estrutural_600C.dat and solve_estrutural_600C.out in versions 1.x).

Results of the report (largest total displacement, mm / largest von Mises stress, MPa):
  600 °C, 97,806 N     1.1151 / 55.858        no temperature, 97,806 N    0.30700 / 34.694
  800 °C, 98,956 N     1.5566 / 73.889        no temperature, 98,956 N    0.31061 / 35.102
  1000 °C, 114,405 N   1.9745 / 92.489        no temperature, 114,405 N   0.35911 / 40.582

The comparison with MEF_Masonry is in examples/paper/benchmark_ansys.py.
