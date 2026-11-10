# Reproducing the paper

Every table and figure of the paper "MEF_Masonry: an open finite element
program for thermal and thermo-mechanical analysis of block masonry in fire"
comes from a script or a ready-made project of this repository. Run the
commands from the repository root; results go to `resultados/` (scripts) or
`resultados_app/` (projects run from the interface). Times are for two cores of
an Intel Xeon at 2.8 GHz.

| Paper | What | How to repeat it | Output | Time |
|---|---|---|---|---|
| Table 3, Section 4 | Verification suite, 240 cases | `python exemplos/00_verificar_tudo.py` | `resultados/resumo_verificacoes.csv` and one CSV per script | 7 min |
| Section 4 | Long verification cases | `python exemplos/00_verificar_tudo.py --completo` | same | 45 min |
| Table 4, Section 5.1 | Convergence of the BCEPQ wall | `python exemplos/artigo/conv_parede.py` | `resultados/artigo/conv_parede.json` | 15 min |
| Table 5, Section 5.2 | Convergence of the dry-stacked prism | `python exemplos/09_estudo_convergencia.py malha_nl --malhas 12x20 10x15 8x12`, and `--elemento H8` / `W15` for the other elements | `resultados/convergencia_malha_nl.csv`, `resultados/comparacao_elementos.csv` | 1 h |
| Section 5.3 | Parallel speed-up | `python exemplos/artigo/speedup.py` | printed | 4 min |
| Tables 7 and 8, Section 6 | Comparison with the ANSYS model of the BCEPQ prism (Sá, 2025) | `python exemplos/artigo/benchmark_ansys.py` (ANSYS input files in `dados/ansys_prisma_bcepq/`) | `resultados/artigo/benchmark_ansys.json` | 6 min |
| Table 9, Section 7.1 | Hot compression of BCEPQ prisms | `python exemplos/07_ensaio_ufv_nao_linear.py --calibrado --temperaturas 30 600 800 1000` | `resultados/nl_resumo_calibrado.csv` | 40 min |
| Table 10, Section 7.3 | Insulation of the walls of Oliveira et al. (2021a) | `python exemplos/artigo/oliveira.py`, `python exemplos/artigo/oliveira_sensibilidade.py` and, for the moisture moved to the unexposed shell, `python exemplos/artigo/oliveira_umidade_migrada.py` | `resultados/artigo/oliveira*.json` | 3 min |
| Fig. 2, Table 11, Section 8.1 | Morris screening of the BCEPQ wall | Interface: ready-made project "Parede BCEPQ 2,94 × 3,00 m — sensibilidade no ISO 834" | `resultados/sensibilidade_parede_bcepq/` | 12 min |
| Fig. 3, Section 8.2 | Illustrative fire-resistance charts | `python exemplos/artigo/abacos.py` | `resultados/artigo/abacos/` | 10 min |

Figures 1 to 3 are drawn from these results by `python exemplos/artigo/figuras.py`,
which writes SVG to `resultados/artigo/figuras/` (and PNG and PDF if `rsvg-convert`
is installed).

The wall analyses of Sections 8.1 and 8.2 use the fast mesh of the sensitivity
studies (10 mm section, 16 beam elements); Section 5.1 shows how far it is from
the recommended mesh.
