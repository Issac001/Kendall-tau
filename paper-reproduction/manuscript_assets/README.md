# Frozen manuscript assets

These are the numerical tables and figure source PDFs used by
`main_scichina_review_combined.pdf` (September 24, 2026; SHA-256
`7474ebd3a6fcd0fe3f0b701e5d81a31453d7a74338c8993578485704159f4b16`).

| Asset | Manuscript location | Regeneration source |
|---|---|---|
| `fig2_tradeoff_theta_1.pdf` | Figure 1 | `code/section5_1/make_geometry_tradeoff_figures.R` |
| `figS_tradeoff_theta_4.pdf` | Figure B1 | `code/section5_1/make_geometry_tradeoff_figures.R` |
| `fig6_pcb_regret.pdf` | Figure 2 | `code/section5_4/make_pcb_bo_figure.R` |
| `table1_model_validation.csv` | Table 1 | `code/section5_2/summarize_paper_results.R` |
| `table2_four_drug.csv` | Table 2 | `code/section5_3/make_paper_table.R` |
| `table3_pcb.csv` | Table 3 | `code/section5_4/make_pcb_table.R` |

The `fig2`, `figS`, and `fig6` prefixes are retained because they are the
actual LaTeX source filenames. They should not be interpreted as manuscript
figure numbers.

Generated working outputs remain under the ignored `outputs/` directory. The
files here are tracked reference assets so that figure drift can be detected by
raster comparison and numerical-table drift by exact CSV comparison.
