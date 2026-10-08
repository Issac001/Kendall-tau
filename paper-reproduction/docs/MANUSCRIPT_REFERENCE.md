# Manuscript reference

This release is audited against the following fixed review manuscript:

- Title: *Space-filling foldover designs for order-of-addition experiments
  under Kendall tau distance criteria*
- Journal template: *Science China Mathematics* manuscript for review
- PDF filename: `main_scichina_review_combined.pdf`
- PDF date: September 24, 2026
- Page count: 24
- Page size: A4
- SHA-256:
  `7474ebd3a6fcd0fe3f0b701e5d81a31453d7a74338c8993578485704159f4b16`

The PDF itself is not redistributed in this repository. Its exact numerical
tables and source figure PDFs are tracked under `manuscript_assets/`, and the
code-to-manuscript map is recorded in `CODE_INDEX.md`.

## Numerical map

| PDF item | Repository source |
|---|---|
| Figure 1 | `code/section5_1/make_geometry_tradeoff_figures.R`; `manuscript_assets/fig2_tradeoff_theta_1.pdf` |
| Table 1 | `code/section5_2/summarize_paper_results.R`; `manuscript_assets/table1_model_validation.csv` |
| Table 2 | `code/section5_3/four_drug/make_paper_table.R`; `manuscript_assets/table2_four_drug.csv` |
| Figure 2 | `code/section5_3/pcb/make_pcb_bo_figure.R`; `manuscript_assets/fig6_pcb_regret.pdf` |
| Table 3 | `code/section5_3/pcb/make_pcb_table.R`; `manuscript_assets/table3_pcb.csv` |
| Figure B1 | `code/section5_1/make_geometry_tradeoff_figures.R`; `manuscript_assets/figS_tradeoff_theta_4.pdf` |

The legacy filename prefixes `fig2`, `fig6`, and `figS` are retained because
they are the names referenced by the final LaTeX source; the table above gives
their actual numbering in the review PDF.

## Scope conventions

- `FSA-KD` in the numerical section means the default `lambda=0.5` method.
- In Sections 5.1--5.2 all four reported methods use strict-foldover designs;
  SRS is a no-search random strict-foldover baseline.
- In the applications, FSA-KD remains strict by construction, while Hamming,
  component-position L2, and SRS retain their native unrestricted classes.
- The four-drug surrogate is the nugget-aware intercept-only Mallows GP.
- The PCB analysis is fixed to `(gamma, omega)=(1,0.5)`.

See `SOURCE_PROVENANCE.md` for the formal-run lineage and the documented
Section 5.2 held-out-domain distinction.
