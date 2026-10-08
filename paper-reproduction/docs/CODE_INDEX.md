# Code-to-manuscript index

All commands are run from `paper-reproduction/`. The manuscript standard is
the PDF fingerprinted in `MANUSCRIPT_REFERENCE.md`.

## Shared implementation

| File | Role |
|---|---|
| `code/common/wcrit_common.R` | permutation, strict-foldover, deterministic-seed, and design-metric utilities |
| `code/common/wcrit_maximin_dist.R` | Hamming and component-position L2 maximin utilities |
| `code/common/sa_core.cpp` | compiled incremental FSA-KD routines |
| `code/common/case_study_common.R` | Mallows-GP, prediction, EI, and application helpers |
| `code/common/paper_plot_style.R` | manuscript plotting theme and method scales |
| `code/validate_release.R` | fast read-only cardinality, method-set, scenario, and published-value audit |

## Section 5.1 and Appendix B.1

| File | PDF output or role |
|---|---|
| `code/section5_1/run_geometry_factorial.R` | 3-by-5, four-method strict-foldover experiment |
| `code/section5_1/section5_1_helpers.R` | construction and diagnostic API |
| `code/section5_1/make_geometry_tradeoff_figures.R` | Figure 1 at `c=1` and Figure B1 at `c=4` |
| `data/section5_1/geometry_four_method_raw_metrics.csv` | 3,000 frozen plotted rows |

SRS is the no-search random strict-foldover design. The removed historical
lambda-sensitivity experiment is not part of the current Appendix B.1.

## Section 5.2 and Appendix B.2

| File | PDF output or role |
|---|---|
| `code/section5_2/run_model_specific.R` | fresh four-method PWO/Mallows-GP experiment |
| `code/section5_2/search_core.R` | paired strict-foldover construction core |
| `code/section5_2/model_core.R` | full-S6 response, fitting, and prediction core |
| `code/section5_2/summarize_paper_results.R` | Table 1 |
| `data/frozen/section5_2/formal_parent_projection/` | exact Table 1 rows, designs, seeds, and parent held-out sets |

The frozen projection reproduces the current PDF numbers. The fresh runner
implements the four-method held-out union written in the PDF; the documented
domain distinction can slightly change nRMSE.

## Section 5.3: four-drug case

| File | PDF output or role |
|---|---|
| `code/section5_3/four_drug/build_initial_design_bank.R` | reconstruct and SHA-check the five reported initial designs |
| `code/section5_3/four_drug/run_four_drug_mallows_gp.R` | nugget-aware intercept-only Mallows-GP prediction and six EI additions |
| `code/section5_3/four_drug/make_paper_table.R` | Table 2 |
| `data/frozen/section5_3/four_drug/experiment21_parent/` | design bank, label maps, held-out folds, seeds, and step-zero summaries |
| `data/frozen/section5_3/four_drug/paper_results/` | compact formal recommendation projection used by Table 2 |

## Section 5.3: PCB case

| File | PDF output or role |
|---|---|
| `code/section5_3/pcb/build_frozen_inputs.R` | reconstruct four reported designs and five-profile exact oracle |
| `code/section5_3/pcb/pcb_common.R` | fixed mixed-physics response and common surrogate |
| `code/section5_3/pcb/run_experiment24_paper.R` | 30 BO blocks at `(gamma,omega)=(1,0.5)` |
| `code/section5_3/pcb/make_pcb_bo_figure.R` | Figure 2, final 4:3 rendering |
| `code/section5_3/pcb/make_pcb_table.R` | Table 3 |
| `data/frozen/section5_3/pcb/` | coordinates, profiles, designs, oracle, pool seeds, curve, and Table 3 source |

## Tracked final assets

`manuscript_assets/` contains the exact three source PDFs and three table CSVs
used by the review manuscript. Generated working files under `outputs/` remain
ignored.

## External case-study data

| File | Source |
|---|---|
| `data/case_studies/four_drug_oofaexp_0.1.0.csv` | `OofAExp::dat.4drug`, package version 0.1.0 |
| `data/case_studies/d493_first10_holes.csv` | fixed depot and first ten selected TSPLIB `d493` holes |
