# Code index

All commands below are run from `paper-reproduction/`.

## Environment setup

| File | Role |
|---|---|
| `config/install_dependencies.R` | installs the R packages used by the four reproduction workflows |

## Shared implementation

| File | Role |
|---|---|
| `code/common/wcrit_common.R` | permutation, strict-foldover, seed, metric, and C++ loading utilities |
| `code/common/wcrit_maximin_dist.R` | Hamming and component-position L2 maximin utilities |
| `code/common/sa_core.cpp` | compiled incremental FSA-KD/Kendall routines |
| `code/common/case_study_common.R` | Mallows-GP, prediction, EI, and application-design helpers |
| `code/common/paper_plot_style.R` | common manuscript plotting theme and method scales |

## Appendix B.1 and Section 5.1

| File | Manuscript output or role |
|---|---|
| `code/section5_1/run_lambda_sensitivity.R` | Appendix B.1 lambda-sensitivity searches |
| `code/section5_1/make_lambda_sensitivity_figures.R` | Figures B1--B2 |
| `code/section5_1/run_geometry_factorial.R` | 3-by-5 strict-foldover four-method experiment |
| `code/section5_1/make_geometry_tradeoff_figures.R` | Section 5.1 `c=1` figure and Appendix B `c=4` figure |
| `code/section5_1/section5_1_helpers.R` | paper-only strict-foldover search and diagnostics |

The plotting inputs under `data/section5_1/` contain only the four published
methods for the formal comparison. SRS is the no-search random strict-foldover
design; it is not described or implemented as an SA method.

## Section 5.2

| File | Role |
|---|---|
| `code/section5_2/run_model_specific.R` | Experiment 29 driver fixed to four methods and four response scenarios |
| `code/section5_2/search_core.R` | paired strict-foldover construction core |
| `code/section5_2/model_core.R` | full-S6 PWO and Mallows-GP evaluation core |
| `code/section5_2/summarize_paper_results.R` | `n=48,60` paper table and paired intervals |

The four response scenarios are PWO at SNR 2 and 5 and Mallows GP at `c=1`
and `c=4`, with the kernel scale estimated by REML. The paper FSA-KD arm is
fixed to `lambda=0.5`.

## Section 5.3

| File | Role |
|---|---|
| `code/section5_3/build_initial_design_bank.R` | reconstructs and SHA-checks only the five published initial designs |
| `code/section5_3/run_four_drug_mallows_gp.R` | held-out prediction and six EI additions under the sole intercept-only Mallows GP |
| `code/section5_3/make_paper_table.R` | held-out nRMSE, cumulative regret, and top-one-at-16 table |

The compact parent input under `data/frozen/section5_3/experiment21_parent/`
contains only the paper methods and intercept-only Mallows step-zero summaries.
`ORIGINAL_PARENT_SHA256.csv` maps every derived subset back to its frozen
Experiment 21 source.

## Section 5.4

| File | Role |
|---|---|
| `code/section5_4/build_frozen_inputs.R` | reconstructs and verifies only the four paper designs and the five-profile `g100_w050` exact oracle |
| `code/section5_4/pcb_common.R` | strong mixed-physics response and common Kendall--adjacency GP |
| `code/section5_4/run_experiment24_paper.R` | 30 paired BO paths, fixed to `g100_w050` and four methods |
| `code/section5_4/make_pcb_bo_figure.R` | manuscript Figure 6 |

The runner reads the paper-only design and exact-oracle projections in
`data/frozen/section5_4/`; its entry point has no alternative gamma, omega,
method, or surrogate branch.

## Data inputs

| File | Source |
|---|---|
| `data/case_studies/four_drug_oofaexp_0.1.0.csv` | `OofAExp::dat.4drug`, package version 0.1.0 |
| `data/case_studies/d493_first10_holes.csv` | fixed depot plus first ten holes of TSPLIB `d493` |

Full attribution and source hashes are in `data/case_studies/README.md`.
