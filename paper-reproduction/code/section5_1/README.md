# Section 5.1 and Appendix B.1 reproduction

This directory contains only the experiments and figures retained in the
paper:

- Appendix B.1: the FSA-KD weight-sensitivity study;
- Section 5.1: the strict-foldover comparison of FSA-KD, Hamming-maximin SA,
  component-position-`L2`-maximin SA, and no-search SRS;
- the Section 5.1 `c = 1` PWO--Mallows criterion-plane figure; and
- its Appendix B `c = 4` counterpart.

There is no executable Kendall-maximin arm and no unrestricted design arm in
this release. `SRS` means a duplicate-free random strict-foldover design drawn
from the same admissible class, with zero search proposals; it is not an SA
method.

## Files

- `section5_1_helpers.R`: minimal strict-foldover search and diagnostic API.
- `run_geometry_factorial.R`: paper-only `3 x 5` four-method experiment.
- `make_geometry_tradeoff_figures.R`: creates the submitted `c = 1` figure
  and the Appendix `c = 4` figure from the archived four-method data.
- `run_lambda_sensitivity.R`: reconstructs historical Experiment 01.
- `make_lambda_sensitivity_figures.R`: creates the two Appendix B.1 figures.

Shared FSA-KD and distance kernels live in `../common/`.

## Reproduce the submitted figures immediately

From `paper-reproduction/` run:

```sh
Rscript code/section5_1/make_geometry_tradeoff_figures.R
Rscript code/section5_1/make_lambda_sensitivity_figures.R
```

The scripts read:

- `data/section5_1/geometry_four_method_raw_metrics.csv`;
- `data/section5_1/lambda_sensitivity_raw.csv`; and
- `data/section5_1/lambda_sensitivity_summary.csv`.

They write to `outputs/section5_1/figures/` by default. The figure filenames
match the manuscript:

- `fig2_geometry_tradeoff_uniform.pdf`;
- `figS_geometry_tradeoff_c4_uniform.pdf`;
- `fig1_distance_components_n5m.png`; and
- `figA1_lambda_sensitivity_n1m_n4m.png`.

## Re-run the Section 5.1 design experiment

Quick smoke test:

```sh
SEC51_SMOKE=true Rscript code/section5_1/run_geometry_factorial.R
```

Formal defaults (15 cells, 50 repetitions, 6,000 logical proposals for each
optimized method):

```sh
SEC51_WORKERS=12 \
SEC51_OUTPUT=outputs/section5_1/geometry_formal_reproduction \
Rscript code/section5_1/run_geometry_factorial.R
```

The formal grid is `m = 6, 10, 20` and `n/m = 1, ..., 5`. The master seed is
`20260831`. FSA-KD uses its native `Phi_0.5` objective, local--global
neighborhood and incremental C++ implementation. Hamming and position-`L2`
use the same general-purpose SA framework, the same initial half-design seed,
the same move seed, and the same logical-proposal budget. This is an
end-to-end method comparison, not a criterion-only claim.

The deterministic seed formulas retain the frozen Experiment-27 namespaces:

```text
FSA-KD:       hash(20260831, "27", "native-FSA-Phi050", m, n, replication)
common init:  hash(20260831, "27", "common-init", m, n, replication, 1)
common moves: hash(20260831, "27", "common-moves", m, n, replication, 1)
```

SRS is exactly the unoptimized strict-foldover design from `common init`.
Checkpoints make the formal run safely resumable.

The frozen formal run was compiled on Linux, whereas a local macOS replay may
compile the C++ sampler against a different standard library.  In a one-cell
full-budget audit, Hamming, position-`L2`, and SRS reproduced their archived
design hashes and all diagnostics exactly.  FSA-KD attained the same archived
PWO and Mallows diagnostics but selected a different tied design because
`std::shuffle` is implementation dependent across those platforms.  Use the
included frozen four-method CSV for exact paper numbers; use the driver for a
scientifically equivalent fresh rerun, or rerun it on the original Linux
toolchain when byte-identical design identities are required.

### Relation to the frozen run

The paper's optimized FSA-KD, Hamming and position-`L2` designs were produced
by frozen Experiment 27 (`protocol hash
38742154048e9e265f747779813d6924ec5acf3bb0cf89fdf5cfb693074a00f3`).
That parent run also evaluated an arm that was later omitted from the paper.
The publication fork in this directory removes that arm completely and adds
the deterministic no-search SRS reconstruction used in the manuscript. The
archived CSV therefore contains only the four published methods. The original
parent configuration is retained in
`config/section5_1/geometry_parent_frozen_config.csv` for provenance, not as an
instruction to run omitted methods.

## Re-run Appendix B.1

Quick smoke test:

```sh
SEC51_LAMBDA_SMOKE=true Rscript code/section5_1/run_lambda_sensitivity.R
```

Formal historical reconstruction:

```sh
SEC51_LAMBDA_WORKERS=8 \
SEC51_LAMBDA_OUTPUT=outputs/section5_1/lambda_sensitivity_reproduction \
Rscript code/section5_1/run_lambda_sensitivity.R
```

Defaults are `m = 5, 10, 20`, `n/m = 1, ..., 5`, lambda from 0 to 1 in steps
of 0.1, 20 repetitions, master seed `20260523`, 6,000 iterations split over
eight restarts, and an 800-step no-improvement cutoff. The archived run has one
missing result at `(m,n,lambda,rep) = (5,10,0,1)`, so the exact-publication
default preserves it and returns 3,299 finite rows out of 3,300. Set
`SEC51_LAMBDA_KEEP_PUBLISHED_MISSING=false` to attempt a fresh complete rerun.
The raw table's `design_path` values were publication-sanitized from the
original server-absolute prefix to project-relative `outputs/...` paths. No
metric, seed, or design filename was changed, and plotting does not use this
column.

The odd-run cases `(5,5)`, `(5,15)`, and `(5,25)` were formed historically by
deleting one run from an even foldover parent. Their displayed A/B values are
exploratory and are never treated as valid strict-foldover scores. The paper's
default-weight argument uses the 12 even-run panels only.

## Software

Required R packages are `Rcpp`, `RcppArmadillo`, `dplyr`, `tidyr`, `ggplot2`,
`gtools`, and `digest`. Cairo support is required to write the two PDF figures.

Archived input hashes (SHA-256):

```text
0a876114075e291b73fc0a1c8265e6db2c3ac6c3aa58928e9996cad2a68c3ba1  geometry_four_method_raw_metrics.csv
daf06c2605f86618c8b47a457976554eae01414bb1b81abb728e77fb9e73c35c  lambda_sensitivity_raw.csv
40bed8a5efcb3c76baaaa347140deeaa47f6f4b7059d0589728d39398e532857  lambda_sensitivity_summary.csv
```
