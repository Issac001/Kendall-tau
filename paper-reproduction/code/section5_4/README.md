# Section 5.4: physics-augmented PCB drilling

This directory is the paper-only reproduction of Experiment 24. Its executable
scope is deliberately fixed to the analysis reported in Section 5.4:

- scenario `g100_w050`, i.e. `(gamma, omega) = (1, 0.5)`;
- 10 holes, 20 initial routes, 40 BO additions and 30 paired replications;
- FSA-KD at `lambda = 0.5` in its intrinsic strict-foldover class;
- unrestricted Hamming maximin, unrestricted component-position (inverse-position)
  `L2` maximin, and unrestricted no-search SRS;
- the same intercept-only additive Kendall--directed-adjacency quotient-kernel GP,
  expected improvement, and 5,000-route raw candidate pool for every method.

There is no executable branch for another response scenario, another FSA weight,
Kendall/adjacency initial-design competitors, or a strict-foldover baseline
sensitivity analysis.

## Run from `paper-reproduction/`

Rebuild the frozen four-method initial-design bank and the five-profile exact
`g100_w050` oracle from their seeds and physical definitions:

```bash
SEC54_BUILD_REPS=30 SEC54_BUILD_WORKERS=4 SEC54_BUILD_EXACT_ORACLE=true \
  Rscript code/section5_4/build_frozen_inputs.R
```

This constructor performs the original searches (6,000 native FSA proposal
iterations; 6,000 complete objective evaluations for each unrestricted SA
competitor), draws the no-search SRS routes, enumerates all
`10!` routes separately for each of the five physics profiles, and verifies the
result against the bundled frozen hashes. For a quick design-construction check,
set `SEC54_BUILD_REPS=1 SEC54_BUILD_EXACT_ORACLE=false`.

Full paired experiment (parallelize replications on a Unix-like host):

```bash
SEC54_REPS=30 SEC54_T=40 SEC54_WORKERS=8 \
  Rscript code/section5_4/run_experiment24_paper.R
```

Short deterministic smoke run (the candidate pool remains the formal 5,000-route
pool, so its first acquisition can be compared with the frozen run):

```bash
SEC54_REPS=1 SEC54_T=1 SEC54_OUT=outputs/sec54_smoke \
  Rscript code/section5_4/run_experiment24_paper.R
```

Recreate the manuscript figure directly from the frozen paper curve:

```bash
Rscript code/section5_4/make_pcb_bo_figure.R
```

This writes the manuscript-matching names `fig6_pcb_gamma1_uniform.pdf` and
`fig6_pcb_gamma1_uniform.png`.

The full runner writes per-replication resumable checkpoints and compact raw
trajectory/acquisition CSVs. When all 30 replications and 40 steps are run, it
also compares the resulting 164 curve means with the frozen manuscript CSV at
a numerical tolerance of `1e-10`.

## Frozen inputs and provenance

`data/frozen/section5_4/initial_designs_paper.rds` is a lossless four-method
projection of the frozen Experiment 24 initial-design bank. It contains no
design from an unreported method. `g100_w050_oracle.rds` is likewise a projection
containing only the five physical-profile oracles and standardization objects
for the reported scenario. The profile table, balanced replication assignment,
candidate-pool seeds, paper curve, paper table, and paired contrasts are included
beside them. `frozen_input_manifest.csv` records their hashes.

The upstream formal run was
`24_pcb_physics_formal_server_20260822_01`, protocol SHA-256
`6cbbce25ee6fdf925c4a63da9b6774f7dda930a6405fc6f7938bc076af537145`.
The PCB coordinate snapshot SHA-256 is
`c375c12d9e67d2791e9a74e2969041ebde41fb3e4c19cfac2aabf62119d3d3b4`.

The original held-out prediction set excluded the union of every design arm in
the larger frozen Experiment 24. Those removed arms are intentionally not
published here. Consequently, the exact paper nRMSE is retained in
`section5_4_pcb_native_core_main_table.csv`; the paper-only runner recomputes the
four BO paths, not that larger-union held-out prediction diagnostic.

## Software

Required R packages are `Rcpp`, `digest`, `gtools`, `dplyr`, `tidyr`, and
`ggplot2`. The runner compiles `code/common/sa_core.cpp` for the Kendall distance.
The frozen formal environment is documented in the repository-level
reproduction notes.
