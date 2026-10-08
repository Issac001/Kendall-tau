# PCB case study (second application in Section 5.3)

This directory is the paper-only reproduction of Experiment 24. Its executable
scope is deliberately fixed to the PCB analysis reported under Section 5.3 of
the September 24, 2026 review PDF:

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
SEC53_PCB_BUILD_REPS=30 SEC53_PCB_BUILD_WORKERS=4 \
  SEC53_PCB_BUILD_EXACT_ORACLE=true \
  Rscript code/section5_3/pcb/build_frozen_inputs.R
```

This constructor performs the original searches (6,000 native FSA proposal
iterations; 6,000 complete objective evaluations for each unrestricted SA
competitor), draws the no-search SRS routes, enumerates all
`10!` routes separately for each of the five physics profiles, and verifies the
result against the bundled frozen hashes. For a quick design-construction check,
set `SEC53_PCB_BUILD_REPS=1 SEC53_PCB_BUILD_EXACT_ORACLE=false`.

Full paired experiment (parallelize replications on a Unix-like host):

```bash
SEC53_PCB_REPS=30 SEC53_PCB_T=40 SEC53_PCB_WORKERS=8 \
  Rscript code/section5_3/pcb/run_experiment24_paper.R
```

Short deterministic smoke run (the candidate pool remains the formal 5,000-route
pool, so its first acquisition can be compared with the frozen run):

```bash
SEC53_PCB_REPS=1 SEC53_PCB_T=1 \
  SEC53_PCB_OUT=outputs/section5_3_pcb_smoke \
  Rscript code/section5_3/pcb/run_experiment24_paper.R
```

Recreate the manuscript figure directly from the frozen paper curve:

```bash
Rscript code/section5_3/pcb/make_pcb_bo_figure.R
Rscript code/section5_3/pcb/make_pcb_table.R
```

This writes the manuscript-source names `fig6_pcb_regret.pdf` and
`fig6_pcb_regret.png` at the final 4:3 aspect ratio (6.4 by 4.8 inches), with
all 41 evaluation checkpoints and the frozen uncertainty ribbons. The second
command writes the compact four-row Table 3 CSV.

The full runner writes per-replication resumable checkpoints and compact raw
trajectory/acquisition CSVs. When all 30 replications and 40 steps are run, it
also compares the resulting 164 curve means with the frozen manuscript CSV at
a numerical tolerance of `1e-10`.

## Frozen inputs and provenance

`data/frozen/section5_3/pcb/initial_designs_paper.rds` is a lossless four-method
projection of the frozen Experiment 24 initial-design bank. It contains no
design from an unreported method. `g100_w050_oracle.rds` is likewise a projection
containing only the five physical-profile oracles and standardization objects
for the reported scenario. The profile table, balanced replication assignment,
candidate-pool seeds, paper curve, and paper table are included beside them.
The previously frozen paired-contrast CSV is retained only as an archival
provenance record and is not used by the current PDF. `frozen_input_manifest.csv`
records all retained hashes.

The upstream formal run was
`24_pcb_physics_formal_server_20260822_01`, protocol SHA-256
`6cbbce25ee6fdf925c4a63da9b6774f7dda930a6405fc6f7938bc076af537145`.
The PCB coordinate snapshot SHA-256 is
`c375c12d9e67d2791e9a74e2969041ebde41fb3e4c19cfac2aabf62119d3d3b4`.

The original held-out prediction set excluded the union of every design arm in
the larger frozen Experiment 24. Those removed arms are intentionally not
published here. That prediction diagnostic is not reported in the current PDF.
The tracked main table and the paper-only runner cover the three reported
sequential quantities: cumulative standardized regret, final standardized
regret, and oracle-defined top-0.1% hit rate.

## Software

Required R packages are `Rcpp`, `digest`, `gtools`, `dplyr`, `tidyr`, and
`ggplot2`. The runner compiles `code/common/sa_core.cpp` for the Kendall distance.
The frozen formal environment is documented in the repository-level
reproduction notes.
