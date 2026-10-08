# Section 5.1 and Appendix B.1 reproduction

This directory contains the design-criterion comparison retained in the
September 24, 2026 review manuscript. It compares four nonreplicated
strict-foldover procedures:

- FSA-KD with the default `lambda = 0.5`;
- Hamming-maximin simulated annealing;
- component-position (inverse-position) L2-maximin simulated annealing; and
- simple random sampling (SRS), a no-search random strict-foldover design.

The formal grid is `m = 6, 10, 20`, `n/m = 1, ..., 5`, with 50 repetitions.
Each optimized method uses 6,000 logical proposals. FSA-KD retains its native
composite objective, neighborhood, and incremental C++ implementation, so this
is an end-to-end comparison within the same admissible class rather than a
criterion-only experiment.

## Reproduce Figure 1 and Figure B1

From `paper-reproduction/`, run:

```sh
Rscript code/section5_1/make_geometry_tradeoff_figures.R
```

The script reads the frozen 3,000-row publication table
`data/section5_1/geometry_four_method_raw_metrics.csv` and writes:

- `outputs/section5_1/figures/fig2_tradeoff_theta_1.pdf`, the main-text
  `c = 1` criterion-plane figure; and
- `outputs/section5_1/figures/figS_tradeoff_theta_4.pdf`, the Appendix B.1
  `c = 4` scale-sensitivity figure.

Both axes are criteria for which larger values are preferred. Every panel has
its own x and y limits, as in the submitted PDF.

## Re-run the design experiment

Smoke test:

```sh
SEC51_SMOKE=true Rscript code/section5_1/run_geometry_factorial.R
```

Formal run:

```sh
SEC51_WORKERS=12 \
SEC51_OUTPUT=outputs/section5_1/geometry_formal_reproduction \
Rscript code/section5_1/run_geometry_factorial.R
```

The master seed is `20260831`. The deterministic seed namespaces are:

```text
FSA-KD:       hash(20260831, "27", "native-FSA-Phi050", m, n, replication)
common init:  hash(20260831, "27", "common-init", m, n, replication, 1)
common moves: hash(20260831, "27", "common-moves", m, n, replication, 1)
```

SRS is the initial strict-foldover design returned without search. The frozen
formal run was compiled on Linux. Because C++ `std::shuffle` is
standard-library dependent, a macOS rerun can select a different design on an
exact FSA-KD tie while reproducing the reported criterion values. The frozen
publication CSV is therefore the source for the exact plotted points.

The earlier lambda-sensitivity Experiment 01 is not part of the current PDF
and is intentionally absent from this publication-only bundle.
