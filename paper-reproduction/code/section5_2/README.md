# Section 5.2: model-specific validation

This directory is a publication-only fork of Experiment 29.  It runs exactly
the four strict-foldover designs in the manuscript:

- FSA-KD at the fixed default `lambda = 0.5`;
- Hamming-maximin SA;
- inverse-position (component-position) L2-maximin SA; and
- random strict foldover (SRS; no search).

The response experiments are also restricted to the four reported settings:
PWO at SNR 2 and 5, and Mallows-GP paths at `c = 1` and `c = 4`, with `c`
estimated by REML.  Extra FSA weights, Kendall-maximin, `c = 8`, and oracle-c
sensitivity fits are not run by this public entry point.

All four designs within a `(n, replication)` block reuse the same admissible
initial half-design.  The three optimized designs use the same pre-generated
proposal tape and 6,000 proposal steps; SRS returns the initial strict-foldover
design without an SA search.  The master seed is `20260831`.

From `paper-reproduction/`, run a smoke test with:

```sh
env WCRIT29_PROFILE=smoke \
  WCRIT29_OUT_SUBDIR=section5_2_smoke \
  Rscript code/section5_2/run_model_specific.R
```

The paper-sized run uses:

```sh
env WCRIT29_PROFILE=formal \
  WCRIT29_OUT_SUBDIR=section5_2_formal \
  WCRIT29_WORKERS=24 \
  Rscript code/section5_2/run_model_specific.R
```

Then create the paper table and paired intervals:

```sh
env SEC52_RUN_DIR=outputs/wcrit/section5_2_formal \
  Rscript code/section5_2/summarize_paper_results.R
```

The formal protocol has `m = 6`, `n = 48, 60`, 50 paired replications, and a
10,000-draw common-replication bootstrap (seed `20260901`).
