# Section 5.2: PWO and Mallows-GP model validation

This directory contains the publication-only Experiment 29 implementation for
the four strict-foldover methods reported in Table 1:

- FSA-KD at the fixed default `lambda = 0.5`;
- Hamming-maximin simulated annealing;
- component-position (inverse-position) L2-maximin simulated annealing; and
- no-search random strict foldover (SRS).

The two prespecified response experiments are PWO at SNR 2 and 5 and
noise-free Mallows-GP paths generated at `c = 1` and `c = 4`. The fitted GP
scale is re-estimated by REML. The formal design settings are `m = 6`,
`n = 48, 60`, 50 paired replications, 6,000 logical proposals for each search
method, and master seed `20260831`.

## Reproduce the current PDF Table 1

The exact publication rows, designs, seed ledger, and held-out memberships are
frozen under `data/frozen/section5_2/formal_parent_projection`. From
`paper-reproduction/`, run:

```sh
Rscript code/section5_2/summarize_paper_results.R
```

This writes `section5_2_main_table.csv` below
`outputs/section5_2_paper_sources/` without requiring a new search.

The frozen parent originally contained nine design arms. Table 1 was formed by
selecting the four reported methods after their prediction metrics had been
evaluated on the nine-arm parent common test set. The current PDF prose instead
defines a four-method union. Both domains and their exact sizes are retained in
`data/frozen/section5_2/formal_parent_projection/config/test_domain_audit.csv`.
See the data-directory README for the full provenance note. This distinction
affects held-out nRMSE but not full-space PWO IPV.

## Fresh four-method rerun

The public driver implements the four-method experiment described in the PDF.
A smoke test is:

```sh
env WCRIT29_PROFILE=smoke \
  WCRIT29_OUT_SUBDIR=section5_2_smoke \
  Rscript code/section5_2/run_model_specific.R
```

A formal rerun is:

```sh
env WCRIT29_PROFILE=formal \
  WCRIT29_OUT_SUBDIR=section5_2_formal_four_method_domain \
  WCRIT29_WORKERS=24 \
  Rscript code/section5_2/run_model_specific.R
```

To summarize a fresh run rather than the frozen PDF projection:

```sh
env SEC52_RUN_DIR=outputs/wcrit/section5_2_formal_four_method_domain \
  SEC52_OUTPUT_DIR=outputs/section5_2_fresh_paper_sources \
  Rscript code/section5_2/summarize_paper_results.R
```

All four methods within an `(n, replication)` block reuse the same admissible
initial half-design. The three optimized methods use the same pregenerated
proposal tape and logical-proposal budget; SRS returns the initial design
without search. FSA-KD retains its own composite objective and incremental
implementation, so the comparison is end-to-end rather than criterion-only.
