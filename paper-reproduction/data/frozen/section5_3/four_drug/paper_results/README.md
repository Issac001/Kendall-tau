# Frozen four-drug paper results

This projection contains the minimum formal Experiment 30 output needed to
rebuild Table 2 in the September 24, 2026 review PDF:

- the five reported design procedures only;
- the nugget-aware intercept-only Mallows GP only;
- seven recommendation checkpoints (`n = 12, ..., 18`); and
- the step-zero held-out prediction summary.

The recommendation projection retains 20 design-seed IDs, all 24 component
relabelings, all three leave-one-replicate-out folds, and the exact
held-out-regret and top-one indicators used by the table. The source formal run
is `30_four_drug_mallows_gp_formal_local_20260902_01`, protocol hash
`fb036365ea878b551e9a14d4816ef06a73ffd347926b719848dbaf4ccd2d6e3b`.

Run `code/section5_3/four_drug/make_paper_table.R` without environment variables to
recreate the table from this projection. A full fresh run can still be supplied
through `SEC53_DRUG_RUN_DIR`.
