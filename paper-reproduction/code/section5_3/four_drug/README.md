# Section 5.3: four-drug experiment

`run_four_drug_mallows_gp.R` is the publication-only fork of Experiment 30.
It retains only the five designs reported in the paper:

- exact strict-foldover FSA-KD with `lambda = 0.5`;
- OofA-OA with its native array structure;
- unrestricted Hamming-maximin SA;
- unrestricted inverse-position (component-position) L2-maximin SA; and
- unrestricted SRS.

Every prediction and expected-improvement decision uses the nugget-aware,
intercept-only Mallows-kernel GP.  The old PWO--Mallows universal-kriging path,
the unrestricted Kendall arm, and the duplicate exact `lambda = 0` and
`lambda = 1` labels are deliberately absent.

The compact frozen parent under
`data/frozen/section5_3/four_drug/experiment21_parent`
contains the audited design-bank subset, label maps, held-out folds,
exact-enumeration records, and seed ledger reused by Experiment 30.  The bank,
manifest, and step-zero summaries have been filtered to the five paper methods
and the intercept-only Mallows GP.  Their public-subset SHA-256 values are
checked before fitting; the provenance notes record the hashes of the original
frozen Experiment 21 inputs from which this subset was derived.
The exact original-to-public mapping is recorded in
`data/frozen/section5_3/four_drug/experiment21_parent/ORIGINAL_PARENT_SHA256.csv`.

## Run

From `paper-reproduction/`:

To reproduce the current PDF Table 2 immediately from the compact frozen
formal projection:

```sh
Rscript code/section5_3/four_drug/make_paper_table.R
```

The result is written to `outputs/section5_3_four_drug_paper_sources/`. The tracked
projection contains the exact five-method recommendation checkpoints used by
the table, so a new 20-by-24-by-3 sequential run is not required merely to
audit the published values.

To reconstruct and SHA-check all five initial-design sequences from the frozen
seeds before running the case study:

```sh
Rscript code/section5_3/four_drug/build_initial_design_bank.R
```

This builder contains no unrestricted Kendall arm and constructs only the
single exact `lambda = 0.5` FSA-KD label used in the paper.

```sh
env WCRIT30_PROFILE=smoke \
  WCRIT30_OUT_SUBDIR=section5_3_four_drug_smoke \
  Rscript code/section5_3/four_drug/run_four_drug_mallows_gp.R
```

For the paper-sized run, use `WCRIT30_PROFILE=formal` and a new output name.
The formal defaults are 20 paired design seeds, all 24 component-label maps,
three held-out response folds, 12 initial runs, and six EI additions.  The
master seed is `20260820`.

After a fresh formal completion, summarize that run with:

```sh
env SEC53_DRUG_RUN_DIR=outputs/wcrit/<formal-run> \
  SEC53_DRUG_OUTPUT_DIR=outputs/section5_3_four_drug_fresh_paper_sources \
  Rscript code/section5_3/four_drug/make_paper_table.R
```

The original formal execution used a broader frozen design bank for auditing.
This publication fork was smoke-tested end to end after method filtering; the
reported arms are pathwise independent of the omitted arms because EI is
computed separately within each design/mapping/fold trajectory.
