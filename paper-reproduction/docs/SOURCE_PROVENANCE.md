# Source provenance and validation

## Manuscript anchor

This bundle is aligned to the 24-page PDF created September 24, 2026 and
fingerprinted in `MANUSCRIPT_REFERENCE.md`. The tracked files under
`manuscript_assets/` are the actual figure-source PDFs and numerical table CSVs
used for the alignment audit.

## Section 5.1 / Experiment 27

- Formal protocol hash:
  `38742154048e9e265f747779813d6924ec5acf3bb0cf89fdf5cfb693074a00f3`
- Master seed: `20260831`
- Grid: `m = 6,10,20`, `n/m = 1,...,5`, 50 repetitions
- Budget: 6,000 logical proposals per optimized method

The archived publication CSV has 3,000 rows: 15 cells by 50 repetitions by
four methods. It passed finite-metric, method-set, and strict-foldover checks.
The current plotting script reproduces the final Figure 1 and Figure B1 raster
pixels exactly at 200 dpi. A cross-platform C++ `std::shuffle` difference can
change FSA-KD design identity on an exact tie while preserving the reported
criterion values; the publication CSV is therefore authoritative for the
plotted points.

The historical lambda-sensitivity Experiment 01 is absent because it is not in
the current PDF.

## Section 5.2 / Experiment 29

- Formal parent protocol hash:
  `db6cb901415724b9f55892557a81f2cbcc3e7f91cf2f712ea45b7b9a2b1eb5a5`
- Formal parent run hash:
  `48a665035bdf41cc33911c2e73f5bd91d631c1fff33b2249a0109cdadfebdf07`
- Master seed: `20260831`
- Formal environment: R 4.5.2, Ubuntu 24.04.3

The parent run contained nine design arms. The current Table 1 selects the four
reported arms, but its nRMSE rows retain the parent common held-out membership
that was formed before the other arms were removed. The tracked formal
projection contains:

- 400 reported design rows and all 400 design matrices;
- 800 PWO, 800 Mallows-GP, and 1,600 cross-model rows;
- the 400-row seed ledger for the retained `c=1,4` scenarios; and
- all 100 exact parent common-test index sets and their hashes.

`test_domain_audit.csv` records the parent and four-method-union sizes side by
side. This makes explicit a discrepancy in the current PDF: the prose defines
a four-method union, whereas its printed nRMSE values use the frozen parent
domain. `summarize_paper_results.R` reproduces the printed table; a fresh
four-method driver run follows the prose definition. Full-space PWO IPV is
independent of this choice.

## Section 5.3 / four-drug Experiment 30

- Formal protocol hash:
  `fb036365ea878b551e9a14d4816ef06a73ffd347926b719848dbaf4ccd2d6e3b`
- Frozen Experiment 21 parent protocol hash:
  `0948b74f6ba13c82d8302950f25212f31b611dbd606e8fa607b1e45f2381cebb`
- Formal environment: R 4.2.3 on x86_64 macOS

The publication subset contains exact strict-foldover FSA-KD, OofA-OA, and
unrestricted Hamming, component-position L2, and SRS. Every reported
prediction and acquisition uses the nugget-aware intercept-only Mallows GP.
The tracked parent subset supplies designs, label maps, folds, exact
enumeration, seeds, and step-zero prediction. The compact formal result
projection retains all 50,400 five-method recommendation checkpoints needed to
rebuild Table 2 across 20 design-seed IDs, 24 relabelings, three held-out folds,
and seven evaluation checkpoints.

## Section 5.3 / PCB Experiment 24

- Formal protocol hash:
  `6cbbce25ee6fdf925c4a63da9b6774f7dda930a6405fc6f7938bc076af537145`
- Formal environment: R 4.5.2, Ubuntu 24.04.3

The publication path is fixed to `g100_w050` and four methods. The bundle
contains the coordinate snapshot, five physics profiles, 30-block assignment,
all four initial-design matrices, the five exact `10!` oracles, candidate-pool
seed ledger, 164 curve rows, and Table 3 source. A one-block/one-acquisition
replay previously matched all formal EI selections with a maximum
standardized-regret difference of `2.8e-14`; all five reconstructed oracles
matched to at most `1.776e-15`. The current plotting script reproduces the
final 4:3 Figure 2 raster pixels exactly at 200 dpi.

## Integrity controls

`MANIFEST.sha256` fingerprints every tracked file below `paper-reproduction/`
except the manifest itself. Hashes establish byte identity; the validation
commands additionally check cardinality, method and scenario sets, design
classes, seed pairing, finite model outputs, table values, and rendered figure
pixels.
