# Source provenance and validation

## Publication-only extraction

The formal research runs were intentionally broader than the final paper so
that sensitivity and audit questions could be answered. This repository is not
a verbatim dump of those omnibus drivers. It is a paper-only fork that removes
unreported executable branches while preserving the retained methods' seeds,
settings, inputs, and pathwise calculations.

The public scripts therefore have their own SHA-256 values. Original formal
source and protocol hashes are recorded below so that the extraction remains
auditable.

## Formal sources

### Appendix B.1 / Experiment 01

- Master seed: `20260523`
- Grid: `m = 5,10,20`, `n/m = 1,...,5`, `lambda = 0,0.1,...,1`
- Planned searches: 3,300; finite archived searches: 3,299
- Search budget: 6,000 iterations over eight restarts; no-improvement cutoff 800

The archived raw and summary CSVs are supplied under `data/section5_1/`.

### Section 5.1 / Experiment 27

- Formal protocol hash: `38742154048e9e265f747779813d6924ec5acf3bb0cf89fdf5cfb693074a00f3`
- Master seed: `20260831`
- Grid: `m = 6,10,20`, `n/m = 1,...,5`; 50 replications
- Budget: 6,000 logical proposals per optimized method

The original parent also evaluated an arm later removed from the manuscript.
The public driver contains only FSA-KD, Hamming, position-L2, and no-search SRS.
The publication CSV contains 3,000 rows (15 cells, 50 paired replications and
four methods), and all rows passed the finite-metric, method-set and
strict-foldover checks.  The `c=1` and `c=4` figure builders were also run from
this CSV; their regenerated PDFs matched the manuscript PDFs pixel for pixel.
The public projection drops the parent's unused `c=16` diagnostic columns;
the unmodified parent setting is retained only in the provenance configuration.

The formal search was compiled on Linux.  A full-budget one-cell macOS replay
reproduced the Hamming, position-L2 and SRS designs and diagnostics exactly.
FSA-KD reproduced the reported PWO and Mallows diagnostics but selected a
different design on an exact tie because C++ `std::shuffle` is standard-library
dependent.  This affects tied design identity, not the reported criterion
values.  The archived publication CSV remains the source for the exact paper
numbers.  In the Appendix-B.1 archive, only the obsolete server prefix in the
unused `design_path` column was replaced by a repository-relative path; metric,
seed and design-filename fields were left unchanged.

### Section 5.2 / Experiment 29

- Formal protocol hash: `db6cb901415724b9f55892557a81f2cbcc3e7f91cf2f712ea45b7b9a2b1eb5a5`
- Formal run hash: `48a665035bdf41cc33911c2e73f5bd91d631c1fff33b2249a0109cdadfebdf07`
- Original source-bundle hash: `b1d5697d50014504fea8964eae8d6ad099291a4087058c07ada854c160fad766`
- Formal environment: R 4.5.2, Ubuntu 24.04.3, Rcpp 1.1.1-1.1,
  dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.3, gtools 3.9.5

The original formal output contains archived sensitivity arms. The public
entry point is fixed to the four methods and four model scenarios in the final
paper. A public smoke run completed with 4 designs, 8 PWO rows, 8 GP rows, and
16 cross-model rows; the paper summarizer independently reproduced the final
`n=48` and `n=60` table values from the formal output.

### Section 5.3 / Experiment 30

- Formal protocol hash: `fb036365ea878b551e9a14d4816ef06a73ffd347926b719848dbaf4ccd2d6e3b`
- Original source-bundle hash: `fe1985554fe93f4acdfed8cadd99890c0748e825ebcf38d90d3ca886e3b7911b`
- Original Experiment 30 driver hash: `d08d9901b5af3526ffa16dc00442d2b7039089b701890730863a89e0db540a17`
- Frozen Experiment 21 parent protocol hash:
  `0948b74f6ba13c82d8302950f25212f31b611dbd606e8fa607b1e45f2381cebb`
- Formal environment: R 4.2.3 on x86_64 macOS

The frozen bank and step-zero summaries were projected to the five manuscript
methods and the intercept-only Mallows GP. Original and public-subset hashes
are recorded in
`data/frozen/section5_3/experiment21_parent/ORIGINAL_PARENT_SHA256.csv`.
The paper-only runner passed an end-to-end smoke run, and the independent bank
builder reproduced all 100 method-by-seed design SHA values.

### Section 5.4 / Experiment 24

- Formal protocol hash: `6cbbce25ee6fdf925c4a63da9b6774f7dda930a6405fc6f7938bc076af537145`
- Parent motion/design protocol hashes:
  `34cd78b073036489839768072767c06ff22f1cb59cf1032074860d912e0d740f`
  and `6182cf56538b7934b271f1fbae0ebd136dedd23f61dcc17146d0bcc12018625c`
- Formal environment: R 4.5.2, Ubuntu 24.04.3, Rcpp 1.1.1-1.1,
  dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.3, gtools 3.9.5

The original formal run covered six physics scenarios and eight design arms.
The public runner contains only the reported `g100_w050` path and four paper
methods. In a one-replication, one-acquisition replay with the formal 5,000-route
candidate pool, all eight incumbent routes and all four EI selections matched
the formal output; the maximum standardized-regret difference was
`2.8e-14`. The paper-only constructor also rebuilt all five full-`10!` exact
oracles: optimal routes matched and the largest oracle/scenario-object numeric
difference was `1.776e-15`. It also reconstructed all 120 paper design matrices
(four methods by 30 replications) from the frozen seeds and matched them
object-for-object. The regenerated manuscript PNG matched the frozen figure
SHA-256.

## Shared formal source hashes

The formal runs used these shared files:

| Component | SHA-256 |
|---|---|
| `wcrit_common.R` | `b189450562c4d0a68e9b0e7f0b94fcde776c48d259aaeb87f6b75d207ecdf670` |
| `wcrit_maximin_dist.R` | `4c42a0a69b66c36b0db5682d571c34a80409f5b74251de3d63f8c1c679638f25` |
| `case_study_common.R` | `41cc655b2f2db37ee020f6720a2847fec7a45af16a53bf8ff316ac119d99ad5c` |
| `sa_core.cpp` | `18be812fb8b34209d18acf3541f09d5627a1861eb1397398ce789a2b5aecb31c` |

The public shared helpers are source-relative and paper-scoped, so their
delivered hashes can differ from these formal originals. `MANIFEST.sha256`
records the exact delivered files.

## What hashes do and do not establish

Hashes verify byte identity and provenance. They are not used as a substitute
for scientific checks. The public validation also verifies method sets,
scenario sets, design classes, row cardinalities, seed pairing, finite metrics,
successful model fits, candidate-pool selections, and numerical agreement with
the retained formal paths.
