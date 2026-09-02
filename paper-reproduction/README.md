# Paper-only reproduction bundle

This directory contains the code and compact inputs for the numerical work
that remains in Section 5 and Appendix B of the final manuscript. It is a
publication-only extraction: discarded experiments and superseded analysis
paths are not shipped as runnable code.

## Included scope

| Manuscript location | Experiment | Published design comparison |
|---|---|---|
| Appendix B.1 | weight sensitivity (Experiment 01) | FSA-KD over the reported lambda grid |
| Section 5.1 | geometry criteria (Experiment 27) | strict-foldover FSA-KD, Hamming, position-L2, and SRS |
| Section 5.2 | model-specific validation (Experiment 29) | strict-foldover FSA-KD, Hamming, position-L2, and random foldover |
| Section 5.3 | four-drug study (Experiment 30 with frozen Experiment 21 inputs) | exact strict-foldover FSA-KD, OofA-OA, unrestricted Hamming, unrestricted position-L2, and unrestricted SRS |
| Section 5.4 | PCB study (paper-only Experiment 24 path) | strict-foldover FSA-KD and unrestricted Hamming, position-L2, and SRS at `g100_w050` |

Throughout this bundle, `L2` means component-position (inverse-position)
L2. In Sections 5.1--5.2 all methods are evaluated within the strict-foldover
class. In Sections 5.3--5.4, strict foldover is intrinsic to FSA-KD; the
comparison methods retain their native design classes.

## Setup

Use R 4.2 or newer with a working C++ toolchain. From this directory:

```sh
Rscript config/install_dependencies.R
```

The formal Linux runs used R 4.5.2 on Ubuntu 24.04; the four-drug formal run
used R 4.2.3 on macOS. Exact package versions and provenance are recorded in
[`docs/SOURCE_PROVENANCE.md`](docs/SOURCE_PROVENANCE.md).

Each section has its own README and commands:

- [`code/section5_1/README.md`](code/section5_1/README.md)
- [`code/section5_2/README.md`](code/section5_2/README.md)
- [`code/section5_3/README.md`](code/section5_3/README.md)
- [`code/section5_4/README.md`](code/section5_4/README.md)

Generated files are written below `outputs/`, which is intentionally ignored
by Git. Formal searches are compute-intensive; run each section's smoke command
before launching a formal job.

## Reproducibility controls

- Master seeds, task-level seeds, label maps, and candidate-pool seeds are
  fixed in code or supplied in frozen ledgers.
- Optimization methods are paired within the replication structure described
  in the manuscript.
- Frozen inputs are checked by SHA-256 before use.
- Formal runners write configuration, seed, checkpoint, and integrity records.
- The repository-level `MANIFEST.sha256` fingerprints the delivered bundle;
  scientific validation also uses cardinality, design-class, path, and numeric
  checks rather than treating hashes as the sole evidence.

See [`docs/CODE_INDEX.md`](docs/CODE_INDEX.md) for the file-to-manuscript map.

## Deliberately excluded

The repository contains no runnable code for the omitted Kendall-maximin arms,
unrestricted Section 5.1/5.2 sensitivity arms, strict-baseline Section 5.3/5.4
sensitivity arms, old four-drug PWO--Mallows universal-kriging BO paths,
duplicate exact-lambda labels, PCB scenarios other than `g100_w050`, directed-
adjacency initial-design competitors, or Experiments 31--33. General Kendall
distance and adjacency-kernel utilities remain because they are part of the
published FSA-KD criterion and PCB surrogate, respectively.
