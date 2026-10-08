# Reproduction bundle for the September 24, 2026 review manuscript

This directory is synchronized to the 24-page review PDF identified in
`docs/MANUSCRIPT_REFERENCE.md`. It contains the code, compact frozen
publication projections, seeds, designs, exact PCB oracle, and final table and
figure assets needed to audit the numerical work in Section 5 and Appendix B.

## Included scope

| PDF location | Published comparison | Primary code |
|---|---|---|
| Section 5.1 and Appendix B.1 | strict-foldover FSA-KD, Hamming, component-position L2, and SRS; Mallows scales `c=1,4` | `code/section5_1/` |
| Section 5.2 | separate PWO and Mallows-GP response experiments for the same four strict-foldover methods | `code/section5_2/` |
| Section 5.3, four-drug case | exact strict-foldover FSA-KD, OofA-OA, unrestricted Hamming, unrestricted component-position L2, and unrestricted SRS | `code/section5_3/` |
| Section 5.3, PCB case | strict-foldover FSA-KD and unrestricted Hamming, component-position L2, and SRS at `(gamma, omega)=(1,0.5)` | `code/section5_4/` |

The internal `section5_4` directory name is retained for path stability; the
final PDF places both applications in Section 5.3. Throughout the bundle, `L2`
means component-position (inverse-position) L2.

The earlier lambda-sensitivity Experiment 01, Kendall-maximin comparison arms,
unrestricted Section 5.1/5.2 sensitivities, strict application baselines,
superseded PWO--Mallows four-drug paths, extra PCB scenarios, and directed-
adjacency initial-design competitors are not part of the current PDF and are
not shipped as runnable publication experiments.

## Immediate reproduction of published assets

From `paper-reproduction/`:

```sh
Rscript code/section5_1/make_geometry_tradeoff_figures.R
Rscript code/section5_2/summarize_paper_results.R
Rscript code/section5_3/make_paper_table.R
Rscript code/section5_4/make_pcb_bo_figure.R
Rscript code/section5_4/make_pcb_table.R
```

These commands use tracked frozen inputs and write working files below
`outputs/`. The exact source PDFs and table CSVs used in the manuscript are
tracked under `manuscript_assets/`.

## Fresh computational runs

Each section README provides smoke and formal commands. The formal searches are
compute-intensive. Use R 4.2 or newer with a working C++ toolchain and install
the required packages with:

```sh
Rscript config/install_dependencies.R
```

The formal Linux runs used R 4.5.2 on Ubuntu 24.04; the four-drug run used
R 4.2.3 on macOS. Frozen session information and source provenance are recorded
under `data/frozen/` and `docs/SOURCE_PROVENANCE.md`.

## Reproducibility controls

- Master and task seeds, label maps, candidate-pool seeds, and reported design
  matrices are fixed in code or tracked frozen ledgers.
- SRS is always a no-search baseline and is never described as sharing an SA
  engine.
- Sections 5.1--5.2 compare complete procedures inside the strict-foldover
  class; the application comparisons preserve each method's native class.
- Frozen inputs are checked by SHA-256 and by cardinality, design-class, path,
  and numerical audits.
- `MANIFEST.sha256` fingerprints every delivered publication file.

Run the fast read-only release audit with:

```sh
Rscript code/validate_release.R
```

## Section 5.2 provenance note

The current PDF Table 1 is reproduced exactly from the tracked four-method
projection of the formal Experiment 29 parent. The parent evaluated prediction
on a common held-out domain formed before five unreported arms were removed,
whereas the current PDF prose defines the union of the four reported methods.
The exact parent membership and both domain sizes are retained under
`data/frozen/section5_2/`; the discrepancy is documented rather than hidden.
A fresh public-driver run implements the four-method definition and can produce
slightly different held-out nRMSE values. Full-space PWO IPV is unaffected.

See `docs/CODE_INDEX.md` for the file-to-PDF map and
`docs/SOURCE_PROVENANCE.md` for formal-run lineage.
