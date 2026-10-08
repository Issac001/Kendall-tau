# Frozen Section 5.2 publication projection

`formal_parent_projection/` is the four-method projection used to print Table 1
in the September 24, 2026 review PDF. It contains the reported methods and
scenarios only, together with the parent seed ledger, paper design matrices,
and exact held-out membership for every `(n, replication)` block.

The numerical projection comes from the frozen nine-method Experiment 29
parent (`protocol hash
db6cb901415724b9f55892557a81f2cbcc3e7f91cf2f712ea45b7b9a2b1eb5a5`).
The four reported methods are:

- FSA-KD with `lambda = 0.5`;
- strict-foldover Hamming-maximin SA;
- strict-foldover component-position L2-maximin SA; and
- no-search random strict foldover.

## Test-domain provenance

The current PDF text defines a held-out set excluding the union of these four
reported designs. The frozen numbers printed in Table 1 were, however,
evaluated before the parent run was projected to four methods, so their
held-out membership excludes the union of all nine parent designs. This is
recorded rather than hidden:

- `config/parent_common_test_indices.rds` stores every exact parent test set;
- `config/test_domain_audit.csv` gives both the parent and four-method test-set
  sizes;
- `config/paper_designs.rds` stores all four reported design matrices; and
- `raw/` contains the exact four-method rows used by the paper table.

Consequently, `summarize_paper_results.R` reproduces the current PDF Table 1
exactly from this frozen projection. A fresh `run_model_specific.R` execution
implements the four-method union described in the PDF and can give slightly
different nRMSE values. Full-space PWO IPV is unaffected by this distinction.

This discrepancy should be resolved in the next manuscript revision either by
recomputing the table on the four-method held-out domain or by revising the
written test-domain definition. It must not be silently erased from the
reproduction record.
