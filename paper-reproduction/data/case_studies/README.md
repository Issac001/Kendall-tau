# Case-study data provenance

## Four-drug order-of-administration data

`four_drug_oofaexp_0.1.0.csv` is a text snapshot of
`OofAExp::dat.4drug` from CRAN package `OofAExp` version `0.1.0`.
The package source archive was downloaded from CRAN on 2026-08-20.

- CRAN source archive SHA-256:
  `2ba82d25ec0e380949b32cb4551452233825e01acb205e0c14153c5582d977f0`
- Original `data/dat.4drug.rda` SHA-256:
  `67ddfd1bbf1e7670f0008e8571f22962ffa3c1439809570855d7ea02fc894347`
- Package reference: Shin-Fu Tsai, `OofAExp: Order-of-Addition
  Experiments`, version 0.1.0, GPL-2.
- Experimental source: Wang, A., Xu, H. and Ding, X. (2020),
  *Simultaneous Optimization of Drug Combination Dose-Ratio Sequence with
  Innovative Design and Active Learning*, Advanced Therapeutics 3:1900135.

Rows `t.1`--`t.24` are the complete four-component permutation space. Row
`t.25` is the simultaneous-administration control and is never treated as a
permutation.

## PCB drilling coordinates

`d493_first10_holes.csv` contains a fixed, rule-based subset of the TSPLIB
`d493` drilling instance: TSPLIB node 1 is the depot, and the first ten
subsequent node IDs (2--11) are the holes. This choice was fixed before any
design or BO result was inspected.

- Source file: `TSP/examples/d493.tsp` from R package `TSP` 1.2.5.
- Full source-file SHA-256:
  `a017d8773a1bb279356929a3dd196df11d3fae7b9d45cc5e92d3385616252820`
- TSPLIB metadata: `COMMENT : Drilling problem (Reinelt)`,
  `EDGE_WEIGHT_TYPE : EUC_2D`.

No machine-motion log or calibration data are available. Consequently, the
PCB study combines real drilling coordinates with a preregistered asymmetric
motion model and must be described as a **semi-synthetic asymmetric PCB
benchmark**, not as a fully measured machine study. The coordinate-to-mm
scale and all motion parameters are recorded in each run configuration.
