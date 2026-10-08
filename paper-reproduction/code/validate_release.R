#!/usr/bin/env Rscript

# Fast, read-only audit of the publication bundle aligned to the September 24,
# 2026 review PDF. Run from any directory with:
#   Rscript paper-reproduction/code/validate_release.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1L]]),
                             winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/",
                      mustWork = TRUE)
stopf <- function(...) stop(sprintf(...), call. = FALSE)
assert <- function(ok, ...) if (!isTRUE(ok)) stopf(...)
read_csv <- function(...) utils::read.csv(file.path(root, ...),
                                          stringsAsFactors = FALSE,
                                          check.names = FALSE)
near <- function(x, y, tolerance = 5e-13) {
  length(x) == length(y) && all(is.finite(x)) &&
    max(abs(as.numeric(x) - as.numeric(y))) <= tolerance
}

paper_methods <- c("FSA-KD", "Hamming", "L2", "SRS")

# Section 5.1: 3 component sizes x 5 run-size ratios x 50 repetitions x 4 methods.
geometry <- read_csv("data", "section5_1",
                     "geometry_four_method_raw_metrics.csv")
assert(nrow(geometry) == 3000L, "Section 5.1 must contain 3,000 rows")
assert(setequal(unique(geometry$paper_label), paper_methods),
       "Section 5.1 method set drifted")
assert(setequal(unique(geometry$m), c(6L, 10L, 20L)),
       "Section 5.1 component grid drifted")
assert(all(geometry$strict_foldover),
       "Section 5.1 contains a non-foldover design")

# Section 5.2: frozen four-method projection of the formal parent run.
sec52 <- file.path("data", "frozen", "section5_2",
                   "formal_parent_projection")
pwo <- read_csv(sec52, "raw", "pwo_prediction_raw.csv")
gp <- read_csv(sec52, "raw", "gp_prediction_raw.csv")
cross <- read_csv(sec52, "raw", "cross_model_loss_raw.csv")
design <- read_csv(sec52, "raw", "design_metrics.csv")
audit <- read_csv(sec52, "config", "test_domain_audit.csv")
assert(nrow(pwo) == 800L && nrow(gp) == 800L && nrow(cross) == 1600L &&
         nrow(design) == 400L,
       "Section 5.2 frozen cardinality drifted")
assert(nrow(audit) == 100L && all(audit$parent_method_count == 9L),
       "Section 5.2 held-out audit drifted")
assert(all(audit$parent_common_test_n < audit$four_method_union_test_n),
       "Section 5.2 domain distinction is no longer present")

table1 <- read_csv("manuscript_assets", "table1_model_validation.csv")
assert(nrow(table1) == 8L && setequal(table1$method, paper_methods),
       "Table 1 method/cardinality audit failed")
expected_ipv <- c(
  0.442735433920410, 0.758540201222552, 0.641725223191939,
  0.822927226777382, 0.308758079804938, 0.443631839045099,
  0.426334448180295, 0.496447650867909
)
assert(near(table1$pwo_ipv_full_space, expected_ipv),
       "Table 1 values drifted")

# Four-drug case: five methods, the primary Mallows GP, and all seven BO checks.
drug_pred <- read_csv("data", "frozen", "section5_3", "paper_results",
                      "results", "initial_prediction_summary.csv")
drug_path <- read_csv("data", "frozen", "section5_3", "paper_results",
                      "raw", "recommendation.csv")
assert(nrow(drug_pred) == 5L &&
         identical(unique(drug_pred$model), "Intercept_Mallows_GP"),
       "Four-drug prediction projection drifted")
assert(nrow(drug_path) == 50400L &&
         setequal(unique(drug_path$bo_step), 0:6),
       "Four-drug BO projection drifted")
table2 <- read_csv("manuscript_assets", "table2_four_drug.csv")
assert(nrow(table2) == 5L && table2$display_method[[1L]] == "FSA-KD",
       "Table 2 method/cardinality audit failed")
assert(near(table2$cumulative_regret[[1L]], 1.65666666666667),
       "Table 2 values drifted")

# PCB case: one preselected scenario, four methods, and 41 checkpoints.
pcb_curve <- read_csv("data", "frozen", "section5_4",
                      "section5_4_pcb_native_core_bo_curve.csv")
assert(nrow(pcb_curve) == 164L &&
         identical(unique(pcb_curve$scenario), "g100_w050") &&
         setequal(unique(pcb_curve$step), 0:40),
       "PCB curve drifted")
table3 <- read_csv("manuscript_assets", "table3_pcb.csv")
assert(nrow(table3) == 4L && setequal(table3$method, paper_methods),
       "Table 3 method/cardinality audit failed")
assert(near(table3$cumulative_regret,
            c(18.9746912825492, 24.9018281036807,
              22.3702697142691, 24.0243960783108)),
       "Table 3 values drifted")

assets <- c(
  "fig2_tradeoff_theta_1.pdf", "figS_tradeoff_theta_4.pdf",
  "fig6_pcb_regret.pdf", "table1_model_validation.csv",
  "table2_four_drug.csv", "table3_pcb.csv"
)
assert(all(file.exists(file.path(root, "manuscript_assets", assets))),
       "A tracked manuscript asset is missing")

message("Release audit passed: September 24, 2026 review manuscript bundle")
