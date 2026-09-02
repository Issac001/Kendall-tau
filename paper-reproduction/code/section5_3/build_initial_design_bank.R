#!/usr/bin/env Rscript

# Rebuild only the five initial designs used in the paper's four-drug study.

options(stringsAsFactors = FALSE, warn = 1)
.args <- commandArgs(trailingOnly = FALSE)
.file <- grep("^--file=", .args, value = TRUE)
.here <- dirname(normalizePath(sub("^--file=", "", .file[[1L]]), mustWork = TRUE))
.root <- normalizePath(file.path(.here, "..", ".."), mustWork = TRUE)
.common <- file.path(.root, "code", "common")

for (pkg in c("Rcpp", "digest", "dplyr", "gtools")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
}
source(file.path(.common, "wcrit_common.R"), local = FALSE)
source(file.path(.common, "wcrit_maximin_dist.R"), local = FALSE)
source(file.path(.common, "case_study_common.R"), local = FALSE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
frozen_dir <- file.path(.root, "data", "frozen", "section5_3", "experiment21_parent")
ledger <- utils::read.csv(file.path(frozen_dir, "config", "seed_ledger.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
frozen <- readRDS(file.path(frozen_dir, "designs", "design_bank.rds"))
output <- path.expand(Sys.getenv(
  "SEC53_BANK_OUTPUT", unset = file.path(.root, "outputs", "section5_3_design_bank")
))
dir.create(output, recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  m = 4L, n0 = 12L, design_seed_ids = 1:20,
  full_rank_retry_max = 50L, rank_tol = 1e-9,
  sa_budget = 6000L, sa_tie_weight = 1e-4,
  sa_replace_prob = 0.25, sa_T0 = 0.05, sa_Tmin = 1e-8
)
methods <- c(
  "Exact_Phi_lambda050", "Unrestricted_Hamming",
  "Unrestricted_Position_L2", "OofA_OA", "SRS"
)

ledger_seed <- function(seed_id, stage, retry_id) {
  hit <- ledger$design_seed_id == seed_id & ledger$stage == stage &
    ledger$retry_id == retry_id
  if (sum(hit) != 1L) stopf("Non-unique ledger lookup")
  as.integer(ledger$seed[hit])
}

pwo_rank <- function(D) {
  qr(wcrit_pwo_matrix(D), tol = cfg$rank_tol, LAPACK = FALSE)$rank
}

entry <- function(D, method, seed_id, stochastic,
                  target_lambda = NA_real_, target_geometry = NA_character_,
                  n_obj_eval = 0L, budget_exhausted = NA,
                  initial_sha256 = NA_character_, initial_seed = NA_integer_,
                  move_seed = NA_integer_, initial_rank_retry = 0L,
                  design_rank_retry = 0L, total_obj_eval = n_obj_eval) {
  D <- case_validate_design(D, cfg$m, cfg$n0, method)
  rank <- pwo_rank(D)
  if (rank != 7L) stopf("Non-full-rank paper design: %s seed %d", method, seed_id)
  list(
    entry_id = sprintf("seed%02d__%s", seed_id, method), method = method,
    design_seed_id = as.integer(seed_id), stochastic = isTRUE(stochastic), D = D,
    design_sha256 = case_matrix_sha256(D),
    initial_design_sha256 = initial_sha256,
    initial_seed = as.integer(initial_seed), move_seed = as.integer(move_seed),
    initial_rank_retry = as.integer(initial_rank_retry),
    design_rank_retry = as.integer(design_rank_retry),
    target_lambda = as.numeric(target_lambda), target_geometry = target_geometry,
    n_obj_eval = as.integer(n_obj_eval), total_obj_eval = as.integer(total_obj_eval),
    budget_exhausted = budget_exhausted, pwo_rank = as.integer(rank)
  )
}

exact <- case_exact_fsa_lambda05_m4n12()
oa <- case_oofa_oa_m4n12()
entries <- list()
cursor <- 0L
for (seed_id in cfg$design_seed_ids) {
  initial <- NULL
  initial_seed <- initial_retry <- NA_integer_
  for (retry_id in 0:cfg$full_rank_retry_max) {
    candidate_seed <- ledger_seed(seed_id, "initial_design_shared", retry_id)
    candidate <- wcrit_sample_unique_permutations(cfg$m, cfg$n0,
                                                   seed = candidate_seed)
    storage.mode(candidate) <- "integer"
    if (pwo_rank(candidate) == 7L) {
      initial <- candidate
      initial_seed <- candidate_seed
      initial_retry <- retry_id
      break
    }
  }
  if (is.null(initial)) stopf("Initial-design retries exhausted: seed %d", seed_id)
  initial_sha <- case_matrix_sha256(initial)

  cursor <- cursor + 1L
  entries[[cursor]] <- entry(
    exact$design, methods[[1L]], seed_id, FALSE, target_lambda = 0.5,
    target_geometry = "kendall_weighted_strict_foldover",
    n_obj_eval = nrow(exact$all_candidates), budget_exhausted = TRUE
  )

  for (j in 2:3) {
    geometry <- c("hamming", "l2_position")[[j - 1L]]
    fit <- NULL
    for (retry_id in 0:cfg$full_rank_retry_max) {
      move_seed <- ledger_seed(seed_id, "sa_move_shared", retry_id)
      candidate_fit <- case_run_unrestricted_sa(
        initial, geometry, cfg$sa_budget, move_seed,
        tie_weight = cfg$sa_tie_weight, replace_prob = cfg$sa_replace_prob,
        T0 = cfg$sa_T0, Tmin = cfg$sa_Tmin
      )
      if (pwo_rank(candidate_fit$D) == 7L) {
        fit <- candidate_fit
        design_retry <- retry_id
        break
      }
    }
    if (is.null(fit)) stopf("SA rank retries exhausted: %s seed %d", geometry, seed_id)
    cursor <- cursor + 1L
    entries[[cursor]] <- entry(
      fit$D, methods[[j]], seed_id, TRUE, target_geometry = geometry,
      n_obj_eval = fit$n_obj_eval, budget_exhausted = fit$budget_exhausted,
      initial_sha256 = initial_sha, initial_seed = initial_seed,
      move_seed = move_seed, initial_rank_retry = initial_retry,
      design_rank_retry = design_retry,
      total_obj_eval = fit$n_obj_eval * (design_retry + 1L)
    )
  }

  cursor <- cursor + 1L
  entries[[cursor]] <- entry(
    oa, methods[[4L]], seed_id, FALSE, target_geometry = "OofA-OA(12,4,2)"
  )
  cursor <- cursor + 1L
  entries[[cursor]] <- entry(
    initial, methods[[5L]], seed_id, TRUE,
    target_geometry = "uniform_without_replacement_conditional_on_full_PWO_rank",
    initial_sha256 = initial_sha, initial_seed = initial_seed,
    initial_rank_retry = initial_retry
  )
}

keys <- function(x) vapply(x, function(z) paste(z$method, z$design_seed_id,
                                                z$design_sha256, sep = "|"), "")
audit <- data.frame(
  check = c("paper_method_set", "entry_count", "design_sha256_sequence"),
  pass = c(
    setequal(unique(vapply(entries, `[[`, "", "method")), methods),
    length(entries) == 100L,
    identical(keys(entries), keys(frozen))
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, file.path(output, "design_bank_audit.csv"), row.names = FALSE)
if (!all(audit$pass)) stopf("Paper design-bank reproduction failed")
saveRDS(entries, file.path(output, "initial_design_bank_paper.rds"), compress = TRUE)
message("Rebuilt and verified the five-method Section 5.3 design bank: ", output)
