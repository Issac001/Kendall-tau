#!/usr/bin/env Rscript

# Paper-only rerun of Section 5.4 (Experiment 24).
# Fixed scope: g100_w050 and the four methods printed in the manuscript.

options(stringsAsFactors = FALSE, warn = 1)

.args <- commandArgs(trailingOnly = FALSE)
.idx <- grep("^--file=", .args)
.this_dir <- if (length(.idx)) {
  dirname(normalizePath(sub("^--file=", "", .args[.idx[[1L]]]), mustWork = TRUE))
} else normalizePath(getwd(), mustWork = TRUE)
.repro_root <- normalizePath(file.path(.this_dir, "..", ".."), mustWork = TRUE)
.common_dir <- file.path(.repro_root, "code", "common")
.frozen_dir <- file.path(.repro_root, "data", "frozen", "section5_4")

for (pkg in c("Rcpp", "digest", "gtools")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R package is missing: ", pkg)
}
source(file.path(.common_dir, "wcrit_common.R"), local = FALSE)
Rcpp::sourceCpp(file.path(.common_dir, "sa_core.cpp"), rebuild = FALSE, showOutput = FALSE)
source(file.path(.this_dir, "pcb_common.R"), local = FALSE)

read_int <- function(name, default, lower, upper) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < lower || value > upper) {
    pcb_stop("%s must be an integer in [%d,%d]", name, lower, upper)
  }
  value
}

reps <- read_int("SEC54_REPS", 30L, 1L, 30L)
n_seq <- read_int("SEC54_T", 40L, 1L, 40L)
workers <- read_int("SEC54_WORKERS", 1L, 1L, 30L)
out_dir <- path.expand(Sys.getenv(
  "SEC54_OUT", unset = file.path(.repro_root, "outputs", "section5_4_experiment24_paper")
))
dir.create(file.path(out_dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  coordinates = file.path(.repro_root, "data", "case_studies", "d493_first10_holes.csv"),
  designs = file.path(.frozen_dir, "initial_designs_paper.rds"),
  oracle = file.path(.frozen_dir, "g100_w050_oracle.rds"),
  profiles = file.path(.frozen_dir, "physics_parameter_profiles.csv"),
  assignment = file.path(.frozen_dir, "replication_profile_assignment.csv"),
  seeds = file.path(.frozen_dir, "seed_ledger_paper.csv")
)
if (any(!file.exists(inputs))) pcb_stop("Missing frozen input: %s", paste(inputs[!file.exists(inputs)], collapse = "; "))

manifest_path <- file.path(.frozen_dir, "frozen_input_manifest.csv")
manifest <- utils::read.csv(manifest_path, check.names = FALSE)
manifest_paths <- file.path(.frozen_dir, manifest$file)
if (any(!file.exists(manifest_paths))) pcb_stop("A manifest-listed frozen input is missing")
manifest_observed <- vapply(
  manifest_paths, digest::digest, character(1L), file = TRUE,
  algo = "sha256", serialize = FALSE
)
if (!identical(unname(manifest_observed), as.character(manifest$sha256))) {
  pcb_stop("Frozen Section 5.4 input manifest verification failed")
}
coordinate_sha <- digest::digest(
  file = inputs[["coordinates"]], algo = "sha256", serialize = FALSE
)
if (!identical(coordinate_sha,
               "c375c12d9e67d2791e9a74e2969041ebde41fb3e4c19cfac2aabf62119d3d3b4")) {
  pcb_stop("PCB coordinate snapshot SHA-256 mismatch")
}

coordinates <- utils::read.csv(inputs[["coordinates"]], check.names = FALSE)
design_bank <- readRDS(inputs[["designs"]])
oracle_bank <- readRDS(inputs[["oracle"]])
profiles <- utils::read.csv(inputs[["profiles"]], check.names = FALSE)
assignment <- utils::read.csv(inputs[["assignment"]], check.names = FALSE)
seeds <- utils::read.csv(inputs[["seeds"]], check.names = FALSE)
C <- pcb_build_cost(coordinates)
kernel_cfg <- pcb_kernel_config()

if (!identical(design_bank$methods, pcb_methods) || length(design_bank$designs) != 30L) {
  pcb_stop("Frozen design bank does not contain exactly the paper method set")
}
if (!identical(oracle_bank$scenario, "g100_w050") ||
    oracle_bank$gamma != 1 || oracle_bank$omega != 0.5) {
  pcb_stop("Frozen oracle is not the paper scenario g100_w050")
}
if (nrow(assignment) != 30L || !identical(assignment$repetition, 1:30)) {
  pcb_stop("Invalid replication-to-profile assignment")
}

# The method-class audit is part of the executable protocol.
class_audit <- do.call(rbind, lapply(pcb_methods, function(method) {
  flags <- vapply(design_bank$designs, function(z) wcrit_is_strict_foldover(z[[method]]), logical(1L))
  data.frame(
    method = method,
    expected_class = if (method == "fsa_lambda05") "strict_foldover" else "unrestricted",
    strict_foldover_repetitions = sum(flags),
    expected_strict_foldover_repetitions = if (method == "fsa_lambda05") 30L else 0L,
    pass = sum(flags) == if (method == "fsa_lambda05") 30L else 0L
  )
}))
if (!all(class_audit$pass)) pcb_stop("Initial-design class audit failed")

pool_for <- function(rep_id, step) {
  row <- seeds[seeds$scope == "bo_candidate_pool" &
                 seeds$repetition == rep_id & seeds$step == step, , drop = FALSE]
  if (nrow(row) != 1L) pcb_stop("No unique pool seed for repetition %d, step %d", rep_id, step)
  expected <- wcrit_hash_seed(20260820L, "22", "bo-raw-pool", rep_id, step)
  if (row$seed[[1L]] != expected) pcb_stop("Frozen pool seed mismatch")
  D <- wcrit_candidate_permutations(10L, 5000L, seed = row$seed[[1L]])
  storage.mode(D) <- "integer"
  list(D = D, seed = row$seed[[1L]], sha256 = wcrit_matrix_sha256(D))
}

run_replication <- function(rep_id) {
  checkpoint <- file.path(out_dir, "checkpoints", sprintf("rep_%03d.rds", rep_id))
  signature <- list(
    experiment = "24-paper-only", scenario = "g100_w050", gamma = 1, omega = 0.5,
    methods = pcb_methods, n_init = 20L, n_seq = n_seq,
    input_sha256 = vapply(inputs, digest::digest, character(1L), file = TRUE,
                          algo = "sha256", serialize = FALSE)
  )
  if (file.exists(checkpoint)) {
    old <- readRDS(checkpoint)
    if (!identical(old$signature, signature)) pcb_stop("Checkpoint signature mismatch: %s", checkpoint)
    return(old)
  }

  profile_id <- assignment$profile_id[assignment$repetition == rep_id]
  profile <- pcb_profile(profiles, profile_id)
  scenario <- oracle_bank$profile_scenarios[[as.character(profile_id)]]
  oracle <- oracle_bank$oracle[oracle_bank$oracle$profile_id == profile_id, , drop = FALSE]
  initial <- design_bank$designs[[rep_id]][pcb_methods]
  states <- lapply(initial, function(D) {
    D <- pcb_validate_routes(D, unique_rows = TRUE)
    list(D = D, y = pcb_eval_response(D, C, profile, scenario))
  })

  trajectory <- acquisition <- list()
  it <- ia <- 0L
  for (step in 0:n_seq) {
    models <- vector("list", length(pcb_methods)); names(models) <- pcb_methods
    for (method in pcb_methods) {
      state <- states[[method]]
      model <- pcb_fit_gp(state$D, state$y, kernel_cfg)
      models[[method]] <- model
      best_id <- which.min(state$y)
      regret <- max(0, state$y[[best_id]] - oracle$optimum[[1L]])
      it <- it + 1L
      trajectory[[it]] <- data.frame(
        repetition = rep_id, physics_profile = profile_id,
        scenario = "g100_w050", gamma = 1, omega = 0.5,
        method = method, method_label = unname(pcb_method_labels[[method]]),
        step = step, n_eval = length(state$y), best_so_far = state$y[[best_id]],
        best_route = paste(state$D[best_id, ], collapse = "-"),
        simple_regret = regret, standardized_regret = regret / oracle$response_sd[[1L]],
        exact_hit = state$y[[best_id]] <= oracle$optimum[[1L]] + 1e-9,
        top_0p1pct_hit = state$y[[best_id]] <= oracle$top_threshold[[1L]] + 1e-9,
        exact_optimum = oracle$optimum[[1L]],
        top_0p1pct_threshold = oracle$top_threshold[[1L]],
        exact_response_sd = oracle$response_sd[[1L]],
        cK = model$cK, cA = model$cA, rho_kendall = model$rho_kendall,
        gp_reml = model$reml
      )
    }
    if (step == n_seq) break

    pool <- pool_for(rep_id, step + 1L)
    pool_keys <- wcrit_row_keys(pool$D)
    pool_y <- pcb_eval_response(pool$D, C, profile, scenario)
    for (method in pcb_methods) {
      state <- states[[method]]
      available <- which(!(pool_keys %in% wcrit_row_keys(state$D)))
      pred <- pcb_predict_gp(models[[method]], pool$D[available, , drop = FALSE])
      incumbent <- min(state$y)
      ei <- pcb_expected_improvement(pred$mean, pred$sd, incumbent)
      local_pick <- pcb_select_ei(ei, available, kernel_cfg$ei_tie_tolerance)
      raw_pick <- available[[local_pick]]
      new_route <- matrix(pool$D[raw_pick, ], nrow = 1L)
      states[[method]]$D <- rbind(state$D, new_route)
      states[[method]]$y <- c(state$y, pool_y[[raw_pick]])
      ia <- ia + 1L
      acquisition[[ia]] <- data.frame(
        repetition = rep_id, physics_profile = profile_id,
        scenario = "g100_w050", gamma = 1, omega = 0.5,
        method = method, method_label = unname(pcb_method_labels[[method]]),
        acquisition_step = step + 1L, candidate_pool_seed = pool$seed,
        candidate_pool_sha256 = pool$sha256, available_n = length(available),
        selected_raw_index = raw_pick,
        selected_route = paste(new_route[1L, ], collapse = "-"),
        selected_y = pool_y[[raw_pick]], expected_improvement = ei[[local_pick]],
        predictive_mean = pred$mean[[local_pick]], predictive_sd = pred$sd[[local_pick]],
        incumbent_before = incumbent
      )
    }
  }
  result <- list(
    signature = signature,
    trajectory = do.call(rbind, trajectory),
    acquisition = do.call(rbind, acquisition)
  )
  saveRDS(result, checkpoint, version = 3)
  result
}

ids <- seq_len(reps)
if (workers > 1L && .Platform$OS.type != "windows") {
  result <- parallel::mclapply(ids, run_replication, mc.cores = workers, mc.preschedule = FALSE)
} else {
  result <- lapply(ids, run_replication)
}
trajectory <- do.call(rbind, lapply(result, `[[`, "trajectory"))
acquisition <- do.call(rbind, lapply(result, `[[`, "acquisition"))

summarize_curve <- function(z) {
  split_data <- split(z, interaction(z$method, z$step, drop = TRUE))
  rows <- lapply(split_data, function(part) {
    n <- nrow(part); avg <- mean(part$standardized_regret); s <- stats::sd(part$standardized_regret)
    se <- s / sqrt(n)
    half <- if (n >= 2L) stats::qt(0.975, n - 1L) * se else NA_real_
    data.frame(
      scenario = "g100_w050", gamma = 1, omega = 0.5,
      method = part$method[[1L]], method_label = part$method_label[[1L]],
      step = part$step[[1L]], metric = "standardized_regret", n = n,
      mean = avg, sd = s, se = se, ci_low = avg - half, ci_high = avg + half,
      total_evaluated = 20L + part$step[[1L]]
    )
  })
  out <- do.call(rbind, rows)
  out[order(match(out$method, pcb_methods), out$step), ]
}
curve <- summarize_curve(trajectory)

utils::write.csv(trajectory, file.path(out_dir, "bo_trajectories.csv"), row.names = FALSE)
utils::write.csv(acquisition, file.path(out_dir, "bo_acquisitions.csv"), row.names = FALSE)
utils::write.csv(curve, file.path(out_dir, "bo_curve_summary.csv"), row.names = FALSE)
utils::write.csv(class_audit, file.path(out_dir, "design_class_audit.csv"), row.names = FALSE)

config <- data.frame(
  key = c("experiment", "scenario", "gamma", "omega", "m", "n_init", "n_seq",
          "replications", "methods", "surrogate", "candidate_pool_size",
          "candidate_pool_seed_rule", "gp_ridge", "observation_noise"),
  value = c("Experiment 24 paper-only", "g100_w050", "1", "0.5", "10", "20",
            n_seq, reps, paste(pcb_methods, collapse = ","),
            "intercept-only additive Kendall--directed-adjacency quotient-kernel GP",
            "5000", "frozen Experiment 24/parent-22 common-pool ledger", "1e-8", "0")
)
utils::write.csv(config, file.path(out_dir, "run_config.csv"), row.names = FALSE)

reference_path <- file.path(.frozen_dir, "section5_4_pcb_native_core_bo_curve.csv")
reference_check <- data.frame(applicable = FALSE, max_abs_mean_difference = NA_real_, pass = NA)
if (reps == 30L && n_seq == 40L && file.exists(reference_path)) {
  reference <- utils::read.csv(reference_path, check.names = FALSE)
  joined <- merge(curve[, c("method", "step", "mean")],
                  reference[, c("method", "step", "mean")],
                  by = c("method", "step"), suffixes = c("_new", "_frozen"))
  delta <- max(abs(joined$mean_new - joined$mean_frozen))
  reference_check <- data.frame(applicable = TRUE, max_abs_mean_difference = delta,
                                pass = nrow(joined) == 164L && delta <= 1e-10)
  if (!reference_check$pass) pcb_stop("Formal paper-curve reproduction check failed")
}
utils::write.csv(reference_check, file.path(out_dir, "frozen_reference_audit.csv"), row.names = FALSE)
message("Section 5.4 paper-only run complete: ", normalizePath(out_dir, mustWork = TRUE))
