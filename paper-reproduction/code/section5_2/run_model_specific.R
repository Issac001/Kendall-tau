# Section 5.2 publication runner: model-specific PWO and Mallows-GP validation.
#
# This driver is intentionally independent of every existing experiment.  It
# sources two function-only cores from the same directory, constructs one
# atomic checkpoint per (n, replication), and writes auditable raw outputs.
# Summaries, bootstrap inference, and paper figures belong to a later
# postprocessor so that the formal raw-data contract remains small and stable.

.wcrit29_driver_args <- commandArgs(trailingOnly = FALSE)
.wcrit29_driver_idx <- grep("^--file=", .wcrit29_driver_args)
.wcrit29_driver_file <- if (length(.wcrit29_driver_idx)) {
  normalizePath(
    sub("^--file=", "", .wcrit29_driver_args[.wcrit29_driver_idx[1L]]),
    winslash = "/", mustWork = FALSE
  )
} else if (!is.null(sys.frames()[[1L]]$ofile)) {
  normalizePath(sys.frames()[[1L]]$ofile, winslash = "/", mustWork = FALSE)
} else {
  normalizePath("29_pwo_gp_model_specific.R", winslash = "/", mustWork = FALSE)
}
.wcrit29_driver_dir <- dirname(.wcrit29_driver_file)

.wcrit29_search_core <- file.path(.wcrit29_driver_dir, "search_core.R")
.wcrit29_model_core <- file.path(.wcrit29_driver_dir, "model_core.R")
if (!file.exists(.wcrit29_search_core) || !file.exists(.wcrit29_model_core)) {
  stop(
    "Section 5.2 requires search_core.R and model_core.R in the driver directory",
    call. = FALSE
  )
}
.wcrit29_public_common <- normalizePath(
  file.path(.wcrit29_driver_dir, "..", "common"), winslash = "/", mustWork = TRUE
)
Sys.setenv(WCRIT29_HELPER_DIR = .wcrit29_public_common)
source(.wcrit29_search_core, chdir = TRUE)
source(.wcrit29_model_core, chdir = TRUE)

wcrit29_driver_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required", call. = FALSE)
  }
  invisible(TRUE)
}

wcrit29_driver_int <- function(name, default, lower = -Inf) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || value < lower) {
    stop(name, " must be an integer >= ", lower, call. = FALSE)
  }
  value
}

wcrit29_driver_num <- function(name, default, lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || !is.finite(value) || value < lower || value > upper) {
    stop(name, " must be finite and lie in [", lower, ", ", upper, "]", call. = FALSE)
  }
  value
}

wcrit29_driver_bool <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, unset = if (isTRUE(default)) "true" else "false"))
  if (!(value %in% c("true", "false", "1", "0", "yes", "no", "y", "n"))) {
    stop(name, " must be boolean", call. = FALSE)
  }
  value %in% c("true", "1", "yes", "y")
}

wcrit29_driver_num_vec <- function(name, default) {
  raw <- Sys.getenv(name, unset = paste(default, collapse = ","))
  value <- suppressWarnings(as.numeric(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])))
  value <- sort(unique(value))
  if (!length(value) || any(!is.finite(value))) stop("Invalid numeric vector in ", name, call. = FALSE)
  value
}

wcrit29_driver_int_vec <- function(name, default) {
  value <- wcrit29_driver_num_vec(name, default)
  if (any(abs(value - round(value)) > 1e-12)) stop(name, " must contain integers", call. = FALSE)
  as.integer(value)
}

wcrit29_driver_char_vec <- function(name, default) {
  raw <- Sys.getenv(name, unset = paste(default, collapse = ","))
  value <- unique(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]]))
  value <- value[nzchar(value)]
  if (!length(value)) stop(name, " is empty", call. = FALSE)
  value
}

wcrit29_driver_sha256_file <- function(path) {
  wcrit29_driver_require("digest")
  if (!file.exists(path)) stop("Cannot hash missing file: ", path, call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

wcrit29_driver_hash_object <- function(object) {
  wcrit29_driver_require("digest")
  digest::digest(object, algo = "sha256", serialize = TRUE)
}

wcrit29_driver_atomic_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = "gzip")
  if (!file.rename(temporary, path)) stop("Atomic checkpoint rename failed: ", path, call. = FALSE)
  invisible(path)
}

wcrit29_driver_atomic_json <- function(object, path, pretty = TRUE) {
  wcrit29_driver_require("jsonlite")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  jsonlite::write_json(object, temporary, auto_unbox = TRUE, pretty = pretty, na = "null")
  if (!file.rename(temporary, path)) stop("Atomic JSON rename failed: ", path, call. = FALSE)
  invisible(path)
}

wcrit29_driver_profile_defaults <- function(profile) {
  profile <- match.arg(profile, c("smoke", "pilot", "formal"))
  switch(
    profile,
    smoke = list(n = 48L, reps = 1L, budget = 120L, workers = 1L),
    pilot = list(n = c(48L, 60L), reps = 5L, budget = 6000L, workers = 4L),
    formal = list(n = c(48L, 60L), reps = 50L, budget = 6000L, workers = 24L)
  )
}

wcrit29_driver_config <- function() {
  profile <- tolower(Sys.getenv("WCRIT29_PROFILE", unset = "smoke"))
  defaults <- wcrit29_driver_profile_defaults(profile)
  output_override <- Sys.getenv("WCRIT29_OUT_DIR", unset = "")
  out_subdir <- Sys.getenv(
    "WCRIT29_OUT_SUBDIR",
    unset = paste0("29_pwo_gp_model_specific_", profile)
  )
  cfg <- list(
    experiment = "29_pwo_gp_model_specific",
    protocol_version = "2026-08-31-v1-model-specific",
    profile = profile,
    m = 6L,
    n_values = wcrit29_driver_int_vec("WCRIT29_N_VALUES", defaults$n),
    reps = wcrit29_driver_int("WCRIT29_REPS", defaults$reps, 1L),
    lambdas = 0.5,
    snr_values = wcrit29_driver_num_vec("WCRIT29_SNR_VALUES", c(2, 5)),
    c_values = c(1, 4),
    fit_modes = "estimated_c",
    c_fit_grid = wcrit29_driver_num_vec("WCRIT29_C_FIT_GRID", wcrit29_default_c_grid()),
    cross_gp_fit_mode = Sys.getenv("WCRIT29_CROSS_GP_FIT_MODE", unset = "estimated_c"),
    primary_loss = Sys.getenv("WCRIT29_PRIMARY_LOSS", unset = "nrmse_common"),
    primary_regret = Sys.getenv("WCRIT29_PRIMARY_REGRET", unset = "absolute"),
    rng_kind = "Mersenne-Twister",
    rng_normal_kind = "Inversion",
    rng_sample_kind = "Rejection",
    master_seed = wcrit29_driver_int("WCRIT29_MASTER_SEED", 20260831L, 1L),
    proposal_budget = wcrit29_driver_int("WCRIT29_PROPOSAL_BUDGET", defaults$budget, 1L),
    T0 = wcrit29_driver_num("WCRIT29_T0", 1, .Machine$double.eps),
    alpha = wcrit29_driver_num("WCRIT29_ALPHA", 0.997, .Machine$double.eps, 1),
    tie_weight = wcrit29_driver_num("WCRIT29_TIE_WEIGHT", 1e-6, 0),
    gp_jitter = wcrit29_driver_num("WCRIT29_GP_JITTER", 1e-10, 0),
    rank_tol = wcrit29_driver_num("WCRIT29_RANK_TOL", 1e-9, .Machine$double.eps),
    confidence = wcrit29_driver_num("WCRIT29_CONFIDENCE", 0.95, .Machine$double.eps, 1),
    workers = wcrit29_driver_int("WCRIT29_WORKERS", defaults$workers, 1L),
    resume = wcrit29_driver_bool("WCRIT29_RESUME", TRUE),
    out_subdir = out_subdir,
    out_dir_override = output_override
  )
  cfg$n_values <- sort(unique(cfg$n_values))
  cfg$fit_modes <- sort(unique(cfg$fit_modes))
  cfg
}

wcrit29_driver_validate_config <- function(cfg) {
  if (!identical(cfg$m, 6L)) stop("Experiment 29 is frozen at m=6", call. = FALSE)
  if (any(cfg$n_values <= 16L) || any(cfg$n_values %% 2L != 0L) || any(cfg$n_values > 720L)) {
    stop("Every n must be even and satisfy 16 < n <= 720", call. = FALSE)
  }
  if (!identical(cfg$lambdas, 0.5)) stop("The paper runner is fixed at lambda=0.5", call. = FALSE)
  if (!identical(cfg$snr_values, c(2, 5))) stop("SNR grid is frozen at {2,5}", call. = FALSE)
  if (!identical(cfg$c_values, c(1, 4))) stop("The paper c grid is {1,4}", call. = FALSE)
  if (!identical(cfg$fit_modes, "estimated_c")) {
    stop("The paper runner uses estimated_c only", call. = FALSE)
  }
  if (!identical(cfg$cross_gp_fit_mode, "estimated_c")) {
    stop("The confirmatory cross-model GP fit mode is estimated_c", call. = FALSE)
  }
  if (!identical(cfg$primary_loss, "nrmse_common")) {
    stop("The frozen primary loss is nrmse_common", call. = FALSE)
  }
  if (!identical(cfg$primary_regret, "absolute")) {
    stop("The frozen primary regret is absolute", call. = FALSE)
  }
  if (cfg$profile == "formal") {
    if (!identical(cfg$n_values, c(48L, 60L)) || cfg$reps != 50L ||
        !identical(cfg$lambdas, 0.5) ||
        cfg$proposal_budget != 6000L || cfg$master_seed != 20260831L) {
      stop("Formal profile differs from the frozen protocol", call. = FALSE)
    }
  }
  invisible(TRUE)
}

wcrit29_driver_method_id <- function(internal) {
  internal <- as.character(internal)
  out <- internal
  out[internal == "Hamming"] <- "SA_hamming_foldover"
  out[internal == "L2"] <- "SA_l2_position_foldover"
  out[internal == "RandomFoldover"] <- "random_foldover"
  out
}

wcrit29_driver_method_dictionary <- function(lambdas) {
  internal <- c(
    vapply(lambdas, function(x) sprintf("FSA_lambda%03d", as.integer(round(100 * x))), character(1L)),
    "Hamming", "L2", "RandomFoldover"
  )
  method_id <- wcrit29_driver_method_id(internal)
  lambda <- rep(NA_real_, length(internal))
  fsa <- grepl("^FSA_lambda", internal)
  lambda[fsa] <- as.numeric(sub("^FSA_lambda", "", internal[fsa])) / 100
  confirmatory <- method_id %in% c(
    "FSA_lambda050", "SA_hamming_foldover",
    "SA_l2_position_foldover", "random_foldover"
  )
  data.frame(
    method_id = method_id,
    internal_method = internal,
    display_label = vapply(internal, wcrit29_method_label, character(1L)),
    family = ifelse(fsa, "FSA_KD", ifelse(internal == "RandomFoldover", "random", "geometry_SA")),
    lambda = lambda,
    criterion = ifelse(
      fsa, "Phi_lambda",
      c(rep(NA_character_, sum(fsa)), "Hamming", "component_position_L2", "none")
    ),
    design_class = "strict_foldover",
    l2_variant = ifelse(internal == "L2", "component_position", NA_character_),
    confirmatory7 = confirmatory,
    extended = FALSE,
    stringsAsFactors = FALSE
  )
}

wcrit29_driver_scenario_dictionary <- function(cfg) {
  pwo <- data.frame(
    model_family = "pwo",
    scenario_id = paste0("pwo_snr", vapply(cfg$snr_values, wcrit29_num_tag, character(1L))),
    fit_mode = "ols",
    snr = cfg$snr_values,
    c_true = NA_real_,
    main_cross_flag = cfg$snr_values %in% c(2, 5),
    analysis_role = "confirmatory",
    stringsAsFactors = FALSE
  )
  gp <- expand.grid(
    c_true = cfg$c_values,
    fit_mode = cfg$fit_modes,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  gp$model_family <- "mallows_gp"
  gp$scenario_id <- paste0("gp_c", vapply(gp$c_true, wcrit29_num_tag, character(1L)))
  gp$snr <- NA_real_
  gp$main_cross_flag <- gp$c_true %in% c(1, 4) & gp$fit_mode == "estimated_c"
  gp$analysis_role <- "confirmatory"
  columns <- c(
    "model_family", "scenario_id", "fit_mode", "snr", "c_true",
    "main_cross_flag", "analysis_role"
  )
  rbind(pwo[, columns], gp[, columns])
}

wcrit29_driver_paths <- function(cfg) {
  if (nzchar(cfg$out_dir_override)) {
    out <- normalizePath(cfg$out_dir_override, winslash = "/", mustWork = FALSE)
  } else {
    out <- file.path(normalizePath(getwd(), winslash = "/", mustWork = FALSE),
                     "outputs", "wcrit", cfg$out_subdir)
  }
  paths <- list(
    out = out,
    config = file.path(out, "config"),
    source_snapshot = file.path(out, "config", "source_snapshot"),
    cache = file.path(out, "cache"),
    checkpoints = file.path(out, "checkpoints"),
    raw = file.path(out, "raw"),
    results = file.path(out, "results"),
    figures = file.path(out, "figures")
  )
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
  paths
}

wcrit29_driver_source_paths <- function() {
  helper_dir <- if (exists(".wcrit29_helper_dir", inherits = TRUE)) {
    get(".wcrit29_helper_dir", inherits = TRUE)
  } else {
    Sys.getenv("WCRIT29_HELPER_DIR", unset = "")
  }
  paths <- c(
    driver = .wcrit29_driver_file,
    search_core = .wcrit29_search_core,
    model_core = .wcrit29_model_core,
    wcrit_common = file.path(helper_dir, "wcrit_common.R"),
    case_study_common = file.path(helper_dir, "case_study_common.R"),
    sa_core = file.path(helper_dir, "sa_core.cpp")
  )
  if (any(!file.exists(paths))) {
    stop("Directed source closure is incomplete: ", paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
  }
  paths
}

wcrit29_driver_scientific_config <- function(cfg) {
  cfg[c(
    "experiment", "protocol_version", "profile", "m", "n_values", "reps",
    "lambdas", "snr_values", "c_values", "fit_modes", "c_fit_grid",
    "cross_gp_fit_mode", "primary_loss", "primary_regret", "master_seed",
    "rng_kind", "rng_normal_kind", "rng_sample_kind",
    "proposal_budget", "T0", "alpha", "tie_weight", "gp_jitter",
    "rank_tol", "confidence"
  )]
}

wcrit29_driver_write_config <- function(cfg, paths, method_dictionary,
                                        scenario_dictionary, source_table,
                                        protocol_hash, source_bundle_hash, run_hash) {
  scientific <- wcrit29_driver_scientific_config(cfg)
  run_manifest_path <- file.path(paths$config, "run_manifest.json")
  if (file.exists(run_manifest_path)) {
    previous <- jsonlite::read_json(run_manifest_path, simplifyVector = TRUE)
    if (!identical(as.character(previous$run_hash), as.character(run_hash))) {
      stop("Output directory contains a different protocol/source bundle", call. = FALSE)
    }
  }
  manifest <- list(
    experiment = cfg$experiment,
    protocol_hash = protocol_hash,
    source_bundle_hash = source_bundle_hash,
    run_hash = run_hash,
    scientific_config = scientific,
    operational = list(
      workers = cfg$workers, resume = cfg$resume,
      out_subdir = cfg$out_subdir, out_dir = paths$out
    )
  )
  wcrit29_driver_atomic_json(manifest, run_manifest_path)
  wcrit29_driver_atomic_json(scientific, file.path(paths$config, "protocol_canonical.json"))
  write.csv(method_dictionary, file.path(paths$config, "method_dictionary.csv"), row.names = FALSE)
  write.csv(scenario_dictionary, file.path(paths$config, "scenario_dictionary.csv"), row.names = FALSE)
  write.csv(source_table, file.path(paths$config, "source_sha256.csv"), row.names = FALSE)

  snapshot_names <- paste0(names(stats::setNames(source_table$source_name, source_table$source_name)),
                           "__", basename(source_table$source_path))
  for (i in seq_len(nrow(source_table))) {
    ok <- file.copy(
      source_table$source_path[[i]],
      file.path(paths$source_snapshot, snapshot_names[[i]]),
      overwrite = TRUE
    )
    if (!isTRUE(ok)) stop("Cannot snapshot source: ", source_table$source_path[[i]], call. = FALSE)
  }
  capture.output(sessionInfo(), file = file.path(paths$config, "session_info.txt"))
  invisible(manifest)
}

wcrit29_driver_design_row <- function(entry, bank, cfg, search_universe, model_universe,
                                      protocol_hash, source_bundle_hash) {
  internal <- entry$method
  method_id <- wcrit29_driver_method_id(internal)
  train_idx <- wcrit29_lookup_design_indices(entry$D, model_universe, method_id)
  universe_idx <- as.integer(unname(search_universe$key_to_index[wcrit29_row_keys(entry$D)]))
  phi05 <- wcrit29_score_indices(entry$H_index, "phi", 0.5, search_universe, cfg$tie_weight)
  kendall <- wcrit29_score_indices(entry$H_index, "kendall", 0.5, search_universe, cfg$tie_weight)
  hamming <- wcrit29_score_indices(entry$H_index, "hamming", 0.5, search_universe, cfg$tie_weight)
  l2 <- wcrit29_score_indices(entry$H_index, "l2", 0.5, search_universe, cfg$tie_weight)
  data.frame(
    experiment = cfg$experiment,
    protocol_hash = protocol_hash,
    source_bundle_hash = source_bundle_hash,
    rep = bank$rep,
    m = bank$m,
    n = bank$n,
    method_id = method_id,
    internal_method = internal,
    method_label = entry$method_label,
    lambda = entry$lambda,
    objective = entry$objective,
    l2_variant = if (internal == "L2") "component_position" else NA_character_,
    strict_foldover = wcrit29_is_strict_foldover(entry$D),
    unique_run_n = length(unique(wcrit29_row_keys(entry$D))),
    reversal_orbit_n = length(unique(search_universe$orbit_id[universe_idx])),
    pwo_rank = qr(model_universe$F[train_idx, , drop = FALSE], tol = cfg$rank_tol)$rank,
    initial_seed = bank$init_seed,
    initial_actual_seed = bank$init_actual_seed,
    initial_sampling_attempt = bank$init_sampling_attempt,
    proposal_tape_seed = bank$tape_seed,
    initial_H_sha256 = entry$initial_H_sha256,
    proposal_tape_sha256 = entry$proposal_tape_sha256,
    design_sha256 = entry$design_sha256,
    design_set_sha256 = entry$design_set_sha256,
    proposals = entry$proposals,
    proposal_budget_requested = bank$proposal_budget,
    proposal_budget_exhausted = if (internal == "RandomFoldover") NA else entry$proposal_budget_exhausted,
    initial_objective_evals = entry$initial_objective_evals,
    candidate_objective_evals = entry$candidate_objective_evals,
    invalid_proposals = entry$invalid_proposals,
    rank_invalid_proposals = entry$rank_invalid_proposals,
    accepted_proposals = entry$accepted_proposals,
    improving_accepts = entry$improving_accepts,
    proposal_partition_pass = if (internal == "RandomFoldover") {
      entry$proposals == 0L
    } else {
      entry$candidate_objective_evals + entry$invalid_proposals == entry$proposals
    },
    budget_pass = if (internal == "RandomFoldover") TRUE else {
      entry$proposals == bank$proposal_budget && isTRUE(entry$proposal_budget_exhausted)
    },
    k_min = phi05$k_min,
    k_m2 = phi05$k_m2,
    A = phi05$A,
    B = phi05$B,
    phi_05 = phi05$phi,
    kendall_min_norm = kendall$primary,
    hamming_min_norm = hamming$primary,
    l2_min_norm = l2$primary,
    stringsAsFactors = FALSE
  )
}

wcrit29_driver_seed_rows <- function(bank, pwo_truth, gp_truths) {
  rows <- list(
    data.frame(
      seed_key = sprintf("initial_half|n=%03d|rep=%03d", bank$n, bank$rep),
      seed_role = "initial_half", rep = bank$rep, n = bank$n,
      c_true = NA_real_, seed = bank$init_actual_seed,
      base_seed = bank$init_seed,
      initial_actual_seed = bank$init_actual_seed,
      initial_sampling_attempt = bank$init_sampling_attempt,
      shared_scope = "all methods within (n,rep)", stringsAsFactors = FALSE
    ),
    data.frame(
      seed_key = sprintf("proposal_tape|n=%03d|rep=%03d", bank$n, bank$rep),
      seed_role = "proposal_tape", rep = bank$rep, n = bank$n,
      c_true = NA_real_, seed = bank$tape_seed,
      base_seed = bank$tape_seed,
      initial_actual_seed = NA_integer_,
      initial_sampling_attempt = NA_integer_,
      shared_scope = "all optimized methods within (n,rep)", stringsAsFactors = FALSE
    ),
    data.frame(
      seed_key = sprintf("pwo_beta|rep=%03d", bank$rep),
      seed_role = "pwo_beta", rep = bank$rep, n = NA_integer_,
      c_true = NA_real_, seed = pwo_truth$beta_seed,
      base_seed = pwo_truth$beta_seed,
      initial_actual_seed = NA_integer_,
      initial_sampling_attempt = NA_integer_,
      shared_scope = "all n, methods, and SNR within rep", stringsAsFactors = FALSE
    ),
    data.frame(
      seed_key = sprintf("pwo_noise|rep=%03d", bank$rep),
      seed_role = "pwo_noise", rep = bank$rep, n = NA_integer_,
      c_true = NA_real_, seed = pwo_truth$noise_seed,
      base_seed = pwo_truth$noise_seed,
      initial_actual_seed = NA_integer_,
      initial_sampling_attempt = NA_integer_,
      shared_scope = "one standard-noise path scaled across SNR and shared over n/method", stringsAsFactors = FALSE
    )
  )
  for (truth in gp_truths) {
    rows[[length(rows) + 1L]] <- data.frame(
      seed_key = sprintf("gp_path|c=%s|rep=%03d", wcrit29_num_tag(truth$c_true), bank$rep),
      seed_role = "gp_path", rep = bank$rep, n = NA_integer_,
      c_true = truth$c_true, seed = truth$path_seed,
      base_seed = truth$path_seed,
      initial_actual_seed = NA_integer_,
      initial_sampling_attempt = NA_integer_,
      shared_scope = "all n, methods, and fit modes within (c,rep)", stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

wcrit29_driver_task <- function(task, cfg, search_universe, model_universe,
                                gp_kernel_cache, protocol_hash, source_bundle_hash) {
  rep_id <- as.integer(task$rep)
  n <- as.integer(task$n)
  started <- proc.time()[["elapsed"]]
  bank <- wcrit29_build_design_bank(
    n = n, rep_id = rep_id, master_seed = cfg$master_seed,
    lambdas = cfg$lambdas, proposal_budget = cfg$proposal_budget,
    universe = search_universe, T0 = cfg$T0, alpha = cfg$alpha,
    tie_weight = cfg$tie_weight
  )
  design_rows <- do.call(rbind, lapply(bank$entries, function(entry) {
    wcrit29_driver_design_row(
      entry, bank, cfg, search_universe, model_universe,
      protocol_hash, source_bundle_hash
    )
  }))
  designs <- lapply(bank$entries, `[[`, "D")
  names(designs) <- wcrit29_driver_method_id(names(bank$entries))
  common_test_idx <- wcrit29_common_test_indices(designs, model_universe)
  common_test_sha256 <- wcrit29_sha256_text(as.character(common_test_idx))

  pwo_truth <- wcrit29_generate_pwo_truth(
    model_universe, rep_id = rep_id, master_seed = cfg$master_seed,
    snr_values = cfg$snr_values
  )
  pwo_rows <- list()
  row_id <- 0L
  for (internal in names(bank$entries)) {
    method_id <- wcrit29_driver_method_id(internal)
    D <- bank$entries[[internal]]$D
    for (snr in cfg$snr_values) {
      row_id <- row_id + 1L
      fit_start <- proc.time()[["elapsed"]]
      metrics <- wcrit29_evaluate_pwo_design(
        D, pwo_truth, snr, model_universe, common_test_idx,
        confidence = cfg$confidence, rank_tol = cfg$rank_tol
      )
      tag <- wcrit29_num_tag(snr)
      metadata <- data.frame(
        experiment = cfg$experiment, protocol_hash = protocol_hash,
        source_bundle_hash = source_bundle_hash,
        rep = rep_id, m = cfg$m, n = n,
        method_id = method_id, internal_method = internal,
        beta_seed = pwo_truth$beta_seed, noise_seed = pwo_truth$noise_seed,
        signal_sd = pwo_truth$signal_sd,
        sigma_noise = as.numeric(pwo_truth$sigma[[tag]]),
        target_snr = snr,
        empirical_standard_noise_sd = stats::sd(pwo_truth$shared_standard_noise),
        common_test_sha256 = common_test_sha256,
        elapsed_fit_predict_sec = proc.time()[["elapsed"]] - fit_start,
        stringsAsFactors = FALSE
      )
      pwo_rows[[row_id]] <- cbind(metadata, metrics)
    }
  }
  pwo_raw <- do.call(rbind, pwo_rows)

  gp_truths <- lapply(cfg$c_values, function(c_value) {
    wcrit29_generate_gp_truth(
      model_universe, c_value, rep_id, cfg$master_seed,
      kernel_cache = gp_kernel_cache[[wcrit29_num_tag(c_value)]]
    )
  })
  names(gp_truths) <- vapply(cfg$c_values, wcrit29_num_tag, character(1L))
  gp_rows <- list()
  row_id <- 0L
  for (truth in gp_truths) {
    for (internal in names(bank$entries)) {
      method_id <- wcrit29_driver_method_id(internal)
      D <- bank$entries[[internal]]$D
      for (fit_mode in cfg$fit_modes) {
        row_id <- row_id + 1L
        fit_start <- proc.time()[["elapsed"]]
        metrics <- wcrit29_evaluate_gp_design(
          D, truth, model_universe, common_test_idx,
          fit_mode = fit_mode, c_grid = cfg$c_fit_grid,
          confidence = cfg$confidence, jitter = cfg$gp_jitter,
          rank_tol = cfg$rank_tol
        )
        metadata <- data.frame(
          experiment = cfg$experiment, protocol_hash = protocol_hash,
          source_bundle_hash = source_bundle_hash,
          rep = rep_id, m = cfg$m, n = n,
          method_id = method_id, internal_method = internal,
          common_test_sha256 = common_test_sha256,
          elapsed_fit_predict_sec = proc.time()[["elapsed"]] - fit_start,
          stringsAsFactors = FALSE
        )
        gp_rows[[row_id]] <- cbind(metadata, metrics)
      }
    }
  }
  gp_raw <- do.call(rbind, gp_rows)

  list(
    schema = "wcrit29-task-checkpoint-v1",
    task_id = sprintf("n%03d_rep%03d", n, rep_id),
    rep = rep_id, n = n,
    protocol_hash = protocol_hash,
    source_bundle_hash = source_bundle_hash,
    design_bank = bank,
    design_raw = design_rows,
    pwo_raw = pwo_raw,
    gp_raw = gp_raw,
    seed_rows = wcrit29_driver_seed_rows(bank, pwo_truth, gp_truths),
    common_test_indices = common_test_idx,
    common_test_sha256 = common_test_sha256,
    elapsed_sec = proc.time()[["elapsed"]] - started
  )
}

wcrit29_driver_checkpoint_path <- function(task, paths) {
  file.path(paths$checkpoints, sprintf("n%03d_rep%03d.rds", task$n, task$rep))
}

wcrit29_driver_load_or_run <- function(task, cfg, paths, search_universe,
                                       model_universe, gp_kernel_cache,
                                       protocol_hash, source_bundle_hash) {
  checkpoint <- wcrit29_driver_checkpoint_path(task, paths)
  if (file.exists(checkpoint)) {
    if (!isTRUE(cfg$resume)) stop("Checkpoint exists but resume is disabled: ", checkpoint, call. = FALSE)
    object <- try(readRDS(checkpoint), silent = TRUE)
    if (inherits(object, "try-error")) stop("Cannot read checkpoint: ", checkpoint, call. = FALSE)
    expected_id <- sprintf("n%03d_rep%03d", task$n, task$rep)
    if (!identical(object$schema, "wcrit29-task-checkpoint-v1") ||
        !identical(object$task_id, expected_id) ||
        !identical(object$protocol_hash, protocol_hash) ||
        !identical(object$source_bundle_hash, source_bundle_hash)) {
      stop("Checkpoint protocol/source mismatch: ", checkpoint, call. = FALSE)
    }
    object$checkpoint_status <- "reused"
    object$checkpoint_path <- checkpoint
    return(object)
  }
  object <- wcrit29_driver_task(
    task, cfg, search_universe, model_universe, gp_kernel_cache,
    protocol_hash, source_bundle_hash
  )
  wcrit29_driver_atomic_rds(object, checkpoint)
  object$checkpoint_status <- "computed"
  object$checkpoint_path <- checkpoint
  object
}

wcrit29_driver_cross_loss <- function(pwo_raw, gp_raw, cfg) {
  pwo <- pwo_raw[pwo_raw$scenario_id %in% c("pwo_snr2", "pwo_snr5"), , drop = FALSE]
  pwo_out <- data.frame(
    experiment = cfg$experiment,
    protocol_hash = pwo$protocol_hash,
    rep = pwo$rep, m = pwo$m, n = pwo$n,
    scenario_id = pwo$scenario_id,
    model_family = "pwo", source_fit_mode = "ols",
    source_c_true = NA_real_, method_id = pwo$method_id,
    loss_name = "nrmse_common", absolute_loss = pwo$nrmse_common,
    stringsAsFactors = FALSE
  )
  gp <- gp_raw[
    gp_raw$fit_mode == cfg$cross_gp_fit_mode & gp_raw$c_true %in% c(1, 4),
    , drop = FALSE
  ]
  gp_out <- data.frame(
    experiment = cfg$experiment,
    protocol_hash = gp$protocol_hash,
    rep = gp$rep, m = gp$m, n = gp$n,
    scenario_id = gp$scenario_id,
    model_family = "mallows_gp", source_fit_mode = gp$fit_mode,
    source_c_true = gp$c_true, method_id = gp$method_id,
    loss_name = "nrmse_common", absolute_loss = gp$nrmse_common,
    stringsAsFactors = FALSE
  )
  rbind(pwo_out, gp_out)
}

wcrit29_driver_unique_key <- function(data, columns) {
  do.call(paste, c(data[columns], sep = "|"))
}

wcrit29_driver_integrity <- function(cfg, method_dictionary, scenario_dictionary,
                                     source_table, task_results, design_raw,
                                     pwo_raw, gp_raw, cross_raw, seed_ledger,
                                     checkpoint_manifest, search_universe,
                                     model_universe, protocol_hash,
                                     source_bundle_hash) {
  gates <- list()
  add_gate <- function(gate_id, pass, observed, expected, scope = "formal_raw",
                       severity = "hard", tolerance = NA_character_, detail = NA_character_) {
    gates[[length(gates) + 1L]] <<- data.frame(
      gate_id = gate_id, severity = severity, scope = scope,
      observed = as.character(observed), expected = as.character(expected),
      tolerance = as.character(tolerance), pass = isTRUE(pass),
      detail = as.character(detail), protocol_hash = protocol_hash,
      stringsAsFactors = FALSE
    )
  }
  methods <- method_dictionary$method_id
  n_tasks <- cfg$reps * length(cfg$n_values)
  n_methods <- length(methods)
  expected_design <- n_tasks * n_methods
  expected_pwo <- expected_design * length(cfg$snr_values)
  expected_gp <- expected_design * length(cfg$c_values) * length(cfg$fit_modes)
  expected_cross <- expected_design * 4L

  add_gate("search_and_model_universe_identity",
           search_universe$size == 720L && model_universe$size == 720L &&
             identical(search_universe$keys, model_universe$keys) &&
             search_universe$q == 15L && model_universe$q == 15L,
           paste(search_universe$size, model_universe$size, search_universe$q, sep = "/"),
           "720/720/15")
  add_gate("frozen_rng_kind", identical(as.character(RNGkind()), c(
             cfg$rng_kind, cfg$rng_normal_kind, cfg$rng_sample_kind
           )), paste(RNGkind(), collapse = "/"),
           paste(c(cfg$rng_kind, cfg$rng_normal_kind, cfg$rng_sample_kind), collapse = "/"))
  add_gate("method_dictionary", nrow(method_dictionary) == n_methods &&
             !anyDuplicated(method_dictionary$method_id) && sum(method_dictionary$confirmatory7) == 4L,
           paste(nrow(method_dictionary), sum(method_dictionary$confirmatory7), sep = "/"),
           paste(n_methods, 4L, sep = "/"))
  add_gate("scenario_main_cross_flags",
           setequal(
             paste(scenario_dictionary$scenario_id[scenario_dictionary$main_cross_flag],
                   scenario_dictionary$fit_mode[scenario_dictionary$main_cross_flag], sep = "/"),
             c("pwo_snr2/ols", "pwo_snr5/ols", "gp_c1/estimated_c", "gp_c4/estimated_c")
           ),
           paste(paste(scenario_dictionary$scenario_id[scenario_dictionary$main_cross_flag],
                       scenario_dictionary$fit_mode[scenario_dictionary$main_cross_flag], sep = "/"), collapse = ";"),
           "pwo_snr2/ols;pwo_snr5/ols;gp_c1/estimated_c;gp_c4/estimated_c")
  add_gate("directed_source_hashes", nrow(source_table) == 6L &&
             all(nchar(source_table$sha256) == 64L) &&
             length(unique(source_table$sha256)) == nrow(source_table),
           nrow(source_table), 6L)
  add_gate("task_checkpoint_count", length(task_results) == n_tasks &&
             nrow(checkpoint_manifest) == n_tasks &&
             all(checkpoint_manifest$protocol_hash == protocol_hash) &&
             all(checkpoint_manifest$source_bundle_hash == source_bundle_hash),
           nrow(checkpoint_manifest), n_tasks)
  add_gate("design_row_count", nrow(design_raw) == expected_design, nrow(design_raw), expected_design)
  add_gate("design_unique_key",
           !anyDuplicated(wcrit29_driver_unique_key(design_raw, c("rep", "n", "method_id"))),
           length(unique(wcrit29_driver_unique_key(design_raw, c("rep", "n", "method_id")))),
           expected_design)
  design_groups <- split(design_raw, interaction(design_raw$rep, design_raw$n, drop = TRUE))
  add_gate("complete_method_bank_per_task",
           all(vapply(design_groups, function(x) setequal(x$method_id, methods), logical(1L))),
           sum(vapply(design_groups, function(x) setequal(x$method_id, methods), logical(1L))),
           n_tasks)
  add_gate("all_designs_strict_foldover", all(design_raw$strict_foldover),
           sum(design_raw$strict_foldover), expected_design)
  add_gate("all_designs_unique_and_complete_orbits",
           all(design_raw$unique_run_n == design_raw$n) &&
             all(design_raw$reversal_orbit_n == design_raw$n / 2),
           sum(design_raw$unique_run_n == design_raw$n &
                 design_raw$reversal_orbit_n == design_raw$n / 2), expected_design)
  add_gate("all_designs_pwo_rank16", all(design_raw$pwo_rank == 16L),
           sum(design_raw$pwo_rank == 16L), expected_design)
  add_gate("shared_initial_half_per_task",
           all(vapply(design_groups, function(x) length(unique(x$initial_H_sha256)) == 1L, logical(1L))),
           sum(vapply(design_groups, function(x) length(unique(x$initial_H_sha256)) == 1L, logical(1L))),
           n_tasks)
  initial_metadata_ok <- vapply(design_groups, function(x) {
    length(unique(x$initial_seed)) == 1L &&
      length(unique(x$initial_actual_seed)) == 1L &&
      length(unique(x$initial_sampling_attempt)) == 1L &&
      all(is.finite(x$initial_seed)) && all(is.finite(x$initial_actual_seed)) &&
      all(x$initial_sampling_attempt >= 1L)
  }, logical(1L))
  add_gate("admissible_initial_resample_metadata",
           all(initial_metadata_ok), sum(initial_metadata_ok), n_tasks)
  add_gate("shared_proposal_tape_per_task",
           all(vapply(design_groups, function(x) length(unique(x$proposal_tape_sha256)) == 1L, logical(1L))),
           sum(vapply(design_groups, function(x) length(unique(x$proposal_tape_sha256)) == 1L, logical(1L))),
           n_tasks)
  add_gate("proposal_budget_and_partition", all(design_raw$budget_pass) &&
             all(design_raw$proposal_partition_pass),
           sum(design_raw$budget_pass & design_raw$proposal_partition_pass), expected_design)
  rank_invalid_ok <- design_raw$rank_invalid_proposals >= 0L &
    design_raw$rank_invalid_proposals <= design_raw$invalid_proposals &
    (design_raw$method_id != "random_foldover" | design_raw$rank_invalid_proposals == 0L)
  add_gate("rank_invalid_proposal_accounting", all(rank_invalid_ok),
           sum(rank_invalid_ok), expected_design)
  add_gate("l2_variant_is_component_position",
           all(design_raw$l2_variant[design_raw$method_id == "SA_l2_position_foldover"] == "component_position"),
           paste(unique(design_raw$l2_variant[design_raw$method_id == "SA_l2_position_foldover"]), collapse = ";"),
           "component_position")

  add_gate("pwo_row_count", nrow(pwo_raw) == expected_pwo, nrow(pwo_raw), expected_pwo)
  add_gate("pwo_unique_key", !anyDuplicated(wcrit29_driver_unique_key(
    pwo_raw, c("rep", "n", "method_id", "scenario_id")
  )), length(unique(wcrit29_driver_unique_key(pwo_raw, c("rep", "n", "method_id", "scenario_id")))), expected_pwo)
  add_gate("pwo_all_fits_success", all(pwo_raw$fit_success), sum(pwo_raw$fit_success), expected_pwo)
  add_gate("pwo_rank16", all(pwo_raw$pwo_rank == 16L), sum(pwo_raw$pwo_rank == 16L), expected_pwo)
  add_gate("pwo_primary_loss_finite", all(is.finite(pwo_raw$nrmse_common)),
           sum(is.finite(pwo_raw$nrmse_common)), expected_pwo)
  pwo_seed_groups <- split(pwo_raw, pwo_raw$rep)
  add_gate("pwo_truth_paired_across_n_method_snr",
           all(vapply(pwo_seed_groups, function(x) {
             length(unique(x$beta_seed)) == 1L && length(unique(x$noise_seed)) == 1L
           }, logical(1L))),
           sum(vapply(pwo_seed_groups, function(x) {
             length(unique(x$beta_seed)) == 1L && length(unique(x$noise_seed)) == 1L
           }, logical(1L))), cfg$reps)
  add_gate("pwo_target_snr_definition",
           all(abs(pwo_raw$signal_sd / pwo_raw$sigma_noise - pwo_raw$target_snr) <= 1e-12),
           max(abs(pwo_raw$signal_sd / pwo_raw$sigma_noise - pwo_raw$target_snr)), 0,
           tolerance = "1e-12")

  add_gate("gp_row_count", nrow(gp_raw) == expected_gp, nrow(gp_raw), expected_gp)
  add_gate("gp_unique_key", !anyDuplicated(wcrit29_driver_unique_key(
    gp_raw, c("rep", "n", "method_id", "scenario_id", "fit_mode")
  )), length(unique(wcrit29_driver_unique_key(
    gp_raw, c("rep", "n", "method_id", "scenario_id", "fit_mode")
  ))), expected_gp)
  add_gate("gp_all_fits_success", all(gp_raw$fit_success), sum(gp_raw$fit_success), expected_gp)
  add_gate("gp_primary_loss_finite", all(is.finite(gp_raw$nrmse_common)),
           sum(is.finite(gp_raw$nrmse_common)), expected_gp)
  add_gate("gp_paths_not_empirically_standardized",
           all(!gp_raw$path_empirically_centered & !gp_raw$path_empirically_standardized),
           sum(!gp_raw$path_empirically_centered & !gp_raw$path_empirically_standardized), expected_gp)
  add_gate("gp_fit_mode_is_estimated_c",
           nrow(gp_raw) > 0L && all(gp_raw$fit_mode == "estimated_c"),
           paste(unique(gp_raw$fit_mode), collapse = ";"), "estimated_c")
  gp_seed_groups <- split(gp_raw, interaction(gp_raw$rep, gp_raw$c_true, drop = TRUE))
  add_gate("gp_truth_paired_across_n_method_mode",
           all(vapply(gp_seed_groups, function(x) length(unique(x$path_seed)) == 1L, logical(1L))),
           sum(vapply(gp_seed_groups, function(x) length(unique(x$path_seed)) == 1L, logical(1L))),
           cfg$reps * length(cfg$c_values))
  add_gate("gp_jitter_bound", all(gp_raw$gp_jitter <= 1e-8),
           max(gp_raw$gp_jitter), "<=1e-8")

  common_groups <- split(
    rbind(
      data.frame(rep = pwo_raw$rep, n = pwo_raw$n, hash = pwo_raw$common_test_sha256,
                 common_n = pwo_raw$common_test_n),
      data.frame(rep = gp_raw$rep, n = gp_raw$n, hash = gp_raw$common_test_sha256,
                 common_n = gp_raw$common_test_n)
    ),
    interaction(c(pwo_raw$rep, gp_raw$rep), c(pwo_raw$n, gp_raw$n), drop = TRUE)
  )
  add_gate("one_positive_common_test_domain_per_task",
           all(vapply(common_groups, function(x) {
             length(unique(x$hash)) == 1L && length(unique(x$common_n)) == 1L && x$common_n[[1L]] > 0L
           }, logical(1L))),
           sum(vapply(common_groups, function(x) {
             length(unique(x$hash)) == 1L && length(unique(x$common_n)) == 1L && x$common_n[[1L]] > 0L
           }, logical(1L))), n_tasks)

  add_gate("cross_model_row_count", nrow(cross_raw) == expected_cross,
           nrow(cross_raw), expected_cross)
  add_gate("cross_model_unique_key", !anyDuplicated(wcrit29_driver_unique_key(
    cross_raw, c("rep", "n", "method_id", "scenario_id")
  )), length(unique(wcrit29_driver_unique_key(
    cross_raw, c("rep", "n", "method_id", "scenario_id")
  ))), expected_cross)
  add_gate("cross_model_exact_scenarios_and_estimated_gp",
           setequal(unique(cross_raw$scenario_id), c("pwo_snr2", "pwo_snr5", "gp_c1", "gp_c4")) &&
             all(cross_raw$source_fit_mode[cross_raw$model_family == "mallows_gp"] == "estimated_c") &&
             !any(cross_raw$source_c_true %in% 8, na.rm = TRUE),
           paste(sort(unique(paste(cross_raw$scenario_id, cross_raw$source_fit_mode, sep = "/"))), collapse = ";"),
           "gp_c1/estimated_c;gp_c4/estimated_c;pwo_snr2/ols;pwo_snr5/ols")
  add_gate("cross_model_loss_finite", all(is.finite(cross_raw$absolute_loss)),
           sum(is.finite(cross_raw$absolute_loss)), expected_cross)

  expected_seed_rows <- n_tasks * 2L + cfg$reps * 2L + cfg$reps * length(cfg$c_values)
  add_gate("seed_ledger_cardinality", nrow(seed_ledger) == expected_seed_rows &&
             !anyDuplicated(seed_ledger$seed_key), nrow(seed_ledger), expected_seed_rows)
  initial_seed_ledger <- seed_ledger[seed_ledger$seed_role == "initial_half", , drop = FALSE]
  add_gate("seed_ledger_initial_resample_metadata",
           nrow(initial_seed_ledger) == n_tasks &&
             all(initial_seed_ledger$seed == initial_seed_ledger$initial_actual_seed) &&
             all(is.finite(initial_seed_ledger$base_seed)) &&
             all(initial_seed_ledger$initial_sampling_attempt >= 1L),
           nrow(initial_seed_ledger), n_tasks)
  endpoint <- wcrit29_endpoint_comparator_selftest()
  add_gate("endpoint_tie_break_selftest", all(endpoint), sum(endpoint), length(endpoint))
  if (cfg$profile == "formal") {
    add_gate("formal_frozen_protocol",
             identical(cfg$n_values, c(48L, 60L)) && cfg$reps == 50L &&
               identical(cfg$lambdas, 0.5) &&
               cfg$proposal_budget == 6000L && cfg$master_seed == 20260831L,
             "see run_manifest", "frozen 2026-08-31 protocol")
  }
  do.call(rbind, gates)
}

wcrit29_driver_main <- function() {
  wcrit29_driver_require("digest")
  wcrit29_driver_require("jsonlite")
  cfg <- wcrit29_driver_config()
  wcrit29_driver_validate_config(cfg)
  RNGkind(
    kind = cfg$rng_kind,
    normal.kind = cfg$rng_normal_kind,
    sample.kind = cfg$rng_sample_kind
  )
  paths <- wcrit29_driver_paths(cfg)
  method_dictionary <- wcrit29_driver_method_dictionary(cfg$lambdas)
  scenario_dictionary <- wcrit29_driver_scenario_dictionary(cfg)
  source_paths <- wcrit29_driver_source_paths()
  source_table <- data.frame(
    source_name = names(source_paths),
    source_path = unname(source_paths),
    sha256 = vapply(source_paths, wcrit29_driver_sha256_file, character(1L)),
    stringsAsFactors = FALSE
  )
  protocol_hash <- wcrit29_driver_hash_object(wcrit29_driver_scientific_config(cfg))
  source_bundle_hash <- wcrit29_driver_hash_object(
    stats::setNames(source_table$sha256, source_table$source_name)
  )
  run_hash <- wcrit29_driver_hash_object(list(
    protocol_hash = protocol_hash, source_bundle_hash = source_bundle_hash
  ))
  wcrit29_driver_write_config(
    cfg, paths, method_dictionary, scenario_dictionary, source_table,
    protocol_hash, source_bundle_hash, run_hash
  )

  message("[29] building S_6 universes and GP kernel cache")
  search_universe <- wcrit29_build_universe(cfg$m)
  model_universe <- wcrit29_model_build_universe(cfg$m)
  if (!identical(search_universe$keys, model_universe$keys)) {
    stop("Search and model universes use different permutation orders", call. = FALSE)
  }
  gp_kernel_cache <- lapply(cfg$c_values, function(c_value) {
    wcrit29_prepare_gp_kernel(model_universe, c_value)
  })
  names(gp_kernel_cache) <- vapply(cfg$c_values, wcrit29_num_tag, character(1L))
  saveRDS(
    list(
      universe_sha256 = search_universe$universe_sha256,
      m = cfg$m, size = search_universe$size, q = search_universe$q,
      c_values = cfg$c_values,
      gp_truth_kernel_jitter = vapply(gp_kernel_cache, `[[`, numeric(1L), "jitter")
    ),
    file.path(paths$cache, "universe_kernel_manifest.rds")
  )

  tasks <- expand.grid(
    n = cfg$n_values, rep = seq_len(cfg$reps),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  tasks <- tasks[order(tasks$rep, tasks$n), , drop = FALSE]
  tasks$task_id <- sprintf("n%03d_rep%03d", tasks$n, tasks$rep)
  workers <- min(cfg$workers, nrow(tasks))
  message(sprintf(
    "[29] profile=%s; tasks=%d; methods=%d; workers=%d; proposals=%d",
    cfg$profile, nrow(tasks), nrow(method_dictionary), workers, cfg$proposal_budget
  ))

  run_one <- function(i) {
    task <- tasks[i, , drop = FALSE]
    tryCatch({
      value <- wcrit29_driver_load_or_run(
        task, cfg, paths, search_universe, model_universe,
        gp_kernel_cache, protocol_hash, source_bundle_hash
      )
      list(ok = TRUE, task_id = task$task_id, value = value, error = NA_character_)
    }, error = function(e) {
      list(ok = FALSE, task_id = task$task_id, value = NULL, error = conditionMessage(e))
    })
  }
  if (workers > 1L && .Platform$OS.type != "windows") {
    task_wrapped <- parallel::mclapply(
      seq_len(nrow(tasks)), run_one, mc.cores = workers,
      mc.preschedule = FALSE, mc.set.seed = FALSE
    )
  } else {
    task_wrapped <- lapply(seq_len(nrow(tasks)), run_one)
  }
  task_status <- data.frame(
    task_id = vapply(task_wrapped, `[[`, character(1L), "task_id"),
    ok = vapply(task_wrapped, `[[`, logical(1L), "ok"),
    error = vapply(task_wrapped, `[[`, character(1L), "error"),
    stringsAsFactors = FALSE
  )
  write.csv(task_status, file.path(paths$results, "task_status.csv"), row.names = FALSE)
  if (any(!task_status$ok)) {
    stop("Experiment 29 task failed: ", task_status$error[which(!task_status$ok)[1L]], call. = FALSE)
  }
  task_results <- lapply(task_wrapped, `[[`, "value")

  design_raw <- do.call(rbind, lapply(task_results, `[[`, "design_raw"))
  pwo_raw <- do.call(rbind, lapply(task_results, `[[`, "pwo_raw"))
  gp_raw <- do.call(rbind, lapply(task_results, `[[`, "gp_raw"))
  cross_raw <- wcrit29_driver_cross_loss(pwo_raw, gp_raw, cfg)
  seed_all <- do.call(rbind, lapply(task_results, `[[`, "seed_rows"))
  seed_split <- split(seed_all, seed_all$seed_key)
  seed_consistent <- vapply(seed_split, function(x) length(unique(x$seed)) == 1L, logical(1L))
  if (!all(seed_consistent)) stop("A seed key maps to multiple seed values", call. = FALSE)
  seed_ledger <- do.call(rbind, lapply(seed_split, function(x) x[1L, , drop = FALSE]))
  seed_ledger <- seed_ledger[order(seed_ledger$seed_role, seed_ledger$seed_key), , drop = FALSE]
  checkpoint_manifest <- data.frame(
    task_id = vapply(task_results, `[[`, character(1L), "task_id"),
    rep = vapply(task_results, `[[`, integer(1L), "rep"),
    n = vapply(task_results, `[[`, integer(1L), "n"),
    status = vapply(task_results, `[[`, character(1L), "checkpoint_status"),
    checkpoint = vapply(task_results, `[[`, character(1L), "checkpoint_path"),
    checkpoint_sha256 = vapply(
      task_results, function(x) wcrit29_driver_sha256_file(x$checkpoint_path), character(1L)
    ),
    elapsed_sec = vapply(task_results, `[[`, numeric(1L), "elapsed_sec"),
    protocol_hash = protocol_hash,
    source_bundle_hash = source_bundle_hash,
    stringsAsFactors = FALSE
  )

  write.csv(design_raw, file.path(paths$raw, "design_metrics.csv"), row.names = FALSE)
  write.csv(pwo_raw, file.path(paths$raw, "pwo_prediction_raw.csv"), row.names = FALSE)
  write.csv(gp_raw, file.path(paths$raw, "gp_prediction_raw.csv"), row.names = FALSE)
  write.csv(cross_raw, file.path(paths$raw, "cross_model_loss_raw.csv"), row.names = FALSE)
  write.csv(seed_ledger, file.path(paths$config, "seed_ledger.csv"), row.names = FALSE)
  write.csv(checkpoint_manifest, file.path(paths$results, "checkpoint_manifest.csv"), row.names = FALSE)

  integrity <- wcrit29_driver_integrity(
    cfg, method_dictionary, scenario_dictionary, source_table,
    task_results, design_raw, pwo_raw, gp_raw, cross_raw, seed_ledger,
    checkpoint_manifest, search_universe, model_universe,
    protocol_hash, source_bundle_hash
  )
  write.csv(integrity, file.path(paths$results, "integrity_audit.csv"), row.names = FALSE)
  hard_pass <- all(integrity$pass[integrity$severity == "hard"])
  completion <- list(
    status = if (hard_pass) "complete" else "integrity_failed",
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    profile = cfg$profile,
    protocol_hash = protocol_hash,
    source_bundle_hash = source_bundle_hash,
    run_hash = run_hash,
    tasks = nrow(tasks), methods = nrow(method_dictionary),
    design_rows = nrow(design_raw), pwo_rows = nrow(pwo_raw),
    gp_rows = nrow(gp_raw), cross_model_rows = nrow(cross_raw),
    checkpoint_computed = sum(checkpoint_manifest$status == "computed"),
    checkpoint_reused = sum(checkpoint_manifest$status == "reused"),
    hard_integrity_pass = hard_pass,
    failed_hard_gates = integrity$gate_id[integrity$severity == "hard" & !integrity$pass],
    output_dir = paths$out
  )
  wcrit29_driver_atomic_json(completion, file.path(paths$results, "completion.json"))
  if (!hard_pass) {
    stop(
      "Experiment 29 hard integrity gates failed: ",
      paste(completion$failed_hard_gates, collapse = ", "),
      call. = FALSE
    )
  }
  message(sprintf(
    "[29] complete: design=%d, PWO=%d, GP=%d, cross=%d; %s",
    nrow(design_raw), nrow(pwo_raw), nrow(gp_raw), nrow(cross_raw), paths$out
  ))
  invisible(list(
    config = cfg, completion = completion, integrity = integrity,
    design = design_raw, pwo = pwo_raw, gp = gp_raw, cross = cross_raw
  ))
}

if (sys.nframe() == 0L) {
  wcrit29_driver_main()
}
