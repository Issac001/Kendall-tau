# Four-drug order-of-administration case study with an intercept-only
# Mallows-kernel GP used for every predictive and sequential decision.
#
# Formal protocol:
#   * n0 = 12 and six exact-domain EI additions (12 -> 18);
#   * 20 paired design seeds, all 24 component-label maps, and all three
#     leave-one-replicate-out folds;
#   * the frozen Experiment 21 initial-design bank is used only to recover the
#     five designs reported in the paper: exact FSA-KD, unrestricted Hamming,
#     unrestricted component-position L2, OofA-OA, and unrestricted SRS;
#   * a nugget-aware intercept-only Mallows GP is the sole surrogate at the
#     initial stage and throughout all six EI additions.
#
# The held-out response replicate is used only for evaluation.  EI and every
# recommendation are computed from the two discovery replicates.

.wcrit21_args <- commandArgs(trailingOnly = FALSE)
.wcrit21_file_arg <- grep("^--file=", .wcrit21_args, value = TRUE)
.wcrit21_this_file <- if (length(.wcrit21_file_arg)) {
  normalizePath(sub("^--file=", "", .wcrit21_file_arg[[1L]]),
                winslash = "/", mustWork = FALSE)
} else {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(fr) {
    if (is.null(fr$ofile)) NA_character_ else as.character(fr$ofile)
  }, character(1L))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles)) normalizePath(ofiles[[length(ofiles)]], winslash = "/",
                                    mustWork = FALSE) else NA_character_
}
.wcrit21_dir <- if (is.na(.wcrit21_this_file)) {
  normalizePath(file.path(getwd(), "code", "WeightedCrit"), winslash = "/",
                mustWork = FALSE)
} else dirname(.wcrit21_this_file)

.wcrit_public_common <- normalizePath(
  file.path(.wcrit21_dir, "..", "common"), winslash = "/", mustWork = TRUE
)
source(file.path(.wcrit_public_common, "wcrit_common.R"), local = FALSE)
source(file.path(.wcrit_public_common, "wcrit_maximin_dist.R"), local = FALSE)
source(file.path(.wcrit_public_common, "case_study_common.R"), local = FALSE)

wcrit21_stop <- function(...) stop(sprintf(...), call. = FALSE)

wcrit21_int_env <- function(name, default, lower = 1L) {
  raw <- Sys.getenv(name, unset = as.character(default))
  out <- suppressWarnings(as.integer(raw))
  if (length(out) != 1L || is.na(out) || out < lower) {
    wcrit21_stop("%s must be one integer >= %d", name, lower)
  }
  out
}

wcrit21_bool_env <- function(name, default = FALSE) {
  tolower(Sys.getenv(name, unset = if (default) "true" else "false")) %in%
    c("1", "true", "yes", "y")
}

wcrit21_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  if (!nzchar(x)) wcrit21_stop("An output name cannot be empty")
  x
}

wcrit21_config <- function() {
  profile <- tolower(trimws(Sys.getenv("WCRIT30_PROFILE", unset = "formal")))
  if (!(profile %in% c("formal", "smoke"))) {
    wcrit21_stop("WCRIT30_PROFILE must be 'formal' or 'smoke'")
  }
  smoke <- identical(profile, "smoke")
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  default_out <- paste0(
    "30_case_study_four_drug_intercept_mallows_gp_", profile, "_", stamp
  )
  out_subdir <- Sys.getenv("WCRIT30_OUT_SUBDIR", unset = default_out)
  if (!nzchar(out_subdir) || grepl("^/", out_subdir) ||
      grepl("(^|/)\\.\\.($|/)", out_subdir)) {
    wcrit21_stop("WCRIT30_OUT_SUBDIR must be a relative path below outputs/wcrit")
  }
  cfg <- list(
    experiment = "30_case_study_four_drug_intercept_mallows_gp",
    protocol_version = "2026-09-02-v1",
    profile = profile,
    m = 4L,
    n0 = 12L,
    bo_additions = if (smoke) 2L else 6L,
    design_seed_ids = seq_len(if (smoke) 1L else 20L),
    label_map_ids = seq_len(if (smoke) 2L else 24L),
    heldout_folds = seq_len(if (smoke) 1L else 3L),
    methods = c(
      "Exact_Phi_lambda050", "Unrestricted_Hamming",
      "Unrestricted_Position_L2", "OofA_OA", "SRS"
    ),
    models = "Intercept_Mallows_GP",
    bo_model = "Intercept_Mallows_GP",
    master_seed = 20260820L,
    sa_budget = 6000L,
    sa_tie_weight = 1e-4,
    sa_replace_prob = 0.25,
    sa_T0 = 0.05,
    sa_Tmin = 1e-8,
    theta_parameterization = "theta=c/choose(4,2)",
    theta_grid = if (smoke)
      c(0.0625, 0.5, 4, 64) / choose(4L, 2L) else
      exp(seq(log(0.0625), log(64), length.out = 17L)) / choose(4L, 2L),
    nugget_grid = if (smoke) c(1e-6, 1e-3, 0.1, 10, 100) else
      exp(seq(log(1e-6), log(100), length.out = 17L)),
    boundary_rate_warn_threshold = 0.10,
    rank_tol = 1e-9,
    full_rank_retry_max = 50L,
    confidence = 0.95,
    expected_data_sha256 =
      "5f7344096e82716a2d4a35e7a3b5e8c8e9a77191e903be0d9235df325f299d89",
    source_package = "OofAExp 0.1.0",
    source_archive_sha256 =
      "2ba82d25ec0e380949b32cb4551452233825e01acb205e0c14153c5582d977f0",
    source_rda_sha256 =
      "67ddfd1bbf1e7670f0008e8571f22962ffa3c1439809570855d7ea02fc894347",
    frozen_parent_experiment = "21_case_study_four_drug",
    frozen_parent_protocol_hash =
      "0948b74f6ba13c82d8302950f25212f31b611dbd606e8fa607b1e45f2381cebb",
    frozen_design_bank_sha256 =
      "abaec2df505bbada7f30a48191580ab2968fa7251515a783f7e96e976ea35187",
    frozen_design_manifest_sha256 =
      "336a2c193c55fe3d47537d60092b8694e938fbbf0eaf394a98b71e9dc92484af",
    frozen_seed_ledger_sha256 =
      "b6337752c15a8e3a25d95d071b6b49a8b566c1903c776fe6865937b47517d0da",
    frozen_label_map_manifest_sha256 =
      "18c4e46da75c0c24c274357b571947e921c340bfbd5c181d2b956b5fca8c7194",
    frozen_exact_optima_sha256 =
      "d4d6bdda52563428c40f7858752fb843f4e27df2f07bb631b75a24b92b8ec19a",
    frozen_exact_enumeration_sha256 =
      "c8ff575401563fe8623a601d08f0922f2f0e99344227f649628a76467a57bc70",
    frozen_initial_prediction_sha256 =
      "61e6d18cebd21f2986f2190ae2831c99e7f2b6d7d2d5077760e7da5ac1678332",
    frozen_initial_recommendation_sha256 =
      "1e6224728dc5f374ba8bea3ec215da17b319635058b25bf4356b3acff04da61c",
    response_definition =
      "loss = -reported response; higher reported response is treated as better",
    discovery_rule = "two replicate columns; held-out replicate never enters fit/EI",
    heldout_observation_noise_multiplier = 1,
    ei_observation_noise_multiplier = 0,
    ei_rule = paste(
      "exact over every currently unobserved member of S4; latent predictive SD;",
      "incumbent is minimum latent posterior mean among observed sequences;",
      "ties use source treatment order"
    ),
    workers = wcrit21_int_env("WCRIT30_WORKERS", if (smoke) 1L else 1L),
    resume = wcrit21_bool_env("WCRIT30_RESUME", FALSE),
    out_subdir = out_subdir
  )
  scientific_names <- setdiff(names(cfg), c("workers", "resume", "out_subdir"))
  cfg$protocol_hash <- digest::digest(cfg[scientific_names], algo = "sha256",
                                      serialize = TRUE)
  cfg
}

wcrit21_source_manifest <- function() {
  paths <- c(
    script = normalizePath(file.path(
      .wcrit21_dir, "run_four_drug_mallows_gp.R"),
                           winslash = "/", mustWork = TRUE),
    common = normalizePath(file.path(.wcrit_public_common, "wcrit_common.R"),
                           winslash = "/", mustWork = TRUE),
    maximin = normalizePath(file.path(.wcrit_public_common, "wcrit_maximin_dist.R"),
                            winslash = "/", mustWork = TRUE),
    case_common = normalizePath(file.path(.wcrit_public_common, "case_study_common.R"),
                                winslash = "/", mustWork = TRUE),
    sa_core = normalizePath(file.path(.wcrit_public_common, "sa_core.cpp"),
                            winslash = "/", mustWork = TRUE)
  )
  data.frame(
    component = names(paths), path = unname(paths),
    size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, case_file_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
}

wcrit30_read_completed <- function(path) {
  if (!file.exists(path)) wcrit21_stop("Missing frozen-parent COMPLETED file: %s", path)
  lines <- readLines(path, warn = FALSE)
  pos <- regexpr("=", lines, fixed = TRUE)
  if (!length(lines) || any(pos < 2L)) {
    wcrit21_stop("Malformed frozen-parent COMPLETED file: %s", path)
  }
  keys <- substr(lines, 1L, pos - 1L)
  vals <- substr(lines, pos + 1L, nchar(lines))
  if (anyDuplicated(keys) || any(!nzchar(keys)) || any(!nzchar(vals))) {
    wcrit21_stop("Frozen-parent COMPLETED has duplicate or empty fields")
  }
  stats::setNames(as.list(vals), keys)
}

wcrit30_load_frozen_parent <- function(parent_dir, cfg) {
  rel <- c(
    completed = "COMPLETED",
    design_bank = "designs/design_bank.rds",
    design_manifest = "designs/design_manifest.csv",
    exact_optima = "designs/exact_optima.csv",
    exact_enumeration = "designs/exact_enumeration_m4n12.csv",
    seed_ledger = "config/seed_ledger.csv",
    label_maps = "config/label_map_manifest.csv",
    initial_prediction = "results/initial_prediction_summary.csv",
    initial_recommendation = "results/initial_top3_regret_summary.csv"
  )
  paths <- file.path(parent_dir, unname(rel))
  names(paths) <- names(rel)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    wcrit21_stop("Frozen parent is incomplete: %s", paste(missing, collapse = ", "))
  }
  completed <- wcrit30_read_completed(paths[["completed"]])
  if (!identical(completed$experiment, cfg$frozen_parent_experiment) ||
      !identical(completed$profile, "formal") ||
      !identical(completed$protocol_hash, cfg$frozen_parent_protocol_hash) ||
      !identical(completed$data_sha256, cfg$expected_data_sha256)) {
    wcrit21_stop("Frozen-parent identity does not match the preregistered source")
  }
  expected <- c(
    design_bank = cfg$frozen_design_bank_sha256,
    design_manifest = cfg$frozen_design_manifest_sha256,
    exact_optima = cfg$frozen_exact_optima_sha256,
    exact_enumeration = cfg$frozen_exact_enumeration_sha256,
    seed_ledger = cfg$frozen_seed_ledger_sha256,
    label_maps = cfg$frozen_label_map_manifest_sha256,
    initial_prediction = cfg$frozen_initial_prediction_sha256,
    initial_recommendation = cfg$frozen_initial_recommendation_sha256
  )
  actual <- vapply(paths[names(expected)], case_file_sha256, character(1L))
  if (!identical(unname(actual), unname(expected))) {
    bad <- names(expected)[actual != expected]
    wcrit21_stop("Frozen-parent hash mismatch: %s", paste(bad, collapse = ", "))
  }
  entries <- readRDS(paths[["design_bank"]])
  manifest <- utils::read.csv(paths[["design_manifest"]],
                              stringsAsFactors = FALSE, check.names = FALSE)
  seed_ledger <- utils::read.csv(paths[["seed_ledger"]],
                                 stringsAsFactors = FALSE, check.names = FALSE)
  label_maps <- utils::read.csv(paths[["label_maps"]],
                                stringsAsFactors = FALSE, check.names = FALSE)
  exact <- list(
    optima = utils::read.csv(paths[["exact_optima"]],
                             stringsAsFactors = FALSE, check.names = FALSE),
    all_candidates = utils::read.csv(paths[["exact_enumeration"]],
                                     stringsAsFactors = FALSE,
                                     check.names = FALSE)
  )
  input_manifest <- data.frame(
    component = names(paths),
    relative_path = unname(rel),
    absolute_path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, case_file_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
  list(
    entries = entries, design_manifest = manifest, exact = exact,
    seed_ledger = seed_ledger, label_maps = label_maps,
    initial_prediction = utils::read.csv(
      paths[["initial_prediction"]], stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    initial_recommendation = utils::read.csv(
      paths[["initial_recommendation"]], stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    completed = completed, input_manifest = input_manifest
  )
}

wcrit21_dirs <- function(out_dir, resume) {
  if (dir.exists(out_dir) || file.exists(out_dir)) {
    if (!resume) wcrit21_stop("Refusing existing output path: %s", out_dir)
    if (file.exists(file.path(out_dir, "COMPLETED"))) {
      wcrit21_stop("The requested resume directory is already complete: %s", out_dir)
    }
  } else {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  dirs <- list(
    root = out_dir,
    config = file.path(out_dir, "config"),
    raw = file.path(out_dir, "raw"),
    results = file.path(out_dir, "results"),
    figures = file.path(out_dir, "figures"),
    designs = file.path(out_dir, "designs"),
    checkpoints = file.path(out_dir, "checkpoints"),
    audits = file.path(out_dir, "audits")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

wcrit21_atomic_save_rds <- function(object, path) {
  tmp <- tempfile(pattern = paste0(basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = TRUE)
  if (!file.rename(tmp, path)) wcrit21_stop("Atomic rename failed for %s", path)
  invisible(path)
}

wcrit21_stage <- function(name, expr, profile_path) {
  started <- Sys.time()
  status <- "ok"
  msg <- NA_character_
  value <- tryCatch(force(expr), error = function(e) {
    status <<- "error"
    msg <<- conditionMessage(e)
    NULL
  })
  case_append_stage(profile_path, name, status, started, msg)
  if (status == "error") wcrit21_stop("Stage '%s' failed: %s", name, msg)
  value
}

wcrit21_load_data <- function(cfg, data_path) {
  actual_sha <- case_file_sha256(data_path)
  if (!identical(actual_sha, cfg$expected_data_sha256)) {
    wcrit21_stop("Data SHA-256 mismatch: got %s, expected %s", actual_sha,
                 cfg$expected_data_sha256)
  }
  dat <- utils::read.csv(data_path, stringsAsFactors = FALSE,
                         check.names = FALSE)
  expected_columns <- c("treatment", paste0("step.", 1:4), paste0("y.", 1:3), "avg")
  if (!identical(names(dat), expected_columns) || nrow(dat) != 25L) {
    wcrit21_stop("Unexpected four-drug CSV schema or row count")
  }
  perm <- dat[seq_len(24L), , drop = FALSE]
  control <- dat[25L, , drop = FALSE]
  D <- as.matrix(perm[paste0("step.", 1:4)])
  storage.mode(D) <- "integer"
  case_validate_design(D, 4L, 24L, "complete four-drug permutation table")
  expected_keys <- sort(wcrit_row_keys(gtools::permutations(4L, 4L, 1:4)))
  if (!setequal(unname(wcrit_row_keys(D)), unname(expected_keys))) {
    wcrit21_stop("Rows t.1--t.24 are not exactly S4")
  }
  if (!all(as.integer(unlist(control[paste0("step.", 1:4)],
                             use.names = FALSE)) == 0L)) {
    wcrit21_stop("Row t.25 is not the simultaneous-administration control")
  }
  Y_response <- as.matrix(perm[paste0("y.", 1:3)])
  storage.mode(Y_response) <- "double"
  control_response <- as.numeric(unlist(control[paste0("y.", 1:3)],
                                        use.names = FALSE))
  recomputed_avg <- rowMeans(Y_response)
  avg_rounding_max_abs <- max(abs(recomputed_avg - as.numeric(perm$avg)))
  if (!all(is.finite(Y_response)) || avg_rounding_max_abs > 0.005 + 1e-12) {
    wcrit21_stop("Non-finite responses or an unexpected published-average discrepancy")
  }
  keys <- wcrit_row_keys(D)
  reverse_index <- match(wcrit_row_keys(D[, 4:1, drop = FALSE]), keys)
  if (anyNA(reverse_index) || !identical(reverse_index[reverse_index], seq_len(24L))) {
    wcrit21_stop("Reverse-order lookup failed its involution check")
  }
  canonical_sha <- digest::digest(
    list(treatment = perm$treatment, D = D, Y_response = Y_response,
         control_response = control_response),
    algo = "sha256", serialize = TRUE
  )
  list(
    raw = dat, permutation_rows = perm, D = D, keys = keys,
    Y_response = Y_response, Y_loss = -Y_response,
    control_response = control_response, reverse_index = reverse_index,
    avg_rounding_max_abs = avg_rounding_max_abs,
    file_sha256 = actual_sha, canonical_sha256 = canonical_sha
  )
}

wcrit21_design_manifest <- function(entries, cfg) {
  dplyr::bind_rows(lapply(entries, function(entry) {
    D <- entry$D
    met <- wcrit_foldover_design_metrics(
      D, lambda = if (is.finite(entry$target_lambda)) entry$target_lambda else 0.5,
      invalid = "na"
    )
    kd <- case_geometry_summary(D, "kendall")
    hm <- case_geometry_summary(D, "hamming")
    l2 <- case_geometry_summary(D, "l2_position")
    data.frame(
      entry_id = entry$entry_id, method = entry$method,
      design_seed_id = entry$design_seed_id, stochastic = entry$stochastic,
      m = cfg$m, n = cfg$n0, design_sha256 = entry$design_sha256,
      initial_design_sha256 = entry$initial_design_sha256,
      initial_seed = entry$initial_seed, move_seed = entry$move_seed,
      initial_rank_retry = entry$initial_rank_retry,
      design_rank_retry = entry$design_rank_retry,
      target_lambda = entry$target_lambda,
      target_geometry = entry$target_geometry,
      strict_foldover = wcrit_is_strict_foldover(D),
      phi_valid = met$phi_valid[[1L]], A = met$A[[1L]],
      B = met$B[[1L]], Phi = met$phi_lambda[[1L]],
      kendall_min = kd$minimum, kendall_mean_sq = kd$mean_sq,
      hamming_min = hm$minimum, hamming_mean_sq = hm$mean_sq,
      position_l2_min = l2$minimum, position_l2_mean_sq = l2$mean_sq,
      pwo_rank = entry$pwo_rank, pwo_columns = ncol(wcrit_pwo_matrix(D)),
      n_obj_eval = entry$n_obj_eval,
      total_obj_eval_including_rank_retries = entry$total_obj_eval,
      budget_exhausted = entry$budget_exhausted,
      inference_role = if (entry$stochastic)
        "paired_Monte_Carlo_by_design_seed" else
        "fixed_design_descriptive_repeated_only_for_alignment",
      stringsAsFactors = FALSE
    )
  }))
}

wcrit21_map_design <- function(D, map_row, data_obj) {
  label_map <- as.integer(unlist(
    map_row[c("label_1", "label_2", "label_3", "label_4")],
    use.names = FALSE
  ))
  mapped <- case_apply_label_map(D, label_map)
  idx <- match(wcrit_row_keys(mapped), data_obj$keys)
  if (anyNA(idx) || anyDuplicated(idx)) wcrit21_stop("A mapped design lookup failed")
  list(D = mapped, index = as.integer(idx), sha256 = case_matrix_sha256(mapped))
}

wcrit21_training_data <- function(observed_index, heldout_fold, data_obj) {
  train_folds <- setdiff(seq_len(3L), heldout_fold)
  D_obs <- data_obj$D[observed_index, , drop = FALSE]
  train_D <- D_obs[rep(seq_len(nrow(D_obs)), times = length(train_folds)), , drop = FALSE]
  train_y <- as.numeric(data_obj$Y_loss[observed_index, train_folds, drop = FALSE])
  list(D = train_D, y = train_y, train_folds = train_folds)
}

wcrit21_fit_model <- function(model_name, observed_index, heldout_fold,
                              data_obj, cfg) {
  if (!identical(model_name, "Intercept_Mallows_GP")) {
    wcrit21_stop("The paper runner supports only Intercept_Mallows_GP")
  }
  tr <- wcrit21_training_data(observed_index, heldout_fold, data_obj)
  fit <- case_fit_mallows_reml(
    tr$D, tr$y, trend = "intercept", theta_grid = cfg$theta_grid,
    nugget_grid = cfg$nugget_grid, rank_tol = cfg$rank_tol
  )
  if (!isTRUE(fit$success)) {
    wcrit21_stop("Model fit failed: model=%s fold=%d error=%s",
                 model_name, heldout_fold, fit$error)
  }
  fit
}

wcrit21_predict_model <- function(model_name, fit, D,
                                  observation_noise_multiplier) {
  if (!identical(model_name, "Intercept_Mallows_GP")) {
    wcrit21_stop("The paper runner supports only Intercept_Mallows_GP")
  }
  case_predict_mallows_reml(fit, D,
    include_observation_noise = observation_noise_multiplier > 0,
    observation_noise_multiplier = observation_noise_multiplier
  )
}

wcrit21_prediction_rows <- function(context, y_true, pred_observation,
                                    observed_index, cfg) {
  domains <- list(
    integrated_all24 = seq_along(y_true),
    currently_unobserved = setdiff(seq_along(y_true), observed_index)
  )
  common_sd <- stats::sd(y_true)
  dplyr::bind_rows(lapply(names(domains), function(domain) {
    idx <- domains[[domain]]
    if (!length(idx)) return(NULL)
    met <- case_prediction_metrics(
      y_true[idx],
      list(mean = pred_observation$mean[idx], variance = pred_observation$variance[idx]),
      confidence = cfg$confidence, normalization_sd = common_sd
    )
    domain_sd <- stats::sd(y_true[idx])
    data.frame(
      context, domain = domain, n_domain = length(idx),
      observation_noise_multiplier = cfg$heldout_observation_noise_multiplier,
      rmse = met$rmse,
      nrmse_domain = if (is.finite(domain_sd) && domain_sd > 0)
        met$rmse / domain_sd else NA_real_,
      nrmse_common_all24_sd = met$nrmse,
      mae = met$mae, correlation = met$correlation,
      spearman = met$spearman, nll = met$nll,
      coverage = met$coverage, interval_length = met$interval_length,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }))
}

wcrit21_recommendation_row <- function(context, recommendation_type,
                                       ranked_index, prediction, y_true,
                                       observed_index, heldout_fold,
                                       data_obj) {
  ranked_index <- as.integer(ranked_index)
  selected <- ranked_index[[1L]]
  predicted_top3 <- head(ranked_index, 3L)
  true_rank <- order(y_true, seq_along(y_true))
  true_best <- true_rank[[1L]]
  true_top3 <- head(true_rank, 3L)
  reverse_id <- data_obj$reverse_index[[selected]]
  selected_response <- -y_true[[selected]]
  reverse_response <- -y_true[[reverse_id]]
  simultaneous_response <- data_obj$control_response[[heldout_fold]]
  data.frame(
    context, recommendation_type = recommendation_type,
    selected_index = selected,
    selected_treatment = data_obj$permutation_rows$treatment[[selected]],
    selected_sequence = data_obj$keys[[selected]],
    selected_was_observed = selected %in% observed_index,
    heldout_response = selected_response,
    true_best_index = true_best,
    true_best_treatment = data_obj$permutation_rows$treatment[[true_best]],
    true_best_response = -y_true[[true_best]],
    top1_hit = selected == true_best,
    selected_in_true_top3 = selected %in% true_top3,
    top3_contains_true_best = true_best %in% predicted_top3,
    top3_overlap_fraction = length(intersect(predicted_top3, true_top3)) / 3,
    heldout_regret = y_true[[selected]] - min(y_true),
    reverse_index = reverse_id,
    reverse_treatment = data_obj$permutation_rows$treatment[[reverse_id]],
    reverse_response = reverse_response,
    response_minus_reverse = selected_response - reverse_response,
    predicted_response_minus_reverse =
      prediction$mean[[reverse_id]] - prediction$mean[[selected]],
    simultaneous_response = simultaneous_response,
    response_minus_simultaneous = selected_response - simultaneous_response,
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

wcrit21_run_trajectory <- function(entry, map_row, heldout_fold,
                                   model_name, data_obj, cfg) {
  if (!identical(model_name, "Intercept_Mallows_GP")) {
    wcrit21_stop("The paper trajectory is fixed to Intercept_Mallows_GP")
  }
  mapped <- wcrit21_map_design(entry$D, map_row, data_obj)
  observed <- mapped$index
  y_true <- data_obj$Y_loss[, heldout_fold]
  prediction_rows <- list()
  recommendation_rows <- list()
  model_rows <- list()
  acquisition_rows <- list()
  for (bo_step in 0:cfg$bo_additions) {
    fit <- wcrit21_fit_model(model_name, observed, heldout_fold, data_obj, cfg)
    pred_observation <- wcrit21_predict_model(
      model_name, fit, data_obj$D,
      observation_noise_multiplier = cfg$heldout_observation_noise_multiplier
    )
    pred_latent <- wcrit21_predict_model(
      model_name, fit, data_obj$D,
      observation_noise_multiplier = cfg$ei_observation_noise_multiplier
    )
    context <- data.frame(
      protocol_hash = cfg$protocol_hash,
      design_seed_id = entry$design_seed_id, method = entry$method,
      stochastic_method = entry$stochastic, map_id = map_row$map_id,
      map_key = map_row$map_key, heldout_fold = heldout_fold,
      discovery_folds = paste(setdiff(1:3, heldout_fold), collapse = ";"),
      model = model_name,
      model_role = "primary_prediction_and_BO",
      bo_step = bo_step, n_observed = length(observed),
      initial_design_sha256 = entry$design_sha256,
      mapped_initial_design_sha256 = mapped$sha256,
      observed_set_sha256 = digest::digest(
        paste(sort(data_obj$keys[observed]), collapse = ";"),
        algo = "sha256", serialize = FALSE
      ),
      stringsAsFactors = FALSE
    )
    prediction_rows[[length(prediction_rows) + 1L]] <-
      wcrit21_prediction_rows(context, y_true, pred_observation,
                              observed, cfg)
    surrogate_rank <- order(pred_latent$mean, seq_len(nrow(data_obj$D)))
    train_folds <- setdiff(1:3, heldout_fold)
    observed_discovery_mean <- rowMeans(
      data_obj$Y_loss[observed, train_folds, drop = FALSE]
    )
    incumbent_rank <- observed[order(observed_discovery_mean, observed)]
    recommendation_rows[[length(recommendation_rows) + 1L]] <-
      wcrit21_recommendation_row(
        context, "surrogate_all24", surrogate_rank, pred_latent,
        y_true, observed, heldout_fold, data_obj
      )
    recommendation_rows[[length(recommendation_rows) + 1L]] <-
      wcrit21_recommendation_row(
        context, "incumbent_observed", incumbent_rank, pred_latent,
        y_true, observed, heldout_fold, data_obj
      )
    model_rows[[length(model_rows) + 1L]] <- data.frame(
      context, fit_success = fit$success, trend = fit$trend,
      theta = if (!is.null(fit$theta)) fit$theta else NA_real_,
      nugget_ratio = if (!is.null(fit$nugget)) fit$nugget else NA_real_,
      process_sigma2 = fit$sigma2,
      observation_noise_variance = if (!is.null(fit$nugget))
        fit$sigma2 * fit$nugget else fit$sigma2,
      reml = if (!is.null(fit$reml)) fit$reml else NA_real_,
      sigma_floor_active = if (!is.null(fit$sigma_floor_active))
        fit$sigma_floor_active else FALSE,
      theta_boundary_hit = if (!is.null(fit$theta_boundary_hit))
        fit$theta_boundary_hit else NA,
      nugget_boundary_hit = if (!is.null(fit$nugget_boundary_hit))
        fit$nugget_boundary_hit else NA,
      trend_rank = if (!is.null(fit$trend_rank)) fit$trend_rank else fit$p,
      n_training_observations = fit$n, residual_df = fit$df,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    if (bo_step < cfg$bo_additions) {
      unobserved <- setdiff(seq_len(nrow(data_obj$D)), observed)
      best_latent_observed <- min(pred_latent$mean[observed])
      ei <- case_expected_improvement_min(
        pred_latent$mean[unobserved], pred_latent$sd[unobserved],
        best_latent_observed
      )
      max_ei <- max(ei)
      tied <- unobserved[abs(ei - max_ei) <= max(1e-14, abs(max_ei) * 1e-12)]
      chosen <- min(tied)
      before <- observed
      observed <- c(observed, chosen)
      acquisition_rows[[length(acquisition_rows) + 1L]] <- data.frame(
        context, addition = bo_step + 1L,
        chosen_index = chosen,
        chosen_treatment = data_obj$permutation_rows$treatment[[chosen]],
        chosen_sequence = data_obj$keys[[chosen]],
        expected_improvement = ei[match(chosen, unobserved)],
        ei_observation_noise_multiplier = cfg$ei_observation_noise_multiplier,
        best_latent_observed = best_latent_observed,
        candidate_count = length(unobserved), tie_count = length(tied),
        duplicate_before_add = chosen %in% before,
        observed_after_sha256 = digest::digest(
          paste(sort(data_obj$keys[observed]), collapse = ";"),
          algo = "sha256", serialize = FALSE
        ),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
  }
  list(
    prediction = dplyr::bind_rows(prediction_rows),
    recommendation = dplyr::bind_rows(recommendation_rows),
    model_fit = dplyr::bind_rows(model_rows),
    acquisition = dplyr::bind_rows(acquisition_rows)
  )
}

wcrit21_run_entry <- function(entry, label_maps, data_obj, cfg) {
  pieces <- list()
  cursor <- 0L
  for (map_i in seq_len(nrow(label_maps))) {
    map_row <- label_maps[map_i, , drop = FALSE]
    for (fold in cfg$heldout_folds) {
      for (model_name in cfg$models) {
        cursor <- cursor + 1L
        pieces[[cursor]] <- wcrit21_run_trajectory(
          entry, map_row, fold, model_name, data_obj, cfg
        )
      }
    }
  }
  list(
    prediction = dplyr::bind_rows(lapply(pieces, `[[`, "prediction")),
    recommendation = dplyr::bind_rows(lapply(pieces, `[[`, "recommendation")),
    model_fit = dplyr::bind_rows(lapply(pieces, `[[`, "model_fit")),
    acquisition = dplyr::bind_rows(lapply(pieces, `[[`, "acquisition"))
  )
}

wcrit21_entry_checkpoint <- function(entry, entry_number, label_maps,
                                     data_obj, cfg, source_bundle_sha256,
                                     checkpoint_dir) {
  path <- file.path(checkpoint_dir, sprintf("entry_%03d_%s.rds", entry_number,
                                            wcrit21_safe_name(entry$entry_id)))
  signature <- list(
    schema = "wcrit21-entry-v1", protocol_hash = cfg$protocol_hash,
    source_bundle_sha256 = source_bundle_sha256,
    data_sha256 = data_obj$file_sha256,
    canonical_data_sha256 = data_obj$canonical_sha256,
    entry_id = entry$entry_id, design_sha256 = entry$design_sha256
  )
  if (file.exists(path)) {
    old <- readRDS(path)
    if (!identical(old$signature, signature)) {
      wcrit21_stop("Checkpoint signature mismatch: %s", path)
    }
    return(old)
  }
  result <- tryCatch(
    list(ok = TRUE, error = NA_character_,
         result = wcrit21_run_entry(entry, label_maps, data_obj, cfg)),
    error = function(e) list(ok = FALSE, error = conditionMessage(e), result = NULL)
  )
  payload <- c(list(signature = signature), result)
  wcrit21_atomic_save_rds(payload, path)
  payload
}

wcrit21_label_invariance_audit <- function(entries, label_maps) {
  rows <- list()
  cursor <- 0L
  for (entry in entries) {
    base <- vapply(c("kendall", "hamming", "l2_position"), function(g) {
      case_geometry_summary(entry$D, g)$minimum
    }, numeric(1L))
    for (i in seq_len(nrow(label_maps))) {
      map_row <- label_maps[i, , drop = FALSE]
      map <- as.integer(unlist(
        map_row[c("label_1", "label_2", "label_3", "label_4")],
        use.names = FALSE
      ))
      mapped <- case_apply_label_map(entry$D, map)
      for (g in names(base)) {
        cursor <- cursor + 1L
        value <- case_geometry_summary(mapped, g)$minimum
        rows[[cursor]] <- data.frame(
          entry_id = entry$entry_id, method = entry$method,
          design_seed_id = entry$design_seed_id, map_id = map_row$map_id,
          geometry = g, original_minimum = base[[g]],
          relabeled_minimum = value, delta = value - base[[g]],
          pass = abs(value - base[[g]]) <= 1e-10,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}

wcrit21_paired_t_summary <- function(df, group_cols, value_col,
                                     confidence = 0.95) {
  alpha <- 1 - confidence
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(.x, .y) {
      x <- as.numeric(.x[[value_col]])
      x <- x[is.finite(x)]
      n <- length(x)
      avg <- if (n) mean(x) else NA_real_
      sd <- if (n >= 2L) stats::sd(x) else NA_real_
      se <- if (n >= 2L) sd / sqrt(n) else NA_real_
      critical <- if (n >= 2L) stats::qt(1 - alpha / 2, df = n - 1L) else NA_real_
      data.frame(
        n = n, mean = avg, sd = sd, se = se,
        ci_low = if (n >= 2L) avg - critical * se else NA_real_,
        ci_high = if (n >= 2L) avg + critical * se else NA_real_,
        wins = sum(x > 1e-12), ties = sum(abs(x) <= 1e-12),
        losses = sum(x < -1e-12), stringsAsFactors = FALSE
      )
    }) |>
    dplyr::ungroup()
}

wcrit21_make_summaries <- function(results, design_manifest, cfg) {
  prediction_summary <- results$prediction |>
    dplyr::group_by(method, stochastic_method, model, bo_step, n_observed, domain) |>
    dplyr::summarise(
      n_rows = dplyr::n(), rmse_mean = mean(rmse), rmse_sd = stats::sd(rmse),
      nrmse_common_mean = mean(nrmse_common_all24_sd),
      nrmse_common_sd = stats::sd(nrmse_common_all24_sd),
      mae_mean = mean(mae), coverage_mean = mean(coverage),
      interval_length_mean = mean(interval_length), .groups = "drop"
    ) |>
    dplyr::mutate(inference_note = ifelse(
      stochastic_method,
      "design-seed Monte Carlo; maps/folds are fixed exhaustive factors",
      "descriptive only; fixed design repeated across seed IDs for alignment"
    ))
  recommendation_summary <- results$recommendation |>
    dplyr::group_by(method, stochastic_method, model, bo_step, n_observed,
                    recommendation_type) |>
    dplyr::summarise(
      n_rows = dplyr::n(), top1_rate = mean(top1_hit),
      selected_in_true_top3_rate = mean(selected_in_true_top3),
      top3_contains_best_rate = mean(top3_contains_true_best),
      top3_overlap_mean = mean(top3_overlap_fraction),
      regret_mean = mean(heldout_regret), regret_sd = stats::sd(heldout_regret),
      reversal_difference_mean = mean(response_minus_reverse),
      simultaneous_difference_mean = mean(response_minus_simultaneous),
      .groups = "drop"
    ) |>
    dplyr::mutate(inference_note = ifelse(
      stochastic_method,
      "design-seed Monte Carlo; maps/folds are fixed exhaustive factors",
      "descriptive only; fixed design repeated across seed IDs for alignment"
    ))
  perf_prediction <- results$prediction |>
    dplyr::filter(domain == "integrated_all24",
                  bo_step %in% c(0L, cfg$bo_additions)) |>
    dplyr::transmute(
      design_seed_id, method, model, bo_step,
      metric = "nrmse_common_all24_sd", value = nrmse_common_all24_sd,
      direction = "minimize"
    )
  perf_rec <- results$recommendation |>
    dplyr::filter(recommendation_type == "surrogate_all24",
                  bo_step %in% c(0L, cfg$bo_additions)) |>
    dplyr::select(design_seed_id, method, model, bo_step,
                  heldout_regret, top1_hit, top3_contains_true_best,
                  selected_in_true_top3, response_minus_simultaneous) |>
    tidyr::pivot_longer(
      cols = c(heldout_regret, top1_hit, selected_in_true_top3,
               top3_contains_true_best,
               response_minus_simultaneous),
      names_to = "metric", values_to = "value"
    ) |>
    dplyr::mutate(
      direction = ifelse(metric == "heldout_regret", "minimize", "maximize")
    )
  perf <- dplyr::bind_rows(perf_prediction, perf_rec) |>
    dplyr::group_by(design_seed_id, method, model, bo_step, metric, direction) |>
    dplyr::summarise(value = mean(as.numeric(value)), .groups = "drop")
  reference <- perf |>
    dplyr::filter(method == "Exact_Phi_lambda050") |>
    dplyr::rename(reference_value = value) |>
    dplyr::select(-method)
  paired <- perf |>
    dplyr::filter(method != "Exact_Phi_lambda050") |>
    dplyr::rename(comparator = method, comparator_value = value) |>
    dplyr::inner_join(
      reference,
      by = c("design_seed_id", "model", "bo_step", "metric", "direction")
    ) |>
    dplyr::mutate(
      reference = "Exact_Phi_lambda050",
      oriented_reference_advantage = ifelse(
        direction == "minimize",
        comparator_value - reference_value,
        reference_value - comparator_value
      )
    )
  stochastic_map <- design_manifest |>
    dplyr::distinct(method, stochastic) |>
    dplyr::rename(comparator = method, comparator_stochastic = stochastic)
  paired <- paired |>
    dplyr::left_join(stochastic_map, by = "comparator") |>
    dplyr::mutate(
      inferential_note = ifelse(
        comparator_stochastic,
        "CI unit is design seed after averaging fixed 24 maps and 3 held-out folds",
        "descriptive only: deterministic comparator is repeated across seed IDs"
      )
    )
  paired_summary <- wcrit21_paired_t_summary(
    paired, c("reference", "comparator", "comparator_stochastic", "model",
              "bo_step", "metric", "direction", "inferential_note"),
    value_col = "oriented_reference_advantage", confidence = cfg$confidence
  )
  hyperparameter_boundary_summary <- results$model_fit |>
    dplyr::group_by(model, model_role, bo_step, n_observed) |>
    dplyr::summarise(
      n_fits = dplyr::n(),
      n_grid_fits = sum(!is.na(theta_boundary_hit)),
      theta_boundary_hits = sum(theta_boundary_hit %in% TRUE),
      nugget_boundary_hits = sum(nugget_boundary_hit %in% TRUE),
      either_boundary_hits = sum(
        (theta_boundary_hit %in% TRUE) | (nugget_boundary_hit %in% TRUE)
      ),
      either_boundary_rate = ifelse(
        n_grid_fits > 0, either_boundary_hits / n_grid_fits, NA_real_
      ),
      sigma_floor_hits = sum(sigma_floor_active %in% TRUE),
      warning_threshold = cfg$boundary_rate_warn_threshold,
      diagnostic_pass = is.na(either_boundary_rate) |
        either_boundary_rate <= warning_threshold,
      action = ifelse(
        diagnostic_pass,
        "none",
        "inspect boundary choices before interpreting uncertainty quantification"
      ),
      .groups = "drop"
    )
  list(
    prediction_summary = prediction_summary,
    recommendation_summary = recommendation_summary,
    paired_raw = paired,
    paired_summary = paired_summary,
    hyperparameter_boundary_summary = hyperparameter_boundary_summary
  )
}

wcrit21_reported_response_tables <- function(data_obj) {
  mean_response <- rowMeans(data_obj$Y_response)
  sd_response <- apply(data_obj$Y_response, 1L, stats::sd)
  se_response <- sd_response / sqrt(3)
  critical <- stats::qt(0.975, df = 2L)
  consensus_order <- order(-mean_response, seq_along(mean_response))
  consensus_rank <- integer(length(mean_response))
  consensus_rank[consensus_order] <- seq_along(consensus_order)
  control_mean <- mean(data_obj$control_response)
  sequence_summary <- data.frame(
    index = seq_len(24L),
    treatment = data_obj$permutation_rows$treatment,
    sequence = data_obj$keys,
    y1 = data_obj$Y_response[, 1L], y2 = data_obj$Y_response[, 2L],
    y3 = data_obj$Y_response[, 3L],
    mean_reported_response = mean_response,
    sd_across_three_replicates = sd_response,
    se_across_three_replicates = se_response,
    ci95_low_t_df2 = mean_response - critical * se_response,
    ci95_high_t_df2 = mean_response + critical * se_response,
    consensus_rank = consensus_rank,
    consensus_top1 = consensus_rank == 1L,
    consensus_top3 = consensus_rank <= 3L,
    consensus_regret = max(mean_response) - mean_response,
    simultaneous_mean_response = control_mean,
    mean_response_minus_simultaneous = mean_response - control_mean,
    interpretation =
      "descriptive full-three-replicate consensus; not held-out CV",
    stringsAsFactors = FALSE
  )
  pair_left <- which(seq_len(24L) < data_obj$reverse_index)
  reverse_rows <- lapply(pair_left, function(i) {
    j <- data_obj$reverse_index[[i]]
    diffs <- data_obj$Y_response[i, ] - data_obj$Y_response[j, ]
    se <- stats::sd(diffs) / sqrt(3)
    data.frame(
      pair_id = sprintf("reverse_pair_%02d", match(i, pair_left)),
      orientation_rule = "lower source-row index minus its reverse",
      left_index = i, left_treatment = data_obj$permutation_rows$treatment[[i]],
      left_sequence = data_obj$keys[[i]],
      right_index = j, right_treatment = data_obj$permutation_rows$treatment[[j]],
      right_sequence = data_obj$keys[[j]],
      left_mean_response = mean_response[[i]],
      right_mean_response = mean_response[[j]],
      difference_y1 = diffs[[1L]], difference_y2 = diffs[[2L]],
      difference_y3 = diffs[[3L]],
      mean_paired_difference = mean(diffs),
      sd_paired_difference = stats::sd(diffs),
      se_paired_difference = se,
      ci95_low_t_df2 = mean(diffs) - critical * se,
      ci95_high_t_df2 = mean(diffs) + critical * se,
      interpretation = "descriptive reversal contrast; n=3 paired replicates",
      stringsAsFactors = FALSE
    )
  })
  control_summary <- data.frame(
    treatment = data_obj$raw$treatment[[25L]],
    y1 = data_obj$control_response[[1L]],
    y2 = data_obj$control_response[[2L]],
    y3 = data_obj$control_response[[3L]],
    mean_reported_response = control_mean,
    sd_across_three_replicates = stats::sd(data_obj$control_response),
    se_across_three_replicates = stats::sd(data_obj$control_response) / sqrt(3),
    interpretation = "simultaneous-administration descriptive control",
    stringsAsFactors = FALSE
  )
  list(
    reported_response_summary = sequence_summary,
    reverse_contrasts_full_three = dplyr::bind_rows(reverse_rows),
    simultaneous_control_summary = control_summary
  )
}

wcrit21_consensus_recommendations <- function(recommendation, data_obj) {
  mean_response <- rowMeans(data_obj$Y_response)
  control_mean <- mean(data_obj$control_response)
  consensus_order <- order(-mean_response, seq_along(mean_response))
  true_best <- consensus_order[[1L]]
  true_top3 <- head(consensus_order, 3L)
  selected <- recommendation$selected_index
  reverse_id <- data_obj$reverse_index[selected]
  raw <- recommendation |>
    dplyr::transmute(
      protocol_hash, design_seed_id, method, stochastic_method,
      map_id, map_key, heldout_fold, model, model_role,
      bo_step, n_observed, recommendation_type,
      selected_index, selected_treatment, selected_sequence,
      consensus_target = "mean of all three reported-response replicates",
      consensus_selected_response = mean_response[selected],
      consensus_best_index = true_best,
      consensus_best_treatment = data_obj$permutation_rows$treatment[[true_best]],
      consensus_best_response = mean_response[[true_best]],
      consensus_top1_hit = selected == true_best,
      consensus_selected_in_true_top3 = selected %in% true_top3,
      consensus_regret = mean_response[[true_best]] - mean_response[selected],
      consensus_reverse_index = reverse_id,
      consensus_response_minus_reverse = mean_response[selected] - mean_response[reverse_id],
      simultaneous_consensus_response = control_mean,
      consensus_response_minus_simultaneous = mean_response[selected] - control_mean,
      inferential_role =
        "post-hoc descriptive consensus only; held-out-fold tables are the CV evidence"
    )
  summary <- raw |>
    dplyr::group_by(method, stochastic_method, model, model_role,
                    bo_step, n_observed, recommendation_type) |>
    dplyr::summarise(
      n_rows = dplyr::n(),
      consensus_top1_rate = mean(consensus_top1_hit),
      consensus_true_top3_rate = mean(consensus_selected_in_true_top3),
      consensus_regret_mean = mean(consensus_regret),
      consensus_regret_sd = stats::sd(consensus_regret),
      consensus_response_minus_control_mean =
        mean(consensus_response_minus_simultaneous),
      consensus_response_minus_reverse_mean =
        mean(consensus_response_minus_reverse),
      .groups = "drop"
    ) |>
    dplyr::mutate(inferential_role = ifelse(
      stochastic_method,
      paste("descriptive consensus; design-seed variation only;",
            "not held-out evaluation"),
      paste("descriptive consensus; fixed comparator repeated for alignment;",
            "not held-out evaluation")
    ))
  list(raw = raw, summary = summary)
}

wcrit21_make_figures <- function(results, response_tables, cfg, dirs) {
  primary_rec <- results$recommendation |>
    dplyr::filter(model == cfg$bo_model,
                  recommendation_type == "incumbent_observed") |>
    dplyr::group_by(design_seed_id, method, stochastic_method,
                    bo_step, n_observed) |>
    dplyr::summarise(
      regret = mean(heldout_regret), top1 = mean(top1_hit),
      true_top3 = mean(selected_in_true_top3), .groups = "drop"
    )
  regret_curve <- primary_rec |>
    dplyr::group_by(method, stochastic_method, bo_step, n_observed) |>
    dplyr::summarise(
      mean = mean(regret), se = case_safe_se(regret),
      n_design_seeds = dplyr::n(), .groups = "drop"
    )
  top_curve <- primary_rec |>
    tidyr::pivot_longer(c(top1, true_top3), names_to = "metric",
                        values_to = "value") |>
    dplyr::group_by(method, stochastic_method, bo_step, n_observed, metric) |>
    dplyr::summarise(
      mean = mean(value), se = case_safe_se(value),
      n_design_seeds = dplyr::n(), .groups = "drop"
    )
  pred_seed <- results$prediction |>
    dplyr::filter(model == cfg$bo_model, domain == "integrated_all24") |>
    dplyr::group_by(design_seed_id, method, stochastic_method,
                    bo_step, n_observed) |>
    dplyr::summarise(
      nrmse = mean(nrmse_common_all24_sd), .groups = "drop"
    )
  prediction_curve <- pred_seed |>
    dplyr::group_by(method, stochastic_method, bo_step, n_observed) |>
    dplyr::summarise(
      mean = mean(nrmse), se = case_safe_se(nrmse),
      n_design_seeds = dplyr::n(), .groups = "drop"
    )
  utils::write.csv(regret_curve,
                   file.path(dirs$results, "mallows_gp_bo_regret_curve_data.csv"),
                   row.names = FALSE)
  utils::write.csv(top_curve,
                   file.path(dirs$results, "mallows_gp_bo_top1_top3_curve_data.csv"),
                   row.names = FALSE)
  utils::write.csv(prediction_curve,
                   file.path(dirs$results, "mallows_gp_bo_prediction_curve_data.csv"),
                   row.names = FALSE)
  base_theme <- ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  p_regret <- ggplot2::ggplot(
    regret_curve,
    ggplot2::aes(n_observed, mean, color = method, group = method)
  ) + ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.3) + base_theme +
    ggplot2::labs(
      x = "Observed sequences", y = "Held-out regret",
      color = "Initial design",
      title = "Intercept-only Mallows-GP BO: held-out incumbent regret",
      subtitle = "Fixed exact/OA curves are descriptive; stochastic CIs use design seed"
    )
  ggplot2::ggsave(file.path(dirs$figures, "mallows_gp_bo_heldout_regret.png"),
                  p_regret, width = 10, height = 6, dpi = 180)
  p_top <- ggplot2::ggplot(
    top_curve,
    ggplot2::aes(n_observed, mean, color = method, group = method)
  ) + ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.3) +
    ggplot2::facet_wrap(~metric, labeller = ggplot2::as_labeller(c(
      top1 = "Selected true best", true_top3 = "Selected in true top 3"
    ))) + base_theme +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Observed sequences", y = "Held-out success rate",
      color = "Initial design",
      title = "Intercept-only Mallows-GP BO: top-order recovery"
    )
  ggplot2::ggsave(file.path(dirs$figures, "mallows_gp_bo_top1_top3.png"),
                  p_top, width = 10, height = 6, dpi = 180)
  p_prediction <- ggplot2::ggplot(
    prediction_curve,
    ggplot2::aes(n_observed, mean, color = method, group = method)
  ) + ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.3) + base_theme +
    ggplot2::labs(
      x = "Observed sequences", y = "Integrated held-out nRMSE",
      color = "Initial design",
      title = "Intercept-only Mallows-GP BO: prediction over all 24 sequences",
      subtitle = "Every fold uses the SD of all 24 held-out responses"
    )
  ggplot2::ggsave(file.path(dirs$figures, "mallows_gp_bo_prediction_nrmse.png"),
                  p_prediction, width = 10, height = 6, dpi = 180)
  resp <- response_tables$reported_response_summary
  p_response <- ggplot2::ggplot(
    resp, ggplot2::aes(stats::reorder(treatment, mean_reported_response),
                       mean_reported_response)
  ) + ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_reported_response - se_across_three_replicates,
                   ymax = mean_reported_response + se_across_three_replicates),
      width = 0.2
    ) + ggplot2::coord_flip() + base_theme +
    ggplot2::labs(
      x = "Order", y = "Mean reported response +/- SE",
      title = "Four-drug reported responses across three replicates"
    )
  ggplot2::ggsave(file.path(dirs$figures, "reported_response_mean_se.png"),
                  p_response, width = 8, height = 7, dpi = 180)
  reverse <- response_tables$reverse_contrasts_full_three
  p_reverse <- ggplot2::ggplot(
    reverse,
    ggplot2::aes(stats::reorder(pair_id, mean_paired_difference),
                 mean_paired_difference)
  ) + ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci95_low_t_df2, ymax = ci95_high_t_df2),
      width = 0.2
    ) + ggplot2::coord_flip() + base_theme +
    ggplot2::labs(
      x = "Reversal pair", y = "Mean paired response difference (95% t interval)",
      title = "Order-versus-reverse contrasts across three paired replicates",
      subtitle = "Orientation is the lower source-row index minus its reverse"
    )
  ggplot2::ggsave(file.path(dirs$figures, "reverse_pair_contrasts.png"),
                  p_reverse, width = 8, height = 6, dpi = 180)
  invisible(list(regret = regret_curve, top = top_curve,
                 prediction = prediction_curve))
}

wcrit21_integrity_audit <- function(cfg, data_obj, entries, design_manifest,
                                    exact, label_maps, label_audit, results,
                                    source_manifest, seed_ledger,
                                    response_tables, consensus, dirs) {
  expected_entries <- length(cfg$design_seed_ids) * length(cfg$methods)
  expected_base_trajectories <- expected_entries * nrow(label_maps) *
    length(cfg$heldout_folds)
  expected_model_stage_rows <- expected_base_trajectories *
    (cfg$bo_additions + 1L)
  expected_prediction <- expected_model_stage_rows * 2L
  expected_recommendation <- expected_model_stage_rows * 2L
  expected_acquisition <- expected_base_trajectories * cfg$bo_additions
  exact_sha <- design_manifest |>
    dplyr::filter(method == "Exact_Phi_lambda050") |>
    dplyr::distinct(method, design_sha256)
  exact_distinct <- dplyr::n_distinct(exact_sha$design_sha256)
  sa_rows <- design_manifest |>
    dplyr::filter(grepl("^Unrestricted_", method))
  checks <- list(
    c("data_file_sha256", data_obj$file_sha256 == cfg$expected_data_sha256,
      data_obj$file_sha256, cfg$expected_data_sha256, "P0"),
    c("complete_S4_plus_control", nrow(data_obj$raw) == 25L,
      nrow(data_obj$raw), 25L, "P0"),
    c("reported_response_table_cardinality",
      nrow(response_tables$reported_response_summary) == 24L,
      nrow(response_tables$reported_response_summary), 24L, "P0"),
    c("reverse_pair_cardinality",
      nrow(response_tables$reverse_contrasts_full_three) == 12L,
      nrow(response_tables$reverse_contrasts_full_three), 12L, "P0"),
    c("published_avg_is_only_rounded", data_obj$avg_rounding_max_abs <= 0.005,
      data_obj$avg_rounding_max_abs, "<=0.005; analysis uses raw y.1:y.3", "P0"),
    c("design_bank_cardinality", length(entries) == expected_entries,
      length(entries), expected_entries, "P0"),
    c("all_initial_designs_unique", all(vapply(entries, function(e) !anyDuplicated(wcrit_row_keys(e$D)), logical(1L))),
      sum(vapply(entries, function(e) anyDuplicated(wcrit_row_keys(e$D)) > 0L, logical(1L))), 0L, "P0"),
    c("all_initial_designs_full_PWO_rank", all(design_manifest$pwo_rank == 7L),
      paste(range(design_manifest$pwo_rank), collapse = ".."), 7L, "P0"),
    c("unrestricted_SA_exact_budget", nrow(sa_rows) > 0L &&
        all(sa_rows$n_obj_eval == cfg$sa_budget) && all(sa_rows$budget_exhausted),
      paste(range(sa_rows$n_obj_eval), collapse = ".."), cfg$sa_budget, "P0"),
    c("rank_retries_within_preregistered_cap", all(
        design_manifest$initial_rank_retry <= cfg$full_rank_retry_max &
          design_manifest$design_rank_retry <= cfg$full_rank_retry_max),
      paste(max(design_manifest$initial_rank_retry),
            max(design_manifest$design_rank_retry), sep = "/"),
      cfg$full_rank_retry_max, "P0"),
    c("exact_enumeration_cardinality", nrow(exact$all_candidates) == choose(12L, 6L),
      nrow(exact$all_candidates), choose(12L, 6L), "P0"),
    c("label_map_cardinality", nrow(label_maps) == length(cfg$label_map_ids),
      nrow(label_maps), length(cfg$label_map_ids), "P0"),
    c("label_invariant_geometries", all(label_audit$pass),
      sum(!label_audit$pass), 0L, "P0"),
    c("prediction_row_cardinality", nrow(results$prediction) == expected_prediction,
      nrow(results$prediction), expected_prediction, "P0"),
    c("recommendation_row_cardinality", nrow(results$recommendation) == expected_recommendation,
      nrow(results$recommendation), expected_recommendation, "P0"),
    c("consensus_row_cardinality", nrow(consensus$raw) == expected_recommendation,
      nrow(consensus$raw), expected_recommendation, "P0"),
    c("acquisition_row_cardinality", nrow(results$acquisition) == expected_acquisition,
      nrow(results$acquisition), expected_acquisition, "P0"),
    c("all_model_fits_successful", nrow(results$model_fit) == expected_model_stage_rows &&
        all(results$model_fit$fit_success),
      sum(!results$model_fit$fit_success), 0L, "P0"),
    c("only_intercept_mallows_gp_is_fitted", all(
        results$model_fit$model == "Intercept_Mallows_GP") &&
        all(results$acquisition$model == "Intercept_Mallows_GP"),
      paste(sort(unique(c(results$model_fit$model,
                          results$acquisition$model))), collapse = ","),
      "Intercept_Mallows_GP", "P0"),
    c("two_discovery_replicates_stacked", all(results$model_fit$n_training_observations ==
        2L * results$model_fit$n_observed),
      paste(range(results$model_fit$n_training_observations -
                    2L * results$model_fit$n_observed), collapse = ".."),
      0L, "P0"),
    c("heldout_prediction_noise_multiplier", all(
        results$prediction$observation_noise_multiplier == 1),
      paste(unique(results$prediction$observation_noise_multiplier), collapse = ","),
      1L, "P0"),
    c("EI_latent_noise_multiplier", all(
        results$acquisition$ei_observation_noise_multiplier == 0),
      paste(unique(results$acquisition$ei_observation_noise_multiplier), collapse = ","),
      0L, "P0"),
    c("EI_never_readds_observed", !any(results$acquisition$duplicate_before_add),
      sum(results$acquisition$duplicate_before_add), 0L, "P0"),
    c("EI_candidate_budget_sequence", all(results$acquisition$candidate_count ==
        24L - results$acquisition$n_observed),
      paste(sort(unique(results$acquisition$candidate_count)), collapse = ","),
      paste((24L - cfg$n0):((24L - cfg$n0) - cfg$bo_additions + 1L), collapse = ","), "P0"),
    c("seed_collision_count", !anyDuplicated(seed_ledger$seed),
      anyDuplicated(seed_ledger$seed), 0L, "P0"),
    c("source_manifest_complete", nrow(source_manifest) == 5L &&
        all(nzchar(source_manifest$sha256)),
      nrow(source_manifest), 5L, "P0"),
    c("primary_mallows_gp_figure_set", all(file.exists(file.path(
        dirs$figures,
        c("mallows_gp_bo_heldout_regret.png", "mallows_gp_bo_top1_top3.png",
          "mallows_gp_bo_prediction_nrmse.png")
      ))),
      sum(file.exists(file.path(
        dirs$figures,
        c("mallows_gp_bo_heldout_regret.png", "mallows_gp_bo_top1_top3.png",
          "mallows_gp_bo_prediction_nrmse.png")
      ))), 3L, "P0"),
    c("reverse_contrast_figure", file.exists(file.path(
        dirs$figures, "reverse_pair_contrasts.png")),
      file.exists(file.path(dirs$figures, "reverse_pair_contrasts.png")),
      TRUE, "P0"),
    c("hyperparameter_boundary_rate_diagnostic", {
      grid_fit <- !is.na(results$model_fit$theta_boundary_hit)
      boundary <- (results$model_fit$theta_boundary_hit %in% TRUE) |
        (results$model_fit$nugget_boundary_hit %in% TRUE)
      mean(boundary[grid_fit]) <= cfg$boundary_rate_warn_threshold
    }, {
      grid_fit <- !is.na(results$model_fit$theta_boundary_hit)
      boundary <- (results$model_fit$theta_boundary_hit %in% TRUE) |
        (results$model_fit$nugget_boundary_hit %in% TRUE)
      sprintf("%.6f", mean(boundary[grid_fit]))
    }, paste0("<=", cfg$boundary_rate_warn_threshold,
              "; WARN only; inspect hyperparameter_boundary_summary.csv"),
      "WARN"),
    c("single_exact_fsa_label", exact_distinct == 1L && nrow(exact_sha) == 1L,
      paste(nrow(exact_sha), exact_distinct, sep = "/"),
      "one exact lambda=0.5 label / one design SHA", "P0")
  )
  dplyr::bind_rows(lapply(checks, function(z) data.frame(
    check = as.character(z[[1L]]), pass = as.logical(z[[2L]]),
    observed = as.character(z[[3L]]), expected = as.character(z[[4L]]),
    severity = as.character(z[[5L]]), stringsAsFactors = FALSE
  )))
}

wcrit21_git_audit <- function(root) {
  run <- function(args) {
    tryCatch(system2("git", c("-C", root, args), stdout = TRUE, stderr = FALSE),
             error = function(e) NA_character_)
  }
  list(
    commit = paste(run(c("rev-parse", "HEAD")), collapse = "\n"),
    status = paste(run(c("status", "--short")), collapse = "\n")
  )
}

wcrit21_main <- function() {
  for (pkg in c("digest", "dplyr", "tidyr", "gtools", "Rcpp", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      wcrit21_stop("Required package is unavailable: %s", pkg)
    }
  }
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
          sample.kind = "Rejection")
  cfg <- wcrit21_config()
  project_root <- wcrit_project_root()
  parent_default <- file.path(
    project_root, "data", "frozen", "section5_3", "experiment21_parent"
  )
  parent_raw <- Sys.getenv("WCRIT30_PARENT_DIR", unset = parent_default)
  if (!grepl("^/", parent_raw)) parent_raw <- file.path(project_root, parent_raw)
  parent_dir <- normalizePath(parent_raw, winslash = "/", mustWork = TRUE)
  frozen <- wcrit30_load_frozen_parent(parent_dir, cfg)
  data_path <- normalizePath(
    file.path(project_root, "data", "case_studies",
              "four_drug_oofaexp_0.1.0.csv"),
    winslash = "/", mustWork = TRUE
  )
  source_manifest <- wcrit21_source_manifest()
  source_bundle_sha256 <- digest::digest(
    source_manifest[c("component", "sha256")], algo = "sha256", serialize = TRUE
  )
  data_obj <- wcrit21_load_data(cfg, data_path)
  out_dir <- file.path(project_root, "outputs", "wcrit", cfg$out_subdir)
  dirs <- wcrit21_dirs(out_dir, cfg$resume)
  profile_path <- file.path(dirs$audits, "stage_profile.csv")
  guard <- list(
    schema = "wcrit30-resume-v1", protocol_hash = cfg$protocol_hash,
    source_bundle_sha256 = source_bundle_sha256,
    data_sha256 = data_obj$file_sha256,
    canonical_data_sha256 = data_obj$canonical_sha256,
    frozen_parent_protocol_hash = frozen$completed$protocol_hash,
    frozen_design_bank_sha256 = cfg$frozen_design_bank_sha256
  )
  guard_path <- file.path(dirs$config, "resume_guard.rds")
  if (cfg$resume) {
    if (!file.exists(guard_path) || !identical(readRDS(guard_path), guard)) {
      wcrit21_stop("Resume guard is missing or differs in protocol/source/data hash")
    }
  } else {
    wcrit21_atomic_save_rds(guard, guard_path)
  }
  utils::write.csv(source_manifest, file.path(dirs$config, "source_manifest.csv"),
                   row.names = FALSE)
  utils::write.csv(frozen$input_manifest,
                   file.path(dirs$config, "frozen_parent_input_manifest.csv"),
                   row.names = FALSE)
  seed_ledger <- frozen$seed_ledger |>
    dplyr::filter(design_seed_id %in% cfg$design_seed_ids)
  utils::write.csv(seed_ledger, file.path(dirs$config, "seed_ledger.csv"),
                   row.names = FALSE)
  label_maps <- frozen$label_maps |>
    dplyr::filter(map_id %in% cfg$label_map_ids)
  utils::write.csv(label_maps, file.path(dirs$config, "label_map_manifest.csv"),
                   row.names = FALSE)
  fold_manifest <- data.frame(
    heldout_fold = cfg$heldout_folds,
    heldout_column = paste0("y.", cfg$heldout_folds),
    discovery_columns = vapply(cfg$heldout_folds, function(f) {
      paste(paste0("y.", setdiff(1:3, f)), collapse = ";")
    }, character(1L)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(fold_manifest, file.path(dirs$config, "fold_manifest.csv"),
                   row.names = FALSE)
  data_manifest <- data.frame(
    path = data_path, file_sha256 = data_obj$file_sha256,
    expected_file_sha256 = cfg$expected_data_sha256,
    canonical_analysis_sha256 = data_obj$canonical_sha256,
    source_package = cfg$source_package,
    source_archive_sha256 = cfg$source_archive_sha256,
    source_rda_sha256 = cfg$source_rda_sha256,
    permutation_rows = 24L, control_rows = 1L,
    avg_rounding_max_abs = data_obj$avg_rounding_max_abs,
    response_used = "raw y.1, y.2, y.3 (not rounded avg)",
    stringsAsFactors = FALSE
  )
  utils::write.csv(data_manifest, file.path(dirs$config, "data_manifest.csv"),
                   row.names = FALSE)
  response_tables <- wcrit21_reported_response_tables(data_obj)
  for (nm in names(response_tables)) {
    utils::write.csv(
      response_tables[[nm]], file.path(dirs$results, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }
  config_flat <- cfg
  config_flat$source_bundle_sha256 <- source_bundle_sha256
  config_flat$data_sha256 <- data_obj$file_sha256
  config_flat$canonical_data_sha256 <- data_obj$canonical_sha256
  wcrit_write_config(file.path(dirs$config, "config.csv"), config_flat)
  jsonlite::write_json(
    config_flat, file.path(dirs$config, "config.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"
  )
  wcrit21_atomic_save_rds(cfg, file.path(dirs$config, "config.rds"))
  audit_manifest <- list(
    config = cfg, source_bundle_sha256 = source_bundle_sha256,
    data_manifest = as.list(data_manifest[1L, ]),
    rng_kind = RNGkind(), R = R.version.string,
    platform = R.version$platform,
    packages = lapply(c("digest", "dplyr", "tidyr", "gtools", "Rcpp", "jsonlite"),
                      function(p) paste0(p, " ", as.character(utils::packageVersion(p)))),
    git = wcrit21_git_audit(project_root),
    inferential_scope = paste(
      "The 24 label maps and three folds are exhaustively averaged fixed",
      "randomization/evaluation factors. Design seed is the paired Monte Carlo",
      "unit only for stochastic methods; this is not a biological-population CI."
    )
  )
  jsonlite::write_json(
    audit_manifest, file.path(dirs$config, "audit_manifest.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"
  )
  writeLines(capture.output(sessionInfo()), file.path(dirs$config, "session_info.txt"))

  entries <- frozen$entries[vapply(
    frozen$entries,
    function(z) z$design_seed_id %in% cfg$design_seed_ids &&
      z$method %in% cfg$methods,
    logical(1L)
  )]
  exact <- frozen$exact
  wcrit21_stage("load_frozen_design_bank", {
    expected_n <- length(cfg$design_seed_ids) * length(cfg$methods)
    if (length(entries) != expected_n) {
      wcrit21_stop("Frozen design subset has %d entries; expected %d",
                   length(entries), expected_n)
    }
    TRUE
  }, profile_path)
  design_manifest <- wcrit21_design_manifest(entries, cfg)
  frozen_manifest_subset <- frozen$design_manifest |>
    dplyr::filter(design_seed_id %in% cfg$design_seed_ids,
                  method %in% cfg$methods) |>
    dplyr::arrange(design_seed_id, match(method, cfg$methods))
  current_design_key <- design_manifest |>
    dplyr::select(entry_id, method, design_seed_id, design_sha256)
  frozen_design_key <- frozen_manifest_subset |>
    dplyr::select(entry_id, method, design_seed_id, design_sha256)
  if (!identical(current_design_key, frozen_design_key)) {
    wcrit21_stop("Reconstructed design manifest does not match frozen Experiment 21")
  }
  utils::write.csv(design_manifest, file.path(dirs$designs, "design_manifest.csv"),
                   row.names = FALSE)
  utils::write.csv(exact$optima, file.path(dirs$designs, "exact_optima.csv"),
                   row.names = FALSE)
  utils::write.csv(exact$all_candidates,
                   file.path(dirs$designs, "exact_enumeration_m4n12.csv"),
                   row.names = FALSE)
  wcrit21_atomic_save_rds(entries, file.path(dirs$designs, "design_bank.rds"))
  label_audit <- wcrit21_stage(
    "label_invariance", wcrit21_label_invariance_audit(entries, label_maps),
    profile_path
  )
  utils::write.csv(label_audit,
                   file.path(dirs$audits, "label_invariance_audit.csv"),
                   row.names = FALSE)

  entry_ids <- seq_along(entries)
  run_one <- function(i) {
    wcrit21_entry_checkpoint(
      entries[[i]], i, label_maps, data_obj, cfg, source_bundle_sha256,
      dirs$checkpoints
    )
  }
  checkpoint_payloads <- wcrit21_stage("prediction_and_bo", {
    if (cfg$workers > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(entry_ids, run_one, mc.cores = cfg$workers,
                         mc.preschedule = FALSE, mc.set.seed = FALSE)
    } else lapply(entry_ids, run_one)
  }, profile_path)
  failures <- dplyr::bind_rows(lapply(seq_along(checkpoint_payloads), function(i) {
    x <- checkpoint_payloads[[i]]
    if (isTRUE(x$ok)) return(NULL)
    data.frame(entry_number = i, entry_id = entries[[i]]$entry_id,
               error = x$error, stringsAsFactors = FALSE)
  }))
  if (nrow(failures)) {
    utils::write.csv(failures, file.path(dirs$audits, "task_failures.csv"),
                     row.names = FALSE)
    wcrit21_stop("%d entry tasks failed; see task_failures.csv", nrow(failures))
  }
  entry_results <- lapply(checkpoint_payloads, function(x) x$result)
  results <- list(
    prediction = dplyr::bind_rows(lapply(entry_results, `[[`, "prediction")),
    recommendation = dplyr::bind_rows(lapply(entry_results, `[[`, "recommendation")),
    model_fit = dplyr::bind_rows(lapply(entry_results, `[[`, "model_fit")),
    acquisition = dplyr::bind_rows(lapply(entry_results, `[[`, "acquisition"))
  )
  for (nm in names(results)) {
    results[[nm]]$data_sha256 <- data_obj$file_sha256
    results[[nm]]$source_bundle_sha256 <- source_bundle_sha256
    utils::write.csv(results[[nm]], file.path(dirs$raw, paste0(nm, ".csv")),
                     row.names = FALSE)
  }
  wcrit21_atomic_save_rds(
    list(config = cfg, data_manifest = data_manifest,
         design_manifest = design_manifest, results = results),
    file.path(dirs$raw, "four_drug_case_results.rds")
  )
  summaries <- wcrit21_stage(
    "summaries", wcrit21_make_summaries(results, design_manifest, cfg),
    profile_path
  )
  for (nm in names(summaries)) {
    utils::write.csv(summaries[[nm]], file.path(dirs$results, paste0(nm, ".csv")),
                     row.names = FALSE)
  }
  utils::write.csv(
    summaries$prediction_summary |>
      dplyr::filter(bo_step == 0L),
    file.path(dirs$results, "initial_prediction_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$recommendation_summary |>
      dplyr::filter(bo_step == 0L),
    file.path(dirs$results, "initial_top3_regret_summary.csv"),
    row.names = FALSE
  )
  reused_initial_prediction <- frozen$initial_prediction |>
    dplyr::filter(
      model == "Intercept_Mallows_GP", bo_step == 0L,
      domain == "currently_unobserved", method %in% cfg$methods
    ) |>
    dplyr::arrange(match(method, cfg$methods))
  reused_initial_recommendation <- frozen$initial_recommendation |>
    dplyr::filter(
      model == "Intercept_Mallows_GP", bo_step == 0L,
      recommendation_type == "incumbent_observed", method %in% cfg$methods
    ) |>
    dplyr::arrange(match(method, cfg$methods))
  utils::write.csv(
    reused_initial_prediction,
    file.path(
      dirs$results,
      "initial_prediction_reused_from_experiment21_intercept_mallows_gp.csv"
    ),
    row.names = FALSE
  )
  utils::write.csv(
    reused_initial_recommendation,
    file.path(
      dirs$results,
      "model_independent_initial_incumbent_reused_from_experiment21.csv"
    ),
    row.names = FALSE
  )
  reuse_checks <- data.frame(
    check = c(
      "frozen_parent_identity",
      "frozen_design_manifest_identity",
      "formal_step0_intercept_prediction_reproduces_experiment21",
      "formal_step0_model_independent_incumbent_reproduces_experiment21"
    ),
    pass = c(TRUE, TRUE, TRUE, TRUE),
    observed = c(
      frozen$completed$protocol_hash,
      cfg$frozen_design_manifest_sha256,
      if (cfg$profile == "formal") "pending comparison" else "smoke: not compared",
      if (cfg$profile == "formal") "pending comparison" else "smoke: not compared"
    ),
    expected = c(
      cfg$frozen_parent_protocol_hash,
      cfg$frozen_design_manifest_sha256,
      "exact equality for formal 20x24x3 summary",
      "exact equality for formal 20x24x3 summary"
    ),
    severity = c("P0", "P0", "P0", "P0"),
    stringsAsFactors = FALSE
  )
  if (cfg$profile == "formal") {
    current_initial_prediction <- summaries$prediction_summary |>
      dplyr::filter(
        model == "Intercept_Mallows_GP", bo_step == 0L,
        domain == "currently_unobserved", method %in% cfg$methods
      ) |>
      dplyr::arrange(match(method, cfg$methods))
    current_initial_recommendation <- summaries$recommendation_summary |>
      dplyr::filter(
        model == "Intercept_Mallows_GP", bo_step == 0L,
        recommendation_type == "incumbent_observed", method %in% cfg$methods
      ) |>
      dplyr::arrange(match(method, cfg$methods))
    pred_cols <- intersect(names(current_initial_prediction),
                           names(reused_initial_prediction))
    rec_cols <- intersect(names(current_initial_recommendation),
                          names(reused_initial_recommendation))
    pred_equal <- isTRUE(all.equal(
      current_initial_prediction[pred_cols], reused_initial_prediction[pred_cols],
      tolerance = 0, check.attributes = FALSE
    ))
    rec_equal <- isTRUE(all.equal(
      current_initial_recommendation[rec_cols],
      reused_initial_recommendation[rec_cols],
      tolerance = 0, check.attributes = FALSE
    ))
    reuse_checks$pass[3:4] <- c(pred_equal, rec_equal)
    reuse_checks$observed[3:4] <- c(pred_equal, rec_equal)
  }
  utils::write.csv(reuse_checks,
                   file.path(dirs$audits, "frozen_reuse_audit.csv"),
                   row.names = FALSE)
  consensus <- wcrit21_consensus_recommendations(
    results$recommendation, data_obj
  )
  utils::write.csv(
    consensus$raw,
    file.path(dirs$raw, "full_three_consensus_recommendation_raw.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    consensus$summary,
    file.path(dirs$results, "full_three_consensus_top_regret_control.csv"),
    row.names = FALSE
  )
  wcrit21_atomic_save_rds(
    list(
      config = cfg, data_manifest = data_manifest,
      response_tables = response_tables, design_manifest = design_manifest,
      results = results, summaries = summaries, consensus = consensus
    ),
    file.path(dirs$results, "four_drug_case_complete.rds")
  )
  wcrit21_stage(
    "figures", wcrit21_make_figures(results, response_tables, cfg, dirs),
    profile_path
  )
  unique_designs <- design_manifest |>
    dplyr::group_by(method, stochastic) |>
    dplyr::summarise(
      seed_rows = dplyr::n(), unique_design_sha256 = dplyr::n_distinct(design_sha256),
      .groups = "drop"
    ) |>
    dplyr::mutate(inference_note = ifelse(
      stochastic,
      "paired Monte Carlo unit is design seed",
      "fixed comparator: repeated seed rows are alignment only, not independent evidence"
    ))
  utils::write.csv(unique_designs,
                   file.path(dirs$results, "unique_designs_by_method.csv"),
                   row.names = FALSE)
  audit <- wcrit21_integrity_audit(
    cfg, data_obj, entries, design_manifest, exact, label_maps,
    label_audit, results, source_manifest, seed_ledger,
    response_tables, consensus, dirs
  )
  audit <- dplyr::bind_rows(audit, reuse_checks)
  utils::write.csv(audit, file.path(dirs$audits, "integrity_audit.csv"),
                   row.names = FALSE)
  if (any(audit$severity == "P0" & !audit$pass)) {
    wcrit21_stop("A P0 integrity gate failed; see integrity_audit.csv")
  }
  completion <- c(
    paste0("experiment=", cfg$experiment),
    paste0("profile=", cfg$profile),
    paste0("protocol_hash=", cfg$protocol_hash),
    paste0("source_bundle_sha256=", source_bundle_sha256),
    paste0("data_sha256=", data_obj$file_sha256),
    paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  )
  writeLines(completion, file.path(dirs$root, "COMPLETED"))
  case_write_file_manifest(dirs$root,
                           file.path(dirs$audits, "file_manifest.csv"))
  message("Completed four-drug case study: ", dirs$root)
  invisible(list(output_dir = dirs$root, config = cfg, audit = audit))
}

if (identical(environment(), globalenv()) &&
    length(.wcrit21_file_arg) &&
    normalizePath(sub("^--file=", "", .wcrit21_file_arg[[1L]]),
                  winslash = "/", mustWork = FALSE) == .wcrit21_this_file) {
  wcrit21_main()
}
