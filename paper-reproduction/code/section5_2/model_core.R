# Experiment 29 model-evaluation helpers.
#
# This module is intentionally independent of Experiment 16.  It sources only
# wcrit_common.R and case_study_common.R, and implements the PWO-only and
# Mallows-GP-only response generation/fitting protocol used by Experiment 29.

wcrit29_locate_helper_dir <- function() {
  candidates <- unique(c(
    Sys.getenv("WCRIT29_HELPER_DIR", unset = ""),
    getwd(),
    file.path(getwd(), "code", "WeightedCrit"),
    file.path(dirname(getwd()), "code", "common")
  ))
  candidates <- candidates[nzchar(candidates)]
  ok <- vapply(candidates, function(path) {
    file.exists(file.path(path, "wcrit_common.R")) &&
      file.exists(file.path(path, "case_study_common.R"))
  }, logical(1L))
  if (!any(ok)) {
    stop(
      "Cannot locate wcrit_common.R and case_study_common.R; set WCRIT29_HELPER_DIR",
      call. = FALSE
    )
  }
  normalizePath(candidates[which(ok)[1L]], winslash = "/", mustWork = TRUE)
}

.wcrit29_helper_dir <- wcrit29_locate_helper_dir()
if (!exists("wcrit_hash_seed", mode = "function", inherits = TRUE)) {
  source(file.path(.wcrit29_helper_dir, "wcrit_common.R"), chdir = TRUE)
}
if (!exists("case_solve_chol", mode = "function", inherits = TRUE)) {
  source(file.path(.wcrit29_helper_dir, "case_study_common.R"), chdir = TRUE)
}

# When this temporary module is run directly, commandArgs() points to
# /private/tmp rather than to the helper directory.  Keep the common helper's
# C++ loader pointed at the actual source directory.
if (!file.exists(file.path(wcrit_script_dir(), "sa_core.cpp"))) {
  wcrit29_helper_dir_fixed <- .wcrit29_helper_dir
  wcrit_script_dir <- function() wcrit29_helper_dir_fixed
}

wcrit29_default_c_grid <- function() {
  c(0.125, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32)
}

wcrit29_num_tag <- function(x) {
  sub(
    "\\.$", "",
    gsub("\\.", "p", trimws(formatC(as.numeric(x), format = "fg", digits = 8)))
  )
}

wcrit29_model_build_universe <- function(m = 6L) {
  m <- as.integer(m)
  if (!is.finite(m) || m < 2L) stop("m must be an integer >= 2", call. = FALSE)
  P <- gtools::permutations(m, m, seq_len(m))
  storage.mode(P) <- "integer"
  keys <- wcrit_row_keys(P)
  F <- wcrit_pwo_matrix(P)
  D <- wcrit_kendall_dmat(P)
  list(
    m = m,
    q = as.integer(choose(m, 2)),
    size = nrow(P),
    P = P,
    keys = keys,
    F = F,
    D = D
  )
}

wcrit29_lookup_design_indices <- function(D, universe, context = "design") {
  D <- case_validate_design(D, universe$m, nrow(D), context)
  idx <- match(wcrit_row_keys(D), universe$keys)
  if (anyNA(idx)) stop(context, " contains a row outside the universe", call. = FALSE)
  as.integer(idx)
}

wcrit29_common_test_indices <- function(designs, universe) {
  if (!is.list(designs) || length(designs) < 1L) {
    stop("designs must be a nonempty named list", call. = FALSE)
  }
  used <- sort(unique(unlist(lapply(seq_along(designs), function(i) {
    label <- names(designs)[i]
    if (is.null(label) || is.na(label) || !nzchar(label)) label <- paste0("design_", i)
    wcrit29_lookup_design_indices(designs[[i]], universe, label)
  }), use.names = FALSE)))
  test <- setdiff(seq_len(universe$size), used)
  if (!length(test)) stop("The common test set is empty", call. = FALSE)
  as.integer(test)
}

wcrit29_chol_with_jitter <- function(K, jitter_grid = c(0, 1e-12, 1e-10)) {
  K <- as.matrix(K)
  if (nrow(K) != ncol(K)) stop("K must be square", call. = FALSE)
  jitter_grid <- sort(unique(as.numeric(jitter_grid)))
  jitter_grid <- jitter_grid[is.finite(jitter_grid) & jitter_grid >= 0]
  for (jitter in jitter_grid) {
    R <- K
    diag(R) <- diag(R) + jitter
    U <- try(chol(R), silent = TRUE)
    if (!inherits(U, "try-error")) {
      return(list(success = TRUE, chol = U, jitter = jitter))
    }
  }
  list(success = FALSE, error = "cholesky_failed")
}

# ---- Dense PWO response -----------------------------------------------------

wcrit29_generate_pwo_truth <- function(universe, rep_id, master_seed,
                                       snr_values = c(2, 5)) {
  snr_values <- sort(unique(as.numeric(snr_values)))
  if (!length(snr_values) || any(!is.finite(snr_values)) || any(snr_values <= 0)) {
    stop("snr_values must be finite and positive", call. = FALSE)
  }
  beta_seed <- wcrit_hash_seed(master_seed, "29", "pwo_beta", as.integer(rep_id))
  noise_seed <- wcrit_hash_seed(master_seed, "29", "pwo_noise", as.integer(rep_id))
  set.seed(beta_seed)
  beta_effects <- stats::rnorm(universe$q, sd = 1 / sqrt(universe$q))
  beta_full <- c(0, beta_effects)
  f <- as.numeric(universe$F %*% beta_full)
  signal_sd <- stats::sd(f)
  if (!is.finite(signal_sd) || signal_sd <= 1e-12) {
    stop("Degenerate dense PWO signal", call. = FALSE)
  }
  set.seed(noise_seed)
  shared_standard_noise <- stats::rnorm(universe$size)
  sigma <- signal_sd / snr_values
  names(sigma) <- vapply(snr_values, wcrit29_num_tag, character(1L))
  y <- lapply(sigma, function(s) f + s * shared_standard_noise)
  names(y) <- paste0("pwo_snr", names(sigma))
  list(
    family = "pwo",
    rep = as.integer(rep_id),
    beta_seed = beta_seed,
    noise_seed = noise_seed,
    beta_effects = beta_effects,
    beta_full = beta_full,
    f = f,
    signal_sd = signal_sd,
    shared_standard_noise = shared_standard_noise,
    snr_values = snr_values,
    sigma = sigma,
    y = y
  )
}

wcrit29_fit_pwo_ols <- function(train_D, y, rank_tol = 1e-9,
                                sigma_floor = .Machine$double.eps) {
  train_D <- as.matrix(train_D)
  storage.mode(train_D) <- "integer"
  y <- as.numeric(y)
  X <- wcrit_pwo_matrix(train_D)
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n || n <= p) {
    return(list(success = FALSE, error = "dimension_or_df"))
  }
  qx <- qr(X, tol = rank_tol, LAPACK = FALSE)
  if (qx$rank < p) {
    return(list(success = FALSE, error = "rank_deficient", rank = qx$rank, p = p))
  }
  beta <- qr.coef(qx, y)
  if (any(!is.finite(beta))) return(list(success = FALSE, error = "nonfinite_beta"))
  residual <- as.numeric(y - X %*% beta)
  df <- n - p
  sse <- sum(residual^2)
  sigma2 <- sse / df
  XtX <- crossprod(X)
  U <- try(chol(XtX), silent = TRUE)
  if (inherits(U, "try-error")) return(list(success = FALSE, error = "xtx_cholesky"))
  list(
    success = TRUE,
    error = NA_character_,
    beta = as.numeric(beta),
    sigma2 = as.numeric(max(sigma2, sigma_floor)),
    sigma2_raw = as.numeric(sigma2),
    sigma_floor_active = sigma2 < sigma_floor,
    residual = residual,
    sse = sse,
    df = df,
    rank = qx$rank,
    p = p,
    n = n,
    chol_XtX = U,
    condition_XtX = kappa(XtX, exact = TRUE),
    train_D = train_D
  )
}

wcrit29_predict_pwo_ols <- function(model, test_D) {
  if (!isTRUE(model$success)) stop("Cannot predict from failed PWO fit", call. = FALSE)
  Xtest <- wcrit_pwo_matrix(test_D)
  mu <- as.numeric(Xtest %*% model$beta)
  inv_times_xt <- case_solve_chol(model$chol_XtX, t(Xtest))
  factor <- colSums(t(Xtest) * inv_times_xt)
  if (any(factor < -1e-10)) stop("Negative PWO variance factor", call. = FALSE)
  factor <- pmax(0, factor)
  list(
    mean = mu,
    variance_factor = factor,
    variance = model$sigma2 * factor
  )
}

wcrit29_prediction_metrics <- function(y_true, prediction, normalization_sd,
                                       confidence = 0.95,
                                       critical = c("normal", "t"), df = Inf,
                                       variance_floor_multiplier = 1e-12) {
  critical <- match.arg(critical)
  y_true <- as.numeric(y_true)
  mu <- as.numeric(prediction$mean)
  variance_raw <- as.numeric(prediction$variance)
  if (length(y_true) != length(mu) || length(mu) != length(variance_raw)) {
    stop("Prediction metric dimensions do not agree", call. = FALSE)
  }
  normalization_sd <- as.numeric(normalization_sd)
  if (length(normalization_sd) != 1L || !is.finite(normalization_sd) || normalization_sd <= 0) {
    stop("normalization_sd must be one positive finite number", call. = FALSE)
  }
  floor_value <- max(
    .Machine$double.eps,
    variance_floor_multiplier * max(stats::var(y_true), normalization_sd^2, 1)
  )
  variance <- pmax(variance_raw, floor_value)
  err <- y_true - mu
  alpha <- 1 - confidence
  crit <- if (critical == "t") {
    stats::qt(1 - alpha / 2, df = df)
  } else {
    stats::qnorm(1 - alpha / 2)
  }
  half <- crit * sqrt(variance)
  list(
    rmse = sqrt(mean(err^2)),
    nrmse = sqrt(mean(err^2)) / normalization_sd,
    mae = mean(abs(err)),
    nll = mean(0.5 * log(2 * pi * variance) + err^2 / (2 * variance)),
    coverage = mean(y_true >= mu - half & y_true <= mu + half),
    interval_length = mean(2 * half),
    ipv = mean(variance_raw),
    ipv_factor = if (!is.null(prediction$variance_factor)) {
      mean(as.numeric(prediction$variance_factor))
    } else {
      NA_real_
    },
    variance_floor = floor_value,
    variance_clipped_n = sum(variance_raw < floor_value)
  )
}

wcrit29_evaluate_pwo_design <- function(D, truth, snr, universe, test_idx,
                                        confidence = 0.95, rank_tol = 1e-9,
                                        return_fit = FALSE) {
  snr <- as.numeric(snr)
  tag <- wcrit29_num_tag(snr)
  scenario_id <- paste0("pwo_snr", tag)
  if (!(scenario_id %in% names(truth$y))) {
    stop("The requested SNR is absent from the PWO truth bank", call. = FALSE)
  }
  train_idx <- wcrit29_lookup_design_indices(D, universe, "PWO training design")
  test_idx <- as.integer(test_idx)
  fit <- wcrit29_fit_pwo_ols(D, truth$y[[scenario_id]][train_idx], rank_tol = rank_tol)
  base <- data.frame(
    model_family = "pwo",
    scenario_id = scenario_id,
    fit_mode = "ols",
    snr = snr,
    fit_success = isTRUE(fit$success),
    fit_error = if (isTRUE(fit$success)) NA_character_ else fit$error,
    train_n = nrow(D),
    common_test_n = length(test_idx),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(fit$success)) {
    metrics <- cbind(base, data.frame(
      pwo_rank = if (!is.null(fit$rank)) fit$rank else NA_integer_,
      pwo_df = NA_integer_, sigma2_hat = NA_real_, sigma2_true = truth$sigma[[tag]]^2,
      condition_XtX = NA_real_, beta_rmse_all = NA_real_, beta_rmse_effects = NA_real_,
      nrmse_full = NA_real_, nrmse_common = NA_real_, rmse_full = NA_real_, rmse_common = NA_real_,
      nll_full = NA_real_, nll_common = NA_real_, coverage_full = NA_real_, coverage_common = NA_real_,
      interval_length_full = NA_real_, interval_length_common = NA_real_,
      ipv_factor_full = NA_real_, ipv_factor_common = NA_real_,
      ipv_hat_full = NA_real_, ipv_hat_common = NA_real_,
      ipv_true_full = NA_real_, ipv_true_common = NA_real_
    ))
    return(if (return_fit) list(metrics = metrics, fit = fit) else metrics)
  }
  pred_full <- wcrit29_predict_pwo_ols(fit, universe$P)
  pred_test <- wcrit29_predict_pwo_ols(fit, universe$P[test_idx, , drop = FALSE])
  score_full <- wcrit29_prediction_metrics(
    truth$f, pred_full, truth$signal_sd, confidence, "t", fit$df
  )
  score_test <- wcrit29_prediction_metrics(
    truth$f[test_idx], pred_test, truth$signal_sd, confidence, "t", fit$df
  )
  sigma2_true <- as.numeric(truth$sigma[[tag]]^2)
  metrics <- cbind(base, data.frame(
    pwo_rank = fit$rank,
    pwo_df = fit$df,
    sigma2_hat = fit$sigma2,
    sigma2_true = sigma2_true,
    condition_XtX = fit$condition_XtX,
    beta_rmse_all = sqrt(mean((fit$beta - truth$beta_full)^2)),
    beta_rmse_effects = sqrt(mean((fit$beta[-1L] - truth$beta_effects)^2)),
    nrmse_full = score_full$nrmse,
    nrmse_common = score_test$nrmse,
    rmse_full = score_full$rmse,
    rmse_common = score_test$rmse,
    nll_full = score_full$nll,
    nll_common = score_test$nll,
    coverage_full = score_full$coverage,
    coverage_common = score_test$coverage,
    interval_length_full = score_full$interval_length,
    interval_length_common = score_test$interval_length,
    ipv_factor_full = score_full$ipv_factor,
    ipv_factor_common = score_test$ipv_factor,
    ipv_hat_full = score_full$ipv,
    ipv_hat_common = score_test$ipv,
    ipv_true_full = sigma2_true * score_full$ipv_factor,
    ipv_true_common = sigma2_true * score_test$ipv_factor,
    stringsAsFactors = FALSE
  ))
  if (return_fit) list(metrics = metrics, fit = fit) else metrics
}

# ---- Raw Mallows-GP response -----------------------------------------------

wcrit29_mallows_kernel <- function(distance_matrix, c_value, q) {
  c_value <- as.numeric(c_value)
  q <- as.numeric(q)
  if (length(c_value) != 1L || !is.finite(c_value) || c_value <= 0 ||
      length(q) != 1L || !is.finite(q) || q <= 0) {
    stop("c_value and q must be positive finite scalars", call. = FALSE)
  }
  exp(-(c_value / q) * as.matrix(distance_matrix))
}

wcrit29_prepare_gp_kernel <- function(universe, c_value,
                                      truth_jitter_grid = c(0, 1e-12, 1e-10)) {
  K <- wcrit29_mallows_kernel(universe$D, c_value, universe$q)
  cf <- wcrit29_chol_with_jitter(K, truth_jitter_grid)
  if (!isTRUE(cf$success)) stop("Full-space Mallows kernel Cholesky failed", call. = FALSE)
  list(c = as.numeric(c_value), q = universe$q, K = K, chol = cf$chol, jitter = cf$jitter)
}

wcrit29_generate_gp_truth <- function(universe, c_value, rep_id, master_seed,
                                      kernel_cache = NULL) {
  if (is.null(kernel_cache)) kernel_cache <- wcrit29_prepare_gp_kernel(universe, c_value)
  if (abs(kernel_cache$c - c_value) > 1e-12 || kernel_cache$q != universe$q) {
    stop("kernel_cache does not match c_value/universe", call. = FALSE)
  }
  seed <- wcrit_hash_seed(
    master_seed, "29", "gp_path", wcrit29_num_tag(c_value), as.integer(rep_id)
  )
  set.seed(seed)
  xi <- stats::rnorm(universe$size)
  g <- as.numeric(t(kernel_cache$chol) %*% xi)
  if (length(g) != universe$size || any(!is.finite(g))) {
    stop("Invalid Mallows-GP realization", call. = FALSE)
  }
  list(
    family = "mallows_gp",
    scenario_id = paste0("gp_c", wcrit29_num_tag(c_value)),
    rep = as.integer(rep_id),
    c_true = as.numeric(c_value),
    theta_true = as.numeric(c_value / universe$q),
    path_seed = seed,
    truth_kernel_jitter = kernel_cache$jitter,
    g = g,
    path_mean = mean(g),
    path_sd = stats::sd(g),
    empirically_centered = FALSE,
    empirically_standardized = FALSE
  )
}

wcrit29_fit_ok_at_c <- function(train_D, y, c_value, q,
                                jitter = 1e-10, rank_tol = 1e-9,
                                sigma_floor = .Machine$double.eps,
                                compute_condition = TRUE) {
  train_D <- as.matrix(train_D)
  storage.mode(train_D) <- "integer"
  y <- as.numeric(y)
  n <- nrow(train_D)
  F <- matrix(1, nrow = n, ncol = 1L)
  p <- 1L
  df <- n - p
  if (length(y) != n || df <= 0L) return(list(success = FALSE, error = "dimension_or_df"))
  dtrain <- wcrit_kendall_dmat(train_D)
  R <- wcrit29_mallows_kernel(dtrain, c_value, q)
  diag(R) <- diag(R) + jitter
  U <- try(chol(R), silent = TRUE)
  if (inherits(U, "try-error")) return(list(success = FALSE, error = "correlation_cholesky"))
  yw <- try(forwardsolve(t(U), y), silent = TRUE)
  Fw <- try(forwardsolve(t(U), F), silent = TRUE)
  if (inherits(yw, "try-error") || inherits(Fw, "try-error")) {
    return(list(success = FALSE, error = "whitening_solve"))
  }
  qF <- qr(Fw, tol = rank_tol, LAPACK = FALSE)
  if (qF$rank < p) return(list(success = FALSE, error = "trend_rank"))
  mu <- qr.coef(qF, yw)
  residual_w <- as.numeric(yw - Fw %*% mu)
  sse <- sum(residual_w^2)
  sigma2_raw <- sse / df
  sigma2 <- max(sigma2_raw, sigma_floor)
  diag_rf <- abs(diag(qr.R(qF))[seq_len(p)])
  if (any(!is.finite(diag_rf)) || any(diag_rf <= 0)) {
    return(list(success = FALSE, error = "trend_logdet"))
  }
  logdet_R <- 2 * sum(log(diag(U)))
  logdet_C <- 2 * sum(log(diag_rf))
  reml <- -0.5 * (
    df * (log(2 * pi) + 1 + log(sigma2)) + logdet_R + logdet_C
  )
  residual <- as.numeric(y - F %*% mu)
  alpha <- try(case_solve_chol(U, residual), silent = TRUE)
  chol_C <- try(chol(crossprod(Fw)), silent = TRUE)
  if (inherits(alpha, "try-error") || inherits(chol_C, "try-error")) {
    return(list(success = FALSE, error = "final_solve"))
  }
  list(
    success = TRUE,
    error = NA_character_,
    c = as.numeric(c_value),
    theta = as.numeric(c_value / q),
    q = as.integer(q),
    mu = as.numeric(mu),
    beta = as.numeric(mu),
    sigma2 = as.numeric(sigma2),
    sigma2_raw = as.numeric(sigma2_raw),
    sigma_floor_active = sigma2_raw < sigma_floor,
    reml = as.numeric(reml),
    sse = sse,
    df = df,
    n = n,
    p = p,
    jitter = jitter,
    chol_R = U,
    chol_C = chol_C,
    F_train = F,
    alpha = as.numeric(alpha),
    train_D = train_D,
    condition_R = if (isTRUE(compute_condition)) kappa(R, exact = TRUE) else NA_real_
  )
}

wcrit29_fit_ok_estimated_c <- function(train_D, y, q,
                                       c_grid = wcrit29_default_c_grid(),
                                       jitter = 1e-10, rank_tol = 1e-9,
                                       sigma_floor = .Machine$double.eps,
                                       optimize_tol = 1e-6) {
  c_grid <- sort(unique(as.numeric(c_grid)))
  c_grid <- c_grid[is.finite(c_grid) & c_grid > 0]
  if (length(c_grid) < 2L) stop("c_grid must contain at least two positive values", call. = FALSE)
  fit_one <- function(c_value) {
    wcrit29_fit_ok_at_c(
      train_D, y, c_value, q, jitter = jitter,
      rank_tol = rank_tol, sigma_floor = sigma_floor,
      compute_condition = FALSE
    )
  }
  grid_fits <- lapply(c_grid, fit_one)
  grid_reml <- vapply(grid_fits, function(z) {
    if (isTRUE(z$success) && is.finite(z$reml)) z$reml else -Inf
  }, numeric(1L))
  if (!any(is.finite(grid_reml))) return(list(success = FALSE, error = "all_grid_candidates_failed"))
  best_i <- which.max(grid_reml)
  lower_i <- max(1L, best_i - 1L)
  upper_i <- min(length(c_grid), best_i + 1L)
  log_lower <- log(c_grid[lower_i])
  log_upper <- log(c_grid[upper_i])
  continuous <- NULL
  if (log_upper > log_lower + .Machine$double.eps) {
    objective <- function(log_c) {
      z <- fit_one(exp(log_c))
      if (isTRUE(z$success) && is.finite(z$reml)) z$reml else -1e300
    }
    opt <- try(stats::optimize(
      objective, interval = c(log_lower, log_upper), maximum = TRUE,
      tol = optimize_tol
    ), silent = TRUE)
    if (!inherits(opt, "try-error") && is.finite(opt$objective)) {
      continuous <- fit_one(exp(opt$maximum))
    }
  }
  candidates <- c(grid_fits, if (!is.null(continuous)) list(continuous) else list())
  candidate_reml <- vapply(candidates, function(z) {
    if (isTRUE(z$success) && is.finite(z$reml)) z$reml else -Inf
  }, numeric(1L))
  best_profile <- candidates[[which.max(candidate_reml)]]
  best <- wcrit29_fit_ok_at_c(
    train_D, y, best_profile$c, q, jitter = jitter,
    rank_tol = rank_tol, sigma_floor = sigma_floor,
    compute_condition = TRUE
  )
  if (!isTRUE(best$success)) {
    return(list(success = FALSE, error = "final_continuous_candidate_failed"))
  }
  c_min <- min(c_grid)
  c_max <- max(c_grid)
  boundary_tol <- max(1e-8, optimize_tol * 10)
  best$fit_mode <- "estimated_c"
  best$c_grid_min <- c_min
  best$c_grid_max <- c_max
  best$c_boundary_hit <-
    abs(log(best$c) - log(c_min)) <= boundary_tol ||
    abs(log(best$c) - log(c_max)) <= boundary_tol
  best$c_grid_best <- c_grid[best_i]
  best$c_continuous_refined <- !is.null(continuous) &&
    abs(best$c - continuous$c) <= max(1e-12, optimize_tol * continuous$c)
  best$profile_grid <- data.frame(
    c = c_grid,
    reml = grid_reml,
    success = vapply(grid_fits, function(z) isTRUE(z$success), logical(1L)),
    stringsAsFactors = FALSE
  )
  best
}

wcrit29_predict_ok <- function(model, test_D) {
  if (!isTRUE(model$success)) stop("Cannot predict from failed GP fit", call. = FALSE)
  test_D <- as.matrix(test_D)
  storage.mode(test_D) <- "integer"
  dtest <- wcrit_kendall_dmat(test_D, model$train_D)
  Rstar <- wcrit29_mallows_kernel(dtest, model$c, model$q)
  Ftest <- matrix(1, nrow = nrow(test_D), ncol = 1L)
  mean_pred <- as.numeric(Ftest %*% model$mu + Rstar %*% model$alpha)
  B <- case_solve_chol(model$chol_R, t(Rstar))
  base <- 1 - colSums(t(Rstar) * B)
  u <- t(Ftest) - crossprod(model$F_train, B)
  Cinv_u <- case_solve_chol(model$chol_C, u)
  trend <- colSums(u * Cinv_u)
  raw_factor <- base + trend
  if (any(raw_factor < -1e-7)) {
    stop("Materially negative ordinary-kriging variance factor", call. = FALSE)
  }
  factor <- pmax(0, raw_factor)
  list(
    mean = mean_pred,
    variance_factor = factor,
    variance = model$sigma2 * factor
  )
}

wcrit29_evaluate_gp_design <- function(D, truth, universe, test_idx,
                                       fit_mode = "estimated_c",
                                       c_grid = wcrit29_default_c_grid(),
                                       confidence = 0.95, jitter = 1e-10,
                                       rank_tol = 1e-9,
                                       return_fit = FALSE) {
  if (!identical(as.character(fit_mode), "estimated_c")) {
    stop("The paper GP evaluation uses estimated_c only", call. = FALSE)
  }
  train_idx <- wcrit29_lookup_design_indices(D, universe, "GP training design")
  test_idx <- as.integer(test_idx)
  y_train <- truth$g[train_idx]
  fit <- wcrit29_fit_ok_estimated_c(
    D, y_train, universe$q, c_grid = c_grid,
    jitter = jitter, rank_tol = rank_tol
  )
  base <- data.frame(
    model_family = "mallows_gp",
    scenario_id = truth$scenario_id,
    fit_mode = fit_mode,
    c_true = truth$c_true,
    theta_true = truth$theta_true,
    fit_success = isTRUE(fit$success),
    fit_error = if (isTRUE(fit$success)) NA_character_ else fit$error,
    train_n = nrow(D),
    common_test_n = length(test_idx),
    path_seed = truth$path_seed,
    path_mean = truth$path_mean,
    path_sd = truth$path_sd,
    path_empirically_centered = truth$empirically_centered,
    path_empirically_standardized = truth$empirically_standardized,
    truth_kernel_jitter = truth$truth_kernel_jitter,
    stringsAsFactors = FALSE
  )
  if (!isTRUE(fit$success)) {
    metrics <- cbind(base, data.frame(
      c_hat = NA_real_, theta_hat = NA_real_, c_boundary_hit = NA,
      c_grid_best = NA_real_, c_continuous_refined = NA,
      sigma2_hat = NA_real_, reml = NA_real_, gp_df = NA_integer_,
      gp_jitter = jitter, condition_R = NA_real_,
      nrmse_full = NA_real_, nrmse_common = NA_real_, rmse_full = NA_real_, rmse_common = NA_real_,
      nll_full = NA_real_, nll_common = NA_real_, coverage_full = NA_real_, coverage_common = NA_real_,
      interval_length_full = NA_real_, interval_length_common = NA_real_,
      ipv_factor_full = NA_real_, ipv_factor_common = NA_real_,
      ipv_hat_full = NA_real_, ipv_hat_common = NA_real_
    ))
    return(if (return_fit) list(metrics = metrics, fit = fit) else metrics)
  }
  pred_full <- wcrit29_predict_ok(fit, universe$P)
  pred_test <- wcrit29_predict_ok(fit, universe$P[test_idx, , drop = FALSE])
  score_full <- wcrit29_prediction_metrics(
    truth$g, pred_full, truth$path_sd, confidence, "normal"
  )
  score_test <- wcrit29_prediction_metrics(
    truth$g[test_idx], pred_test, truth$path_sd, confidence, "normal"
  )
  metrics <- cbind(base, data.frame(
    c_hat = fit$c,
    theta_hat = fit$theta,
    c_boundary_hit = fit$c_boundary_hit,
    c_grid_best = fit$c_grid_best,
    c_continuous_refined = fit$c_continuous_refined,
    sigma2_hat = fit$sigma2,
    reml = fit$reml,
    gp_df = fit$df,
    gp_jitter = fit$jitter,
    condition_R = fit$condition_R,
    nrmse_full = score_full$nrmse,
    nrmse_common = score_test$nrmse,
    rmse_full = score_full$rmse,
    rmse_common = score_test$rmse,
    nll_full = score_full$nll,
    nll_common = score_test$nll,
    coverage_full = score_full$coverage,
    coverage_common = score_test$coverage,
    interval_length_full = score_full$interval_length,
    interval_length_common = score_test$interval_length,
    ipv_factor_full = score_full$ipv_factor,
    ipv_factor_common = score_test$ipv_factor,
    ipv_hat_full = score_full$ipv,
    ipv_hat_common = score_test$ipv,
    stringsAsFactors = FALSE
  ))
  if (return_fit) list(metrics = metrics, fit = fit) else metrics
}

# ---- Numerical unit/smoke test ---------------------------------------------

wcrit29_model_eval_self_test <- function(verbose = TRUE) {
  universe <- wcrit29_model_build_universe(4L)
  H <- wcrit_sample_foldover_base_permutations(4L, 6L, seed = 290041L)
  D <- wcrit_foldover_design(H)
  stopifnot(wcrit_is_strict_foldover(D), qr(wcrit_pwo_matrix(D))$rank == 7L)
  test_idx <- setdiff(seq_len(universe$size), wcrit29_lookup_design_indices(D, universe))

  pwo <- wcrit29_generate_pwo_truth(
    universe, rep_id = 1L, master_seed = 20260831L, snr_values = c(2, 5)
  )
  z2 <- (pwo$y[["pwo_snr2"]] - pwo$f) / pwo$sigma[["2"]]
  z5 <- (pwo$y[["pwo_snr5"]] - pwo$f) / pwo$sigma[["5"]]
  stopifnot(max(abs(z2 - z5)) < 1e-12)
  pwo_eval <- wcrit29_evaluate_pwo_design(D, pwo, 2, universe, test_idx, return_fit = TRUE)
  stopifnot(
    isTRUE(pwo_eval$fit$success),
    pwo_eval$fit$rank == 7L,
    all(is.finite(unlist(pwo_eval$metrics[c(
      "nrmse_full", "nrmse_common", "beta_rmse_all",
      "coverage_common", "ipv_factor_full"
    )])))
  )
  noiseless <- wcrit29_fit_pwo_ols(D, pwo$f[wcrit29_lookup_design_indices(D, universe)])
  noiseless_pred <- wcrit29_predict_pwo_ols(noiseless, universe$P)
  stopifnot(max(abs(noiseless_pred$mean - pwo$f)) < 1e-9)

  cache <- wcrit29_prepare_gp_kernel(universe, 1)
  gp <- wcrit29_generate_gp_truth(
    universe, c_value = 1, rep_id = 1L, master_seed = 20260831L,
    kernel_cache = cache
  )
  stopifnot(
    identical(gp$empirically_centered, FALSE),
    identical(gp$empirically_standardized, FALSE),
    abs(gp$path_mean) > 1e-10,
    abs(gp$path_sd - 1) > 1e-10
  )
  gp_est <- wcrit29_evaluate_gp_design(
    D, gp, universe, test_idx, fit_mode = "estimated_c",
    c_grid = c(0.125, 0.25, 0.5, 1, 2, 4, 8), return_fit = TRUE
  )
  stopifnot(
    isTRUE(gp_est$fit$success),
    gp_est$fit$c >= 0.125 - 1e-10,
    gp_est$fit$c <= 8 + 1e-10,
    all(is.finite(unlist(gp_est$metrics[c(
      "c_hat", "nrmse_common", "nll_common", "ipv_factor_common"
    )]))),
    !anyDuplicated(names(pwo_eval$metrics)),
    !anyDuplicated(names(gp_est$metrics))
  )
  train_pred <- wcrit29_predict_ok(gp_est$fit, D)
  stopifnot(max(abs(train_pred$mean - gp$g[wcrit29_lookup_design_indices(D, universe)])) < 1e-6)

  out <- data.frame(
    test = c(
      "pwo_shared_noise", "pwo_noiseless_recovery", "pwo_metrics",
      "gp_raw_path", "gp_estimated_c", "gp_interpolation"
    ),
    pass = TRUE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(verbose)) {
    print(out, row.names = FALSE)
    message(sprintf(
      "[wcrit29 self-test] estimated c=%.6f; PWO test nRMSE=%.6f; GP test nRMSE=%.6f",
      gp_est$fit$c, pwo_eval$metrics$nrmse_common, gp_est$metrics$nrmse_common
    ))
  }
  invisible(list(
    audit = out,
    pwo = pwo_eval$metrics,
    gp_estimated = gp_est$metrics
  ))
}

if (sys.nframe() == 0L) {
  wcrit29_model_eval_self_test(verbose = TRUE)
}
