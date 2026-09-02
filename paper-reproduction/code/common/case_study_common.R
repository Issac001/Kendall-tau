# Shared, auditable helpers for the two applied case studies.

case_file_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

case_matrix_sha256 <- function(D) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  payload <- paste(
    paste(dim(D), collapse = "x"),
    paste(apply(D, 1L, paste, collapse = ","), collapse = ";"),
    sep = "|"
  )
  digest::digest(payload, algo = "sha256", serialize = FALSE)
}

case_safe_mean <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) return(NA_real_)
  mean(x[is.finite(x)])
}

case_safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

case_safe_se <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

case_validate_design <- function(D, m, n, context = "design") {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  if (!identical(dim(D), c(as.integer(n), as.integer(m)))) {
    stop(context, " has dimension ", paste(dim(D), collapse = "x"),
         "; expected ", n, "x", m)
  }
  target <- seq_len(m)
  valid <- apply(D, 1L, function(z) identical(sort(as.integer(z)), target))
  if (!all(valid)) stop(context, " contains a non-permutation row")
  if (anyDuplicated(wcrit_row_keys(D))) stop(context, " contains duplicate rows")
  D
}

case_position_matrix <- function(D) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  P <- matrix(NA_integer_, nrow(D), ncol(D))
  for (i in seq_len(nrow(D))) P[i, D[i, ]] <- seq_len(ncol(D))
  P
}

case_geometry_max <- function(m, geometry) {
  switch(
    geometry,
    kendall = choose(m, 2),
    hamming = m,
    l2_position = sqrt(m * (m^2 - 1) / 3),
    stop("Unknown geometry: ", geometry)
  )
}

case_geometry_dmat <- function(D, geometry) {
  switch(
    geometry,
    kendall = wcrit_kendall_dmat(D),
    hamming = wcrit_distance_matrix(D, criterion = "hamming"),
    l2_position = as.matrix(stats::dist(case_position_matrix(D))),
    stop("Unknown geometry: ", geometry)
  )
}

case_geometry_summary <- function(D, geometry) {
  dm <- case_geometry_dmat(D, geometry)
  d <- as.numeric(dm[upper.tri(dm)])
  if (!length(d)) stop("A design needs at least two rows")
  scale <- case_geometry_max(ncol(D), geometry)
  dmin <- min(d)
  tol <- max(1e-12, abs(dmin) * 1e-12)
  list(
    minimum = dmin,
    minimum_norm = dmin / scale,
    minimum_multiplicity = sum(abs(d - dmin) <= tol),
    mean = mean(d),
    mean_sq = mean(d^2),
    mean_sq_norm = mean(d^2) / scale^2,
    dvec = d
  )
}

case_maximin_score <- function(D, geometry, tie_weight = 1e-4) {
  out <- case_geometry_summary(D, geometry)
  out$score <- out$minimum_norm + tie_weight * out$mean_sq_norm
  out
}

case_mutate_permutation <- function(row, replace_prob = 0.25) {
  if (stats::runif(1L) < replace_prob) return(sample.int(length(row), length(row)))
  ij <- sample.int(length(row), 2L)
  out <- row
  out[ij] <- rev(out[ij])
  out
}

# The application sections use this unrestricted SA engine for Hamming and
# component-position L2.  The initial design and proposal seed are shared
# across the two geometries; invalid duplicate proposals count as attempts but
# not as objective evaluations.
case_run_unrestricted_sa <- function(D_init, geometry, budget, move_seed,
                                     tie_weight = 1e-4, replace_prob = 0.25,
                                     T0 = 1, Tmin = 1e-8,
                                     proposal_multiplier = 30L) {
  if (!(geometry %in% c("hamming", "l2_position"))) {
    stop("The paper-only unrestricted SA supports hamming or l2_position")
  }
  D <- as.matrix(D_init)
  storage.mode(D) <- "integer"
  m <- ncol(D)
  n <- nrow(D)
  case_validate_design(D, m, n, paste0("initial ", geometry, " design"))
  budget <- as.integer(budget)
  if (budget < 1L) stop("budget must be positive")
  set.seed(wcrit_safe_seed(move_seed))
  current <- case_maximin_score(D, geometry, tie_weight)
  best <- list(D = D, score = current)
  n_obj_eval <- 1L
  proposals <- 0L
  invalid <- 0L
  noop <- 0L
  accepted <- 0L
  proposal_cap <- max(100L, as.integer(proposal_multiplier * budget))
  while (n_obj_eval < budget && proposals < proposal_cap) {
    proposals <- proposals + 1L
    row_id <- sample.int(n, 1L)
    candidate <- case_mutate_permutation(D[row_id, ], replace_prob)
    if (identical(as.integer(candidate), as.integer(D[row_id, ]))) {
      noop <- noop + 1L
      next
    }
    other_keys <- wcrit_row_keys(D[-row_id, , drop = FALSE])
    if (paste(candidate, collapse = ",") %in% other_keys) {
      invalid <- invalid + 1L
      next
    }
    D_new <- D
    D_new[row_id, ] <- candidate
    score_new <- case_maximin_score(D_new, geometry, tie_weight)
    n_obj_eval <- n_obj_eval + 1L
    frac <- if (budget <= 2L) 1 else (n_obj_eval - 2) / (budget - 2)
    temperature <- T0 * (Tmin / T0)^min(1, max(0, frac))
    delta <- score_new$score - current$score
    accept <- is.finite(delta) &&
      (delta >= 0 || stats::runif(1L) < exp(delta / max(temperature, 1e-14)))
    if (accept) {
      D <- D_new
      current <- score_new
      accepted <- accepted + 1L
      if (current$score > best$score$score) best <- list(D = D, score = current)
    }
  }
  list(
    D = best$D,
    score = best$score,
    n_obj_eval = n_obj_eval,
    proposals = proposals,
    invalid_proposals = invalid,
    noop_proposals = noop,
    accepted_proposals = accepted,
    budget_exhausted = n_obj_eval == budget
  )
}

case_reversal_pair_universe <- function(m) {
  allp <- gtools::permutations(m, m, seq_len(m))
  storage.mode(allp) <- "integer"
  keys <- wcrit_row_keys(allp)
  reverse_keys <- wcrit_row_keys(allp[, m:1L, drop = FALSE])
  canonical <- pmin(keys, reverse_keys)
  H <- allp[!duplicated(canonical), , drop = FALSE]
  storage.mode(H) <- "integer"
  H
}

# Exact FSA-KD design at the paper default lambda=0.5 for m=4,n=12.
case_exact_fsa_lambda05_m4n12 <- function() {
  m <- 4L
  n <- 12L
  H_universe <- case_reversal_pair_universe(m)
  combos <- utils::combn(nrow(H_universe), n / 2L)
  rows <- vector("list", ncol(combos))
  designs <- vector("list", ncol(combos))
  for (j in seq_len(ncol(combos))) {
    D <- wcrit_foldover_design(H_universe[combos[, j], , drop = FALSE])
    storage.mode(D) <- "integer"
    met <- wcrit_design_metrics(D, lambda = 0.5)
    kd <- case_geometry_summary(D, "kendall")
    designs[[j]] <- D
    rows[[j]] <- data.frame(
      candidate = j,
      k_min = met$k_min[[1L]],
      k_m2 = met$k_m2[[1L]],
      A = met$A[[1L]],
      B = met$B[[1L]],
      Phi050 = met$phi_lambda[[1L]],
      nearest_pair_multiplicity = kd$minimum_multiplicity,
      pair_indices = paste(combos[, j], collapse = ";"),
      design_sha256 = case_matrix_sha256(D),
      stringsAsFactors = FALSE
    )
  }
  tab <- dplyr::bind_rows(rows)
  pick05 <- with(tab, order(-Phi050, design_sha256))[1L]
  list(
    design = designs[[pick05]],
    optimum = cbind(lambda_label = "lambda050",
                    tab[pick05, , drop = FALSE]),
    all_candidates = tab
  )
}

case_oofa_oa_m4n12 <- function() {
  full <- gtools::permutations(4L, 4L, 1:4)
  D <- full[c(2, 3, 5, 7, 10, 12, 14, 15, 17, 20, 21, 24), , drop = FALSE]
  storage.mode(D) <- "integer"
  case_validate_design(D, 4L, 12L, "OofA-OA(12,4,2)")
}

case_apply_label_map <- function(D, label_map) {
  D <- as.matrix(D)
  out <- matrix(as.integer(label_map[D]), nrow(D), ncol(D))
  storage.mode(out) <- "integer"
  out
}

case_solve_chol <- function(U, B) {
  backsolve(U, forwardsolve(t(U), B))
}

case_trend_matrix <- function(D, trend = c("pwo", "intercept")) {
  trend <- match.arg(trend)
  if (trend == "pwo") return(wcrit_pwo_matrix(D))
  matrix(1, nrow(D), 1L, dimnames = list(NULL, "Intercept"))
}

# REML fit for a universal Mallows GP with an estimable nugget ratio.  Repeated
# permutations are allowed and provide the information needed to estimate
# experimental noise in the four-drug case.
case_fit_mallows_reml <- function(train_D, y, trend = c("pwo", "intercept"),
                                  theta_grid, nugget_grid,
                                  rank_tol = 1e-9, jitter = 1e-10,
                                  sigma_floor = 1e-12) {
  trend <- match.arg(trend)
  train_D <- as.matrix(train_D)
  storage.mode(train_D) <- "integer"
  y <- as.numeric(y)
  F <- case_trend_matrix(train_D, trend)
  n <- length(y)
  p <- ncol(F)
  if (nrow(train_D) != n || nrow(F) != n || n <= p) {
    return(list(success = FALSE, error = "dimension_or_df"))
  }
  dmat <- wcrit_kendall_dmat(train_D)
  best <- NULL
  for (theta in sort(unique(as.numeric(theta_grid)))) {
    if (!is.finite(theta) || theta <= 0) next
    K0 <- exp(-theta * dmat)
    for (nugget in sort(unique(as.numeric(nugget_grid)))) {
      if (!is.finite(nugget) || nugget < 0) next
      R <- K0
      diag(R) <- diag(R) + nugget + jitter
      U <- try(chol(R), silent = TRUE)
      if (inherits(U, "try-error")) next
      yw <- try(forwardsolve(t(U), y), silent = TRUE)
      Fw <- try(forwardsolve(t(U), F), silent = TRUE)
      if (inherits(yw, "try-error") || inherits(Fw, "try-error")) next
      qF <- qr(Fw, tol = rank_tol, LAPACK = FALSE)
      if (qF$rank < p) next
      beta <- qr.coef(qF, yw)
      if (any(!is.finite(beta))) next
      residual_w <- as.numeric(yw - Fw %*% beta)
      df <- n - p
      sse <- sum(residual_w^2)
      sigma2 <- max(sse / df, sigma_floor)
      diag_rf <- abs(diag(qr.R(qF))[seq_len(p)])
      if (any(!is.finite(diag_rf)) || any(diag_rf <= 0)) next
      logdet_R <- 2 * sum(log(diag(U)))
      logdet_F <- 2 * sum(log(diag_rf))
      reml <- -0.5 * (df * (log(2 * pi) + log(sigma2)) +
                        sse / sigma2 + logdet_R + logdet_F)
      if (!is.finite(reml)) next
      residual <- as.numeric(y - F %*% beta)
      alpha <- try(case_solve_chol(U, residual), silent = TRUE)
      chol_C <- try(chol(crossprod(Fw)), silent = TRUE)
      if (inherits(alpha, "try-error") || inherits(chol_C, "try-error")) next
      candidate <- list(
        success = TRUE,
        error = NA_character_,
        trend = trend,
        theta = theta,
        nugget = nugget,
        sigma2 = sigma2,
        reml = reml,
        beta = as.numeric(beta),
        alpha = as.numeric(alpha),
        chol_R = U,
        chol_C = chol_C,
        F_train = F,
        train_D = train_D,
        trend_rank = qF$rank,
        n = n,
        p = p,
        df = df,
        jitter = jitter
      )
      if (is.null(best) || candidate$reml > best$reml) best <- candidate
    }
  }
  if (is.null(best)) return(list(success = FALSE, error = "all_candidates_failed"))
  theta_values <- sort(unique(as.numeric(theta_grid)))
  theta_values <- theta_values[is.finite(theta_values) & theta_values > 0]
  nugget_values <- sort(unique(as.numeric(nugget_grid)))
  nugget_values <- nugget_values[is.finite(nugget_values) & nugget_values >= 0]
  best$sigma_floor_active <- isTRUE(best$sigma2 <= sigma_floor * (1 + 1e-12))
  best$theta_grid_min <- min(theta_values)
  best$theta_grid_max <- max(theta_values)
  best$nugget_grid_min <- min(nugget_values)
  best$nugget_grid_max <- max(nugget_values)
  best$theta_boundary_hit <- best$theta %in% range(theta_values)
  best$nugget_boundary_hit <- best$nugget %in% range(nugget_values)
  best
}

case_predict_mallows_reml <- function(model, test_D,
                                      include_observation_noise = TRUE,
                                      observation_noise_multiplier = NULL) {
  if (!isTRUE(model$success)) stop("Cannot predict from a failed model")
  test_D <- as.matrix(test_D)
  storage.mode(test_D) <- "integer"
  F_test <- case_trend_matrix(test_D, model$trend)
  dtest <- wcrit_kendall_dmat(test_D, model$train_D)
  Rstar <- exp(-model$theta * dtest)
  mean_pred <- as.numeric(F_test %*% model$beta + Rstar %*% model$alpha)
  B <- case_solve_chol(model$chol_R, t(Rstar))
  base <- 1 - colSums(t(Rstar) * B)
  u <- t(F_test) - crossprod(model$F_train, B)
  Cinv_u <- case_solve_chol(model$chol_C, u)
  trend_term <- colSums(u * Cinv_u)
  raw_factor <- base + trend_term
  if (any(raw_factor < -1e-8)) {
    stop("Materially negative kriging variance factor")
  }
  factor <- pmax(0, raw_factor)
  if (is.null(observation_noise_multiplier)) {
    observation_noise_multiplier <- if (isTRUE(include_observation_noise)) 1 else 0
  }
  observation_noise_multiplier <- as.numeric(observation_noise_multiplier)
  if (length(observation_noise_multiplier) != 1L ||
      !is.finite(observation_noise_multiplier) || observation_noise_multiplier < 0) {
    stop("observation_noise_multiplier must be one finite nonnegative number")
  }
  factor <- factor + observation_noise_multiplier * model$nugget
  variance <- pmax(0, model$sigma2 * factor)
  list(mean = mean_pred, variance = variance, sd = sqrt(variance))
}

case_fit_pwo_linear <- function(train_D, y, rank_tol = 1e-9,
                                sigma_floor = 1e-12) {
  X <- wcrit_pwo_matrix(train_D)
  y <- as.numeric(y)
  qx <- qr(X, tol = rank_tol, LAPACK = FALSE)
  if (qx$rank < ncol(X) || length(y) <= ncol(X)) {
    return(list(success = FALSE, error = "rank_or_df"))
  }
  beta <- qr.coef(qx, y)
  residual <- as.numeric(y - X %*% beta)
  df <- length(y) - ncol(X)
  sigma2 <- max(sum(residual^2) / df, sigma_floor)
  XtX_chol <- try(chol(crossprod(X)), silent = TRUE)
  if (inherits(XtX_chol, "try-error")) {
    return(list(success = FALSE, error = "xtx_cholesky"))
  }
  list(
    success = TRUE, error = NA_character_, trend = "pwo_only",
    beta = as.numeric(beta), sigma2 = sigma2, df = df,
    chol_XtX = XtX_chol, p = ncol(X), n = length(y)
  )
}

case_predict_pwo_linear <- function(model, test_D,
                                    include_observation_noise = TRUE,
                                    observation_noise_multiplier = NULL) {
  if (!isTRUE(model$success)) stop("Cannot predict from a failed PWO model")
  Xtest <- wcrit_pwo_matrix(test_D)
  mu <- as.numeric(Xtest %*% model$beta)
  Cinv_xt <- case_solve_chol(model$chol_XtX, t(Xtest))
  leverage <- colSums(t(Xtest) * Cinv_xt)
  if (is.null(observation_noise_multiplier)) {
    observation_noise_multiplier <- if (isTRUE(include_observation_noise)) 1 else 0
  }
  observation_noise_multiplier <- as.numeric(observation_noise_multiplier)
  if (length(observation_noise_multiplier) != 1L ||
      !is.finite(observation_noise_multiplier) || observation_noise_multiplier < 0) {
    stop("observation_noise_multiplier must be one finite nonnegative number")
  }
  factor <- leverage + observation_noise_multiplier
  variance <- model$sigma2 * pmax(0, factor)
  list(mean = mu, variance = variance, sd = sqrt(variance))
}

case_prediction_metrics <- function(y_true, prediction, confidence = 0.95,
                                    variance_floor = 1e-12,
                                    normalization_sd = NULL) {
  y_true <- as.numeric(y_true)
  mu <- as.numeric(prediction$mean)
  variance <- pmax(as.numeric(prediction$variance), variance_floor)
  err <- y_true - mu
  rmse <- sqrt(mean(err^2))
  ysd <- if (is.null(normalization_sd)) stats::sd(y_true) else as.numeric(normalization_sd)
  if (length(ysd) != 1L || !is.finite(ysd) || ysd <= 0) ysd <- NA_real_
  z <- stats::qnorm(0.5 + confidence / 2)
  sd_pred <- sqrt(variance)
  lower <- mu - z * sd_pred
  upper <- mu + z * sd_pred
  data.frame(
    rmse = rmse,
    nrmse = if (is.finite(ysd) && ysd > 1e-12) rmse / ysd else NA_real_,
    mae = mean(abs(err)),
    correlation = if (stats::sd(mu) > 0 && ysd > 0) stats::cor(mu, y_true) else NA_real_,
    spearman = suppressWarnings(stats::cor(mu, y_true, method = "spearman")),
    nll = mean(0.5 * log(2 * pi * variance) + err^2 / (2 * variance)),
    coverage = mean(y_true >= lower & y_true <= upper),
    interval_length = mean(upper - lower),
    stringsAsFactors = FALSE
  )
}

case_expected_improvement_min <- function(mu, sigma, best_y) {
  sigma <- pmax(as.numeric(sigma), 1e-12)
  z <- (best_y - mu) / sigma
  (best_y - mu) * stats::pnorm(z) + sigma * stats::dnorm(z)
}

case_paired_t_summary <- function(df, group_cols, value_col = "difference",
                                  confidence = 0.95) {
  alpha <- 1 - confidence
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n = sum(is.finite(.data[[value_col]])),
      mean = case_safe_mean(.data[[value_col]]),
      sd = case_safe_sd(.data[[value_col]]),
      se = case_safe_se(.data[[value_col]]),
      ci_low = dplyr::if_else(
        n >= 2L,
        mean - stats::qt(1 - alpha / 2, df = n - 1L) * se,
        NA_real_
      ),
      ci_high = dplyr::if_else(
        n >= 2L,
        mean + stats::qt(1 - alpha / 2, df = n - 1L) * se,
        NA_real_
      ),
      wins = sum(.data[[value_col]] > 0, na.rm = TRUE),
      ties = sum(abs(.data[[value_col]]) <= 1e-12, na.rm = TRUE),
      losses = sum(.data[[value_col]] < 0, na.rm = TRUE),
      .groups = "drop"
    )
}

case_prepare_output_dirs <- function(out_dir) {
  if (dir.exists(out_dir) || file.exists(out_dir)) {
    stop("Refusing to overwrite an existing output path: ", out_dir)
  }
  dirs <- list(
    root = out_dir,
    config = file.path(out_dir, "config"),
    raw = file.path(out_dir, "raw"),
    results = file.path(out_dir, "results"),
    figures = file.path(out_dir, "figures"),
    designs = file.path(out_dir, "designs"),
    checkpoints = file.path(out_dir, "checkpoints")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

case_append_stage <- function(path, stage, status, started, message = NA_character_) {
  row <- data.frame(
    stage = stage,
    status = status,
    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
    message = as.character(message),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    row, path, sep = ",", row.names = FALSE,
    col.names = !file.exists(path), append = file.exists(path), quote = TRUE
  )
  invisible(row)
}

case_write_file_manifest <- function(root, path) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[normalizePath(files, mustWork = FALSE) != normalizePath(path, mustWork = FALSE)]
  rel <- substring(files, nchar(normalizePath(root)) + 2L)
  manifest <- data.frame(
    relative_path = rel,
    size_bytes = as.numeric(file.info(files)$size),
    sha256 = vapply(files, case_file_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$relative_path), , drop = FALSE]
  utils::write.csv(manifest, path, row.names = FALSE)
  invisible(manifest)
}
