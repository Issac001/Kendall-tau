suppressPackageStartupMessages({
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gtools)
})

.wcrit_sourced_dir <- local({
  frames <- sys.frames()
  ofiles <- vapply(frames, function(fr) {
    if (!is.null(fr$ofile)) as.character(fr$ofile) else NA_character_
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles) == 0L) return(NA_character_)
  dirname(normalizePath(ofiles[[length(ofiles)]], winslash = "/", mustWork = FALSE))
})

wcrit_script_dir <- function() {
  # In the public paper-reproduction layout this helper is sourced from
  # code/common by section-specific drivers.  Resolve companion C++ and R
  # files relative to this file before considering the calling script.
  if (is.character(.wcrit_sourced_dir) && nzchar(.wcrit_sourced_dir) &&
      file.exists(file.path(.wcrit_sourced_dir, "wcrit_common.R"))) {
    return(.wcrit_sourced_dir)
  }
  args <- commandArgs(trailingOnly = FALSE)
  idx <- grep("^--file=", args)
  if (length(idx) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", args[idx[1]]), winslash = "/", mustWork = FALSE)))
  }
  if (file.exists(file.path(getwd(), "code", "WeightedCrit", "wcrit_common.R"))) {
    return(normalizePath(file.path(getwd(), "code", "WeightedCrit"), winslash = "/", mustWork = FALSE))
  }
  if (file.exists(file.path(getwd(), "R", "WeightedCrit", "wcrit_common.R"))) {
    return(normalizePath(file.path(getwd(), "R", "WeightedCrit"), winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

wcrit_project_root <- function() {
  normalizePath(file.path(wcrit_script_dir(), "..", ".."), winslash = "/", mustWork = FALSE)
}

wcrit_path <- function(...) file.path(wcrit_project_root(), ...)

wcrit_file_in_script_dir <- function(file) {
  normalizePath(file.path(wcrit_script_dir(), file), winslash = "/", mustWork = FALSE)
}

wcrit_is_script_main <- function(file) {
  args <- commandArgs(trailingOnly = FALSE)
  idx <- grep("^--file=", args)
  if (length(idx) == 0) return(FALSE)
  called <- normalizePath(sub("^--file=", "", args[idx[1]]), winslash = "/", mustWork = FALSE)
  target <- normalizePath(file.path(wcrit_script_dir(), file), winslash = "/", mustWork = FALSE)
  identical(called, target)
}

wcrit_output_dir <- function(...) {
  out <- wcrit_path("outputs", "wcrit", ...)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

wcrit_parse_int_vec <- function(x, default) {
  if (!nzchar(trimws(x))) return(default)
  as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

wcrit_parse_num_vec <- function(x, default) {
  if (!nzchar(trimws(x))) return(default)
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

wcrit_parse_char_vec <- function(x, default) {
  if (!nzchar(trimws(x))) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

wcrit_bool_env <- function(name, default = FALSE) {
  x <- tolower(Sys.getenv(name, unset = if (isTRUE(default)) "1" else "0"))
  x %in% c("1", "true", "yes", "y")
}

wcrit_safe_seed <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x) || is.na(x)) return(1L)
  y <- as.integer(abs(floor(x)))
  if (is.na(y) || y < 1L) y <- 1L
  if (y > .Machine$integer.max) y <- as.integer((y %% .Machine$integer.max) + 1L)
  y
}

wcrit_hash_seed <- function(...) {
  key <- paste(..., sep = "|")
  vals <- utf8ToInt(enc2utf8(key))
  h <- 17
  for (v in vals) h <- (h * 131 + as.numeric(v)) %% 2147483629
  wcrit_safe_seed(h + 1)
}

wcrit_load_cpp <- local({
  loaded <- FALSE
  function() {
    if (isTRUE(loaded) &&
        exists("sa_optimize_H_weighted_foldover_cpp", mode = "function", inherits = TRUE) &&
        exists("kendall_dmat_between_cpp", mode = "function", inherits = TRUE) &&
        exists("pwo_matrix_cpp", mode = "function", inherits = TRUE)) {
      return(invisible(TRUE))
    }
    cpp_path <- file.path(wcrit_script_dir(), "sa_core.cpp")
    if (!file.exists(cpp_path)) stop("Cannot find sa_core.cpp at ", cpp_path)
    Rcpp::sourceCpp(cpp_path)
    loaded <<- TRUE
    invisible(TRUE)
  }
})

wcrit_constants <- function(m, n) {
  C1 <- floor(m * (m - 1) / 4)
  C2 <- n * m * (9 * m^3 - 14 * m^2 + 15 * m - 10) / (144 * (n - 1))
  U2 <- (n * m^2 * (m - 1)^2 - 4 * (n - 2) * (m * (m - 1) - 2)) / (8 * (n - 1))
  list(C1 = C1, C2 = C2, U2 = U2)
}

wcrit_kendall_dmat <- function(D1, D2 = D1) {
  wcrit_load_cpp()
  D1 <- as.matrix(D1)
  D2 <- as.matrix(D2)
  storage.mode(D1) <- "integer"
  storage.mode(D2) <- "integer"
  kendall_dmat_between_cpp(D1, D2)
}

wcrit_is_strict_foldover <- function(D) {
  D <- as.matrix(D)
  if (length(D) == 0L || nrow(D) < 2L || nrow(D) %% 2L != 0L) return(FALSE)
  storage.mode(D) <- "integer"
  if (anyNA(D)) return(FALSE)
  target <- seq_len(ncol(D))
  valid_rows <- apply(D, 1L, function(x) length(x) == length(target) && all(sort(x) == target))
  if (!all(valid_rows)) return(FALSE)
  keys <- apply(D, 1L, paste, collapse = ",")
  if (anyDuplicated(keys)) return(FALSE)
  reverse_keys <- apply(D, 1L, function(x) paste(rev(x), collapse = ","))
  all(reverse_keys %in% keys)
}

wcrit_design_metrics <- function(D, lambda = NA_real_) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  n <- nrow(D)
  m <- ncol(D)
  const <- wcrit_constants(m, n)
  if (n <= 1L) {
    k_min <- choose(m, 2)
    k_ave <- 0
    k_m2 <- 0
  } else {
    dm <- wcrit_kendall_dmat(D)
    dvals <- dm[upper.tri(dm)]
    k_min <- min(dvals)
    k_ave <- mean(dvals)
    k_m2 <- 2 * sum(dvals^2) / (n * (n - 1))
  }
  A <- if (const$C1 > 1) (k_min - 1) / (const$C1 - 1) else 1
  B <- (const$U2 - k_m2) / (const$U2 - const$C2)
  phi <- if (is.finite(lambda)) lambda * A + (1 - lambda) * B else NA_real_
  phi_valid <- wcrit_is_strict_foldover(D)
  if (!phi_valid) {
    A <- NA_real_
    B <- NA_real_
    phi <- NA_real_
  }
  data.frame(
    k_min = as.numeric(k_min),
    k_ave = as.numeric(k_ave),
    k_m2 = as.numeric(k_m2),
    A = as.numeric(A),
    B = as.numeric(B),
    phi_lambda = as.numeric(phi),
    C1 = as.numeric(const$C1),
    C2 = as.numeric(const$C2),
    U2 = as.numeric(const$U2),
    phi_valid = phi_valid,
    design_class = if (phi_valid) "strict_foldover" else "nonfoldover_or_reduced",
    stringsAsFactors = FALSE
  )
}

wcrit_sample_unique_permutations <- function(m, n_target, seed = NULL, exclude = NULL, max_tries = 1000000L) {
  if (log(n_target) > lfactorial(m) + 1e-12) {
    stop("n_target exceeds m!; cannot sample unique permutations")
  }
  if (!is.null(seed)) set.seed(wcrit_safe_seed(seed))
  seen <- new.env(parent = emptyenv())
  if (!is.null(exclude) && length(exclude) > 0) {
    exclude <- as.matrix(exclude)
    for (i in seq_len(nrow(exclude))) assign(paste(exclude[i, ], collapse = ","), TRUE, envir = seen)
  }
  out <- matrix(NA_integer_, nrow = n_target, ncol = m)
  i <- 1L
  tries <- 0L
  while (i <= n_target) {
    tries <- tries + 1L
    if (tries > max_tries) stop("Exceeded max_tries while sampling unique permutations")
    perm <- sample.int(m, m)
    key <- paste(perm, collapse = ",")
    if (!exists(key, envir = seen, inherits = FALSE)) {
      assign(key, TRUE, envir = seen)
      out[i, ] <- perm
      i <- i + 1L
    }
  }
  out
}

wcrit_candidate_permutations <- function(m, n_target, seed = NULL, exclude = NULL, enumerate_limit = 5040L) {
  if (!is.null(seed)) set.seed(wcrit_safe_seed(seed))
  exclude_keys <- character(0)
  if (!is.null(exclude) && length(exclude) > 0) {
    exclude <- as.matrix(exclude)
    exclude_keys <- apply(exclude, 1, paste, collapse = ",")
  }
  if (factorial(m) <= enumerate_limit) {
    allp <- gtools::permutations(m, m, seq_len(m))
    keys <- apply(allp, 1, paste, collapse = ",")
    keep <- !(keys %in% exclude_keys)
    allp <- allp[keep, , drop = FALSE]
    if (nrow(allp) <= n_target) return(allp)
    return(allp[sample.int(nrow(allp), n_target), , drop = FALSE])
  }
  wcrit_sample_unique_permutations(m, n_target, seed = seed, exclude = exclude)
}

wcrit_foldover_perm <- function(perm, m = length(perm)) {
  as.integer(rev(perm))
}

wcrit_foldover_design <- function(H) {
  H <- as.matrix(H)
  storage.mode(H) <- "integer"
  H_fold <- t(apply(H, 1, wcrit_foldover_perm, m = ncol(H)))
  storage.mode(H_fold) <- "integer"
  rbind(H, H_fold)
}

# The theoretical C1/C2/U2 normalization is valid only for a duplicate-free,
# even-run design that contains the reversal of every run.  Keep this check in
# one place so non-foldover baselines and odd-run reductions are never ranked
# by a foldover-only Phi value.

wcrit_design_class <- function(D) {
  if (wcrit_is_strict_foldover(D)) "strict_foldover" else "nonfoldover_or_reduced"
}

wcrit_foldover_design_metrics <- function(D, lambda = NA_real_, invalid = c("na", "error")) {
  invalid <- match.arg(invalid)
  out <- wcrit_design_metrics(D, lambda = lambda)
  valid <- wcrit_is_strict_foldover(D)
  if (!valid && identical(invalid, "error")) {
    stop("Foldover-normalized A, B, and Phi require a strict even foldover design")
  }
  if (!valid) {
    out$A <- NA_real_
    out$B <- NA_real_
    out$phi_lambda <- NA_real_
  }
  out$phi_valid <- valid
  out$design_class <- wcrit_design_class(D)
  out
}

wcrit_row_keys <- function(D) {
  apply(as.matrix(D), 1, paste, collapse = ",")
}

wcrit_matrix_sha256 <- function(D) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required to compute auditable SHA-256 fingerprints")
  }
  payload <- paste(wcrit_row_keys(as.matrix(D)), collapse = "\n")
  digest::digest(payload, algo = "sha256", serialize = FALSE)
}

wcrit_sample_foldover_base_permutations <- function(m, K, seed = NULL, max_tries = 1000000L) {
  if (!is.null(seed)) set.seed(wcrit_safe_seed(seed))
  if (log(2 * K) > lfactorial(m) + 1e-12) {
    stop("2*K exceeds m!; cannot build a duplicate-free foldover design")
  }
  seen <- new.env(parent = emptyenv())
  out <- matrix(NA_integer_, nrow = K, ncol = m)
  i <- 1L
  tries <- 0L
  while (i <= K) {
    tries <- tries + 1L
    if (tries > max_tries) stop("Exceeded max_tries while sampling foldover base permutations")
    perm <- sample.int(m, m)
    fold <- wcrit_foldover_perm(perm, m)
    key <- paste(perm, collapse = ",")
    fkey <- paste(fold, collapse = ",")
    if (!exists(key, envir = seen, inherits = FALSE) &&
        !exists(fkey, envir = seen, inherits = FALSE)) {
      assign(key, TRUE, envir = seen)
      assign(fkey, TRUE, envir = seen)
      out[i, ] <- perm
      i <- i + 1L
    }
  }
  out
}

wcrit_reduce_design_by_one <- function(D_full, method = c("min_impact", "random")) {
  method <- match.arg(method)
  D_full <- as.matrix(D_full)
  storage.mode(D_full) <- "integer"
  if (nrow(D_full) <= 1L) stop("Cannot reduce a design with <= 1 rows")
  idx_remove <- if (identical(method, "random")) {
    sample.int(nrow(D_full), 1L)
  } else {
    dm <- wcrit_kendall_dmat(D_full)
    which.min(rowSums(dm))
  }
  list(D = D_full[-idx_remove, , drop = FALSE], removed = as.integer(idx_remove))
}

wcrit_finalize_foldover_design <- function(H, n_target, reduce_method = "min_impact") {
  D_full <- wcrit_foldover_design(H)
  removed <- NA_integer_
  if (nrow(D_full) > n_target) {
    if (nrow(D_full) != n_target + 1L) {
      stop("Foldover design row count is incompatible with n_target")
    }
    red <- wcrit_reduce_design_by_one(D_full, method = reduce_method)
    D_full <- red$D
    removed <- red$removed
  }
  storage.mode(D_full) <- "integer"
  list(D = D_full, removed = removed)
}

wcrit_build_weighted_design <- function(m, n, lambda, seed,
                                        max_iter = 6000L,
                                        restarts = 8L,
                                        T0 = 1.0,
                                        alpha = 0.997,
                                        no_improve_stop = 800L,
                                        max_seconds = 0,
                                        trace = FALSE,
                                        foldover = TRUE,
                                        foldover_reduce_method = "min_impact") {
  wcrit_load_cpp()
  if (!isTRUE(foldover)) {
    stop("The paper-reproduction build supports the weighted criterion only in the foldover class")
  }
  max_iter <- as.integer(max_iter)
  restarts <- max(1L, as.integer(restarts))
  base_iter <- max(1L, max_iter %/% restarts)
  rem <- max_iter %% restarts

  best <- NULL
  total_iter <- 0L
  start_time <- proc.time()[["elapsed"]]
  for (r in seq_len(restarts)) {
    iter_r <- base_iter + as.integer(r <= rem)
    seed_r <- wcrit_hash_seed(seed, "weighted-sa", m, n, lambda, r)
    tr_vec <- if (isTRUE(trace)) numeric(iter_r) else NULL
    n_even <- ifelse(n %% 2L == 0L, as.integer(n), as.integer(n + 1L))
    K <- as.integer(n_even / 2L)
    H_init <- wcrit_sample_foldover_base_permutations(m = m, K = K, seed = seed_r)
    fit <- sa_optimize_H_weighted_foldover_cpp(
      H_init = H_init,
      m = as.integer(m),
      lambda = as.numeric(lambda),
      T0 = as.numeric(T0),
      alpha = as.numeric(alpha),
      max_iter = as.integer(iter_r),
      no_improve_stop = as.integer(no_improve_stop),
      duplicate_free = TRUE,
      verbose = FALSE,
      incremental = TRUE,
      max_seconds = as.numeric(max_seconds),
      rng_seed = as.integer(seed_r),
      trace = tr_vec
    )
    fin <- wcrit_finalize_foldover_design(
      fit$H_best,
      n_target = as.integer(n),
      reduce_method = foldover_reduce_method
    )
    D_fit <- fin$D
    removed_row <- fin$removed
    total_iter <- total_iter + as.integer(fit$n_iter)
    met_fit <- wcrit_foldover_design_metrics(D_fit, lambda = lambda, invalid = "na")
    obj_fit <- if (isTRUE(met_fit$phi_valid[[1L]])) {
      as.numeric(met_fit$phi_lambda[[1L]])
    } else {
      NA_real_
    }
    # For an odd target n, the search itself was performed on the valid n+1
    # parent foldover design.  Use that search score only to select a restart;
    # the reduced final design deliberately receives no Phi value.
    selection_score <- if (is.finite(obj_fit)) obj_fit else as.numeric(fit$best_obj)
    cand <- c(
      fit,
      list(
        D_final = D_fit,
        metrics_final = met_fit,
        obj_final = obj_fit,
        selection_score = selection_score,
        removed_row = removed_row
      )
    )
    if (is.null(best) || selection_score > best$selection_score) best <- cand
  }
  D <- as.matrix(best$D_final)
  storage.mode(D) <- "integer"
  met <- best$metrics_final
  list(
    D = D,
    H = if (!is.null(best$H_best)) as.matrix(best$H_best) else NULL,
    metrics = met,
    objective = as.numeric(best$obj_final),
    search_objective = as.numeric(best$best_obj),
    n_iter = total_iter,
    elapsed_sec = as.numeric(proc.time()[["elapsed"]] - start_time),
    foldover = wcrit_is_strict_foldover(D),
    requested_foldover = TRUE,
    removed_row = best$removed_row
  )
}

wcrit_pwo_matrix <- function(D) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  has_cpp <- tryCatch({
    wcrit_load_cpp()
    exists("pwo_matrix_cpp", mode = "function", inherits = TRUE)
  }, error = function(e) FALSE)
  if (isTRUE(has_cpp)) {
    return(pwo_matrix_cpp(D, intercept = TRUE))
  }
  n <- nrow(D)
  m <- ncol(D)
  pairs <- utils::combn(seq_len(m), 2)
  Z <- matrix(0, nrow = n, ncol = ncol(pairs))
  for (i in seq_len(n)) {
    pos <- integer(m)
    pos[D[i, ]] <- seq_len(m)
    Z[i, ] <- ifelse(pos[pairs[1, ]] < pos[pairs[2, ]], 1, -1)
  }
  cbind(Intercept = 1, Z)
}

wcrit_pwo_logdet <- function(D) {
  X <- wcrit_pwo_matrix(D)
  if (nrow(X) < ncol(X)) return(-Inf)
  XtX <- crossprod(X)
  det_obj <- determinant(XtX, logarithm = TRUE)
  if (as.numeric(det_obj$sign) <= 0) return(-Inf)
  as.numeric(det_obj$modulus)
}

wcrit_pwo_ms <- function(D) {
  X <- wcrit_pwo_matrix(D)
  M <- crossprod(X)
  as.numeric(sum(diag(M %*% M)))
}

# For the intercept + {-1,+1} PWO coding, this is
# tr{E(x x')^2} under the uniform distribution on all m! permutations.
wcrit_pwo_full_ms_normalized <- function(m) {
  m <- as.integer(m)
  if (!is.finite(m) || m < 2L) stop("m must be an integer >= 2")
  1 + choose(m, 2) + m * (m - 1) * (m - 2) / 9
}

# Absolute MS efficiency, valid for any permutation design.  Unlike the older
# sample-relative helper, this does not normalize by the best observed method.
wcrit_pwo_ms_efficiency <- function(D) {
  D <- as.matrix(D)
  if (nrow(D) < 1L || ncol(D) < 2L) return(NA_real_)
  denom <- wcrit_pwo_ms(D) / (nrow(D)^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  as.numeric(wcrit_pwo_full_ms_normalized(ncol(D)) / denom)
}

wcrit_read_design_matrix <- function(design_path) {
  if (is.na(design_path) || !nzchar(design_path) || !file.exists(design_path)) {
    return(NULL)
  }
  obj <- readRDS(design_path)
  if (is.list(obj) && !is.null(obj$D)) obj$D else obj
}

wcrit_attach_pwo_efficiencies <- function(raw_df, recompute_ms_from_designs = FALSE) {
  need_ms <- recompute_ms_from_designs ||
    !("pwo_ms" %in% names(raw_df)) ||
    all(!is.finite(raw_df$pwo_ms))
  if (need_ms) {
    if (!"design_path" %in% names(raw_df)) {
      stop("raw_df has no design_path column; cannot compute pwo_ms")
    }
    message(sprintf("Computing PWO MS for %d designs...", nrow(raw_df)))
    raw_df$pwo_ms <- vapply(seq_len(nrow(raw_df)), function(i) {
      D <- wcrit_read_design_matrix(raw_df$design_path[[i]])
      if (is.null(D)) return(NA_real_)
      wcrit_pwo_ms(D)
    }, numeric(1))
  }

  raw_df %>%
    group_by(m, n) %>%
    mutate(
      max_pwo_logdet = if (any(is.finite(pwo_logdet))) {
        max(pwo_logdet[is.finite(pwo_logdet)])
      } else {
        NA_real_
      },
      pwo_efficiency = ifelse(
        is.finite(pwo_logdet) & is.finite(max_pwo_logdet),
        exp((pwo_logdet - max_pwo_logdet) / (choose(m, 2) + 1)),
        0
      ),
      max_pwo_ms = if (any(is.finite(pwo_ms))) max(pwo_ms[is.finite(pwo_ms)]) else NA_real_,
      pwo_ms_efficiency = ifelse(
        is.finite(pwo_ms) & pwo_ms > 0,
        (1 + choose(m, 2) + m * (m - 1) * (m - 2) / 9) / (pwo_ms / n^2),
        NA_real_
      )
    ) %>%
    ungroup()
}

wcrit_write_pwo_derived_outputs <- function(raw_df, out_dir) {
  raw_path <- file.path(out_dir, "pwo_raw.csv")
  raw_no_eff_path <- file.path(out_dir, "pwo_raw_no_efficiency.csv")
  write.csv(raw_df, raw_path, row.names = FALSE)
  drop_eff <- c("max_pwo_logdet", "pwo_efficiency", "max_pwo_ms", "pwo_ms_efficiency")
  raw_no <- raw_df[, setdiff(names(raw_df), drop_eff), drop = FALSE]
  write.csv(raw_no, raw_no_eff_path, row.names = FALSE)

  summary_df <- wcrit_summarise_mean_sd(
    raw_df,
    group_cols = c("m", "n", "lambda"),
    metric_cols = c(
      "pwo_logdet", "pwo_efficiency", "pwo_ms", "pwo_ms_efficiency",
      "k_min", "k_m2", "A", "B", "elapsed_sec"
    )
  )
  write.csv(summary_df, file.path(out_dir, "pwo_summary.csv"), row.names = FALSE)

  by_lambda <- raw_df %>%
    group_by(lambda) %>%
    summarise(
      mean_pwo_efficiency = mean(pwo_efficiency, na.rm = TRUE),
      se_pwo_efficiency = stats::sd(pwo_efficiency, na.rm = TRUE) / sqrt(dplyr::n()),
      mean_pwo_ms_efficiency = mean(pwo_ms_efficiency, na.rm = TRUE),
      se_pwo_ms_efficiency = stats::sd(pwo_ms_efficiency, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )
  write.csv(by_lambda, file.path(out_dir, "pwo_by_lambda.csv"), row.names = FALSE)

  p_d <- summary_df %>%
    ggplot(aes(x = lambda, y = pwo_efficiency_mean, color = factor(m), group = interaction(m, n))) +
    geom_line(alpha = 0.55, linewidth = 0.6) +
    geom_point(size = 1.6) +
    theme_bw(base_size = 10) +
    labs(x = expression(lambda), y = "Empirical PWO D-efficiency", color = "m")
  ggsave(file.path(out_dir, "pwo_efficiency_by_lambda.png"), p_d, width = 9, height = 5, dpi = 250)

  p_ms <- summary_df %>%
    ggplot(aes(x = lambda, y = pwo_ms_efficiency_mean, color = factor(m), group = interaction(m, n))) +
    geom_line(alpha = 0.55, linewidth = 0.6) +
    geom_point(size = 1.6) +
    theme_bw(base_size = 10) +
    labs(x = expression(lambda), y = "Empirical PWO MS-efficiency", color = "m")
  ggsave(file.path(out_dir, "pwo_ms_efficiency_by_lambda.png"), p_ms, width = 9, height = 5, dpi = 250)

  invisible(list(raw = raw_df, summary = summary_df, out_dir = out_dir))
}

wcrit_mallows_logdet <- function(D, theta, jitter = 1e-10) {
  dm <- wcrit_kendall_dmat(D)
  K <- exp(-theta * dm) + diag(jitter, nrow(dm))
  cholK <- try(chol(K), silent = TRUE)
  if (inherits(cholK, "try-error")) return(-Inf)
  as.numeric(2 * sum(log(diag(cholK))))
}

wcrit_fit_gp_mallows <- function(train_D, y, theta_grid = c(0.01, 0.05, 0.1, 0.2),
                                 noise = 1e-6, jitter = 1e-10) {
  y <- as.numeric(y)
  y_mean <- mean(y)
  y_sd <- stats::sd(y)
  if (!is.finite(y_sd) || y_sd < 1e-12) y_sd <- 1
  yz <- (y - y_mean) / y_sd
  dm <- wcrit_kendall_dmat(train_D)
  best <- list(theta = NA_real_, lml = -Inf, cholK = NULL, alpha = NULL)
  for (theta in theta_grid) {
    K <- exp(-theta * dm) + diag(noise + jitter, nrow(dm))
    cholK <- try(chol(K), silent = TRUE)
    if (inherits(cholK, "try-error")) next
    alpha <- backsolve(cholK, forwardsolve(t(cholK), yz))
    logdet <- 2 * sum(log(diag(cholK)))
    lml <- as.numeric(-0.5 * (crossprod(yz, alpha) + logdet + length(yz) * log(2 * pi)))
    if (is.finite(lml) && lml > best$lml) {
      best <- list(theta = theta, lml = lml, cholK = cholK, alpha = alpha)
    }
  }
  if (!is.finite(best$lml)) stop("Mallows GP fitting failed")
  list(
    theta = best$theta,
    lml = best$lml,
    cholK = best$cholK,
    alpha = best$alpha,
    train_D = as.matrix(train_D),
    y_mean = y_mean,
    y_sd = y_sd,
    noise = noise
  )
}

wcrit_predict_gp_mallows <- function(model, test_D, return_var = FALSE) {
  dtest <- wcrit_kendall_dmat(as.matrix(test_D), model$train_D)
  Kstar <- exp(-model$theta * dtest)
  mu_z <- as.numeric(Kstar %*% model$alpha)
  mu <- model$y_mean + model$y_sd * mu_z
  if (!isTRUE(return_var)) return(mu)
  v <- forwardsolve(t(model$cholK), t(Kstar))
  var_z <- pmax(1e-12, 1 - colSums(v^2))
  data.frame(mean = mu, sd = sqrt(var_z) * model$y_sd)
}

wcrit_expected_improvement_min <- function(mu, sigma, best_y) {
  sigma <- pmax(as.numeric(sigma), 1e-12)
  z <- (best_y - mu) / sigma
  (best_y - mu) * stats::pnorm(z) + sigma * stats::dnorm(z)
}

wcrit_append_csv <- function(row, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    row,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(path),
    append = file.exists(path),
    quote = TRUE
  )
}

wcrit_write_config <- function(path, values) {
  cfg <- data.frame(
    key = names(values),
    value = vapply(values, function(x) paste(x, collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(cfg, path, row.names = FALSE)
}

wcrit_to_json_value <- function(x) {
  if (is.null(x)) return("null")
  if (is.logical(x) && length(x) == 1L) return(ifelse(isTRUE(x), "true", "false"))
  if (is.numeric(x) && length(x) == 1L) return(ifelse(is.finite(x), as.character(x), "null"))
  if (is.character(x) && length(x) == 1L) {
    esc <- gsub("\\\\", "\\\\\\\\", x)
    esc <- gsub("\"", "\\\\\"", esc)
    return(sprintf("\"%s\"", esc))
  }

  if (is.atomic(x) && length(x) > 1L) {
    vals <- vapply(as.list(x), wcrit_to_json_value, character(1))
    return(sprintf("[%s]", paste(vals, collapse = ",")))
  }

  if (is.list(x)) {
    if (is.null(names(x))) {
      vals <- vapply(x, wcrit_to_json_value, character(1))
      return(sprintf("[%s]", paste(vals, collapse = ",")))
    }
    kv <- vapply(names(x), function(k) {
      esc_key <- gsub("\\\\", "\\\\\\\\", k)
      esc_key <- gsub("\"", "\\\\\"", esc_key)
      sprintf("\"%s\":%s", esc_key, wcrit_to_json_value(x[[k]]))
    }, character(1))
    return(sprintf("{%s}", paste(kv, collapse = ",")))
  }

  "\"\""
}

wcrit_write_json <- function(path, values) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  json_txt <- wcrit_to_json_value(values)
  writeLines(json_txt, con = path, useBytes = TRUE)
}

wcrit_lambda_key <- function(x) sprintf("%.10g", as.numeric(x))

wcrit_row_key <- function(...) paste(..., sep = "|")

wcrit_summarise_mean_sd <- function(df, group_cols, metric_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      across(
        all_of(metric_cols),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~ifelse(dplyr::n() > 1, stats::sd(.x, na.rm = TRUE), 0),
          se = ~ifelse(dplyr::n() > 1, stats::sd(.x, na.rm = TRUE) / sqrt(sum(is.finite(.x))), 0)
        ),
        .names = "{.col}_{.fn}"
      ),
      n_runs = dplyr::n(),
      .groups = "drop"
    )
}
