# Section 5.1 paper-only helpers.
#
# These functions are the minimal API used by the published four-method
# strict-foldover comparison.  They were extracted from the frozen
# Experiment-27 implementation; no unrestricted or Kendall-maximin arm is
# implemented here.

sec51_stop <- function(...) stop(sprintf(...), call. = FALSE)

sec51_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args)
  if (length(hit)) {
    return(normalizePath(sub("^--file=", "", args[hit[[1L]]]),
                         winslash = "/", mustWork = FALSE))
  }
  frames <- sys.frames()
  ofiles <- vapply(frames, function(fr) {
    if (is.null(fr$ofile)) NA_character_ else as.character(fr$ofile)
  }, character(1L))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles)) {
    return(normalizePath(ofiles[[length(ofiles)]], winslash = "/",
                         mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

sec51_code_dir <- dirname(sec51_script_path())
sec51_common_dir <- normalizePath(file.path(sec51_code_dir, "..", "common"),
                                  winslash = "/", mustWork = FALSE)
source(file.path(sec51_common_dir, "wcrit_common.R"))
source(file.path(sec51_common_dir, "wcrit_maximin_dist.R"))

sec51_is_main <- function(filename) {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args)
  if (!length(hit)) return(FALSE)
  called <- normalizePath(sub("^--file=", "", args[hit[[1L]]]),
                          winslash = "/", mustWork = FALSE)
  target <- normalizePath(file.path(sec51_code_dir, filename),
                          winslash = "/", mustWork = FALSE)
  identical(called, target)
}

sec51_project_root <- function() {
  normalizePath(file.path(sec51_code_dir, "..", ".."),
                winslash = "/", mustWork = FALSE)
}

sec51_int_env <- function(name, default, lower = -Inf) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || value < lower) {
    sec51_stop("%s must be an integer >= %s", name, format(lower))
  }
  value
}

sec51_num_env <- function(name, default, lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper) {
    sec51_stop("%s must lie in [%s, %s]", name, format(lower), format(upper))
  }
  value
}

sec51_parse_settings <- function(text) {
  tokens <- trimws(strsplit(text, ",", fixed = TRUE)[[1L]])
  tokens <- tokens[nzchar(tokens)]
  if (!length(tokens)) sec51_stop("The settings list is empty")
  rows <- lapply(tokens, function(token) {
    bits <- strsplit(tolower(token), "x", fixed = TRUE)[[1L]]
    values <- suppressWarnings(as.integer(bits))
    if (length(values) != 2L || anyNA(values) ||
        values[[1L]] < 3L || values[[2L]] < 2L) {
      sec51_stop("Bad setting '%s'; use m x n, for example 6x12", token)
    }
    data.frame(
      m = values[[1L]], n = values[[2L]],
      ratio = values[[2L]] / values[[1L]],
      setting = sprintf("m%d_n%d", values[[1L]], values[[2L]]),
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  if (anyDuplicated(answer$setting)) sec51_stop("The settings list has duplicates")
  if (any(answer$n %% 2L != 0L)) sec51_stop("Strict foldover requires even n")
  if (any(abs(answer$ratio - round(answer$ratio)) > 1e-12)) {
    sec51_stop("Every n/m ratio must be an integer")
  }
  answer$ratio <- as.integer(round(answer$ratio))
  answer
}

sec51_file_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (!requireNamespace("digest", quietly = TRUE)) sec51_stop("digest is required")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sec51_parallel_lapply <- function(X, FUN, ..., workers = 1L) {
  workers <- max(1L, as.integer(workers))
  if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      X, FUN, ..., mc.cores = workers, mc.preschedule = FALSE
    )
  } else {
    lapply(X, FUN, ...)
  }
}

sec51_validate_design <- function(D, m, n, context = "design") {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  if (!identical(dim(D), c(as.integer(n), as.integer(m)))) {
    sec51_stop("%s has dimension %s; expected %dx%d", context,
               paste(dim(D), collapse = "x"), n, m)
  }
  target <- seq_len(m)
  valid <- apply(D, 1L, function(row) identical(sort(as.integer(row)), target))
  if (!all(valid)) sec51_stop("%s contains a non-permutation row", context)
  if (anyDuplicated(wcrit_row_keys(D))) sec51_stop("%s contains duplicate rows", context)
  if (!wcrit_is_strict_foldover(D)) sec51_stop("%s is not strict foldover", context)
  D
}

sec51_position_matrix <- function(D) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  P <- matrix(NA_integer_, nrow = nrow(D), ncol = ncol(D))
  for (i in seq_len(nrow(D))) P[i, D[i, ]] <- seq_len(ncol(D))
  P
}

sec51_geometry_max <- function(m, geometry) {
  switch(
    geometry,
    hamming = as.numeric(m),
    position_l2 = sqrt(m * (m^2 - 1) / 3),
    sec51_stop("Unknown paper geometry '%s'", geometry)
  )
}

sec51_distance_matrix <- function(D, geometry) {
  switch(
    geometry,
    hamming = wcrit_distance_matrix(D, criterion = "hamming"),
    position_l2 = as.matrix(stats::dist(sec51_position_matrix(D))),
    sec51_stop("Unknown paper geometry '%s'", geometry)
  )
}

sec51_distance_summary <- function(D, geometry) {
  dm <- sec51_distance_matrix(D, geometry)
  distances <- dm[upper.tri(dm)]
  if (!length(distances)) sec51_stop("At least two runs are required")
  scale <- sec51_geometry_max(ncol(D), geometry)
  minimum <- min(distances)
  tolerance <- max(1e-12, abs(minimum) * 1e-12)
  list(
    minimum = as.numeric(minimum),
    minimum_multiplicity = as.integer(sum(abs(distances - minimum) <= tolerance)),
    mean = as.numeric(mean(distances)),
    mean_sq = as.numeric(mean(distances^2)),
    minimum_norm = as.numeric(minimum / scale),
    mean_sq_norm = as.numeric(mean(distances^2) / scale^2)
  )
}

sec51_maximin_score <- function(D, geometry, tie_weight) {
  result <- sec51_distance_summary(D, geometry)
  result$score <- result$minimum_norm + tie_weight * result$mean_sq_norm
  if (length(result$score) != 1L || !is.finite(result$score)) {
    sec51_stop("Non-finite %s maximin score", geometry)
  }
  result
}

sec51_mutate <- function(row, replace_prob) {
  if (stats::runif(1L) < replace_prob) return(sample.int(length(row), length(row)))
  positions <- sample.int(length(row), 2L)
  candidate <- row
  candidate[positions] <- rev(candidate[positions])
  candidate
}

sec51_valid_foldover_candidate <- function(H, row_id, candidate) {
  other <- H[-row_id, , drop = FALSE]
  if (!nrow(other)) return(TRUE)
  keys <- wcrit_row_keys(other)
  key <- paste(candidate, collapse = ",")
  reverse_key <- paste(rev(candidate), collapse = ",")
  !(key %in% keys || reverse_key %in% keys)
}

# General-purpose strict-foldover maximin SA used for Hamming and position-L2.
# Invalid duplicate/reversal proposals consume one logical proposal and cool
# the chain, exactly as in frozen Experiment 27.
sec51_sa_once <- function(H_initial, geometry, proposal_budget, move_seed, cfg) {
  set.seed(wcrit_safe_seed(move_seed))
  H <- as.matrix(H_initial)
  D <- wcrit_foldover_design(H)
  current <- sec51_maximin_score(D, geometry, cfg$tie_weight)
  best <- list(H = H, D = D, score = current)
  candidate_evaluations <- 0L
  invalid <- 0L
  for (proposal in seq_len(proposal_budget)) {
    temperature <- cfg$sa_T0 * cfg$sa_alpha^(proposal - 1L)
    row_id <- sample.int(nrow(H), 1L)
    candidate <- sec51_mutate(H[row_id, ], cfg$replace_prob)
    if (!sec51_valid_foldover_candidate(H, row_id, candidate)) {
      invalid <- invalid + 1L
      next
    }
    H_new <- H
    H_new[row_id, ] <- candidate
    D_new <- wcrit_foldover_design(H_new)
    score_new <- sec51_maximin_score(D_new, geometry, cfg$tie_weight)
    candidate_evaluations <- candidate_evaluations + 1L
    delta <- score_new$score - current$score
    accept <- is.finite(delta) &&
      (delta >= 0 || stats::runif(1L) < exp(delta / max(temperature, 1e-14)))
    if (accept) {
      H <- H_new
      D <- D_new
      current <- score_new
      if (current$score > best$score$score) best <- list(H = H, D = D, score = current)
    }
  }
  list(
    D = best$D, H = best$H, score = best$score,
    proposals = as.integer(proposal_budget),
    candidate_evaluations = candidate_evaluations,
    invalid_proposals = invalid
  )
}

sec51_mallows_metrics <- function(D, c_values = c(1, 4), jitter = 1e-10) {
  m <- ncol(D)
  n <- nrow(D)
  q <- choose(m, 2)
  dm <- wcrit_kendall_dmat(D)
  out <- list()
  for (c_value in c_values) {
    tag <- format(c_value, scientific = FALSE, trim = TRUE)
    theta <- c_value / q
    K <- exp(-theta * dm) + diag(jitter, n)
    chol_K <- try(chol(K), silent = TRUE)
    logdet <- if (inherits(chol_K, "try-error")) {
      -Inf
    } else {
      as.numeric(2 * sum(log(diag(chol_K))))
    }
    out[[paste0("mallows_theta_c", tag)]] <- theta
    out[[paste0("mallows_logdet_c", tag)]] <- logdet
    out[[paste0("mallows_det_root_c", tag)]] <-
      if (is.finite(logdet)) exp(logdet / n) else 0
  }
  out
}

sec51_design_metrics <- function(D, entry, cfg) {
  D <- sec51_validate_design(D, entry$m, entry$n, entry$method)
  fold <- wcrit_foldover_design_metrics(D, lambda = 0.5, invalid = "error")
  hamming <- sec51_distance_summary(D, "hamming")
  position_l2 <- sec51_distance_summary(D, "position_l2")
  X <- wcrit_pwo_matrix(D)
  rank <- qr(X, tol = cfg$rank_tol)$rank
  mallows <- sec51_mallows_metrics(D, cfg$mallows_c, cfg$mallows_jitter)
  row <- data.frame(
    protocol_hash = cfg$protocol_hash,
    method = entry$method,
    paper_label = entry$paper_label,
    design_class = "strict_foldover",
    m = entry$m, n = entry$n, ratio = entry$ratio, rep = entry$rep,
    seed = entry$seed, initial_seed = entry$initial_seed,
    move_seed = entry$move_seed,
    design_sha256 = wcrit_matrix_sha256(D),
    initial_design_sha256 = entry$initial_design_sha256,
    search_geometry = entry$search_geometry,
    search_objective = entry$search_objective,
    proposals = entry$proposals,
    candidate_objective_evaluations = entry$candidate_objective_evaluations,
    invalid_proposals = entry$invalid_proposals,
    elapsed_sec = entry$elapsed_sec,
    A = fold$A[[1L]], B = fold$B[[1L]], Phi = fold$phi_lambda[[1L]],
    hamming_min_norm = hamming$minimum_norm,
    position_l2_min_norm = position_l2$minimum_norm,
    pwo_p = ncol(X), pwo_rank = rank, pwo_full_rank = rank == ncol(X),
    pwo_ms_efficiency = wcrit_pwo_ms_efficiency(D),
    stringsAsFactors = FALSE
  )
  for (name in names(mallows)) row[[name]] <- mallows[[name]]
  row
}
