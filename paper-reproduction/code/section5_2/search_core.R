# Experiment 29: strict-foldover search engine with a shared proposal tape.
#
# This temporary module is deliberately independent of every existing
# WeightedCrit experiment.  It defines functions only; sourcing the file does
# not run an experiment or write any output.  The formal domain is m = 6 and
# even n <= 720.  All optimized methods in one (replication, n) bank receive
# exactly the same ordered half-design and exactly the same proposal tape.
# Invalid/no-op/reversal-orbit proposals consume one proposal and one cooling
# step, but do not trigger an objective evaluation.

wcrit29_stop <- function(...) stop(sprintf(...), call. = FALSE)

wcrit29_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    wcrit29_stop("Package '%s' is required", package)
  }
  invisible(TRUE)
}

wcrit29_safe_seed <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || !is.finite(x)) return(1L)
  x <- as.integer(abs(floor(x)) %% 2147483629)
  if (is.na(x) || x < 1L) 1L else x
}

wcrit29_hash_seed <- function(...) {
  bytes <- utf8ToInt(enc2utf8(paste(..., sep = "|")))
  h <- 17
  for (b in bytes) h <- (h * 131 + as.numeric(b)) %% 2147483629
  wcrit29_safe_seed(h + 1)
}

wcrit29_with_seed <- function(seed, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(wcrit29_safe_seed(seed))
  force(code)
}

wcrit29_sha256_text <- function(text) {
  wcrit29_require("digest")
  digest::digest(paste(text, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

wcrit29_row_keys <- function(D) {
  apply(as.matrix(D), 1L, paste, collapse = ",")
}

wcrit29_matrix_sha256 <- function(D, sort_rows = FALSE) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  keys <- wcrit29_row_keys(D)
  if (isTRUE(sort_rows)) keys <- sort(keys)
  wcrit29_sha256_text(c(paste(dim(D), collapse = "x"), keys))
}

wcrit29_canonical_set_key <- function(D) {
  paste(sort(wcrit29_row_keys(D)), collapse = "|")
}

wcrit29_validate_permutation_matrix <- function(D, m = 6L, n = nrow(D), context = "design") {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  if (!identical(dim(D), c(as.integer(n), as.integer(m)))) {
    wcrit29_stop("%s has dimension %s; expected %dx%d",
                 context, paste(dim(D), collapse = "x"), n, m)
  }
  target <- seq_len(m)
  valid <- apply(D, 1L, function(x) identical(sort(as.integer(x)), target))
  if (!all(valid)) wcrit29_stop("%s contains a non-permutation row", context)
  if (anyDuplicated(wcrit29_row_keys(D))) wcrit29_stop("%s contains duplicate rows", context)
  D
}

wcrit29_is_strict_foldover <- function(D) {
  D <- as.matrix(D)
  if (!length(D) || nrow(D) %% 2L != 0L) return(FALSE)
  ok <- try(wcrit29_validate_permutation_matrix(D, ncol(D), nrow(D)), silent = TRUE)
  if (inherits(ok, "try-error")) return(FALSE)
  keys <- wcrit29_row_keys(D)
  reverse_keys <- wcrit29_row_keys(D[, ncol(D):1L, drop = FALSE])
  all(reverse_keys %in% keys)
}

wcrit29_position_matrix <- function(D) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  P <- matrix(NA_integer_, nrow(D), ncol(D))
  for (i in seq_len(nrow(D))) P[i, D[i, ]] <- seq_len(ncol(D))
  P
}

# Build the complete S_6 universe and its 360 reversal orbits.  Distance
# matrices are precomputed once, so the R search loop only performs indexed
# subsetting rather than recomputing distances at every proposal.
wcrit29_build_universe <- function(m = 6L) {
  m <- as.integer(m)
  if (!identical(m, 6L)) wcrit29_stop("Experiment 29 search is frozen at m=6")
  wcrit29_require("gtools")
  P <- gtools::permutations(m, m, seq_len(m))
  storage.mode(P) <- "integer"
  keys <- wcrit29_row_keys(P)
  if (nrow(P) != factorial(m) || anyDuplicated(keys)) {
    wcrit29_stop("S_6 enumeration invariant failed")
  }
  key_to_index <- stats::setNames(seq_len(nrow(P)), keys)
  reverse_keys <- wcrit29_row_keys(P[, m:1L, drop = FALSE])
  reverse_index <- as.integer(unname(key_to_index[reverse_keys]))
  if (anyNA(reverse_index) || any(reverse_index[reverse_index] != seq_len(nrow(P)))) {
    wcrit29_stop("Reversal-index involution invariant failed")
  }

  # For m=6 every textual row key contains only single-digit labels, so the
  # lexicographic minimum of the two orientations is a canonical orbit key.
  orbit_key <- ifelse(keys <= reverse_keys, keys, reverse_keys)
  orbit_levels <- sort(unique(orbit_key))
  orbit_id <- match(orbit_key, orbit_levels)
  orbit_representative_index <- as.integer(unname(key_to_index[orbit_levels]))
  if (length(orbit_levels) != factorial(m) / 2L || anyNA(orbit_representative_index)) {
    wcrit29_stop("Reversal-orbit enumeration invariant failed")
  }

  # Pairwise-precedence signatures give Kendall distance through
  # d_K(x,y) = {q - z(x)'z(y)}/2.
  positions <- wcrit29_position_matrix(P)
  pairs <- utils::combn(seq_len(m), 2L)
  Z <- matrix(NA_integer_, nrow(P), ncol(pairs))
  for (j in seq_len(ncol(pairs))) {
    Z[, j] <- ifelse(
      positions[, pairs[1L, j]] < positions[, pairs[2L, j]], 1L, -1L
    )
  }
  q <- ncol(Z)
  kendall <- round((q - tcrossprod(Z)) / 2)
  storage.mode(kendall) <- "integer"

  hamming <- matrix(0L, nrow(P), nrow(P))
  for (j in seq_len(m)) {
    hamming <- hamming + outer(P[, j], P[, j], FUN = "!=")
  }
  storage.mode(hamming) <- "integer"

  pos_norm <- rowSums(positions^2)
  position_l2_sq <- outer(pos_norm, pos_norm, "+") - 2 * tcrossprod(positions)
  # Put the matrix first so pmax preserves its dimensions.
  position_l2 <- sqrt(pmax(position_l2_sq, 0))
  diag(position_l2) <- 0

  universe <- list(
    m = m,
    size = nrow(P),
    q = q,
    P = P,
    keys = keys,
    key_to_index = key_to_index,
    reverse_index = reverse_index,
    orbit_key = orbit_key,
    orbit_id = as.integer(orbit_id),
    orbit_levels = orbit_levels,
    orbit_representative_index = orbit_representative_index,
    positions = positions,
    precedence = Z,
    distance = list(
      kendall = kendall,
      hamming = hamming,
      l2 = position_l2
    ),
    distance_definition = c(
      kendall = "pairwise-precedence disagreement count",
      hamming = "positionwise component mismatch count",
      l2 = "Euclidean distance between inverse-permutation component-position vectors"
    )
  )
  universe$universe_sha256 <- wcrit29_sha256_text(c(
    paste0("m=", m),
    paste0("permutations=", paste(keys, collapse = ";")),
    paste0("orbits=", paste(orbit_levels, collapse = ";"))
  ))
  class(universe) <- c("wcrit29_universe", "list")
  universe
}

wcrit29_validate_universe <- function(universe) {
  if (!inherits(universe, "wcrit29_universe") || universe$m != 6L ||
      universe$size != 720L || length(universe$orbit_levels) != 360L) {
    wcrit29_stop("Invalid Experiment-29 universe")
  }
  invisible(TRUE)
}

wcrit29_design_indices <- function(H_index, universe) {
  c(as.integer(H_index), universe$reverse_index[as.integer(H_index)])
}

wcrit29_design_from_indices <- function(H_index, universe) {
  universe$P[wcrit29_design_indices(H_index, universe), , drop = FALSE]
}

wcrit29_pwo_rank_indices <- function(H_index, universe, rank_tol = 1e-9) {
  D_index <- wcrit29_design_indices(H_index, universe)
  X <- cbind(Intercept = 1, universe$precedence[D_index, , drop = FALSE])
  as.integer(qr(X, tol = rank_tol, LAPACK = FALSE)$rank)
}

wcrit29_sample_initial_half <- function(n, seed, universe,
                                        require_pwo_full_rank = TRUE,
                                        rank_tol = 1e-9,
                                        max_attempts = 1000L) {
  wcrit29_validate_universe(universe)
  n <- as.integer(n)
  if (n < 4L || n %% 2L != 0L || n > universe$size) {
    wcrit29_stop("n must be even and satisfy 4 <= n <= 720")
  }
  h <- n %/% 2L
  base_seed <- wcrit29_safe_seed(seed)
  max_attempts <- as.integer(max_attempts)
  if (max_attempts < 1L) wcrit29_stop("max_attempts must be positive")
  selected <- NULL
  actual_seed <- NA_integer_
  pwo_rank <- NA_integer_
  attempt <- 0L
  for (attempt_id in seq_len(max_attempts)) {
    attempt_seed <- if (attempt_id == 1L) {
      base_seed
    } else {
      wcrit29_hash_seed(base_seed, "estimable-initial-attempt", attempt_id)
    }
    candidate <- wcrit29_with_seed(attempt_seed, {
      sample(universe$orbit_representative_index, h, replace = FALSE)
    })
    candidate_rank <- wcrit29_pwo_rank_indices(candidate, universe, rank_tol)
    if (!isTRUE(require_pwo_full_rank) || candidate_rank == universe$q + 1L) {
      selected <- candidate
      actual_seed <- attempt_seed
      pwo_rank <- candidate_rank
      attempt <- attempt_id
      break
    }
  }
  if (is.null(selected)) {
    wcrit29_stop("Could not sample a full-rank initial foldover design in %d attempts", max_attempts)
  }
  D <- wcrit29_design_from_indices(selected, universe)
  if (!wcrit29_is_strict_foldover(D)) wcrit29_stop("Initial half-design is not strict foldover")
  list(
    H_index = as.integer(selected),
    H = universe$P[selected, , drop = FALSE],
    D = D,
    base_seed = base_seed,
    seed = actual_seed,
    sampling_attempt = as.integer(attempt),
    pwo_rank = as.integer(pwo_rank),
    rank_required = isTRUE(require_pwo_full_rank),
    H_sha256 = wcrit29_matrix_sha256(universe$P[selected, , drop = FALSE]),
    D_sha256 = wcrit29_matrix_sha256(D),
    D_set_sha256 = wcrit29_matrix_sha256(D, sort_rows = TRUE)
  )
}

wcrit29_tape_sha256 <- function(tape) {
  required <- c(
    "proposal", "row_id", "temperature", "global_move", "global_index",
    "swap_i", "swap_j", "accept_u"
  )
  if (!all(required %in% names(tape))) wcrit29_stop("Proposal tape has an invalid schema")
  lines <- c(
    sprintf("schema=%s", attr(tape, "schema")),
    sprintf("seed=%d", attr(tape, "seed")),
    sprintf("m=%d", attr(tape, "m")),
    sprintf("h=%d", attr(tape, "h")),
    sprintf("T0=%.17g", attr(tape, "T0")),
    sprintf("alpha=%.17g", attr(tape, "alpha")),
    sprintf("budget=%d", nrow(tape)),
    vapply(seq_len(nrow(tape)), function(i) {
      paste(
        as.integer(tape$proposal[[i]]), as.integer(tape$row_id[[i]]),
        sprintf("%.17g", as.numeric(tape$temperature[[i]])),
        as.integer(tape$global_move[[i]]), as.integer(tape$global_index[[i]]),
        as.integer(tape$swap_i[[i]]), as.integer(tape$swap_j[[i]]),
        sprintf("%.17g", as.numeric(tape$accept_u[[i]])),
        sep = "|"
      )
    }, character(1L))
  )
  wcrit29_sha256_text(lines)
}

# The same tape is replayed by every optimized method.  A global proposal is a
# fixed S_6 permutation; a local proposal swaps the same two positions of the
# method's current target row.  The acceptance uniform is also shared.
wcrit29_make_proposal_tape <- function(h, proposal_budget = 6000L, seed,
                                       universe, T0 = 1, alpha = 0.997) {
  wcrit29_validate_universe(universe)
  h <- as.integer(h)
  proposal_budget <- as.integer(proposal_budget)
  if (h < 2L || h > length(universe$orbit_levels)) wcrit29_stop("Invalid half-design size")
  if (proposal_budget < 1L) wcrit29_stop("proposal_budget must be positive")
  if (!is.finite(T0) || !is.finite(alpha) || T0 <= 0 || alpha <= 0 || alpha > 1) {
    wcrit29_stop("Require T0 > 0 and 0 < alpha <= 1")
  }
  tape <- wcrit29_with_seed(seed, {
    temperature <- T0 * alpha^(seq_len(proposal_budget) - 1L)
    global_u <- stats::runif(proposal_budget)
    swap_i <- sample.int(universe$m, proposal_budget, replace = TRUE)
    swap_j0 <- sample.int(universe$m - 1L, proposal_budget, replace = TRUE)
    swap_j <- swap_j0 + as.integer(swap_j0 >= swap_i)
    data.frame(
      proposal = seq_len(proposal_budget),
      row_id = sample.int(h, proposal_budget, replace = TRUE),
      temperature = as.numeric(temperature),
      global_move = global_u < (temperature / T0),
      global_index = sample.int(universe$size, proposal_budget, replace = TRUE),
      swap_i = as.integer(swap_i),
      swap_j = as.integer(swap_j),
      accept_u = stats::runif(proposal_budget),
      stringsAsFactors = FALSE
    )
  })
  attr(tape, "schema") <- "wcrit29-proposal-tape-v1"
  attr(tape, "seed") <- wcrit29_safe_seed(seed)
  attr(tape, "m") <- universe$m
  attr(tape, "h") <- h
  attr(tape, "T0") <- as.numeric(T0)
  attr(tape, "alpha") <- as.numeric(alpha)
  attr(tape, "sha256") <- wcrit29_tape_sha256(tape)
  class(tape) <- c("wcrit29_tape", "data.frame")
  tape
}

wcrit29_validate_tape <- function(tape, h, proposal_budget, universe) {
  if (!inherits(tape, "wcrit29_tape") || nrow(tape) != as.integer(proposal_budget) ||
      attr(tape, "m") != universe$m || attr(tape, "h") != as.integer(h) ||
      !is.finite(attr(tape, "T0")) || attr(tape, "T0") <= 0 ||
      !is.finite(attr(tape, "alpha")) || attr(tape, "alpha") <= 0 || attr(tape, "alpha") > 1 ||
      any(tape$row_id < 1L | tape$row_id > h) ||
      any(tape$global_index < 1L | tape$global_index > universe$size) ||
      any(tape$swap_i == tape$swap_j) ||
      any(tape$accept_u < 0 | tape$accept_u > 1) ||
      !identical(attr(tape, "sha256"), wcrit29_tape_sha256(tape))) {
    wcrit29_stop("Proposal-tape invariant failed")
  }
  invisible(TRUE)
}

wcrit29_phi_constants <- function(m, n) {
  C1 <- floor(m * (m - 1) / 4)
  C2 <- n * m * (9 * m^3 - 14 * m^2 + 15 * m - 10) / (144 * (n - 1))
  U2 <- (n * m^2 * (m - 1)^2 - 4 * (n - 2) * (m * (m - 1) - 2)) / (8 * (n - 1))
  list(C1 = C1, C2 = C2, U2 = U2)
}

wcrit29_score_indices <- function(H_index, objective = c("phi", "hamming", "l2"),
                                  lambda = 0.5, universe, tie_weight = 1e-6) {
  # "kendall" is accepted only as a report-only geometry diagnostic used by
  # the frozen driver.  It is deliberately not a searchable method.
  diagnostic_kendall <- identical(objective, "kendall")
  if (!diagnostic_kendall) objective <- match.arg(objective)
  H_index <- as.integer(H_index)
  D_index <- wcrit29_design_indices(H_index, universe)
  n <- length(D_index)
  geometry <- if (objective == "phi" || diagnostic_kendall) "kendall" else objective
  dm <- universe$distance[[geometry]][D_index, D_index, drop = FALSE]
  d <- as.numeric(dm[upper.tri(dm)])
  dmin <- min(d)
  nearest <- sum(abs(d - dmin) <= 1e-12)
  distance_m2 <- mean(d^2)
  dmax <- switch(
    geometry,
    kendall = universe$q,
    hamming = universe$m,
    l2 = sqrt(universe$m * (universe$m^2 - 1) / 3)
  )
  min_norm <- dmin / dmax
  mean_sq_norm <- distance_m2 / dmax^2
  D <- universe$P[D_index, , drop = FALSE]
  canonical_key <- wcrit29_canonical_set_key(D)

  if (objective == "phi") {
    lambda <- as.numeric(lambda)
    if (length(lambda) != 1L || !is.finite(lambda) || !isTRUE(all.equal(lambda, 0.5))) {
      wcrit29_stop("The paper search is fixed at lambda=0.5")
    }
    constants <- wcrit29_phi_constants(universe$m, n)
    A <- (dmin - 1) / (constants$C1 - 1)
    B <- (constants$U2 - distance_m2) / (constants$U2 - constants$C2)
    phi <- lambda * A + (1 - lambda) * B
    return(list(
      objective = objective,
      lambda = lambda,
      accept_value = as.numeric(phi),
      primary = as.numeric(phi),
      phi = as.numeric(phi), A = as.numeric(A), B = as.numeric(B),
      k_min = as.numeric(dmin), k_m2 = as.numeric(distance_m2),
      nearest_pair_multiplicity = as.integer(nearest),
      distance_m2 = as.numeric(distance_m2),
      min_norm = as.numeric(min_norm), mean_sq_norm = as.numeric(mean_sq_norm),
      canonical_key = canonical_key
    ))
  }

  list(
    objective = objective,
    lambda = NA_real_,
    accept_value = as.numeric(min_norm + tie_weight * mean_sq_norm),
    primary = as.numeric(min_norm),
    phi = NA_real_, A = NA_real_, B = NA_real_,
    k_min = if (diagnostic_kendall) as.numeric(dmin) else NA_real_,
    k_m2 = if (diagnostic_kendall) as.numeric(distance_m2) else NA_real_,
    nearest_pair_multiplicity = as.integer(nearest),
    distance_m2 = as.numeric(distance_m2),
    min_norm = as.numeric(min_norm), mean_sq_norm = as.numeric(mean_sq_norm),
    canonical_key = canonical_key
  )
}

wcrit29_num_prefer <- function(a, b, direction = c("high", "low"), tolerance = 1e-12) {
  direction <- match.arg(direction)
  if (!is.finite(a) || !is.finite(b)) return(0L)
  if (abs(a - b) <= tolerance) return(0L)
  if (direction == "high") if (a > b) 1L else -1L else if (a < b) 1L else -1L
}

# Deterministic best-state comparator for the three paper search objectives.
wcrit29_is_better <- function(a, b, objective, lambda = 0.5, tolerance = 1e-12) {
  if (is.null(b)) return(TRUE)
  objective <- match.arg(objective, c("phi", "hamming", "l2"))
  if (objective == "phi" && !isTRUE(all.equal(as.numeric(lambda), 0.5))) {
    wcrit29_stop("The paper search is fixed at lambda=0.5")
  }
  cmp <- wcrit29_num_prefer(a$primary, b$primary, "high", tolerance)
  if (cmp != 0L) return(cmp > 0L)

  if (objective == "phi") {
    cmp <- wcrit29_num_prefer(a$k_min, b$k_min, "high", tolerance)
    if (cmp != 0L) return(cmp > 0L)
    cmp <- wcrit29_num_prefer(a$k_m2, b$k_m2, "low", tolerance)
    if (cmp != 0L) return(cmp > 0L)
    cmp <- wcrit29_num_prefer(
      a$nearest_pair_multiplicity, b$nearest_pair_multiplicity, "low", tolerance = 0
    )
    if (cmp != 0L) return(cmp > 0L)
  } else {
    cmp <- wcrit29_num_prefer(
      a$nearest_pair_multiplicity, b$nearest_pair_multiplicity, "low", tolerance = 0
    )
    if (cmp != 0L) return(cmp > 0L)
    cmp <- wcrit29_num_prefer(a$distance_m2, b$distance_m2, "high", tolerance)
    if (cmp != 0L) return(cmp > 0L)
  }
  isTRUE(a$canonical_key < b$canonical_key)
}

wcrit29_method_id <- function(objective, lambda = NA_real_) {
  if (objective == "phi") {
    if (!isTRUE(all.equal(as.numeric(lambda), 0.5))) {
      wcrit29_stop("The paper search is fixed at lambda=0.5")
    }
    return("FSA_lambda050")
  }
  switch(objective,
         hamming = "Hamming",
         l2 = "L2",
         wcrit29_stop("Unknown objective '%s'", objective))
}

wcrit29_method_label <- function(method) {
  if (identical(method, "FSA_lambda050")) return("FSA-KD (lambda=0.5)")
  switch(method,
         Hamming = "Hamming",
         L2 = "L2",
         RandomFoldover = "Random foldover",
         method)
}

wcrit29_search_one <- function(initial, tape, objective = c("phi", "hamming", "l2"),
                               lambda = 0.5, universe, tie_weight = 1e-6,
                               tolerance = 1e-12) {
  objective <- match.arg(objective)
  H_index <- as.integer(initial$H_index)
  h <- length(H_index)
  wcrit29_validate_tape(tape, h, nrow(tape), universe)
  if (anyDuplicated(universe$orbit_id[H_index])) {
    wcrit29_stop("Initial half-design contains duplicate reversal orbits")
  }
  if (wcrit29_pwo_rank_indices(H_index, universe) != universe$q + 1L) {
    wcrit29_stop("Initial half-design is not PWO estimable")
  }

  current <- wcrit29_score_indices(H_index, objective, lambda, universe, tie_weight)
  best <- current
  best_H_index <- H_index
  valid_evals <- 0L
  invalid <- 0L
  rank_invalid <- 0L
  accepted <- 0L
  improving_accepts <- 0L

  for (p in seq_len(nrow(tape))) {
    row_id <- tape$row_id[[p]]
    if (isTRUE(tape$global_move[[p]])) {
      candidate_index <- tape$global_index[[p]]
    } else {
      candidate <- universe$P[H_index[[row_id]], ]
      i <- tape$swap_i[[p]]
      j <- tape$swap_j[[p]]
      candidate[c(i, j)] <- candidate[c(j, i)]
      candidate_index <- as.integer(unname(
        universe$key_to_index[[paste(candidate, collapse = ",")]]
      ))
    }

    candidate_orbit <- universe$orbit_id[[candidate_index]]
    other_orbits <- universe$orbit_id[H_index[-row_id]]
    current_orbit <- universe$orbit_id[[H_index[[row_id]]]]
    # No-op, reversal of the current representative, or collision with any
    # other occupied reversal orbit: the proposal is invalid but still cools.
    if (candidate_orbit == current_orbit || candidate_orbit %in% other_orbits) {
      invalid <- invalid + 1L
      next
    }

    H_new <- H_index
    H_new[[row_id]] <- candidate_index
    if (wcrit29_pwo_rank_indices(H_new, universe) != universe$q + 1L) {
      invalid <- invalid + 1L
      rank_invalid <- rank_invalid + 1L
      next
    }
    candidate_score <- wcrit29_score_indices(
      H_new, objective, lambda, universe, tie_weight
    )
    valid_evals <- valid_evals + 1L
    delta <- candidate_score$accept_value - current$accept_value
    accept <- is.finite(delta) && (
      delta >= 0 || tape$accept_u[[p]] < exp(delta / max(tape$temperature[[p]], 1e-14))
    )
    if (accept) {
      accepted <- accepted + 1L
      if (candidate_score$accept_value > current$accept_value + tolerance) {
        improving_accepts <- improving_accepts + 1L
      }
      H_index <- H_new
      current <- candidate_score
      if (wcrit29_is_better(current, best, objective, lambda, tolerance)) {
        best <- current
        best_H_index <- H_index
      }
    }
  }

  D <- wcrit29_design_from_indices(best_H_index, universe)
  if (!wcrit29_is_strict_foldover(D)) wcrit29_stop("Search returned a non-foldover design")
  method <- wcrit29_method_id(objective, lambda)
  list(
    method = method,
    method_label = wcrit29_method_label(method),
    objective = objective,
    lambda = if (objective == "phi") as.numeric(lambda) else NA_real_,
    H_index = as.integer(best_H_index),
    H = universe$P[best_H_index, , drop = FALSE],
    D = D,
    score = best,
    proposals = as.integer(nrow(tape)),
    initial_objective_evals = 1L,
    candidate_objective_evals = as.integer(valid_evals),
    n_obj_eval = as.integer(1L + valid_evals),
    invalid_proposals = as.integer(invalid),
    rank_invalid_proposals = as.integer(rank_invalid),
    accepted_proposals = as.integer(accepted),
    improving_accepts = as.integer(improving_accepts),
    proposal_budget_exhausted = TRUE,
    budget_unit = "logical_proposal_attempts",
    initial_H_sha256 = initial$H_sha256,
    initial_D_sha256 = initial$D_sha256,
    proposal_tape_sha256 = attr(tape, "sha256"),
    design_sha256 = wcrit29_matrix_sha256(D),
    design_set_sha256 = wcrit29_matrix_sha256(D, sort_rows = TRUE),
    strict_foldover = TRUE,
    pwo_rank = wcrit29_pwo_rank_indices(best_H_index, universe)
  )
}

wcrit29_random_baseline <- function(initial, tape) {
  D <- initial$D
  list(
    method = "RandomFoldover",
    method_label = wcrit29_method_label("RandomFoldover"),
    objective = "none",
    lambda = NA_real_,
    H_index = initial$H_index,
    H = initial$H,
    D = D,
    score = NULL,
    proposals = 0L,
    initial_objective_evals = 0L,
    candidate_objective_evals = 0L,
    n_obj_eval = 0L,
    invalid_proposals = 0L,
    rank_invalid_proposals = 0L,
    accepted_proposals = 0L,
    improving_accepts = 0L,
    proposal_budget_exhausted = NA,
    budget_unit = "one_unoptimized_strict_foldover_design",
    initial_H_sha256 = initial$H_sha256,
    initial_D_sha256 = initial$D_sha256,
    proposal_tape_sha256 = attr(tape, "sha256"),
    design_sha256 = wcrit29_matrix_sha256(D),
    design_set_sha256 = wcrit29_matrix_sha256(D, sort_rows = TRUE),
    strict_foldover = wcrit29_is_strict_foldover(D),
    pwo_rank = initial$pwo_rank
  )
}

# Generate the fully paired publication design bank.
wcrit29_build_design_bank <- function(n = 48L, rep_id = 1L, master_seed = 20260901L,
                                      lambdas = 0.5,
                                      proposal_budget = 6000L,
                                      universe = NULL,
                                      T0 = 1, alpha = 0.997,
                                      tie_weight = 1e-6) {
  if (is.null(universe)) universe <- wcrit29_build_universe(6L)
  wcrit29_validate_universe(universe)
  n <- as.integer(n)
  rep_id <- as.integer(rep_id)
  lambdas <- as.numeric(lambdas)
  if (!identical(lambdas, 0.5)) {
    wcrit29_stop("The paper design bank is fixed at lambda=0.5")
  }
  init_seed <- wcrit29_hash_seed(master_seed, "29", "shared-initial-half", n, rep_id)
  tape_seed <- wcrit29_hash_seed(master_seed, "29", "shared-proposal-tape", n, rep_id)
  initial <- wcrit29_sample_initial_half(n, init_seed, universe)
  tape <- wcrit29_make_proposal_tape(
    n %/% 2L, proposal_budget, tape_seed, universe, T0 = T0, alpha = alpha
  )

  started <- proc.time()[["elapsed"]]
  entries <- list()
  id <- wcrit29_method_id("phi", 0.5)
  entries[[id]] <- wcrit29_search_one(
    initial, tape, "phi", 0.5, universe, tie_weight
  )
  for (objective in c("hamming", "l2")) {
    id <- wcrit29_method_id(objective)
    entries[[id]] <- wcrit29_search_one(
      initial, tape, objective, 0.5, universe, tie_weight
    )
  }
  entries[["RandomFoldover"]] <- wcrit29_random_baseline(initial, tape)

  elapsed <- proc.time()[["elapsed"]] - started
  if (length(unique(vapply(entries, `[[`, character(1L), "initial_H_sha256"))) != 1L ||
      length(unique(vapply(entries, `[[`, character(1L), "proposal_tape_sha256"))) != 1L) {
    wcrit29_stop("Shared-initial/shared-tape invariant failed")
  }
  if (any(vapply(entries[names(entries) != "RandomFoldover"], `[[`, integer(1L), "proposals") != proposal_budget)) {
    wcrit29_stop("Proposal-budget invariant failed")
  }

  list(
    schema = "wcrit29-design-bank-v1",
    m = universe$m,
    n = n,
    rep = rep_id,
    master_seed = as.integer(master_seed),
    init_seed = init_seed,
    init_actual_seed = initial$seed,
    init_sampling_attempt = initial$sampling_attempt,
    tape_seed = tape_seed,
    lambdas = lambdas,
    proposal_budget = as.integer(proposal_budget),
    T0 = T0,
    alpha = alpha,
    tie_weight = tie_weight,
    universe_sha256 = universe$universe_sha256,
    initial = initial,
    tape = tape,
    entries = entries,
    elapsed_sec = as.numeric(elapsed)
  )
}

wcrit29_bank_audit <- function(bank) {
  entries <- bank$entries
  methods <- names(entries)
  optimized <- methods != "RandomFoldover"
  data.frame(
    method = methods,
    method_label = vapply(entries, `[[`, character(1L), "method_label"),
    strict_foldover = vapply(entries, `[[`, logical(1L), "strict_foldover"),
    initial_H_sha256 = vapply(entries, `[[`, character(1L), "initial_H_sha256"),
    proposal_tape_sha256 = vapply(entries, `[[`, character(1L), "proposal_tape_sha256"),
    T0 = rep(as.numeric(attr(bank$tape, "T0")), length(entries)),
    alpha = rep(as.numeric(attr(bank$tape, "alpha")), length(entries)),
    proposals = vapply(entries, `[[`, integer(1L), "proposals"),
    candidate_objective_evals = vapply(entries, `[[`, integer(1L), "candidate_objective_evals"),
    invalid_proposals = vapply(entries, `[[`, integer(1L), "invalid_proposals"),
    rank_invalid_proposals = vapply(entries, `[[`, integer(1L), "rank_invalid_proposals"),
    pwo_rank = vapply(entries, `[[`, integer(1L), "pwo_rank"),
    proposal_partition_pass = vapply(entries, function(x) {
      if (x$method == "RandomFoldover") return(x$proposals == 0L)
      x$candidate_objective_evals + x$invalid_proposals == x$proposals
    }, logical(1L)),
    budget_pass = vapply(entries, function(x) {
      if (x$method == "RandomFoldover") return(TRUE)
      x$proposals == bank$proposal_budget && isTRUE(x$proposal_budget_exhausted)
    }, logical(1L)),
    optimized = optimized,
    design_sha256 = vapply(entries, `[[`, character(1L), "design_sha256"),
    design_set_sha256 = vapply(entries, `[[`, character(1L), "design_set_sha256"),
    stringsAsFactors = FALSE
  )
}

# Legacy function name retained because the frozen driver calls this integrity
# hook.  The public paper version tests only the lambda=0.5 deterministic
# canonical tie break; no endpoint-weight search is implemented.
wcrit29_endpoint_comparator_selftest <- function() {
  base <- list(
    primary = 0.5, k_m2 = 10, k_min = 3,
    nearest_pair_multiplicity = 7L, distance_m2 = 10,
    canonical_key = "b"
  )
  canonical <- base
  canonical$canonical_key <- "a"
  c(phi05_final_canonical_tie =
      wcrit29_is_better(canonical, base, "phi", 0.5))
}

# Full source-level smoke test.  The default deliberately replays two complete
# 6000-proposal banks, rather than weakening the budget invariant for speed.
wcrit29_selftest <- function(n = 48L, proposal_budget = 6000L,
                             master_seed = 20260901L,
                             lambdas = 0.5, verbose = TRUE) {
  total_start <- proc.time()[["elapsed"]]
  universe_start <- proc.time()[["elapsed"]]
  universe <- wcrit29_build_universe(6L)
  universe_elapsed <- proc.time()[["elapsed"]] - universe_start

  first_start <- proc.time()[["elapsed"]]
  bank1 <- wcrit29_build_design_bank(
    n = n, rep_id = 1L, master_seed = master_seed,
    lambdas = lambdas, proposal_budget = proposal_budget,
    universe = universe
  )
  first_elapsed <- proc.time()[["elapsed"]] - first_start
  second_start <- proc.time()[["elapsed"]]
  bank2 <- wcrit29_build_design_bank(
    n = n, rep_id = 1L, master_seed = master_seed,
    lambdas = lambdas, proposal_budget = proposal_budget,
    universe = universe
  )
  second_elapsed <- proc.time()[["elapsed"]] - second_start

  audit1 <- wcrit29_bank_audit(bank1)
  audit2 <- wcrit29_bank_audit(bank2)
  checks <- data.frame(
    check = c(
      "universe_has_720_permutations",
      "universe_has_360_reversal_orbits",
      "l2_is_inverse_component_position_not_native_row",
      "initial_and_final_designs_are_pwo_estimable",
      "all_designs_strict_foldover",
      "one_shared_initial_half_hash",
      "one_shared_proposal_tape_hash",
      "optimized_methods_use_exact_proposal_budget",
      "invalid_plus_valid_partitions_every_proposal",
      "same_seed_reproduces_initial_hash",
      "same_seed_reproduces_tape_hash",
      "same_seed_reproduces_every_design",
      "random_baseline_is_initial_design"
    ),
    pass = c(
      universe$size == 720L,
      length(universe$orbit_levels) == 360L,
      abs(universe$distance$l2[2L, 4L] - sqrt(8)) <= 1e-12 &&
        abs(sqrt(sum((universe$P[2L, ] - universe$P[4L, ])^2)) - sqrt(2)) <= 1e-12,
      bank1$initial$pwo_rank == universe$q + 1L && all(audit1$pwo_rank == universe$q + 1L),
      all(audit1$strict_foldover),
      length(unique(audit1$initial_H_sha256)) == 1L,
      length(unique(audit1$proposal_tape_sha256)) == 1L,
      all(audit1$budget_pass),
      all(audit1$proposal_partition_pass),
      identical(bank1$initial$H_sha256, bank2$initial$H_sha256),
      identical(attr(bank1$tape, "sha256"), attr(bank2$tape, "sha256")),
      identical(audit1$design_sha256, audit2$design_sha256),
      identical(
        bank1$entries$RandomFoldover$design_sha256,
        bank1$initial$D_sha256
      )
    ),
    stringsAsFactors = FALSE
  )
  if (any(!checks$pass)) {
    wcrit29_stop("Self-test failed: %s", paste(checks$check[!checks$pass], collapse = ", "))
  }
  timing <- data.frame(
    stage = c("universe", "first_full_bank", "deterministic_replay", "total"),
    elapsed_sec = c(
      universe_elapsed, first_elapsed, second_elapsed,
      proc.time()[["elapsed"]] - total_start
    ),
    stringsAsFactors = FALSE
  )
  if (isTRUE(verbose)) {
    message(sprintf(
      "wcrit29 self-test passed: n=%d, methods=%d, proposals/optimized-method=%d, total=%.2fs",
      n, length(bank1$entries), proposal_budget, timing$elapsed_sec[timing$stage == "total"]
    ))
  }
  invisible(list(
    pass = TRUE,
    checks = checks,
    timing = timing,
    audit = audit1,
    bank = bank1,
    replay_bank = bank2,
    universe = universe
  ))
}
