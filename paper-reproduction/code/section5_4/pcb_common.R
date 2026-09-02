# Core functions for the paper-only Experiment 24 reproduction.
#
# This file intentionally contains only the response and surrogate machinery
# used by Section 5.4: one strong mixed-physics response and one additive
# Kendall--directed-adjacency GP.  Initial-design construction is not repeated;
# the four frozen paper designs are supplied in data/frozen/section5_4.

pcb_stop <- function(...) stop(sprintf(...), call. = FALSE)

pcb_methods <- c(
  "fsa_lambda05", "unrestricted_hamming",
  "unrestricted_position_l2", "srs"
)

pcb_method_labels <- c(
  fsa_lambda05 = "FSA-KD",
  unrestricted_hamming = "Hamming",
  unrestricted_position_l2 = "L2",
  srs = "SRS"
)

pcb_validate_routes <- function(D, m = 10L, unique_rows = FALSE) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  if (nrow(D) < 1L || ncol(D) != m) pcb_stop("Invalid route dimensions")
  ok <- apply(D, 1L, function(x) identical(sort(as.integer(x)), seq_len(m)))
  if (!all(ok)) pcb_stop("A route is not a permutation of 1:m")
  if (unique_rows && anyDuplicated(wcrit_row_keys(D))) pcb_stop("Duplicate routes")
  D
}

pcb_axis_time <- function(distance, vmax, acceleration) {
  distance <- abs(as.numeric(distance))
  threshold <- vmax^2 / acceleration
  ifelse(distance <= threshold,
         2 * sqrt(distance / acceleration),
         distance / vmax + vmax / acceleration)
}

pcb_build_cost <- function(coordinates) {
  required <- c("role", "local_id", "tsplib_node_id", "x_tsplib", "y_tsplib")
  if (!identical(names(coordinates), required) || nrow(coordinates) != 11L) {
    pcb_stop("Unexpected PCB coordinate snapshot")
  }
  x <- coordinates$x_tsplib * 0.001
  y <- coordinates$y_tsplib * 0.001
  C <- matrix(0, 11L, 11L)
  for (i in seq_len(11L)) for (j in seq_len(11L)) {
    if (i == j) next
    dx <- x[[j]] - x[[i]]
    dy <- y[[j]] - y[[i]]
    tx <- if (dx >= 0) pcb_axis_time(dx, 1.00, 4.00) else pcb_axis_time(dx, 0.85, 3.40)
    ty <- if (dy >= 0) pcb_axis_time(dy, 0.80, 3.20) else pcb_axis_time(dy, 0.70, 2.80)
    C[i, j] <- max(tx, ty) + 0.05 + 0.04 * (dx < 0) + 0.03 * (dy < 0)
  }
  C
}

pcb_profile <- function(profile_table, profile_id) {
  z <- profile_table[profile_table$profile_id == profile_id, , drop = FALSE]
  z <- z[order(z$hole), , drop = FALSE]
  if (nrow(z) != 10L || !identical(as.integer(z$hole), 1:10)) {
    pcb_stop("Invalid physics profile %s", profile_id)
  }
  list(
    uW = z$uW, sW = z$sW, uH = z$uH, sH = z$sH,
    tau_H = unique(z$tau_H), profile_id = as.integer(profile_id)
  )
}

pcb_eval_components <- function(D, C, profile, dwell = 0.20) {
  D <- pcb_validate_routes(D)
  move <- wear <- thermal <- numeric(nrow(D))
  for (i in seq_len(nrow(D))) {
    route <- D[i, ]
    idx <- route + 1L
    move[[i]] <- C[1L, idx[[1L]]] +
      sum(C[cbind(idx[-length(idx)], idx[-1L])]) +
      C[idx[[length(idx)]], 1L] + ncol(D) * dwell
    W <- H <- 0
    previous <- 0L
    for (r in seq_len(ncol(D))) {
      hole <- route[[r]]
      wear[[i]] <- wear[[i]] + profile$sW[[hole]] * W
      W <- W + profile$uW[[hole]]
      gap <- C[previous + 1L, hole + 1L] + dwell
      Hbar <- exp(-gap / profile$tau_H) * H
      thermal[[i]] <- thermal[[i]] + profile$sH[[hole]] * Hbar
      H <- Hbar + profile$uH[[hole]]
      previous <- hole
    }
  }
  data.frame(move = move, wear = wear, thermal = thermal)
}

pcb_eval_response <- function(D, C, profile, scenario) {
  raw <- pcb_eval_components(D, C, profile)
  zW <- (raw$wear - scenario$wear_mean) / scenario$wear_sd
  zH <- (raw$thermal - scenario$thermal_mean) / scenario$thermal_sd
  mixed <- scenario$omega * zW + (1 - scenario$omega) * zH
  mixed <- (mixed - scenario$mixed_z_mean) / scenario$mixed_z_sd
  raw$move + scenario$gamma * scenario$motion_sd * mixed
}

pcb_make_c_grid <- function(cmax = 1024, log2_step = 0.5) {
  small <- c(0, 2^-12, 2^-10, 2^-8)
  exponents <- seq(-6, floor(log2(cmax) / log2_step) * log2_step, by = log2_step)
  sort(unique(c(small, 2^exponents, cmax)))
}

pcb_kernel_config <- function() {
  list(
    m = 10L, c_min = 0, c_max = 1024,
    c_grid = pcb_make_c_grid(),
    pure_local_maximum_starts = 10L,
    pure_continuous_multistarts = 5L,
    additive_main_rho_starts = c(0.2, 0.5, 0.8),
    additive_rho_epsilon = 1e-6,
    additive_maxit = 100L,
    optimizer_tie_tolerance = 1e-8,
    gp_ridge = 1e-8, gp_jitter = 0,
    variance_floor = 1e-12,
    ei_tie_tolerance = 1e-12,
    extended_start_t_offsets = c(-0.25, 0, 0.25)
  )
}

pcb_arc_keys <- function(D) {
  D <- as.matrix(D)
  q <- ncol(D) + 1L
  out <- matrix(NA_integer_, nrow(D), q)
  for (i in seq_len(nrow(D))) {
    from <- c(0L, D[i, ])
    to <- c(D[i, ], 0L)
    out[i, ] <- from * q + to + 1L
  }
  out
}

pcb_adjacency_dmat <- function(D1, D2 = D1) {
  e1 <- pcb_arc_keys(D1)
  e2 <- pcb_arc_keys(D2)
  q <- ncol(e1)
  out <- matrix(NA_real_, nrow(e1), nrow(e2))
  for (j in seq_len(nrow(e2))) {
    out[, j] <- (q - rowSums(matrix(e1 %in% e2[j, ], nrow(e1), q))) / q
  }
  out
}

pcb_distance_matrices <- function(D1, D2 = D1) {
  m <- ncol(as.matrix(D1))
  list(
    kendall = wcrit_kendall_dmat(as.matrix(D1), as.matrix(D2)) / choose(m, 2),
    adjacency = pcb_adjacency_dmat(D1, D2)
  )
}

pcb_quotient_kernel <- function(S, c) {
  S <- as.matrix(S)
  S[S < 0] <- 0
  S[S > 1] <- 1
  if (c == 0) return(S)
  if (is.infinite(c)) return(matrix(as.numeric(S >= 1 - 1e-15), nrow(S), ncol(S)))
  if (c < 1e-5) {
    poly <- function(x) x * (1 + x * (1 / 2 + x * (1 / 6 + x / 24)))
    return(poly(c * S) / poly(c))
  }
  if (c > 20) return(exp(c * (S - 1)) * (-expm1(-c * S)) / (-expm1(-c)))
  expm1(c * S) / expm1(c)
}

pcb_component_kernel <- function(d, geometry, c) {
  pcb_quotient_kernel(1 - d[[geometry]], c)
}

pcb_combined_kernel <- function(d, cK, cA, rho) {
  if (rho == 1) return(pcb_component_kernel(d, "kendall", cK))
  if (rho == 0) return(pcb_component_kernel(d, "adjacency", cA))
  rho * pcb_component_kernel(d, "kendall", cK) +
    (1 - rho) * pcb_component_kernel(d, "adjacency", cA)
}

pcb_reml_matrix <- function(Rlatent, yz, ridge, variance_floor, keep = FALSE) {
  R <- as.matrix(Rlatent)
  diag(R) <- diag(R) + ridge
  U <- try(chol(R), silent = TRUE)
  if (inherits(U, "try-error")) return(list(reml = -Inf))
  solve_R <- function(B) backsolve(U, forwardsolve(t(U), B))
  one <- rep(1, length(yz))
  Ri_one <- solve_R(one)
  denom <- sum(Ri_one)
  if (!is.finite(denom) || denom <= 0) return(list(reml = -Inf))
  beta <- sum(solve_R(yz)) / denom
  residual <- yz - beta
  alpha <- solve_R(residual)
  df <- length(yz) - 1L
  sse <- sum(residual * alpha)
  if (!is.finite(sse) || sse < -1e-8) return(list(reml = -Inf))
  sigma2 <- max(max(0, sse) / df, variance_floor)
  reml <- -0.5 * (df * (log(2 * pi) + log(sigma2)) + sse / sigma2 +
                    2 * sum(log(diag(U))) + log(denom))
  out <- list(reml = reml)
  if (keep) out <- c(out, list(
    beta = beta, sigma2 = sigma2, chol_R = U, alpha = alpha,
    trend_information = denom, latent_kernel = Rlatent, ridge = ridge
  ))
  out
}

pcb_reml_at <- function(d, yz, cfg, cK, cA, rho, keep = FALSE) {
  valid_c <- function(x) length(x) == 1L && !is.na(x) && x >= cfg$c_min &&
    ((is.finite(x) && x <= cfg$c_max) || (is.infinite(x) && x > 0))
  if (!is.finite(rho) || rho < 0 || rho > 1 ||
      (rho > 0 && !valid_c(cK)) || (rho < 1 && !valid_c(cA))) return(list(reml = -Inf))
  R <- pcb_combined_kernel(d, cK, cA, rho)
  out <- pcb_reml_matrix(R, yz, cfg$gp_ridge + cfg$gp_jitter, cfg$variance_floor, keep)
  if (keep) {
    out$cK <- if (rho > 0) cK else NA_real_
    out$cA <- if (rho < 1) cA else NA_real_
    out$rho_kendall <- rho
  }
  out
}

pcb_pure_optimum <- function(d, yz, geometry, cfg) {
  rho <- if (geometry == "kendall") 1 else 0
  score <- function(tval) {
    cval <- expm1(tval)
    pcb_reml_at(d, yz, cfg,
                if (rho == 1) cval else NA_real_,
                if (rho == 0) cval else NA_real_, rho)$reml
  }
  t_grid <- log1p(cfg$c_grid)
  grid_reml <- vapply(t_grid, score, numeric(1L))
  local <- which(vapply(seq_along(grid_reml), function(i) {
    left <- if (i == 1L) -Inf else grid_reml[[i - 1L]]
    right <- if (i == length(grid_reml)) -Inf else grid_reml[[i + 1L]]
    is.finite(grid_reml[[i]]) && grid_reml[[i]] >= left && grid_reml[[i]] >= right
  }, logical(1L)))
  local <- head(local[order(grid_reml[local], decreasing = TRUE)],
                cfg$pure_local_maximum_starts)
  global <- head(order(grid_reml, decreasing = TRUE, na.last = NA),
                 cfg$pure_continuous_multistarts)
  starts <- unique(c(local, global))
  continuous <- list()
  for (k in seq_along(starts)) {
    i <- starts[[k]]
    lower <- if (i == 1L) t_grid[[1L]] else (t_grid[[i - 1L]] + t_grid[[i]]) / 2
    upper <- if (i == length(t_grid)) tail(t_grid, 1L) else (t_grid[[i]] + t_grid[[i + 1L]]) / 2
    opt <- try(stats::optimize(score, c(lower, upper), maximum = TRUE,
                               tol = .Machine$double.eps^0.25), silent = TRUE)
    if (!inherits(opt, "try-error") && is.finite(opt$objective)) {
      continuous[[length(continuous) + 1L]] <- data.frame(
        t = opt$maximum, reml = opt$objective, source = sprintf("continuous_%d", k),
        complexity = 1L, priority = 3L
      )
    }
  }
  candidates <- rbind(
    data.frame(
      t = c(0, log1p(cfg$c_max), Inf),
      reml = c(score(0), score(log1p(cfg$c_max)), score(Inf)),
      source = c("exact_linear_limit", "finite_upper", "exact_independence"),
      complexity = c(0L, 1L, 0L), priority = c(0L, 2L, 1L)
    ),
    data.frame(t = t_grid, reml = grid_reml, source = "grid",
               complexity = 1L, priority = 4L),
    if (length(continuous)) do.call(rbind, continuous) else NULL
  )
  candidates$c <- ifelse(is.infinite(candidates$t), Inf, expm1(candidates$t))
  eligible <- which(candidates$reml >= max(candidates$reml) - cfg$optimizer_tie_tolerance)
  ord <- order(candidates$complexity[eligible], candidates$priority[eligible],
               candidates$c[eligible], candidates$source[eligible])
  best <- candidates[eligible[ord[[1L]]], ]
  model <- pcb_reml_at(d, yz, cfg,
                       if (rho == 1) best$c else NA_real_,
                       if (rho == 0) best$c else NA_real_, rho, keep = TRUE)
  model$active_c <- best$c
  model
}

pcb_three_starts <- function(single_c, cfg) {
  tmax <- log1p(cfg$c_max)
  single_t <- if (is.finite(single_c)) log1p(single_c) else tmax
  target <- pmin(tmax, pmax(0, single_t + cfg$extended_start_t_offsets * tmax))
  pool <- unique(c(target, 0, tmax / 2, tmax, log1p(cfg$c_grid)))
  chosen <- numeric(0)
  for (x in target) {
    available <- pool[!vapply(pool, function(v) any(abs(v - chosen) <= 1e-12), logical(1L))]
    if (length(available)) chosen <- c(chosen, available[[which.min(abs(available - x))]])
  }
  if (length(chosen) < 3L) {
    available <- pool[!vapply(pool, function(v) any(abs(v - chosen) <= 1e-12), logical(1L))]
    chosen <- c(chosen, head(available, 3L - length(chosen)))
  }
  expm1(chosen[1:3])
}

pcb_additive_optimum <- function(d, yz, pureK, pureA, cfg) {
  candidates <- list(
    c(pcb_reml_at(d, yz, cfg, pureK$active_c, NA_real_, 1, TRUE),
      list(source = "exact_rho_1", complexity = 0L)),
    c(pcb_reml_at(d, yz, cfg, NA_real_, pureA$active_c, 0, TRUE),
      list(source = "exact_rho_0", complexity = 0L))
  )
  starts <- expand.grid(
    cK = pcb_three_starts(pureK$active_c, cfg),
    cA = pcb_three_starts(pureA$active_c, cfg),
    rho = cfg$additive_main_rho_starts
  )
  objective <- function(par) {
    value <- pcb_reml_at(d, yz, cfg, expm1(par[[1L]]), expm1(par[[2L]]), par[[3L]])$reml
    if (is.finite(value)) -value else .Machine$double.xmax^0.25
  }
  for (i in seq_len(nrow(starts))) {
    opt <- try(stats::optim(
      c(log1p(starts$cK[[i]]), log1p(starts$cA[[i]]), starts$rho[[i]]),
      objective, method = "L-BFGS-B",
      lower = c(0, 0, cfg$additive_rho_epsilon),
      upper = c(log1p(cfg$c_max), log1p(cfg$c_max), 1 - cfg$additive_rho_epsilon),
      control = list(maxit = cfg$additive_maxit, factr = 1e7)
    ), silent = TRUE)
    if (!inherits(opt, "try-error") && is.finite(opt$value) && opt$convergence == 0L) {
      fit <- pcb_reml_at(d, yz, cfg, expm1(opt$par[[1L]]), expm1(opt$par[[2L]]),
                         opt$par[[3L]], TRUE)
      fit$source <- sprintf("interior_%d", i)
      fit$complexity <- 1L
      candidates[[length(candidates) + 1L]] <- fit
    }
  }
  tab <- data.frame(
    id = seq_along(candidates),
    reml = vapply(candidates, `[[`, numeric(1L), "reml"),
    complexity = vapply(candidates, `[[`, integer(1L), "complexity"),
    cK = vapply(candidates, `[[`, numeric(1L), "cK"),
    cA = vapply(candidates, `[[`, numeric(1L), "cA"),
    rho = vapply(candidates, `[[`, numeric(1L), "rho_kendall")
  )
  eligible <- which(tab$reml >= max(tab$reml) - cfg$optimizer_tie_tolerance)
  ord <- order(tab$complexity[eligible], ifelse(is.na(tab$cK[eligible]), Inf, tab$cK[eligible]),
               ifelse(is.na(tab$cA[eligible]), Inf, tab$cA[eligible]), tab$rho[eligible])
  candidates[[eligible[ord[[1L]]]]]
}

pcb_fit_gp <- function(D, y, cfg = pcb_kernel_config()) {
  D <- pcb_validate_routes(D, cfg$m, TRUE)
  y_center <- mean(y)
  y_scale <- stats::sd(y)
  if (!is.finite(y_scale) || y_scale < 1e-12) y_scale <- 1
  yz <- (y - y_center) / y_scale
  d <- pcb_distance_matrices(D)
  pureK <- pcb_pure_optimum(d, yz, "kendall", cfg)
  pureA <- pcb_pure_optimum(d, yz, "adjacency", cfg)
  model <- pcb_additive_optimum(d, yz, pureK, pureA, cfg)
  model$train_D <- D
  model$y_center <- y_center
  model$y_scale <- y_scale
  model
}

pcb_predict_gp <- function(model, Dnew) {
  Dnew <- pcb_validate_routes(Dnew)
  d <- pcb_distance_matrices(Dnew, model$train_D)
  Rstar <- pcb_combined_kernel(d, model$cK, model$cA, model$rho_kendall)
  mean_z <- as.numeric(model$beta + Rstar %*% model$alpha)
  B <- backsolve(model$chol_R, forwardsolve(t(model$chol_R), t(Rstar)))
  factor <- pmax(0, 1 - colSums(t(Rstar) * B) +
                   (1 - colSums(B))^2 / model$trend_information)
  variance <- pmax(0, model$y_scale^2 * model$sigma2 * factor)
  list(mean = model$y_center + model$y_scale * mean_z,
       variance = variance, sd = sqrt(variance))
}

pcb_expected_improvement <- function(mu, sigma, incumbent) {
  sigma <- pmax(as.numeric(sigma), 1e-12)
  z <- (incumbent - mu) / sigma
  (incumbent - mu) * stats::pnorm(z) + sigma * stats::dnorm(z)
}

pcb_select_ei <- function(ei, raw_index, tolerance = 1e-12) {
  maximum <- max(ei)
  eligible <- which(maximum - ei <= tolerance * max(1, abs(maximum)))
  eligible[[which.min(raw_index[eligible])]]
}
