#!/usr/bin/env Rscript

# Rebuild the two compact inputs used by the paper-only Experiment 24 runner.
# The design constructor has exactly four arms and the oracle constructor has
# exactly one response scenario: g100_w050, (gamma, omega) = (1, 0.5).

options(stringsAsFactors = FALSE, warn = 1)
.args <- commandArgs(trailingOnly = FALSE)
.idx <- grep("^--file=", .args)
.this_dir <- if (length(.idx)) {
  dirname(normalizePath(sub("^--file=", "", .args[.idx[[1L]]]), mustWork = TRUE))
} else normalizePath(getwd(), mustWork = TRUE)
.repro_root <- normalizePath(file.path(.this_dir, "..", ".."), mustWork = TRUE)
.common_dir <- file.path(.repro_root, "code", "common")
.reference_dir <- file.path(.repro_root, "data", "frozen", "section5_4")

for (pkg in c("Rcpp", "digest", "gtools", "dplyr", "tidyr", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R package is missing: ", pkg)
}
source(file.path(.common_dir, "wcrit_common.R"), local = FALSE)
source(file.path(.common_dir, "wcrit_maximin_dist.R"), local = FALSE)
source(file.path(.common_dir, "case_study_common.R"), local = FALSE)
Rcpp::sourceCpp(file.path(.common_dir, "sa_core.cpp"), rebuild = FALSE, showOutput = FALSE)
source(file.path(.this_dir, "pcb_common.R"), local = FALSE)

read_int <- function(name, default, lower, upper) {
  x <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(x) || x < lower || x > upper) stop(name, " must be in [", lower, ",", upper, "]")
  x
}
read_bool <- function(name, default) {
  x <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
  if (!x %in% c("true", "false", "1", "0", "yes", "no")) stop("Invalid Boolean: ", name)
  x %in% c("true", "1", "yes")
}

reps <- read_int("SEC54_BUILD_REPS", 30L, 1L, 30L)
build_workers <- read_int("SEC54_BUILD_WORKERS", 1L, 1L, 30L)
build_exact_oracle <- read_bool("SEC54_BUILD_EXACT_ORACLE", TRUE)
out_dir <- path.expand(Sys.getenv(
  "SEC54_BUILD_OUT", unset = file.path(.repro_root, "outputs", "section5_4_frozen_rebuild")
))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coordinate_path <- file.path(.repro_root, "data", "case_studies", "d493_first10_holes.csv")
coordinates <- utils::read.csv(coordinate_path, check.names = FALSE)
C <- pcb_build_cost(coordinates)
dwell <- 0.20

build_profiles <- function() {
  cv <- 0.25
  sdlog <- sqrt(log1p(cv^2))
  meanlog <- -0.5 * sdlog^2
  base <- stats::qlnorm((1:10 - 0.5) / 10, meanlog, sdlog)
  base <- base / mean(base)
  typical_gap <- stats::median(C[row(C) != col(C)] + dwell)
  tau <- -typical_gap / log(0.70)
  component_names <- c("uW", "sW", "uH", "sH")
  profiles <- vector("list", 5L)
  rows <- vector("list", 5L)
  for (profile_id in 1:5) {
    values <- list(); component_seeds <- integer(4L)
    for (j in seq_along(component_names)) {
      seed <- wcrit_hash_seed(20260824L, "24", "physics-profile", profile_id,
                              component_names[[j]])
      set.seed(wcrit_safe_seed(seed))
      values[[component_names[[j]]]] <- base[sample.int(10L, 10L)]
      values[[component_names[[j]]]] <-
        values[[component_names[[j]]]] / mean(values[[component_names[[j]]]])
      component_seeds[[j]] <- seed
    }
    profiles[[profile_id]] <- c(values, list(
      profile_id = profile_id, tau_H = tau, typical_gap = typical_gap,
      component_seeds = stats::setNames(component_seeds, component_names)
    ))
    rows[[profile_id]] <- data.frame(
      profile_id = profile_id, hole = 1:10,
      uW = values$uW, sW = values$sW, uH = values$uH, sH = values$sH,
      tau_H = tau, typical_gap = typical_gap,
      seed_uW = component_seeds[[1L]], seed_sW = component_seeds[[2L]],
      seed_uH = component_seeds[[3L]], seed_sH = component_seeds[[4L]]
    )
  }
  assignment_seed <- wcrit_hash_seed(20260824L, "24", "balanced-profile-assignment")
  ids <- rep(1:5, length.out = 30L)
  set.seed(wcrit_safe_seed(assignment_seed))
  assignment <- data.frame(
    repetition = 1:30, profile_id = sample(ids, 30L, replace = FALSE),
    assignment_seed = assignment_seed
  )
  list(profiles = profiles, parameters = do.call(rbind, rows), assignment = assignment)
}

physics <- build_profiles()
reference_profiles <- utils::read.csv(
  file.path(.reference_dir, "physics_parameter_profiles.csv"), check.names = FALSE
)
reference_assignment <- utils::read.csv(
  file.path(.reference_dir, "replication_profile_assignment.csv"), check.names = FALSE
)
profile_pass <- isTRUE(all.equal(physics$parameters, reference_profiles,
                                 tolerance = 1e-14, check.attributes = FALSE))
assignment_pass <- identical(physics$assignment, reference_assignment)
if (!profile_pass || !assignment_pass) stop("Physics profile reconstruction differs from frozen input")

build_one_design_set <- function(rep_id) {
  fsa_seed <- wcrit_hash_seed(20260820L, "22", "initial-design", rep_id, "fsa_lambda05")
  fsa <- wcrit_build_weighted_design(
    m = 10L, n = 20L, lambda = 0.5, seed = fsa_seed,
    max_iter = 6000L, restarts = 1L, T0 = 1.0, alpha = 0.997,
    no_improve_stop = 6001L, max_seconds = 0, foldover = TRUE
  )$D

  common_init_seed <- wcrit_hash_seed(
    20260820L, "22", "unrestricted-common-init", rep_id, 1L
  )
  common_move_seed <- wcrit_hash_seed(
    20260820L, "22", "unrestricted-common-moves", rep_id, 1L
  )
  D_init <- wcrit_sample_unique_permutations(10L, 20L, seed = common_init_seed)
  hamming <- case_run_unrestricted_sa(
    D_init, "hamming", 6000L, common_move_seed,
    tie_weight = 1e-4, replace_prob = 0.25, T0 = 0.05, Tmin = 1e-8,
    proposal_multiplier = 50L
  )$D
  position_l2 <- case_run_unrestricted_sa(
    D_init, "l2_position", 6000L, common_move_seed,
    tie_weight = 1e-4, replace_prob = 0.25, T0 = 0.05, Tmin = 1e-8,
    proposal_multiplier = 50L
  )$D
  srs_seed <- wcrit_hash_seed(20260820L, "22", "initial-design", rep_id, "srs")
  srs <- wcrit_sample_unique_permutations(10L, 20L, seed = srs_seed)

  out <- list(
    fsa_lambda05 = fsa,
    unrestricted_hamming = hamming,
    unrestricted_position_l2 = position_l2,
    srs = srs
  )
  lapply(out, function(D) case_validate_design(D, 10L, 20L))
}

reference_designs <- readRDS(file.path(.reference_dir, "initial_designs_paper.rds"))
rebuild_and_check <- function(rep_id) {
  message(sprintf("Rebuilding paper design set %d/%d", rep_id, reps))
  z <- build_one_design_set(rep_id)
  expected <- reference_designs$designs[[rep_id]]
  observed_hash <- vapply(z, wcrit_matrix_sha256, character(1L))
  expected_hash <- vapply(expected, wcrit_matrix_sha256, character(1L))
  if (!identical(observed_hash, expected_hash)) {
    stop("Initial-design hash mismatch at replication ", rep_id)
  }
  z
}
if (build_workers > 1L && .Platform$OS.type != "windows") {
  designs <- parallel::mclapply(seq_len(reps), rebuild_and_check,
                                mc.cores = build_workers, mc.preschedule = FALSE)
} else {
  designs <- lapply(seq_len(reps), rebuild_and_check)
}
failed <- which(vapply(designs, inherits, logical(1L), what = "try-error"))
if (length(failed)) {
  stop("Design reconstruction failed for replication(s): ", paste(failed, collapse = ","),
       ". Rerun those replications with SEC54_BUILD_WORKERS=1 for full diagnostics.")
}

design_output <- list(
  experiment = "24", scenario = "g100_w050", methods = pcb_methods,
  n_init = 20L, replications = reps, designs = designs,
  reconstruction = "paper-only seed-based reconstruction; 6000 search units"
)
saveRDS(design_output, file.path(out_dir, "initial_designs_paper.rds"), version = 3)

design_audit <- do.call(rbind, lapply(seq_len(reps), function(rep_id) {
  do.call(rbind, lapply(pcb_methods, function(method) data.frame(
    repetition = rep_id, method = method,
    design_sha256 = wcrit_matrix_sha256(designs[[rep_id]][[method]]),
    frozen_sha256 = wcrit_matrix_sha256(reference_designs$designs[[rep_id]][[method]]),
    exact_match = identical(designs[[rep_id]][[method]],
                            reference_designs$designs[[rep_id]][[method]])
  )))
}))
utils::write.csv(design_audit, file.path(out_dir, "initial_design_hash_audit.csv"), row.names = FALSE)

stream_env <- new.env(parent = globalenv())
Rcpp::cppFunction(
  depends = "Rcpp", env = stream_env,
  code = '
  Rcpp::List pcb_paper_stream_components_cpp(
      const Rcpp::NumericMatrix& cost,
      const Rcpp::NumericVector& uW,
      const Rcpp::NumericVector& sW,
      const Rcpp::NumericVector& uH,
      const Rcpp::NumericVector& sH,
      const double tau,
      const double dwell) {
    const int m = uW.size();
    long long total = 1;
    for (int j = 2; j <= m; ++j) total *= j;
    Rcpp::NumericVector move(total), wear(total), thermal(total);
    std::vector<int> p(m);
    for (int j = 0; j < m; ++j) p[j] = j;
    long long idx = 0;
    do {
      double mv = 0.0, wr = 0.0, th = 0.0, W = 0.0, H = 0.0;
      int prev = -1;
      for (int r = 0; r < m; ++r) {
        const int hole = p[r];
        wr += sW[hole] * W;
        W += uW[hole];
        const double edge = cost(prev + 1, hole + 1);
        const double Hbar = std::exp(-(edge + dwell) / tau) * H;
        th += sH[hole] * Hbar;
        H = Hbar + uH[hole];
        mv += edge;
        prev = hole;
      }
      move[idx] = mv + cost(prev + 1, 0) + m * dwell;
      wear[idx] = wr;
      thermal[idx] = th;
      ++idx;
    } while (std::next_permutation(p.begin(), p.end()));
    return Rcpp::List::create(
      Rcpp::Named("move") = move,
      Rcpp::Named("wear") = wear,
      Rcpp::Named("thermal") = thermal
    );
  }'
)

unrank_lex <- function(index, m = 10L) {
  rank <- as.numeric(index) - 1
  available <- seq_len(m); route <- integer(m)
  for (position in seq_len(m)) {
    block <- factorial(m - position)
    digit <- if (block == 0) 0 else floor(rank / block)
    rank <- if (block == 0) 0 else rank %% block
    route[[position]] <- available[[digit + 1L]]
    available <- available[-(digit + 1L)]
  }
  route
}

if (build_exact_oracle) {
  oracle_rows <- vector("list", 5L)
  profile_scenarios <- vector("list", 5L)
  for (profile_id in 1:5) {
    message(sprintf("Enumerating all 10! routes for physics profile %d/5", profile_id))
    p <- physics$profiles[[profile_id]]
    raw <- stream_env$pcb_paper_stream_components_cpp(
      C, p$uW, p$sW, p$uH, p$sH, p$tau_H, dwell
    )
    zW <- (raw$wear - mean(raw$wear)) / stats::sd(raw$wear)
    zH <- (raw$thermal - mean(raw$thermal)) / stats::sd(raw$thermal)
    mixed <- 0.5 * zW + 0.5 * zH
    mixed_sd <- stats::sd(mixed)
    mixed_z <- (mixed - mean(mixed)) / mixed_sd
    motion_sd <- stats::sd(raw$move)
    response <- raw$move + motion_sd * mixed_z
    best_id <- which.min(response)
    optimum <- response[[best_id]]
    threshold <- as.numeric(stats::quantile(response, 0.001, names = FALSE, type = 1))
    route <- unrank_lex(best_id)
    oracle_rows[[profile_id]] <- data.frame(
      profile_id = profile_id, scenario = "g100_w050", gamma = 1, omega = 0.5,
      n_routes = length(response), full_factorial_n = factorial(10), exact_is_full = TRUE,
      optimum = optimum, optimum_route = paste(route, collapse = "-"),
      optimum_tie_count = sum(abs(response - optimum) <= 1e-9),
      top_probability = 0.001, top_threshold = threshold,
      top_count_including_ties = sum(response <= threshold + 1e-9),
      response_mean = mean(response), response_sd = stats::sd(response),
      response_min = min(response), response_median = stats::median(response),
      response_max = max(response), motion_mean = mean(raw$move), motion_sd = motion_sd,
      memory_component_sd = stats::sd(motion_sd * mixed_z),
      memory_to_motion_sd_ratio = stats::sd(motion_sd * mixed_z) / motion_sd
    )
    profile_scenarios[[profile_id]] <- list(
      scenario = "g100_w050", gamma = 1, omega = 0.5,
      wear_mean = mean(raw$wear), wear_sd = stats::sd(raw$wear),
      thermal_mean = mean(raw$thermal), thermal_sd = stats::sd(raw$thermal),
      mixed_z_mean = mean(mixed), mixed_z_sd = mixed_sd,
      motion_mean = mean(raw$move), motion_sd = motion_sd,
      optimum = optimum, top_threshold = threshold, random_mean = mean(response),
      response_sd = stats::sd(response), optimum_route = route
    )
    rm(raw, zW, zH, mixed, mixed_z, response); gc(FALSE)
  }
  names(profile_scenarios) <- as.character(1:5)
  oracle <- do.call(rbind, oracle_rows)
  reference_oracle <- readRDS(file.path(.reference_dir, "g100_w050_oracle.rds"))
  numeric_columns <- names(oracle)[vapply(oracle, is.numeric, logical(1L))]
  max_delta <- max(vapply(numeric_columns, function(column) {
    max(abs(oracle[[column]] - reference_oracle$oracle[[column]]))
  }, numeric(1L)))
  route_pass <- identical(oracle$optimum_route, reference_oracle$oracle$optimum_route)
  object_delta <- max(vapply(1:5, function(profile_id) {
    a <- profile_scenarios[[profile_id]]; b <- reference_oracle$profile_scenarios[[profile_id]]
    fields <- names(a)[vapply(a, is.numeric, logical(1L)) & lengths(a) == 1L]
    max(abs(unlist(a[fields]) - unlist(b[fields])))
  }, numeric(1L)))
  oracle_pass <- max_delta <= 1e-10 && object_delta <= 1e-10 && route_pass
  utils::write.csv(data.frame(
    max_oracle_numeric_delta = max_delta,
    max_scenario_object_numeric_delta = object_delta,
    optimum_routes_identical = route_pass, pass = oracle_pass
  ), file.path(out_dir, "oracle_hash_audit.csv"), row.names = FALSE)
  if (!oracle_pass) stop("Rebuilt g100_w050 oracle differs from frozen reference")
  saveRDS(list(
    experiment = "24", scenario = "g100_w050", gamma = 1, omega = 0.5,
    oracle = oracle, profile_scenarios = profile_scenarios,
    reconstruction = "paper-only full enumeration of 10! routes per profile"
  ), file.path(out_dir, "g100_w050_oracle.rds"), version = 3)
}

utils::write.csv(physics$parameters, file.path(out_dir, "physics_parameter_profiles.csv"), row.names = FALSE)
utils::write.csv(physics$assignment, file.path(out_dir, "replication_profile_assignment.csv"), row.names = FALSE)
utils::write.csv(data.frame(
  check = c("physics_profiles_match", "profile_assignment_matches",
            "all_rebuilt_design_hashes_match", "exact_oracle_requested"),
  pass = c(profile_pass, assignment_pass, all(design_audit$exact_match), build_exact_oracle)
), file.path(out_dir, "reconstruction_audit.csv"), row.names = FALSE)
message("Paper-only frozen-input reconstruction complete: ", normalizePath(out_dir, mustWork = TRUE))
