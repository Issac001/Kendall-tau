#!/usr/bin/env Rscript

# Appendix B.1: FSA-KD weight-sensitivity experiment.
#
# This is a paper-only reconstruction of the historical Experiment 01.  The
# three odd-n panels are deliberately marked exploratory: their A/B values are
# computed after deleting one run from an even foldover parent and are not
# valid instances of the paper's strict-foldover normalization.  Formal claims
# use only the 12 even-n panels.

options(stringsAsFactors = FALSE, warn = 1)

.sec51_args <- commandArgs(trailingOnly = FALSE)
.sec51_hit <- grep("^--file=", .sec51_args)
.sec51_here <- if (length(.sec51_hit)) {
  dirname(normalizePath(sub("^--file=", "", .sec51_args[.sec51_hit[[1L]]]),
                        winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
source(file.path(.sec51_here, "section5_1_helpers.R"))

# Recreate the numerical convention used when the historical odd-run panels
# were generated.  This function is used only to report those exploratory
# reduced designs; it does not make them strict foldover.
sec51_historical_lambda_metrics <- function(D, lambda) {
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  n <- nrow(D)
  m <- ncol(D)
  dm <- wcrit_kendall_dmat(D)
  distances <- dm[upper.tri(dm)]
  k_min <- min(distances)
  k_ave <- mean(distances)
  k_m2 <- mean(distances^2)
  constants <- wcrit_constants(m, n)
  A <- if (constants$C1 > 1) {
    (k_min - 1) / (constants$C1 - 1)
  } else 1
  B <- (constants$U2 - k_m2) / (constants$U2 - constants$C2)
  data.frame(
    k_min = k_min, k_ave = k_ave, k_m2 = k_m2,
    A = A, B = B, phi_lambda = lambda * A + (1 - lambda) * B,
    C1 = constants$C1, C2 = constants$C2, U2 = constants$U2,
    strict_foldover = wcrit_is_strict_foldover(D),
    interpretation = if (wcrit_is_strict_foldover(D)) {
      "strict_foldover"
    } else {
      "historical_exploratory_odd_run_reduction"
    },
    stringsAsFactors = FALSE
  )
}

sec51_lambda_config <- function() {
  smoke <- wcrit_bool_env("SEC51_LAMBDA_SMOKE", FALSE)
  list(
    experiment = "section5_1_lambda_sensitivity",
    protocol_version = "2026-09-03-historical-reconstruction-v1",
    smoke = smoke,
    m_values = wcrit_parse_int_vec(
      Sys.getenv("SEC51_LAMBDA_M", unset = if (smoke) "5" else "5,10,20"),
      c(5L, 10L, 20L)
    ),
    lambda_values = wcrit_parse_num_vec(
      Sys.getenv(
        "SEC51_LAMBDAS",
        unset = if (smoke) "0.5" else "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1"
      ),
      seq(0, 1, by = 0.1)
    ),
    reps = sec51_int_env("SEC51_LAMBDA_REPS", if (smoke) 1L else 20L, 1L),
    master_seed = sec51_int_env("SEC51_LAMBDA_MASTER_SEED", 20260523L, 1L),
    max_iter = sec51_int_env(
      "SEC51_LAMBDA_MAX_ITER", if (smoke) 24L else 6000L, 1L
    ),
    restarts = sec51_int_env("SEC51_LAMBDA_RESTARTS", if (smoke) 1L else 8L, 1L),
    no_improve_stop = sec51_int_env(
      "SEC51_LAMBDA_NO_IMPROVE", if (smoke) 25L else 800L, 1L
    ),
    workers = sec51_int_env("SEC51_LAMBDA_WORKERS", if (smoke) 1L else 8L, 1L),
    keep_published_missing = wcrit_bool_env(
      "SEC51_LAMBDA_KEEP_PUBLISHED_MISSING", !smoke
    ),
    save_designs = wcrit_bool_env("SEC51_LAMBDA_SAVE_DESIGNS", !smoke),
    output = Sys.getenv(
      "SEC51_LAMBDA_OUTPUT",
      unset = file.path(
        sec51_project_root(), "outputs", "section5_1",
        if (smoke) "lambda_sensitivity_smoke" else "lambda_sensitivity"
      )
    )
  )
}

sec51_lambda_failure_row <- function(m, n, lambda, rep_id, seed) {
  data.frame(
    experiment = "lambda_sensitivity", m = m, n = n, lambda = lambda,
    rep = rep_id, seed = seed, design_path = NA_character_,
    objective = NA_real_, n_iter = NA_integer_, elapsed_sec = NA_real_,
    k_min = NA_real_, k_ave = NA_real_, k_m2 = NA_real_,
    A = NA_real_, B = NA_real_, phi_lambda = NA_real_,
    C1 = NA_real_, C2 = NA_real_, U2 = NA_real_,
    strict_foldover = NA, interpretation = "historical_failed_run",
    stringsAsFactors = FALSE
  )
}

sec51_lambda_one <- function(row, cfg, design_dir) {
  m <- as.integer(row$m)
  n <- as.integer(row$n)
  lambda <- as.numeric(row$lambda)
  rep_id <- as.integer(row$rep)
  seed <- wcrit_hash_seed(cfg$master_seed, "distance", m, n, lambda, rep_id)

  # The archived published data contain one failed job.  Keeping that missing
  # row is the default for exact reproduction of the reported 3299/3300 count.
  historical_missing <- m == 5L && n == 10L &&
    abs(lambda) < 1e-12 && rep_id == 1L
  if (cfg$keep_published_missing && historical_missing) {
    return(sec51_lambda_failure_row(m, n, lambda, rep_id, seed))
  }

  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    wcrit_build_weighted_design(
      m = m, n = n, lambda = lambda, seed = seed,
      max_iter = cfg$max_iter, restarts = cfg$restarts,
      no_improve_stop = cfg$no_improve_stop, foldover = TRUE
    ),
    error = function(error) NULL
  )
  if (is.null(fit)) return(sec51_lambda_failure_row(m, n, lambda, rep_id, seed))

  metrics <- sec51_historical_lambda_metrics(fit$D, lambda)
  path <- NA_character_
  if (cfg$save_designs) {
    tag <- sprintf("%03d", as.integer(round(lambda * 100)))
    path <- file.path(
      design_dir,
      sprintf("lambda_m%02d_n%03d_lambda%s_rep%02d.rds", m, n, tag, rep_id)
    )
    saveRDS(list(
      D = fit$D, m = m, n = n, lambda = lambda, rep = rep_id, seed = seed
    ), path, version = 3)
  }
  cbind(
    data.frame(
      experiment = "lambda_sensitivity", m = m, n = n, lambda = lambda,
      rep = rep_id, seed = seed,
      # Store a path relative to this run's output directory so that the
      # resulting CSV can be moved between machines without leaking a host
      # filesystem prefix.
      design_path = if (is.na(path)) NA_character_ else
        file.path("designs", basename(path)),
      objective = metrics$phi_lambda[[1L]], n_iter = fit$n_iter,
      elapsed_sec = as.numeric(proc.time()[["elapsed"]] - started),
      stringsAsFactors = FALSE
    ),
    metrics
  )
}

sec51_lambda_summary <- function(raw) {
  groups <- split(raw, interaction(raw$m, raw$n, raw$lambda,
                                   drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(data) {
    result <- data.frame(
      m = data$m[[1L]], n = data$n[[1L]], lambda = data$lambda[[1L]],
      n_runs = nrow(data), n_finite = sum(is.finite(data$A) & is.finite(data$B)),
      formal_even_n = data$n[[1L]] %% 2L == 0L,
      stringsAsFactors = FALSE
    )
    for (metric in c(
      "k_min", "k_ave", "k_m2", "A", "B", "phi_lambda",
      "elapsed_sec", "n_iter"
    )) {
      x <- as.numeric(data[[metric]])
      finite <- x[is.finite(x)]
      result[[paste0(metric, "_mean")]] <- if (length(finite)) mean(finite) else NA_real_
      result[[paste0(metric, "_sd")]] <- if (length(finite) > 1L) stats::sd(finite) else 0
      result[[paste0(metric, "_se")]] <- if (length(finite) > 1L) {
        stats::sd(finite) / sqrt(length(finite))
      } else 0
    }
    result
  })
  answer <- do.call(rbind, rows)
  answer[order(answer$m, answer$n, answer$lambda), ]
}

sec51_lambda_main <- function() {
  cfg <- sec51_lambda_config()
  if (!cfg$smoke) {
    if (!identical(as.integer(cfg$m_values), c(5L, 10L, 20L)) ||
        !isTRUE(all.equal(cfg$lambda_values, seq(0, 1, by = 0.1))) ||
        cfg$reps != 20L || cfg$max_iter != 6000L || cfg$restarts != 8L ||
        cfg$master_seed != 20260523L || cfg$no_improve_stop != 800L) {
      sec51_stop("Formal mode requires the Appendix B.1 configuration")
    }
  }
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  output <- path.expand(cfg$output)
  design_dir <- file.path(output, "designs")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (cfg$save_designs) dir.create(design_dir, recursive = TRUE, showWarnings = FALSE)
  grid <- expand.grid(
    m = cfg$m_values, n_multiplier = 1:5,
    lambda = cfg$lambda_values, rep = seq_len(cfg$reps),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid$n <- grid$m * grid$n_multiplier
  grid$n_multiplier <- NULL
  wcrit_load_cpp()
  rows <- sec51_parallel_lapply(
    split(grid, seq_len(nrow(grid))),
    sec51_lambda_one, cfg = cfg, design_dir = design_dir,
    workers = cfg$workers
  )
  raw <- dplyr::bind_rows(rows)
  raw <- raw[order(raw$m, raw$n, raw$lambda, raw$rep), ]
  summary <- sec51_lambda_summary(raw)
  wcrit_write_config(file.path(output, "run_config.csv"), cfg)
  utils::write.csv(raw, file.path(output, "lambda_sensitivity_raw.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output, "lambda_sensitivity_summary.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
  message(sprintf(
    "Appendix B.1 lambda sensitivity complete: %d/%d finite searches; %s",
    sum(is.finite(raw$A) & is.finite(raw$B)), nrow(raw), output
  ))
  invisible(list(config = cfg, raw = raw, summary = summary, output = output))
}

if (sec51_is_main("run_lambda_sensitivity.R")) sec51_lambda_main()
