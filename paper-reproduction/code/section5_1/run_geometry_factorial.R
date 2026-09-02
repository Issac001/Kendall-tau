#!/usr/bin/env Rscript

# Section 5.1: paper-only 3 x 5 strict-foldover comparison.
#
# Methods: FSA-KD (lambda = 0.5), Hamming-maximin SA,
# component-position-L2-maximin SA, and no-search random strict foldover (SRS).
# This publication driver intentionally contains no Kendall-maximin or
# unrestricted arm.

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

sec51_default_settings <- function() {
  paste(unlist(lapply(c(6L, 10L, 20L), function(m) {
    sprintf("%dx%d", m, m * 1:5)
  })), collapse = ",")
}

sec51_geometry_config <- function() {
  smoke <- wcrit_bool_env("SEC51_SMOKE", FALSE)
  cfg <- list(
    experiment = "section5_1_geometry_factorial_paper_only",
    protocol_version = "2026-09-03-paper-only-v1",
    smoke = smoke,
    settings_text = Sys.getenv(
      "SEC51_SETTINGS",
      unset = if (smoke) "6x6" else sec51_default_settings()
    ),
    reps = sec51_int_env("SEC51_REPS", if (smoke) 1L else 50L, 1L),
    master_seed = sec51_int_env("SEC51_MASTER_SEED", 20260831L, 1L),
    workers = sec51_int_env("SEC51_WORKERS", if (smoke) 1L else 12L, 1L),
    proposal_budget = sec51_int_env(
      "SEC51_PROPOSAL_BUDGET", if (smoke) 24L else 6000L, 4L
    ),
    tie_weight = sec51_num_env("SEC51_TIE_WEIGHT", 1e-4, 0, 1),
    sa_T0 = sec51_num_env("SEC51_SA_T0", 0.05, 0),
    sa_alpha = sec51_num_env("SEC51_SA_ALPHA", 0.995, 0, 1),
    fsa_T0 = sec51_num_env("SEC51_FSA_T0", 1.0, 0),
    fsa_alpha = sec51_num_env("SEC51_FSA_ALPHA", 0.997, 0, 1),
    replace_prob = sec51_num_env("SEC51_REPLACE_PROB", 0.25, 0, 1),
    mallows_c = c(1, 4),
    mallows_jitter = 1e-10,
    rank_tol = 1e-9,
    resume = wcrit_bool_env("SEC51_RESUME", TRUE),
    save_design_bank = wcrit_bool_env("SEC51_SAVE_DESIGN_BANK", TRUE),
    output = Sys.getenv(
      "SEC51_OUTPUT",
      unset = file.path(
        sec51_project_root(), "outputs", "section5_1",
        if (smoke) "geometry_smoke" else "geometry_formal"
      )
    )
  )
  scientific <- cfg[setdiff(names(cfg), c(
    "workers", "resume", "save_design_bank", "output"
  ))]
  cfg$protocol_hash <- digest::digest(scientific, algo = "sha256", serialize = TRUE)
  cfg
}

sec51_methods <- function() c("FSA_KD", "Hamming", "L2", "SRS")

sec51_tasks <- function(settings, cfg) {
  tasks <- tidyr::crossing(
    setting_id = seq_len(nrow(settings)),
    rep = seq_len(cfg$reps),
    method = sec51_methods()
  )
  tasks <- dplyr::left_join(
    tasks,
    transform(settings, setting_id = seq_len(nrow(settings))),
    by = "setting_id"
  )
  tasks$task_id <- sprintf(
    "%s__rep%03d__%s", tasks$setting, tasks$rep, tasks$method
  )
  tasks
}

sec51_initial_state <- function(m, n, rep_id, cfg) {
  seed <- wcrit_hash_seed(
    cfg$master_seed, "27", "common-init", m, n, rep_id, 1L
  )
  H <- wcrit_sample_foldover_base_permutations(m, n / 2L, seed = seed)
  list(
    H = H,
    seed = seed,
    design_sha256 = wcrit_matrix_sha256(wcrit_foldover_design(H))
  )
}

sec51_run_fsa <- function(m, n, rep_id, ratio, cfg) {
  seed <- wcrit_hash_seed(
    cfg$master_seed, "27", "native-FSA-Phi050", m, n, rep_id
  )
  started <- proc.time()[["elapsed"]]
  fit <- wcrit_build_weighted_design(
    m = m, n = n, lambda = 0.5, seed = seed,
    max_iter = cfg$proposal_budget, restarts = 1L,
    T0 = cfg$fsa_T0, alpha = cfg$fsa_alpha,
    no_improve_stop = cfg$proposal_budget + 1L,
    max_seconds = 0, foldover = TRUE
  )
  D <- sec51_validate_design(fit$D, m, n, "FSA-KD")
  entry <- list(
    D = D, method = "FSA_KD", paper_label = "FSA-KD",
    m = m, n = n, ratio = ratio, rep = rep_id,
    seed = seed, initial_seed = NA_integer_, move_seed = seed,
    initial_design_sha256 = NA_character_,
    search_geometry = "Phi_0.5",
    search_objective = wcrit_design_metrics(D, lambda = 0.5)$phi_lambda[[1L]],
    proposals = as.integer(fit$n_iter),
    candidate_objective_evaluations = NA_integer_,
    invalid_proposals = NA_integer_,
    elapsed_sec = as.numeric(proc.time()[["elapsed"]] - started)
  )
  list(entry = entry, metrics = sec51_design_metrics(D, entry, cfg))
}

sec51_run_geometric_sa <- function(m, n, rep_id, ratio, method, cfg) {
  geometry <- switch(
    method,
    Hamming = "hamming",
    L2 = "position_l2",
    sec51_stop("Unknown geometric method: %s", method)
  )
  initial <- sec51_initial_state(m, n, rep_id, cfg)
  move_seed <- wcrit_hash_seed(
    cfg$master_seed, "27", "common-moves", m, n, rep_id, 1L
  )
  started <- proc.time()[["elapsed"]]
  fit <- sec51_sa_once(
    initial$H, geometry, cfg$proposal_budget, move_seed, cfg
  )
  D <- sec51_validate_design(fit$D, m, n, method)
  entry <- list(
    D = D, method = method, paper_label = method,
    m = m, n = n, ratio = ratio, rep = rep_id,
    seed = move_seed, initial_seed = initial$seed, move_seed = move_seed,
    initial_design_sha256 = initial$design_sha256,
    search_geometry = geometry,
    search_objective = fit$score$score,
    proposals = fit$proposals,
    candidate_objective_evaluations = fit$candidate_evaluations,
    invalid_proposals = fit$invalid_proposals,
    elapsed_sec = as.numeric(proc.time()[["elapsed"]] - started)
  )
  list(entry = entry, metrics = sec51_design_metrics(D, entry, cfg))
}

sec51_run_srs <- function(m, n, rep_id, ratio, cfg) {
  initial <- sec51_initial_state(m, n, rep_id, cfg)
  D <- sec51_validate_design(
    wcrit_foldover_design(initial$H), m, n, "SRS strict foldover"
  )
  entry <- list(
    D = D, method = "SRS", paper_label = "SRS",
    m = m, n = n, ratio = ratio, rep = rep_id,
    seed = initial$seed, initial_seed = initial$seed, move_seed = NA_integer_,
    initial_design_sha256 = initial$design_sha256,
    search_geometry = "none_no_search", search_objective = NA_real_,
    proposals = 0L, candidate_objective_evaluations = 0L,
    invalid_proposals = 0L, elapsed_sec = 0
  )
  list(entry = entry, metrics = sec51_design_metrics(D, entry, cfg))
}

sec51_run_task <- function(task, cfg) {
  method <- task$method[[1L]]
  args <- list(
    m = task$m[[1L]], n = task$n[[1L]], rep_id = task$rep[[1L]],
    ratio = task$ratio[[1L]], cfg = cfg
  )
  if (identical(method, "FSA_KD")) {
    do.call(sec51_run_fsa, args)
  } else if (identical(method, "SRS")) {
    do.call(sec51_run_srs, args)
  } else {
    do.call(sec51_run_geometric_sa, c(args[1:4], list(method = method), args[5]))
  }
}

sec51_checkpoint_path <- function(checkpoint_dir, task_id) {
  file.path(checkpoint_dir, paste0(task_id, ".rds"))
}

sec51_checkpoint_task <- function(task, cfg, checkpoint_dir) {
  path <- sec51_checkpoint_path(checkpoint_dir, task$task_id[[1L]])
  if (file.exists(path)) {
    if (!cfg$resume) sec51_stop("Checkpoint exists and resume is false: %s", path)
    object <- readRDS(path)
    if (!identical(object$protocol_hash, cfg$protocol_hash) ||
        !identical(object$task_id, task$task_id[[1L]])) {
      sec51_stop("Checkpoint metadata mismatch: %s", path)
    }
    return(data.frame(
      task_id = task$task_id[[1L]], status = "resumed", path = path,
      stringsAsFactors = FALSE
    ))
  }
  result <- sec51_run_task(task, cfg)
  object <- list(
    protocol_hash = cfg$protocol_hash,
    task_id = task$task_id[[1L]],
    design_sha256 = wcrit_matrix_sha256(result$entry$D),
    entry = result$entry,
    metrics = result$metrics
  )
  temporary <- sprintf("%s.tmp-pid%d", path, Sys.getpid())
  saveRDS(object, temporary, version = 3)
  if (!file.rename(temporary, path)) sec51_stop("Could not finalize %s", path)
  data.frame(
    task_id = task$task_id[[1L]], status = "computed", path = path,
    stringsAsFactors = FALSE
  )
}

sec51_read_checkpoint <- function(task, cfg, checkpoint_dir) {
  path <- sec51_checkpoint_path(checkpoint_dir, task$task_id[[1L]])
  object <- readRDS(path)
  if (!identical(object$protocol_hash, cfg$protocol_hash) ||
      !identical(object$task_id, task$task_id[[1L]]) ||
      !identical(object$design_sha256, wcrit_matrix_sha256(object$entry$D))) {
    sec51_stop("Checkpoint integrity failure: %s", path)
  }
  object
}

sec51_welch <- function(raw) {
  metrics <- c("pwo_ms_efficiency", "mallows_det_root_c1")
  comparators <- c("Hamming", "L2", "SRS")
  rows <- list()
  cursor <- 0L
  cells <- unique(raw[c("m", "n", "ratio")])
  for (i in seq_len(nrow(cells))) for (metric in metrics) for (comparator in comparators) {
    keep <- raw$m == cells$m[[i]] & raw$n == cells$n[[i]]
    x <- raw[keep & raw$method == "FSA_KD", metric]
    y <- raw[keep & raw$method == comparator, metric]
    nx <- length(x); ny <- length(y)
    vx <- stats::var(x); vy <- stats::var(y)
    se2 <- vx / nx + vy / ny
    se <- sqrt(max(0, se2))
    denominator <- (vx / nx)^2 / (nx - 1L) + (vy / ny)^2 / (ny - 1L)
    df <- if (is.finite(denominator) && denominator > 0) {
      se2^2 / denominator
    } else Inf
    critical <- if (is.finite(df)) stats::qt(0.975, df) else stats::qnorm(0.975)
    difference <- mean(x) - mean(y)
    cursor <- cursor + 1L
    rows[[cursor]] <- data.frame(
      m = cells$m[[i]], n = cells$n[[i]], ratio = cells$ratio[[i]],
      metric = metric, reference = "FSA-KD", comparator = comparator,
      n_reference = nx, n_comparator = ny,
      reference_mean = mean(x), comparator_mean = mean(y),
      fsa_minus_comparator = difference, welch_se = se, welch_df = df,
      ci95_lower = difference - critical * se,
      ci95_upper = difference + critical * se,
      interval_note = "unadjusted two-sided 95% Welch interval",
      estimand = "end-to-end method difference; positive favors FSA-KD",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

sec51_geometry_main <- function() {
  required <- c("digest", "dplyr", "tidyr", "Rcpp", "gtools")
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) sec51_stop("Missing packages: %s", paste(missing, collapse = ", "))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
  )
  cfg <- sec51_geometry_config()
  settings <- sec51_parse_settings(cfg$settings_text)
  if (!cfg$smoke) {
    expected <- expand.grid(m = c(6L, 10L, 20L), ratio = 1:5)
    expected$n <- expected$m * expected$ratio
    if (nrow(settings) != 15L ||
        !setequal(paste(settings$m, settings$n), paste(expected$m, expected$n)) ||
        cfg$reps != 50L || cfg$proposal_budget != 6000L) {
      sec51_stop("Formal mode requires the paper's 3 x 5 grid, 50 reps, and 6000 proposals")
    }
  }
  output <- path.expand(cfg$output)
  results_dir <- file.path(output, "results")
  checkpoint_dir <- file.path(output, "checkpoints")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  protocol_path <- file.path(output, "protocol_hash.txt")
  if (file.exists(protocol_path)) {
    old_hash <- trimws(readLines(protocol_path, n = 1L, warn = FALSE))
    if (!identical(old_hash, cfg$protocol_hash)) {
      sec51_stop("Output directory has a different protocol hash: %s", output)
    }
  } else {
    writeLines(cfg$protocol_hash, protocol_path)
  }

  tasks <- sec51_tasks(settings, cfg)
  wcrit_write_config(file.path(results_dir, "run_config.csv"), cfg)
  utils::write.csv(settings, file.path(results_dir, "settings.csv"), row.names = FALSE)
  utils::write.csv(tasks, file.path(results_dir, "task_ledger.csv"), row.names = FALSE)

  seed_rows <- do.call(rbind, lapply(seq_len(nrow(settings)), function(i) {
    do.call(rbind, lapply(seq_len(cfg$reps), function(rep_id) {
      m <- settings$m[[i]]; n <- settings$n[[i]]
      data.frame(
        m = m, n = n, ratio = settings$ratio[[i]], rep = rep_id,
        role = c("FSA-KD", "common_initial_Hamming_L2_SRS", "common_moves_Hamming_L2"),
        seed = c(
          wcrit_hash_seed(cfg$master_seed, "27", "native-FSA-Phi050", m, n, rep_id),
          wcrit_hash_seed(cfg$master_seed, "27", "common-init", m, n, rep_id, 1L),
          wcrit_hash_seed(cfg$master_seed, "27", "common-moves", m, n, rep_id, 1L)
        ),
        stringsAsFactors = FALSE
      )
    }))
  }))
  seed_rows$protocol_hash <- cfg$protocol_hash
  utils::write.csv(seed_rows, file.path(results_dir, "seed_ledger.csv"), row.names = FALSE)

  wcrit_load_cpp()
  task_list <- split(tasks, seq_len(nrow(tasks)))
  status <- sec51_parallel_lapply(task_list, function(task) {
    tryCatch(
      sec51_checkpoint_task(task, cfg, checkpoint_dir),
      error = function(error) data.frame(
        task_id = task$task_id[[1L]], status = "error",
        path = NA_character_, error = conditionMessage(error),
        stringsAsFactors = FALSE
      )
    )
  }, workers = cfg$workers)
  status <- dplyr::bind_rows(status)
  utils::write.csv(status, file.path(results_dir, "task_status.csv"), row.names = FALSE)
  if (any(status$status == "error")) {
    sec51_stop("%d task(s) failed; inspect task_status.csv", sum(status$status == "error"))
  }

  checkpoints <- lapply(task_list, sec51_read_checkpoint,
                        cfg = cfg, checkpoint_dir = checkpoint_dir)
  raw <- dplyr::bind_rows(lapply(checkpoints, `[[`, "metrics"))
  raw <- raw[order(raw$m, raw$n, raw$rep, match(raw$method, sec51_methods())), ]
  rownames(raw) <- NULL

  expected_rows <- nrow(settings) * cfg$reps * length(sec51_methods())
  counts <- table(raw$method, raw$m, raw$ratio)
  optimized <- raw$method != "SRS"
  integrity <- data.frame(
    check = c(
      "row_count", "method_set", "cell_replication_count", "strict_foldover",
      "unique_design_per_method_cell_rep", "finite_paper_metrics",
      "optimized_proposal_budget", "SRS_has_no_search"
    ),
    pass = c(
      nrow(raw) == expected_rows,
      setequal(unique(raw$method), sec51_methods()),
      all(counts == cfg$reps),
      all(raw$design_class == "strict_foldover"),
      !anyDuplicated(paste(raw$method, raw$m, raw$n, raw$rep, sep = "|")),
      all(is.finite(raw$pwo_ms_efficiency)) &&
        all(is.finite(raw$mallows_det_root_c1)) &&
        all(is.finite(raw$mallows_det_root_c4)),
      all(raw$proposals[optimized] == cfg$proposal_budget),
      all(raw$proposals[!optimized] == 0L)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(integrity, file.path(results_dir, "integrity_audit.csv"), row.names = FALSE)
  if (any(!integrity$pass)) {
    sec51_stop("Integrity failure: %s", paste(integrity$check[!integrity$pass], collapse = ", "))
  }

  utils::write.csv(raw, file.path(results_dir, "four_method_raw_metrics.csv"), row.names = FALSE)
  welch <- sec51_welch(raw)
  utils::write.csv(welch, file.path(results_dir, "fsa_vs_baseline_welch.csv"), row.names = FALSE)
  if (cfg$save_design_bank) {
    bank <- lapply(checkpoints, function(object) object$entry$D)
    names(bank) <- tasks$task_id
    saveRDS(bank, file.path(results_dir, "design_bank.rds"), version = 3)
  }
  writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))
  wcrit_write_json(file.path(results_dir, "completion.json"), list(
    status = "complete", protocol_hash = cfg$protocol_hash,
    rows = nrow(raw), methods = sec51_methods(), settings = nrow(settings),
    reps = cfg$reps, proposal_budget = cfg$proposal_budget
  ))
  message("Section 5.1 geometry experiment complete: ", output)
  invisible(list(config = cfg, raw = raw, output = output))
}

if (sec51_is_main("run_geometry_factorial.R")) sec51_geometry_main()
