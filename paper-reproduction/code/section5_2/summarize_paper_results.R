# Publication-only postprocessing for Section 5.2.
#
# Produces the n=48/n=60 table and paired baseline-minus-FSA-KD intervals for
# exactly the four methods and four response scenarios reported in the paper.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1L]]),
                             winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."),
                              winslash = "/", mustWork = TRUE)
stopf <- function(...) stop(sprintf(...), call. = FALSE)

run_raw <- Sys.getenv("SEC52_RUN_DIR", unset = "")
if (!nzchar(run_raw)) stopf("Set SEC52_RUN_DIR to a completed Section 5.2 run")
run_dir <- if (grepl("^/", run_raw)) run_raw else file.path(project_root, run_raw)
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
out_raw <- Sys.getenv("SEC52_OUTPUT_DIR", unset = file.path(run_dir, "paper_sources"))
out_dir <- if (grepl("^/", out_raw)) out_raw else file.path(project_root, out_raw)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  pwo = file.path(run_dir, "raw", "pwo_prediction_raw.csv"),
  gp = file.path(run_dir, "raw", "gp_prediction_raw.csv"),
  cross = file.path(run_dir, "raw", "cross_model_loss_raw.csv"),
  completion = file.path(run_dir, "results", "completion.json")
)
if (any(!file.exists(files))) {
  stopf("Run is incomplete; missing: %s",
        paste(names(files)[!file.exists(files)], collapse = ", "))
}

methods <- c("FSA_lambda050", "SA_hamming_foldover",
             "SA_l2_position_foldover", "random_foldover")
labels <- c(FSA_lambda050 = "FSA-KD", SA_hamming_foldover = "Hamming",
            SA_l2_position_foldover = "L2", random_foldover = "SRS")
scenarios <- c("pwo_snr2", "pwo_snr5", "gp_c1", "gp_c4")

pwo <- utils::read.csv(files[["pwo"]], stringsAsFactors = FALSE,
                       check.names = FALSE)
gp <- utils::read.csv(files[["gp"]], stringsAsFactors = FALSE,
                      check.names = FALSE)
cross <- utils::read.csv(files[["cross"]], stringsAsFactors = FALSE,
                         check.names = FALSE)

if (!all(methods %in% unique(cross$method_id)) ||
    !all(scenarios %in% unique(cross$scenario_id)) ||
    !all(c(48L, 60L) %in% unique(cross$n))) {
  stopf("Input does not contain the four paper methods/scenarios at n=48,60")
}
cross <- cross[cross$method_id %in% methods &
                 cross$scenario_id %in% scenarios & cross$n %in% c(48L, 60L),
               , drop = FALSE]
pwo <- pwo[pwo$method_id %in% methods &
             pwo$scenario_id %in% c("pwo_snr2", "pwo_snr5") &
             pwo$n %in% c(48L, 60L), , drop = FALSE]
gp <- gp[gp$method_id %in% methods & gp$c_true %in% c(1, 4) &
           gp$fit_mode == "estimated_c" & gp$n %in% c(48L, 60L),
         , drop = FALSE]
if (!nrow(gp) || !setequal(unique(gp$c_true), c(1, 4))) {
  stopf("Mallows-GP input lacks estimated-c fits at c=1 and c=4")
}

cell_mean <- stats::aggregate(
  cross$absolute_loss,
  list(n = cross$n, method_id = cross$method_id,
       scenario_id = cross$scenario_id), mean
)
names(cell_mean)[4L] <- "mean_nrmse"
wide <- reshape(cell_mean, idvar = c("n", "method_id"),
                timevar = "scenario_id", direction = "wide")
names(wide) <- sub("^mean_nrmse\\.", "nrmse_", names(wide))

ipv_source <- pwo[pwo$scenario_id == "pwo_snr2",
                  c("rep", "n", "method_id", "ipv_factor_full")]
if (anyDuplicated(ipv_source[c("rep", "n", "method_id")])) {
  stopf("PWO IPV source is not unique by replication, n, and method")
}
ipv <- stats::aggregate(ipv_source$ipv_factor_full,
                        list(n = ipv_source$n, method_id = ipv_source$method_id),
                        mean)
names(ipv)[3L] <- "pwo_ipv_full_space"
table_out <- merge(wide, ipv, by = c("n", "method_id"), all = TRUE,
                   sort = FALSE)
table_out$method <- unname(labels[table_out$method_id])
table_out <- table_out[order(table_out$n,
                             match(table_out$method_id, methods)),
                       c("n", "method", "method_id", "pwo_ipv_full_space",
                         paste0("nrmse_", scenarios)), drop = FALSE]
if (nrow(table_out) != 8L || anyNA(table_out)) {
  stopf("The Section 5.2 table failed its 2 x 4 cardinality audit")
}
utils::write.csv(table_out,
                 file.path(out_dir, "section5_2_main_table.csv"),
                 row.names = FALSE)

# One frozen common-replication bootstrap tape is shared across all cells.
B <- as.integer(Sys.getenv("SEC52_BOOTSTRAP_B", unset = "10000"))
bootstrap_seed <- as.integer(Sys.getenv("SEC52_BOOTSTRAP_SEED",
                                        unset = "20260901"))
if (!is.finite(B) || B < 100L || !is.finite(bootstrap_seed)) {
  stopf("Invalid bootstrap configuration")
}
rep_ids <- sort(unique(cross$rep))
set.seed(bootstrap_seed)
counts <- t(replicate(B, tabulate(sample(seq_along(rep_ids),
                                     length(rep_ids), replace = TRUE),
                                  nbins = length(rep_ids))))

paired_rows <- list()
cursor <- 0L
for (n_value in sort(unique(cross$n))) {
  for (scenario in scenarios) {
    reference <- cross[cross$n == n_value & cross$scenario_id == scenario &
                         cross$method_id == "FSA_lambda050",
                       c("rep", "absolute_loss")]
    reference <- reference[match(rep_ids, reference$rep), ]
    for (method in methods[-1L]) {
      comparator <- cross[cross$n == n_value &
                            cross$scenario_id == scenario &
                            cross$method_id == method,
                          c("rep", "absolute_loss")]
      comparator <- comparator[match(rep_ids, comparator$rep), ]
      if (anyNA(reference) || anyNA(comparator)) stopf("Incomplete paired cell")
      difference <- comparator$absolute_loss - reference$absolute_loss
      boot <- as.numeric(counts %*% difference) / length(rep_ids)
      ci <- stats::quantile(boot, c(0.025, 0.975), type = 8, names = FALSE)
      cursor <- cursor + 1L
      paired_rows[[cursor]] <- data.frame(
        n = n_value, scenario_id = scenario,
        comparator = unname(labels[[method]]), comparator_id = method,
        reference = "FSA-KD", reference_id = "FSA_lambda050",
        mean_difference_baseline_minus_fsa = mean(difference),
        ci95_low = ci[[1L]], ci95_high = ci[[2L]],
        bootstrap_B = B, bootstrap_seed = bootstrap_seed,
        ci_method = "common-replication percentile bootstrap, quantile type 8",
        stringsAsFactors = FALSE
      )
    }
  }
}
paired <- do.call(rbind, paired_rows)
if (nrow(paired) != 24L) stopf("Expected 2 x 4 x 3 paired contrasts")
utils::write.csv(paired,
                 file.path(out_dir, "section5_2_paired_intervals.csv"),
                 row.names = FALSE)
message("Wrote publication sources under ", out_dir)
