# Build the Section 5.3 table from a completed paper-only Experiment 30 run.
#
# Usage from paper-reproduction/:
#   SEC53_RUN_DIR=outputs/wcrit/<formal-run> \
#     Rscript code/section5_3/make_paper_table.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1L]]),
                             winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_file), "..", ".."),
                              winslash = "/", mustWork = TRUE)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
run_raw <- Sys.getenv("SEC53_RUN_DIR", unset = "")
if (!nzchar(run_raw)) stopf("Set SEC53_RUN_DIR to a completed Section 5.3 run")
run_dir <- if (grepl("^/", run_raw)) run_raw else file.path(project_root, run_raw)
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

required <- c(
  recommendation = file.path(run_dir, "raw", "recommendation.csv"),
  prediction = file.path(run_dir, "results", "initial_prediction_summary.csv"),
  completed = file.path(run_dir, "COMPLETED")
)
if (any(!file.exists(required))) {
  stopf("Run is incomplete; missing: %s",
        paste(names(required)[!file.exists(required)], collapse = ", "))
}

methods <- c(
  "Exact_Phi_lambda050", "OofA_OA", "Unrestricted_Hamming",
  "Unrestricted_Position_L2", "SRS"
)
labels <- c(
  Exact_Phi_lambda050 = "FSA-KD", OofA_OA = "OofA-OA",
  Unrestricted_Hamming = "Hamming",
  Unrestricted_Position_L2 = "L2", SRS = "SRS"
)

prediction <- utils::read.csv(required[["prediction"]],
                              stringsAsFactors = FALSE, check.names = FALSE)
prediction <- prediction[
  prediction$model == "Intercept_Mallows_GP" & prediction$bo_step == 0L &
    prediction$domain == "currently_unobserved" &
    prediction$method %in% methods,
  c("method", "nrmse_common_mean"), drop = FALSE
]
if (nrow(prediction) != length(methods) || anyDuplicated(prediction$method)) {
  stopf("Expected one initial held-out prediction row for each paper method")
}

recommendation <- utils::read.csv(required[["recommendation"]],
                                  stringsAsFactors = FALSE, check.names = FALSE)
recommendation <- recommendation[
  recommendation$model == "Intercept_Mallows_GP" &
    recommendation$model_role == "primary_prediction_and_BO" &
    recommendation$recommendation_type == "incumbent_observed" &
    recommendation$bo_step %in% 0:6 & recommendation$method %in% methods,
  , drop = FALSE
]
key <- c("design_seed_id", "method", "map_id", "heldout_fold")
counts <- stats::aggregate(recommendation$bo_step, recommendation[key], length)
if (!nrow(counts) || any(counts$x != 7L) ||
    any(recommendation$n_observed != 12L + recommendation$bo_step)) {
  stopf("Every paper-method path must contain the seven checkpoints n=12,...,18")
}

auc_path <- stats::aggregate(recommendation$heldout_regret,
                             recommendation[key], sum)
names(auc_path)[ncol(auc_path)] <- "cumulative_regret"
auc <- stats::aggregate(auc_path$cumulative_regret,
                        list(method = auc_path$method), mean)
names(auc)[2L] <- "cumulative_regret"

n16 <- recommendation[recommendation$n_observed == 16L, , drop = FALSE]
top1 <- stats::aggregate(as.numeric(n16$top1_hit),
                         list(method = n16$method), mean)
names(top1)[2L] <- "top_one_at_16"

out <- Reduce(function(x, y) merge(x, y, by = "method", all = TRUE,
                                   sort = FALSE),
              list(prediction, auc, top1))
out <- out[match(methods, out$method), , drop = FALSE]
if (anyNA(out) || !identical(out$method, methods)) {
  stopf("Paper table assembly failed its method/cardinality audit")
}
out$display_method <- unname(labels[out$method])
out <- out[c("display_method", "method", "nrmse_common_mean",
             "cumulative_regret", "top_one_at_16")]

output_dir <- file.path(run_dir, "paper_sources")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "section5_3_main_table.csv")
utils::write.csv(out, output_path, row.names = FALSE)
message("Wrote ", output_path)
