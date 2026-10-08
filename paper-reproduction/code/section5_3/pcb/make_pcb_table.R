#!/usr/bin/env Rscript

# Recreate the compact PCB Table 3 source from the frozen four-method summary.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1L]]),
                             winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_file), "..", "..", ".."),
                      winslash = "/", mustWork = TRUE)

input_raw <- Sys.getenv(
  "SEC53_PCB_TABLE_INPUT",
  unset = file.path(root, "data", "frozen", "section5_3", "pcb",
                    "section5_3_pcb_native_core_main_table.csv")
)
output_raw <- Sys.getenv(
  "SEC53_PCB_TABLE_OUTPUT",
  unset = file.path(root, "outputs", "section5_3_pcb_table", "table3_pcb.csv")
)
input <- path.expand(input_raw)
output <- path.expand(output_raw)

x <- utils::read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
method_ids <- c("fsa_lambda05", "unrestricted_hamming",
                "unrestricted_position_l2", "srs")
x <- x[match(method_ids, x$method), , drop = FALSE]
if (nrow(x) != 4L || anyNA(x$method) || !identical(x$method, method_ids)) {
  stop("Input is not the frozen four-method PCB summary", call. = FALSE)
}

out <- data.frame(
  method = x$method_label,
  cumulative_regret = x$standardized_regret_sum_auc,
  final_regret = x$final_standardized_regret,
  top_0p1_hit_rate = x$final_top_0p1pct_hit_rate,
  stringsAsFactors = FALSE
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, output, row.names = FALSE)
message("Wrote ", normalizePath(output, mustWork = TRUE))
