#!/usr/bin/env Rscript

# Recreate the square Section 5.4 BO figure from the frozen paper-only CSV.

options(stringsAsFactors = FALSE)
.args <- commandArgs(trailingOnly = FALSE)
.idx <- grep("^--file=", .args)
.this_dir <- if (length(.idx)) {
  dirname(normalizePath(sub("^--file=", "", .args[.idx[[1L]]]), mustWork = TRUE))
} else normalizePath(getwd(), mustWork = TRUE)
.repro_root <- normalizePath(file.path(.this_dir, "..", ".."), mustWork = TRUE)
source(file.path(.repro_root, "code", "common", "paper_plot_style.R"), local = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
input <- path.expand(Sys.getenv(
  "SEC54_CURVE",
  unset = file.path(.repro_root, "data", "frozen", "section5_4",
                    "section5_4_pcb_native_core_bo_curve.csv")
))
output <- path.expand(Sys.getenv(
  "SEC54_FIGURE_DIR", unset = file.path(.repro_root, "outputs", "section5_4_figure")
))
dir.create(output, recursive = TRUE, showWarnings = FALSE)

curve <- utils::read.csv(input, check.names = FALSE)
method_ids <- c("fsa_lambda05", "unrestricted_hamming", "unrestricted_position_l2", "srs")
labels <- c(
  fsa_lambda05 = "FSA-KD", unrestricted_hamming = "Hamming",
  unrestricted_position_l2 = "L2", srs = "SRS"
)
if (!identical(unique(curve$scenario), "g100_w050") ||
    !setequal(unique(curve$method), method_ids) || nrow(curve) != 164L ||
    !identical(sort(unique(curve$step)), 0:40)) {
  stop("Input is not the frozen one-scenario/four-method/41-checkpoint paper curve")
}
curve$method_label_paper <- factor(unname(labels[curve$method]), levels = unname(labels))
if (!"total_evaluated" %in% names(curve)) curve$total_evaluated <- 20L + curve$step
curve <- curve[order(match(curve$method, method_ids), curve$step), ]

p <- ggplot2::ggplot(
  curve,
  ggplot2::aes(total_evaluated, mean, color = method_label_paper,
               linetype = method_label_paper, group = method_label_paper)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_low, ymax = ci_high, fill = method_label_paper),
    alpha = 0.06, color = NA, show.legend = FALSE
  ) +
  ggplot2::geom_line(linewidth = 0.80) +
  ggplot2::geom_point(size = 0.80, stroke = 0.18, alpha = 0.90) +
  ggplot2::scale_x_continuous(breaks = seq(20, 60, 10)) +
  ggplot2::labs(x = "Total evaluated routes", y = "Mean standardized regret") +
  paper_theme(8.8) +
  ggplot2::theme(
    legend.text = ggplot2::element_text(size = 7.6),
    legend.key.width = grid::unit(1.15, "lines")
  )
p <- paper_add_method_scales(
  p, c("fsa", "hamming", "position_l2", "srs"), unname(labels),
  aesthetics = c("color", "fill", "linetype")
) + ggplot2::guides(
  color = ggplot2::guide_legend(nrow = 1),
  linetype = ggplot2::guide_legend(nrow = 1)
)

pdf <- file.path(output, "fig6_pcb_gamma1_uniform.pdf")
png <- file.path(output, "fig6_pcb_gamma1_uniform.png")
paper_save_pdf(p, pdf, 5.2, 5.2)
ggplot2::ggsave(png, p, width = 5.2, height = 5.2, units = "in",
                dpi = 400, bg = "white", limitsize = FALSE)
message("Section 5.4 figure written to: ", normalizePath(output, mustWork = TRUE))
