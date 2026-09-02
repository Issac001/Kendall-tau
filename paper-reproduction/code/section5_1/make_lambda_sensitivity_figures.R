#!/usr/bin/env Rscript

# Rebuild the two Appendix B.1 lambda-sensitivity figures from the archived
# Experiment-01 summary.  Odd-run panels are rendered for historical
# transparency but are not used for formal inference.

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

sec51_lambda_figure_main <- function() {
  root <- sec51_project_root()
  input <- path.expand(Sys.getenv(
    "SEC51_LAMBDA_SUMMARY",
    unset = file.path(root, "data", "section5_1", "lambda_sensitivity_summary.csv")
  ))
  raw_input <- path.expand(Sys.getenv(
    "SEC51_LAMBDA_RAW",
    unset = file.path(root, "data", "section5_1", "lambda_sensitivity_raw.csv")
  ))
  output <- path.expand(Sys.getenv(
    "SEC51_LAMBDA_FIGURE_OUTPUT",
    unset = file.path(root, "outputs", "section5_1", "figures")
  ))
  if (!file.exists(input)) sec51_stop("Missing lambda summary: %s", input)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  data <- utils::read.csv(input, check.names = FALSE)
  required <- c("m", "n", "lambda", "A_mean", "B_mean", "n_runs")
  missing <- setdiff(required, names(data))
  if (length(missing)) sec51_stop("Lambda summary is missing: %s", paste(missing, collapse = ", "))
  expected <- expand.grid(
    m = c(5L, 10L, 20L), multiplier = 1:5,
    lambda = seq(0, 1, by = 0.1), KEEP.OUT.ATTRS = FALSE
  )
  expected$n <- expected$m * expected$multiplier
  if (nrow(data) != 165L ||
      !setequal(paste(data$m, data$n, data$lambda),
                paste(expected$m, expected$n, expected$lambda)) ||
      any(data$n_runs != 20L)) {
    sec51_stop("Archived lambda-summary grid audit failed")
  }
  if (file.exists(raw_input)) {
    raw <- utils::read.csv(raw_input, check.names = FALSE)
    finite <- is.finite(raw$A) & is.finite(raw$B)
    if (nrow(raw) != 3300L || sum(finite) != 3299L) {
      sec51_stop("Archived lambda-raw completion audit failed")
    }
  }

  open_png <- function(filename, height) {
    grDevices::png(
      file.path(output, filename), width = 2100L, height = height,
      res = 220L, bg = "white"
    )
    graphics::par(
      family = "sans", las = 1, bty = "o", cex = 0.82,
      mgp = c(2.35, 0.7, 0)
    )
  }

  draw_panel <- function(m_value, multiplier, legend = FALSE,
                         title_style = c("numeric", "multiplier")) {
    title_style <- match.arg(title_style)
    panel <- data[data$m == m_value & data$n == multiplier * m_value, ]
    panel <- panel[order(panel$lambda), ]
    y_range <- range(c(panel$A_mean, panel$B_mean), finite = TRUE)
    panel_title <- if (identical(title_style, "numeric")) {
      paste0("m = ", m_value, ", n = ", multiplier * m_value)
    } else if (multiplier == 1L) {
      paste0("m = ", m_value, ", n = m")
    } else {
      paste0("m = ", m_value, ", n = ", multiplier, "m")
    }
    graphics::plot(
      panel$lambda, panel$A_mean, type = "n",
      ylim = y_range,
      xlab = expression(lambda), ylab = "Mean normalized component",
      main = panel_title
    )
    graphics::lines(panel$lambda, panel$A_mean,
                    col = "#D55E5E", lty = 1, lwd = 1.7)
    graphics::points(panel$lambda, panel$A_mean,
                     col = "#D55E5E", pch = 16, cex = 0.8)
    graphics::lines(panel$lambda, panel$B_mean,
                    col = "#00A7B5", lty = 2, lwd = 1.7)
    graphics::points(panel$lambda, panel$B_mean,
                     col = "#00A7B5", pch = 5, cex = 0.8)
    if (legend) {
      graphics::legend(
        "topright", legend = c(expression(A(D)), expression(B(D))),
        col = c("#D55E5E", "#00A7B5"), lty = c(1, 2),
        pch = c(16, 5), bty = "n", cex = 0.75
      )
    }
  }

  open_png("fig1_distance_components_n5m.png", 700L)
  graphics::par(mfrow = c(1, 3), mar = c(4.0, 4.2, 2.3, 0.8))
  for (m_value in c(5L, 10L, 20L)) {
    draw_panel(m_value, 5L, legend = m_value == 5L,
               title_style = "numeric")
  }
  grDevices::dev.off()

  open_png("figA1_lambda_sensitivity_n1m_n4m.png", 2550L)
  graphics::par(mfrow = c(4, 3), mar = c(3.2, 3.85, 1.95, 0.65))
  for (multiplier in 1:4) for (m_value in c(5L, 10L, 20L)) {
    draw_panel(m_value, multiplier,
               legend = multiplier == 1L && m_value == 5L,
               title_style = "multiplier")
  }
  grDevices::dev.off()

  files <- file.path(output, c(
    "fig1_distance_components_n5m.png",
    "figA1_lambda_sensitivity_n1m_n4m.png"
  ))
  if (any(!file.exists(files)) || any(file.info(files)$size <= 0)) {
    sec51_stop("A lambda-sensitivity figure was not created")
  }
  message("Appendix B.1 figures written to: ", output)
  invisible(files)
}

if (sec51_is_main("make_lambda_sensitivity_figures.R")) {
  sec51_lambda_figure_main()
}
