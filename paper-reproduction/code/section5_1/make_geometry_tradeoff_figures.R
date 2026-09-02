#!/usr/bin/env Rscript

# Rebuild the Section 5.1 c=1 trade-off figure and its Appendix B c=4
# sensitivity counterpart from the four-method paper data.  Every one of the
# 15 panels uses its own x and y limits, matching the submitted figures.

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

sec51_tradeoff_main <- function() {
  root <- sec51_project_root()
  input <- path.expand(Sys.getenv(
    "SEC51_GEOMETRY_METRICS",
    unset = file.path(root, "data", "section5_1", "geometry_four_method_raw_metrics.csv")
  ))
  output <- path.expand(Sys.getenv(
    "SEC51_GEOMETRY_FIGURE_OUTPUT",
    unset = file.path(root, "outputs", "section5_1", "figures")
  ))
  if (!file.exists(input)) sec51_stop("Missing four-method input: %s", input)
  if (!capabilities("cairo")) sec51_stop("Cairo PDF support is required")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)

  raw <- utils::read.csv(input, check.names = FALSE)
  required <- c(
    "paper_label", "method", "m", "n", "ratio", "rep",
    "pwo_ms_efficiency", "mallows_det_root_c1", "mallows_det_root_c4"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing)) sec51_stop("Input is missing: %s", paste(missing, collapse = ", "))
  method_order <- c("FSA-KD", "Hamming", "L2", "SRS")
  m_order <- c(6L, 10L, 20L)
  ratio_order <- 1:5
  raw <- raw[
    raw$paper_label %in% method_order & raw$m %in% m_order &
      raw$ratio %in% ratio_order,
    required,
    drop = FALSE
  ]
  raw$method_order <- match(raw$paper_label, method_order)
  raw <- raw[order(raw$m, raw$ratio, raw$method_order, raw$rep), ]
  key <- paste(raw$paper_label, raw$m, raw$ratio, raw$rep, sep = "|")
  counts <- table(raw$paper_label, raw$m, raw$ratio)
  if (nrow(raw) != 3000L || anyDuplicated(key) || any(counts != 50L) ||
      any(raw$n != raw$m * raw$ratio) ||
      any(!is.finite(raw$pwo_ms_efficiency)) ||
      any(!is.finite(raw$mallows_det_root_c1)) ||
      any(!is.finite(raw$mallows_det_root_c4))) {
    sec51_stop("Four-method paper-data audit failed")
  }

  styles <- data.frame(
    paper_label = method_order,
    color = c("#000000", "#D55E5E", "#AA77CC", "#00A7B5"),
    pch = c(16L, 2L, 4L, 5L),
    stringsAsFactors = FALSE
  )
  draw_order <- c("SRS", "Hamming", "L2", "FSA-KD")

  padded_range <- function(x, fraction = 0.12) {
    range_x <- range(x, finite = TRUE)
    span <- diff(range_x)
    padding <- if (is.finite(span) && span > 0) {
      fraction * span
    } else {
      0.05 * max(abs(range_x), 1)
    }
    range_x + c(-padding, padding)
  }

  compact_ticks <- function(limits) {
    at <- pretty(limits, n = 2L)
    at <- at[at >= limits[[1L]] & at <= limits[[2L]]]
    if (length(at) > 3L) {
      at <- at[unique(round(seq(1, length(at), length.out = 3L)))]
    }
    if (length(at) < 2L) at <- seq(limits[[1L]], limits[[2L]], length.out = 3L)
    list(at = at, labels = formatC(at, format = "f", digits = 2L))
  }

  draw_row_strip <- function(m_value) {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::text(0.50, 0.50, paste0("m = ", m_value),
                   srt = -90, font = 2, cex = 0.68)
  }

  draw_legend <- function() {
    graphics::par(mar = c(0, 0, 0, 0), xpd = NA)
    graphics::plot.new()
    graphics::legend(
      "center", legend = styles$paper_label,
      col = styles$color, pch = styles$pch,
      pt.cex = 0.82, pt.lwd = 1.0, horiz = TRUE,
      bty = "n", cex = 0.66, x.intersp = 0.65
    )
  }

  draw_figure <- function(c_value, means) {
    y_metric <- paste0("mallows_det_root_c", c_value)
    limits <- list()
    for (m_value in m_order) for (ratio_value in ratio_order) {
      panel <- raw[raw$m == m_value & raw$ratio == ratio_value, ]
      limits[[paste(m_value, ratio_value, sep = "|")]] <- list(
        x = padded_range(panel$pwo_ms_efficiency),
        y = padded_range(panel[[y_metric]])
      )
    }
    panels <- matrix(1:15, nrow = 3L, ncol = 5L, byrow = TRUE)
    layout_matrix <- rbind(
      c(panels[1L, ], 16L), c(panels[2L, ], 17L),
      c(panels[3L, ], 18L), rep(19L, 6L)
    )
    graphics::layout(
      layout_matrix, widths = c(rep(1, 5), 0.10),
      heights = c(1, 1, 1, 0.17)
    )
    graphics::par(
      family = "sans", las = 1, bty = "o", cex = 0.60,
      mgp = c(1.8, 0.42, 0), mar = c(1.58, 1.75, 0.72, 0.28),
      oma = c(2.30, 3.15, 0.45, 0.05)
    )
    for (m_index in seq_along(m_order)) {
      m_value <- m_order[[m_index]]
      for (ratio_index in seq_along(ratio_order)) {
        ratio_value <- ratio_order[[ratio_index]]
        panel_key <- paste(m_value, ratio_value, sep = "|")
        panel_limits <- limits[[panel_key]]
        x_ticks <- compact_ticks(panel_limits$x)
        y_ticks <- compact_ticks(panel_limits$y)
        panel_raw <- raw[raw$m == m_value & raw$ratio == ratio_value, ]
        panel_mean <- means[means$m == m_value & means$ratio == ratio_value, ]
        graphics::plot(panel_limits$x, panel_limits$y, type = "n", axes = FALSE,
                       xlab = "", ylab = "")
        graphics::box(lwd = 0.75)
        graphics::axis(1, at = x_ticks$at, labels = x_ticks$labels,
                       cex.axis = 0.47, tck = -0.027, lwd = 0.62,
                       mgp = c(0, 0.34, 0))
        graphics::axis(2, at = y_ticks$at, labels = y_ticks$labels,
                       cex.axis = 0.47, tck = -0.027, lwd = 0.62,
                       mgp = c(0, 0.42, 0))
        if (m_index == 1L) {
          graphics::mtext(paste0("n/m = ", ratio_value), side = 3,
                          line = 0.10, font = 2, cex = 0.66)
        }
        for (label in draw_order) {
          style <- styles[styles$paper_label == label, ]
          points <- panel_raw[panel_raw$paper_label == label, ]
          graphics::points(
            points$pwo_ms_efficiency, points[[y_metric]],
            col = grDevices::adjustcolor(style$color[[1L]], alpha.f = 0.18),
            pch = style$pch[[1L]], cex = 0.35, lwd = 0.45
          )
        }
        for (label in draw_order) {
          style <- styles[styles$paper_label == label, ]
          point <- panel_mean[panel_mean$paper_label == label, ]
          graphics::points(
            point$pwo_ms_efficiency, point[[y_metric]],
            col = style$color[[1L]], pch = style$pch[[1L]],
            cex = 0.90, lwd = 1.12
          )
        }
      }
    }
    for (m_value in m_order) draw_row_strip(m_value)
    draw_legend()
    graphics::mtext("PWO absolute MS efficiency", side = 1, outer = TRUE,
                    line = 1.12, cex = 0.69)
    graphics::mtext(
      sprintf("Mallows determinant-root (c = %d)", c_value),
      side = 2, outer = TRUE, line = 1.68, cex = 0.69, las = 0
    )
  }

  means_rows <- list()
  for (c_value in c(1L, 4L)) {
    y_metric <- paste0("mallows_det_root_c", c_value)
    groups <- split(
      raw,
      interaction(raw$paper_label, raw$method, raw$m, raw$n, raw$ratio,
                  drop = TRUE, lex.order = TRUE)
    )
    means <- do.call(rbind, lapply(groups, function(data) {
      data.frame(
        paper_label = data$paper_label[[1L]], method = data$method[[1L]],
        m = data$m[[1L]], n = data$n[[1L]], ratio = data$ratio[[1L]],
        n_reps = nrow(data),
        pwo_ms_efficiency = mean(data$pwo_ms_efficiency),
        value = mean(data[[y_metric]]),
        stringsAsFactors = FALSE
      )
    }))
    names(means)[names(means) == "value"] <- y_metric
    if (nrow(means) != 60L || any(means$n_reps != 50L)) {
      sec51_stop("Expected 60 means for c=%d", c_value)
    }
    means_rows[[as.character(c_value)]] <- means
    filename <- if (c_value == 1L) {
      "fig2_geometry_tradeoff_uniform.pdf"
    } else {
      "figS_geometry_tradeoff_c4_uniform.pdf"
    }
    path <- file.path(output, filename)
    grDevices::cairo_pdf(path, width = 5.6, height = 4.5,
                         family = "sans", bg = "white")
    ok <- FALSE
    tryCatch({
      draw_figure(c_value, means)
      ok <- TRUE
    }, finally = grDevices::dev.off())
    if (!ok || !file.exists(path) || file.info(path)$size <= 0) {
      sec51_stop("Failed to create %s", path)
    }
  }
  utils::write.csv(
    means_rows[["1"]], file.path(output, "geometry_tradeoff_c1_means.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    means_rows[["4"]], file.path(output, "geometry_tradeoff_c4_means.csv"),
    row.names = FALSE
  )
  message("Section 5.1 trade-off figures written to: ", output)
  invisible(output)
}

if (sec51_is_main("make_geometry_tradeoff_figures.R")) sec51_tradeoff_main()
