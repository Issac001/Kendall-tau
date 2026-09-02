# Shared publication style for the revised manuscript figures.
#
# This deliberately follows the visual language of the original manuscript:
# white panels, black boxes, no grid, sans-serif text, compact legends, and
# redundant color/line/shape encodings.  Figure numbers, long titles,
# subtitles, and captions belong in LaTeX, not inside the graphics.

paper_palette <- c(
  primary = "#000000",
  coral = "#D55E5E",
  blue = "#0072B2",
  green = "#33AA33",
  purple = "#AA77CC",
  cyan = "#00A7B5",
  orange = "#E69F00",
  grey = "#666666"
)

paper_method_style <- data.frame(
  key = c("fsa", "hamming", "position_l2", "srs"),
  label = c("FSA-KD (.5)", "Hamming", "Position-L2", "SRS"),
  color = unname(paper_palette[c("primary", "coral", "purple", "cyan")]),
  linetype = c("solid", "dashed", "dotted", "dashed"),
  shape = c(16, 2, 4, 5),
  stringsAsFactors = FALSE
)

paper_style_values <- function(keys, field) {
  stopifnot(field %in% names(paper_method_style))
  rows <- match(keys, paper_method_style$key)
  if (anyNA(rows)) {
    stop("Unknown paper style key(s): ", paste(keys[is.na(rows)], collapse = ", "))
  }
  values <- paper_method_style[[field]][rows]
  stats::setNames(values, paper_method_style$label[rows])
}

paper_theme <- function(base_size = 9, legend_position = "bottom") {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA, color = "black", linewidth = 0.45
      ),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.35),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.02), hjust = 0
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", color = "black"),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.width = grid::unit(1.35, "lines"),
      legend.spacing.x = grid::unit(0.25, "lines"),
      plot.margin = ggplot2::margin(4, 5, 4, 5)
    )
}

paper_heatmap_theme <- function(base_size = 8.5, legend_position = "right") {
  paper_theme(base_size = base_size, legend_position = legend_position) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA, color = "black", linewidth = 0.45
      ),
      legend.key.height = grid::unit(1.35, "lines")
    )
}

paper_add_method_scales <- function(
  plot, keys, labels = NULL,
  aesthetics = c("color", "fill", "linetype", "shape")
) {
  if (is.null(labels)) labels <- paper_method_style$label[match(keys, paper_method_style$key)]
  colors <- stats::setNames(paper_style_values(keys, "color"), labels)
  linetypes <- stats::setNames(paper_style_values(keys, "linetype"), labels)
  shapes <- stats::setNames(paper_style_values(keys, "shape"), labels)
  unknown <- setdiff(aesthetics, c("color", "fill", "linetype", "shape"))
  if (length(unknown)) stop("Unknown aesthetics: ", paste(unknown, collapse = ", "))
  out <- plot
  if ("color" %in% aesthetics) {
    out <- out + ggplot2::scale_color_manual(values = colors, breaks = labels, drop = FALSE)
  }
  if ("fill" %in% aesthetics) {
    out <- out + ggplot2::scale_fill_manual(values = colors, breaks = labels, drop = FALSE)
  }
  if ("linetype" %in% aesthetics) {
    out <- out + ggplot2::scale_linetype_manual(values = linetypes, breaks = labels, drop = FALSE)
  }
  if ("shape" %in% aesthetics) {
    out <- out + ggplot2::scale_shape_manual(values = shapes, breaks = labels, drop = FALSE)
  }
  out
}

paper_save_pdf <- function(plot, filename, width, height) {
  if (!capabilities("cairo")) stop("R was built without Cairo PDF support")
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
  info <- file.info(filename)
  if (!file.exists(filename) || is.na(info$size) || info$size <= 0) {
    stop("Failed to create non-empty PDF: ", filename)
  }
  invisible(normalizePath(filename, winslash = "/", mustWork = TRUE))
}
