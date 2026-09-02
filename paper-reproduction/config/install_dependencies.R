packages <- c(
  "Rcpp", "RcppArmadillo", "digest", "dplyr", "tidyr", "ggplot2",
  "gtools", "jsonlite"
)

missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All reproduction dependencies are already installed.")
}
