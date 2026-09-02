# Hamming and component-position L2 distances used by the paper comparisons.
wcrit_row_keys <- function(D) {
  apply(as.matrix(D), 1, paste, collapse = ",")
}

wcrit_permutation_hamming <- function(a, b) {
  sum(as.integer(a) != as.integer(b))
}

wcrit_permutation_l2 <- function(a, b) {
  sqrt(sum((as.integer(a) - as.integer(b))^2))
}

wcrit_distance_matrix <- function(D, criterion = c("hamming", "l2")) {
  criterion <- match.arg(criterion)
  D <- as.matrix(D)
  storage.mode(D) <- "integer"
  has_cpp <- tryCatch({
    wcrit_load_cpp()
    exists("permutation_dist_dmat_cpp", mode = "function", inherits = TRUE)
  }, error = function(e) FALSE)
  if (isTRUE(has_cpp)) {
    return(permutation_dist_dmat_cpp(D, criterion = criterion))
  }
  n <- nrow(D)
  out <- matrix(0, nrow = n, ncol = n)
  if (n <= 1L) return(out)
  dist_fn <- switch(
    criterion,
    hamming = wcrit_permutation_hamming,
    l2 = wcrit_permutation_l2
  )
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      dij <- dist_fn(D[i, ], D[j, ])
      out[i, j] <- dij
      out[j, i] <- dij
    }
  }
  out
}
