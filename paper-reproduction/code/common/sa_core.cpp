// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// Kendall distance (raw inversion count)
// ---------------------------------------------------------------------------
long long count_inversions(IntegerVector arr) {
  int n = arr.size();
  if (n <= 1) return 0;

  int mid = n / 2;
  IntegerVector left = head(arr, mid);
  IntegerVector right = tail(arr, n - mid);

  long long inv = count_inversions(left) + count_inversions(right);

  int i = 0, j = 0, k = 0;
  int n_left = left.size(), n_right = right.size();

  while (i < n_left && j < n_right) {
    if (left[i] <= right[j]) {
      arr[k++] = left[i++];
    } else {
      arr[k++] = right[j++];
      inv += (n_left - i);
    }
  }
  while (i < n_left) arr[k++] = left[i++];
  while (j < n_right) arr[k++] = right[j++];

  return inv;
}

int kendall_dist(const IntegerVector& a, const IntegerVector& b) {
  int n = a.size();
  IntegerVector indexes_b(n + 1);
  for (int i = 0; i < n; i++) {
    indexes_b[b[i]] = i + 1;
  }

  IntegerVector rank_a(n);
  for (int i = 0; i < n; i++) {
    rank_a[i] = indexes_b[a[i]];
  }

  return (int)count_inversions(clone(rank_a));
}

// [[Rcpp::export]]
IntegerMatrix kendall_dmat_between_cpp(const IntegerMatrix& D1, const IntegerMatrix& D2) {
  int n1 = D1.nrow();
  int n2 = D2.nrow();
  IntegerMatrix out(n1, n2);

  for (int i = 0; i < n1; i++) {
    IntegerVector a = D1(i, _);
    for (int j = 0; j < n2; j++) {
      IntegerVector b = D2(j, _);
      out(i, j) = kendall_dist(a, b);
    }
  }

  return out;
}

static void rebuild_dmat_H(const IntegerMatrix& H, int K, int m, imat& dmat) {
  for (int i = 0; i < K; i++) {
    for (int j = i + 1; j < K; j++) {
      IntegerVector a = H(i, _);
      IntegerVector b = H(j, _);
      int d = kendall_dist(a, b);
      dmat(i, j) = d;
      dmat(j, i) = d;
    }
  }
}

static bool row_equals_vec(const IntegerMatrix& D, int row, const IntegerVector& x, int m) {
  for (int j = 0; j < m; j++) {
    if (D(row, j) != x[j]) return false;
  }
  return true;
}

static bool row_exists_except(const IntegerMatrix& D, const IntegerVector& x, int skip_row, int m) {
  int n = D.nrow();
  for (int i = 0; i < n; i++) {
    if (i == skip_row) continue;
    if (row_equals_vec(D, i, x, m)) return true;
  }
  return false;
}

static IntegerVector reverse_perm_cpp(const IntegerVector& x) {
  int m = x.size();
  IntegerVector out(m);
  for (int j = 0; j < m; j++) out[j] = x[m - 1 - j];
  return out;
}

static IntegerMatrix build_foldover_D_cpp(const IntegerMatrix& H) {
  int K = H.nrow();
  int m = H.ncol();
  IntegerMatrix D(2 * K, m);
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < m; j++) {
      D(i, j) = H(i, j);
      D(K + i, j) = H(i, m - 1 - j);
    }
  }
  return D;
}

// ---------------------------------------------------------------------------
// FSA-KD on K reversal-pair representatives for the weighted criterion
// Phi_lambda = lambda * k_min / C1 + (1 - lambda) * (U2 - k_m2) / (U2 - C2)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List sa_optimize_H_weighted_foldover_cpp(
    const IntegerMatrix& H_init,
    int m,
    double lambda,
    double T0,
    double alpha,
    int max_iter,
    int no_improve_stop,
    bool duplicate_free,
    bool verbose,
    bool incremental,
    double max_seconds,
    unsigned int rng_seed,
    Nullable<NumericVector> trace
) {
  int K = H_init.nrow();
  int n = 2 * K;
  int p = m * (m - 1) / 2;
  IntegerMatrix H = clone(H_init);
  IntegerMatrix H_best = clone(H_init);

  imat dmat(K, K, fill::zeros);
  rebuild_dmat_H(H, K, m, dmat);

  std::vector<std::vector<int>> pos_maps;
  if (incremental) {
    pos_maps.resize(K, std::vector<int>(m + 1));
    for (int i = 0; i < K; i++) {
      for (int j = 0; j < m; j++) {
        pos_maps[i][H(i, j)] = j;
      }
    }
  }

  const double C1 = std::floor((double)m * (double)(m - 1) / 4.0);
  const double C2 = ((double)n * (double)m *
                     (9.0 * std::pow((double)m, 3.0) -
                      14.0 * std::pow((double)m, 2.0) +
                      15.0 * (double)m - 10.0)) /
                    (144.0 * (double)(n - 1));
  const double U2 = (((double)n * std::pow((double)m, 2.0) *
                      std::pow((double)(m - 1), 2.0)) -
                     (4.0 * (double)(n - 2) *
                      ((double)m * (double)(m - 1) - 2.0))) /
                    (8.0 * (double)(n - 1));
  const double denomB = U2 - C2;

  struct WeightedMetrics {
    double obj;
    double minK;
    double k_m2;
    double A;
    double B;
  };

  auto compute_weighted = [&](const imat& dm) -> WeightedMetrics {
    double minK_H = (K <= 1) ? (double)p : 1e100;
    double maxK_H = -1e100;
    double sum_d = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < K; i++) {
      for (int j = i + 1; j < K; j++) {
        int d = dm(i, j);
        if (d < minK_H) minK_H = (double)d;
        if (d > maxK_H) maxK_H = (double)d;
        sum_d += (double)d;
        sum_sq += (double)d * (double)d;
      }
    }
    double minK;
    if (K <= 1) {
      minK = (double)p;
    } else {
      minK = std::min(minK_H, (double)p - maxK_H);
    }
    double Kd = (double)K;
    double pd = (double)p;
    double k_m2 = (4.0 * sum_sq + Kd * Kd * pd * pd - 4.0 * pd * sum_d) / (Kd * (2.0 * Kd - 1.0));
    double A = (C1 > 1.0) ? (minK - 1.0) / (C1 - 1.0) : 1.0;
    double B = (std::abs(denomB) > 1e-12) ? (U2 - k_m2) / denomB : 0.0;
    double obj = lambda * A + (1.0 - lambda) * B;
    return {obj, minK, k_m2, A, B};
  };

  std::mt19937 gen(rng_seed == 0u ? std::random_device{}() : static_cast<std::mt19937::result_type>(rng_seed));
  std::uniform_real_distribution<> unif(0.0, 1.0);
  std::uniform_int_distribution<> unif_int(0, K - 1);
  std::uniform_int_distribution<> unif_m(0, m - 1);

  WeightedMetrics cur = compute_weighted(dmat);
  WeightedMetrics best = cur;
  double T = T0;
  int no_improve = 0;

  NumericVector tr;
  bool do_trace = trace.isNotNull();
  if (do_trace) {
    tr = NumericVector(max_iter, NA_REAL);
  }

  auto t0 = std::chrono::steady_clock::now();
  int it_done = 0;

  for (int it = 0; it < max_iter; it++) {
    if (max_seconds > 0) {
      auto now = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed = now - t0;
      if (elapsed.count() >= max_seconds) break;
    }

    int idx = unif_int(gen);

    IntegerVector old_row = H(idx, _);
    IntegerVector new_row;
    bool is_local = false;
    int pos1 = 0, pos2 = 0;

    if (unif(gen) < T / T0) {
      new_row = IntegerVector(m);
      for (int j = 0; j < m; j++) new_row[j] = j + 1;
      std::shuffle(new_row.begin(), new_row.end(), gen);
    } else {
      is_local = true;
      pos1 = unif_m(gen);
      pos2 = unif_m(gen);
      while (pos1 == pos2) pos2 = unif_m(gen);
      new_row = clone(old_row);
      std::swap(new_row[pos1], new_row[pos2]);
    }

    if (duplicate_free) {
      if (row_exists_except(H, new_row, idx, m)) {
        no_improve++;
        T *= alpha;
        it_done = it + 1;
        if (do_trace) tr[it] = best.obj;
        if (no_improve > no_improve_stop || T < 1e-8) break;
        continue;
      }
      IntegerVector new_fold = reverse_perm_cpp(new_row);
      if (row_exists_except(H, new_fold, idx, m)) {
        no_improve++;
        T *= alpha;
        it_done = it + 1;
        if (do_trace) tr[it] = best.obj;
        if (no_improve > no_improve_stop || T < 1e-8) break;
        continue;
      }
    }

    ivec old_dists;
    if (incremental) {
      old_dists = dmat.col(idx);
      if (is_local) {
        int s = std::min(pos1, pos2);
        int t = std::max(pos1, pos2);
        int u = old_row[s];
        int v = old_row[t];
        for (int j = 0; j < K; j++) {
          if (j == idx) continue;
          int delta_dist = 0;
          delta_dist += (pos_maps[j][u] < pos_maps[j][v]) ? 1 : -1;
          for (int l = s + 1; l < t; l++) {
            int w = old_row[l];
            delta_dist += (pos_maps[j][u] < pos_maps[j][w]) ? 1 : -1;
            delta_dist -= (pos_maps[j][v] < pos_maps[j][w]) ? 1 : -1;
          }
          dmat(idx, j) = old_dists[j] + delta_dist;
          dmat(j, idx) = dmat(idx, j);
        }
      } else {
        for (int j = 0; j < K; j++) {
          if (j == idx) continue;
          IntegerVector b = H(j, _);
          int d = kendall_dist(new_row, b);
          dmat(idx, j) = d;
          dmat(j, idx) = d;
        }
      }
      H(idx, _) = new_row;
      for (int j = 0; j < m; j++) {
        pos_maps[idx][new_row[j]] = j;
      }
    } else {
      H(idx, _) = new_row;
      rebuild_dmat_H(H, K, m, dmat);
    }

    WeightedMetrics cand = compute_weighted(dmat);
    double delta = cand.obj - cur.obj;

    bool accept = false;
    if (delta >= 0) {
      accept = true;
    } else {
      double prob = std::exp(delta / T);
      if (unif(gen) < prob) accept = true;
    }

    if (accept) {
      cur = cand;
      if (cur.obj > best.obj) {
        best = cur;
        H_best = clone(H);
        no_improve = 0;
      }
    } else {
      if (incremental) {
        H(idx, _) = old_row;
        dmat.col(idx) = old_dists;
        dmat.row(idx) = old_dists.t();
        for (int j = 0; j < m; j++) {
          pos_maps[idx][old_row[j]] = j;
        }
      } else {
        H(idx, _) = old_row;
        rebuild_dmat_H(H, K, m, dmat);
      }
      no_improve++;
    }

    T *= alpha;
    it_done = it + 1;
    if (do_trace) tr[it] = best.obj;
    if (no_improve > no_improve_stop || T < 1e-8) break;
  }

  IntegerMatrix D_best = build_foldover_D_cpp(H_best);

  List out = List::create(
      Named("H_best") = H_best,
      Named("D_best") = D_best,
      Named("best_obj") = best.obj,
      Named("best_minK") = best.minK,
      Named("best_k_m2") = best.k_m2,
      Named("best_A") = best.A,
      Named("best_B") = best.B,
      Named("n_iter") = it_done);

  if (do_trace) {
    out["trace"] = tr;
  }

  return out;
}

// ---------------------------------------------------------------------------
// Fast PWO model matrix and paper-baseline distance matrices
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix pwo_matrix_cpp(const IntegerMatrix& D, bool intercept = true) {
  int n = D.nrow();
  int m = D.ncol();
  int p = m * (m - 1) / 2;
  int offset = intercept ? 1 : 0;
  NumericMatrix X(n, p + offset);

  if (intercept) {
    for (int i = 0; i < n; i++) X(i, 0) = 1.0;
  }

  IntegerVector pos(m + 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      int item = D(i, j);
      if (item >= 1 && item <= m) pos[item] = j;
    }
    int col = offset;
    for (int a = 1; a <= m - 1; a++) {
      for (int b = a + 1; b <= m; b++) {
        X(i, col) = (pos[a] < pos[b]) ? 1.0 : -1.0;
        col++;
      }
    }
  }

  return X;
}

// Pairwise distances needed by the Hamming and component-position L2 figures.
static double permutation_dist_cpp_rows(const IntegerMatrix& D, int row1, int row2,
                                        int criterion_id) {
  int m = D.ncol();
  double d = 0.0;
  for (int j = 0; j < m; j++) {
    if (criterion_id == 1) {
      d += (D(row1, j) != D(row2, j));
    } else {
      double diff = (double)D(row1, j) - (double)D(row2, j);
      d += diff * diff;
    }
  }
  return (criterion_id == 2) ? std::sqrt(d) : d;
}

static int criterion_to_id(const std::string& criterion) {
  if (criterion == "hamming") return 1;
  if (criterion == "l2") return 2;
  stop("Unknown maximin distance criterion: " + criterion);
}

// [[Rcpp::export]]
NumericMatrix permutation_dist_dmat_cpp(const IntegerMatrix& D,
                                        std::string criterion = "hamming") {
  int n = D.nrow();
  int criterion_id = criterion_to_id(criterion);
  NumericMatrix out(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double d = permutation_dist_cpp_rows(D, i, j, criterion_id);
      out(i, j) = d;
      out(j, i) = d;
    }
  }
  return out;
}
