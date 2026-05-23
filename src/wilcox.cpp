// wilcox_test_p_asym -- two-sample Wilcoxon rank-sum test, asymptotic
// approximation with continuity correction, returning the two-sided p-value.
//
// Ties are handled by average ranks; the variance estimate includes the
// standard correction for tied groups. The continuity correction shifts the
// test statistic toward the null mean by 0.5 before standardising. The
// resulting z is mapped to a two-sided p via the standard normal CDF and
// capped at 1.
//
// Useful inside tight loops that compare a large simulated null against a
// small observed sample (the typical TEATIME peak-test workload), where the
// fixed overhead of a generic R wilcox.test wrapper would dominate.
//
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>

// beta_reassign_core(probs, vaf_index, unique_vaf, freq, n_clusters)
// Computes the (vaf, cluster, freq) tidy table that beta_reassign returns:
//   1. group rows of `probs` by `vaf_index` and sum them per group
//   2. divide each group's row by its row sum and multiply by `freq[group]`
//   3. round the result to integer allocations
//   4. expand each non-zero allocation into `alloc` rows with
//      (vaf = unique_vaf[group], cluster = column, freq = alloc)
// Iteration order is row-major (group ascending, then column ascending) so the
// output matches the R loop's row order bit-identically.
//
// [[Rcpp::export]]
Rcpp::List beta_reassign_core(Rcpp::NumericMatrix probs,
                              Rcpp::IntegerVector vaf_index,
                              Rcpp::NumericVector unique_vaf,
                              Rcpp::IntegerVector freq) {
  const int N = probs.nrow();
  const int K = probs.ncol();
  const int G = unique_vaf.size();
  // Group-sum probs by vaf_index (1-based group ids).
  Rcpp::NumericMatrix prob_sums(G, K);
  for (int i = 0; i < N; ++i) {
    const int g = vaf_index[i] - 1;
    for (int j = 0; j < K; ++j) prob_sums(g, j) += probs(i, j);
  }
  // Per-group totals.
  std::vector<double> totals(G, 0.0);
  for (int g = 0; g < G; ++g)
    for (int j = 0; j < K; ++j) totals[g] += prob_sums(g, j);
  // Allocation: round(prob_sums / totals * freq[group]).
  Rcpp::IntegerMatrix alloc(G, K);
  for (int g = 0; g < G; ++g) {
    if (totals[g] <= 0.0 || !std::isfinite(totals[g])) continue;
    const double t = totals[g];
    const double fg = static_cast<double>(freq[g]);
    for (int j = 0; j < K; ++j) {
      const double v = prob_sums(g, j) / t * fg;
      if (!std::isfinite(v)) continue;
      // R's round() is banker's rounding to even on half values; std::round is
      // away-from-zero. Replicate R's by using nearbyint with FE_TONEAREST.
      alloc(g, j) = static_cast<int>(std::nearbyint(v));
    }
  }
  // Total output rows = sum of alloc.
  int total_rows = 0;
  for (int g = 0; g < G; ++g)
    for (int j = 0; j < K; ++j)
      if (alloc(g, j) > 0) total_rows += alloc(g, j);
  Rcpp::NumericVector out_vaf(total_rows);
  Rcpp::IntegerVector out_cluster(total_rows);
  Rcpp::IntegerVector out_freq(total_rows);
  int p = 0;
  for (int g = 0; g < G; ++g) {
    for (int j = 0; j < K; ++j) {
      const int n = alloc(g, j);
      if (n <= 0) continue;
      const double v = unique_vaf[g];
      const int cluster_id = j + 1;  // R is 1-based
      for (int r = 0; r < n; ++r) {
        out_vaf[p]     = v;
        out_cluster[p] = cluster_id;
        out_freq[p]    = n;
        ++p;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::_["vaf"]     = out_vaf,
    Rcpp::_["cluster"] = out_cluster,
    Rcpp::_["freq"]    = out_freq
  );
}

// dbeta_matrix_cpp(x, a, b) returns a length(x) x length(a) matrix with
// column j = stats::dbeta(x, a[j], b[j]). The tight C++ double loop calls
// R's `dbeta` density directly per (i, j), skipping the R-level allocation of
// `rep()` arguments that the vectorised R wrapper needs.
//
// [[Rcpp::export]]
Rcpp::NumericMatrix dbeta_matrix_cpp(Rcpp::NumericVector x,
                                     Rcpp::NumericVector a,
                                     Rcpp::NumericVector b) {
  const int N = x.size();
  const int K = a.size();
  Rcpp::NumericMatrix out(N, K);
  for (int j = 0; j < K; ++j) {
    const double aj = a[j];
    const double bj = b[j];
    for (int i = 0; i < N; ++i) {
      out(i, j) = R::dbeta(x[i], aj, bj, /*give_log=*/0);
    }
  }
  return out;
}

// [[Rcpp::export]]
double wilcox_test_p_asym(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  const int nx = x.size();
  const int ny = y.size();
  const int n  = nx + ny;
  if (nx == 0 || ny == 0) return NA_REAL;

  // Combined vector with origin tracking (0 = from x, 1 = from y).
  std::vector<double> v(n);
  std::vector<int>    src(n);
  for (int i = 0; i < nx; ++i) { v[i] = x[i];      src[i] = 0; }
  for (int j = 0; j < ny; ++j) { v[nx + j] = y[j]; src[nx + j] = 1; }

  // Permutation that sorts v.
  std::vector<int> ord(n);
  for (int i = 0; i < n; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return v[a] < v[b]; });

  // Average ranks for ties; collect tied-group sizes for the variance correction.
  std::vector<double> ranks(n);
  std::vector<int>    ties;
  int i = 0;
  while (i < n) {
    int j = i;
    while (j + 1 < n && v[ord[j + 1]] == v[ord[i]]) ++j;
    // Ranks are 1-based in R; values [i..j] (0-based positions in sorted order)
    // get the average of (i+1)..(j+1) which equals ((i+1)+(j+1))/2 = (i+j+2)/2.
    const double avg = (i + j + 2) / 2.0;
    for (int k = i; k <= j; ++k) ranks[ord[k]] = avg;
    if (j > i) ties.push_back(j - i + 1);
    i = j + 1;
  }

  // STATISTIC = sum(ranks of x) - nx*(nx+1)/2
  double sum_rx = 0.0;
  for (int k = 0; k < nx; ++k) sum_rx += ranks[k];
  const double statistic = sum_rx - static_cast<double>(nx) * (nx + 1.0) / 2.0;

  // SIGMA -- ties correction term in the variance.
  double tie_sum = 0.0;
  for (int t : ties) {
    const double td = static_cast<double>(t);
    tie_sum += td * td * td - td;
  }
  const double denom = static_cast<double>(n) * (n - 1.0);
  const double sigma_sq = (static_cast<double>(nx) * ny / 12.0) *
                          ((nx + ny + 1.0) - tie_sum / denom);
  if (!(sigma_sq > 0.0)) return NA_REAL;  // degenerate (all ties)
  const double sigma = std::sqrt(sigma_sq);

  // Continuity-corrected z (matches stats::wilcox.test alternative="two.sided")
  double z = statistic - static_cast<double>(nx) * ny / 2.0;
  const double correction = (z > 0.0) ? 0.5 : ((z < 0.0) ? -0.5 : 0.0);
  z = (z - correction) / sigma;

  // Two-sided p-value via the standard normal.
  const double p_lower = R::pnorm(z, 0.0, 1.0, 1, 0);
  const double p_upper = R::pnorm(z, 0.0, 1.0, 0, 0);
  double p = 2.0 * std::min(p_lower, p_upper);
  if (p > 1.0) p = 1.0;
  return p;
}
