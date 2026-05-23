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
