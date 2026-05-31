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

// ks_stat_2s_cpp -- two-sample Kolmogorov-Smirnov statistic D (two.sided),
// bit-identical to stats::ks.test's internal computation:
//
//   w <- c(x, y)
//   z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
//   if (ties) z <- z[c(which(diff(sort(w)) != 0), N)]   # collapse to run-ends
//   D <- max(abs(z))
//
// The floating-point accumulation order is preserved exactly: elements are
// visited in R's order(w) order (ascending value, ties broken by original
// index, so x-origin elements -- indices 0..n.x-1 -- precede y-origin ones for
// equal values). The running sum adds +1/n.x or -1/n.y one term at a time, and
// |sum| is recorded only at the end of each distinct-value run (which is what
// the tie-collapse selects; with no ties every position is its own run-end, so
// max(abs(.)) over run-ends equals max(abs(z)) either way).
//
// The (cheap, scalar) asymptotic p-value is left to R's own stats::psmirnov so
// it stays bit-identical without re-deriving C_pkolmogorov_two_limit here.
//
// [[Rcpp::export]]
double ks_stat_2s_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  const int nx = x.size();
  const int ny = y.size();
  const int n  = nx + ny;
  if (nx == 0 || ny == 0) return NA_REAL;

  // Combined vector; original indices 0..nx-1 are x, nx..n-1 are y.
  std::vector<double> v(n);
  for (int i = 0; i < nx; ++i) v[i]      = x[i];
  for (int j = 0; j < ny; ++j) v[nx + j] = y[j];

  // order(w): ascending value, ties by ascending original index. Including the
  // index as the tiebreak reproduces R's stable radix order (x before y).
  std::vector<int> ord(n);
  for (int i = 0; i < n; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return v[a] < v[b] || (v[a] == v[b] && a < b); });

  const double inc_x =  1.0 / nx;   // matches R's 1/n.x (double division)
  const double inc_y = -1.0 / ny;   // matches R's -1/n.y

  double s = 0.0;
  double D = 0.0;
  for (int k = 0; k < n; ++k) {
    s += (ord[k] < nx) ? inc_x : inc_y;   // ord[k] < nx  <=>  order(w) <= n.x
    if (k == n - 1 || v[ord[k + 1]] != v[ord[k]]) {   // end of a distinct-value run
      const double as = std::fabs(s);
      if (as > D) D = as;
    }
  }
  return D;
}

// ks_d_batch_cpp -- the two-sample KS statistic D for many simulated sets that
// share a FIXED first sample `x`. `sims` holds `ntry` consecutive segments,
// each of length ny = length(sims)/ntry, playing the role of `y` in
// ks_stat_2s_cpp. `x` is sorted once and reused; each segment is sorted and
// merged against it (x-before-y on ties, matching c(x,y)'s original-index
// order), accumulating the same +1/nx, -1/ny running sum and recording |sum|
// at distinct-value run-ends. The per-try D is therefore bit-identical to
// ks_stat_2s_cpp(x, segment); only the O(nx log nx) sort of x is amortised
// across the ntry tries (the dominant cost in simulate_from_estimation's loop
// when x = compare_set is large and fixed).
//
// [[Rcpp::export]]
Rcpp::NumericVector ks_d_batch_cpp(Rcpp::NumericVector x, Rcpp::NumericVector sims, int ntry) {
  const int nx = x.size();
  const int total = sims.size();
  if (nx == 0 || ntry <= 0 || total % ntry != 0) return Rcpp::NumericVector(ntry, NA_REAL);
  const int ny = total / ntry;
  Rcpp::NumericVector out(ntry);
  if (ny == 0) { for (int t = 0; t < ntry; ++t) out[t] = NA_REAL; return out; }

  std::vector<double> xs(x.begin(), x.end());
  std::sort(xs.begin(), xs.end());
  const double inc_x =  1.0 / nx;
  const double inc_y = -1.0 / ny;

  std::vector<double> ys(ny);
  for (int t = 0; t < ntry; ++t) {
    const int base = t * ny;
    for (int j = 0; j < ny; ++j) ys[j] = sims[base + j];
    std::sort(ys.begin(), ys.end());

    double s = 0.0, D = 0.0;
    int i = 0, j = 0;
    while (i < nx || j < ny) {
      // next value; ties take x first (x's original indices precede y's)
      bool takeX;
      if (i < nx && j < ny) takeX = !(ys[j] < xs[i]);   // xs<=ys -> x
      else takeX = (i < nx);
      const double cur = takeX ? xs[i] : ys[j];
      if (takeX) { s += inc_x; ++i; } else { s += inc_y; ++j; }
      bool hasNext = (i < nx || j < ny);
      double nextval = 0.0;
      if (hasNext) {
        if (i < nx && j < ny)      nextval = std::min(xs[i], ys[j]);
        else if (i < nx)           nextval = xs[i];
        else                       nextval = ys[j];
      }
      if (!hasNext || nextval != cur) {
        const double as = std::fabs(s);
        if (as > D) D = as;
      }
    }
    out[t] = D;
  }
  return out;
}

// wilcox_p_batch_cpp -- two-sided asymptotic Wilcoxon rank-sum p-value for many
// simulated `y` sets that share a FIXED `x` (= compare_set). Mirrors
// wilcox_test_p_asym (statistic = sum of ranks of x; tie-corrected variance;
// continuity correction; 2*min(pnorm), capped at 1), but sorts x once and
// merges each sorted segment, grouping ties to get average ranks from counts.
// All the rank arithmetic is exact (half-integers, integer m^3 - m sums well
// under 2^53 for these n), so the per-try p is bit-identical to
// wilcox_test_p_asym(x, segment). Amortises x's sort across the ntry tries.
//
// [[Rcpp::export]]
Rcpp::NumericVector wilcox_p_batch_cpp(Rcpp::NumericVector x, Rcpp::NumericVector sims, int ntry) {
  const int nx = x.size();
  const int total = sims.size();
  if (nx == 0 || ntry <= 0 || total % ntry != 0) return Rcpp::NumericVector(ntry, NA_REAL);
  const int ny = total / ntry;
  const int n  = nx + ny;
  Rcpp::NumericVector out(ntry);
  if (ny == 0) { for (int t = 0; t < ntry; ++t) out[t] = NA_REAL; return out; }

  std::vector<double> xs(x.begin(), x.end());
  std::sort(xs.begin(), xs.end());
  const double denom = static_cast<double>(n) * (n - 1.0);

  std::vector<double> ys(ny);
  for (int t = 0; t < ntry; ++t) {
    const int base = t * ny;
    for (int j = 0; j < ny; ++j) ys[j] = sims[base + j];
    std::sort(ys.begin(), ys.end());

    double sum_rx = 0.0, tie_sum = 0.0, cum = 0.0;   // cum = # ranked before group
    int i = 0, j = 0;
    while (i < nx || j < ny) {
      double cur;
      if (i < nx && j < ny) cur = std::min(xs[i], ys[j]);
      else if (i < nx)      cur = xs[i];
      else                  cur = ys[j];
      double cx = 0.0, cy = 0.0;
      while (i < nx && xs[i] == cur) { cx += 1.0; ++i; }
      while (j < ny && ys[j] == cur) { cy += 1.0; ++j; }
      const double m = cx + cy;
      const double avg_rank = cum + (m + 1.0) / 2.0;
      sum_rx  += cx * avg_rank;
      tie_sum += m * m * m - m;
      cum     += m;
    }
    const double statistic = sum_rx - static_cast<double>(nx) * (nx + 1.0) / 2.0;
    const double sigma_sq  = (static_cast<double>(nx) * ny / 12.0) *
                             ((nx + ny + 1.0) - tie_sum / denom);
    if (!(sigma_sq > 0.0)) { out[t] = NA_REAL; continue; }
    double z = statistic - static_cast<double>(nx) * ny / 2.0;
    const double corr = (z > 0.0) ? 0.5 : ((z < 0.0) ? -0.5 : 0.0);
    z = (z - corr) / std::sqrt(sigma_sq);
    double p = 2.0 * std::min(R::pnorm(z, 0.0, 1.0, 1, 0), R::pnorm(z, 0.0, 1.0, 0, 0));
    if (p > 1.0) p = 1.0;
    out[t] = p;
  }
  return out;
}
