.ok <- function(value) {
  list(status = "ok", value = value, message = NA_character_)
}

.err <- function(cond) {
  list(status = "error", value = NULL, message = conditionMessage(cond))
}

.log <- function(ctx, step, msg, timer = NULL) {
  if (!isTRUE(ctx$verbose)) {
    return(invisible(NULL))
  }

  elapsed <- ""
  if (!is.null(timer)) {
    elapsed <- sprintf(" (%.2fs)", unname((proc.time() - timer)[["elapsed"]]))
  }

  cat(sprintf("[%s] %s%s\n", step, msg, elapsed))
  invisible(NULL)
}

# Cap a candidate sweep (cell.div / mu candidate vector) to at most
# getOption("teatime.max_candidates") evenly-spaced entries, preserving order
# and endpoints. Unset (or Inf) -> returned unchanged (no cap), which keeps
# the default/faithful path bit-identical. The fast dispatch sets the option
# (default 2000); because the mu likelihood/peak score is smooth in mu,
# subsampling the candidate grid shifts the picked value <1% while turning the
# O(n_snv) sweep into O(max_candidates).
.cap_candidates <- function(x) {
  n <- getOption("teatime.max_candidates", Inf)
  if (!is.finite(n) || length(x) <= n) {
    return(x)
  }
  idx <- unique(round(seq.int(1L, length(x), length.out = as.integer(n))))
  x[idx]
}

collect_rows <- function(fn, items) {
  n <- length(items)
  if (n == 0L) {
    return(data.frame())
  }

  out <- vector("list", n)
  for (idx in seq_len(n)) {
    out[[idx]] <- fn(items[[idx]])
  }

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0L) {
    return(data.frame())
  }

  do.call(rbind, out)
}

# dbeta_matrix(x, a, b) returns the N x K matrix with column j =
# stats::dbeta(x, a[j], b[j]). In fast mode x is deduplicated before the
# Rcpp dbeta call (depth is bounded so x rarely has more than ~depth+1
# distinct values even when length(x) reaches 10^4); rows are then expanded
# back to the original length(x) order. Output is bit-identical to the
# non-dedup version. Default mode keeps the per-column sapply form.
dbeta_matrix <- function(x, a, b) {
  if (isTRUE(getOption("teatime.fast_version", FALSE))) {
    xn <- as.numeric(x)
    nx <- length(xn)
    ux <- unique(xn)
    # Only dedup when it actually helps; below ~5x compression the match()
    # overhead breaks even, so use the direct path.
    if (nx >= 200L && length(ux) * 5L < nx) {
      small <- dbeta_matrix_cpp(ux, as.numeric(a), as.numeric(b))
      idx   <- match(xn, ux)
      small[idx, , drop = FALSE]
    } else {
      dbeta_matrix_cpp(xn, as.numeric(a), as.numeric(b))
    }
  } else {
    sapply(seq_along(a), function(i) stats::dbeta(x, a[i], b[i]))
  }
}

.as_prob_matrix <- function(df) {
  prob_cols <- grep("^prob", names(df))
  if (length(prob_cols) == 0L) {
    stop("No probability columns found.", call. = FALSE)
  }
  list(
    probs = as.matrix(df[, prob_cols, drop = FALSE]),
    prob_cols = names(df)[prob_cols]
  )
}

beta_reassign <- function(df) {
  if (!("vaf" %in% names(df))) {
    stop("beta_reassign requires a 'vaf' column.", call. = FALSE)
  }

  prob_info <- .as_prob_matrix(df)
  probs <- prob_info$probs
  vaf <- df$vaf
  unique_vaf <- sort(unique(vaf))
  vaf_index <- match(vaf, unique_vaf)
  freq <- tabulate(vaf_index, nbins = length(unique_vaf))
  prob_sums <- rowsum(probs, group = vaf_index, reorder = TRUE)
  totals <- rowSums(prob_sums)
  alloc <- round(prob_sums / totals * freq)
  alloc[!is.finite(alloc)] <- 0

  # Fast mode delegates the group-by + rounding + expansion to a single C++
  # routine (beta_reassign_core) that walks the (vaf x cluster) matrix once and
  # writes the tidy (vaf, cluster, freq) output. Row order is row-major (group
  # ascending then column ascending) to match the per-row R loop output.
  if (isTRUE(getOption("teatime.fast_version", FALSE))) {
    parts <- beta_reassign_core(probs, as.integer(vaf_index),
                                as.numeric(unique_vaf), as.integer(freq))
    final_df <- data.frame(vaf = parts$vaf, cluster = parts$cluster, freq = parts$freq)
  } else {
    expanded_vaf <- c()
    expanded_cluster <- c()
    expanded_freq <- c()
    for (i in seq_len(nrow(alloc))) {
      keep <- alloc[i, ] > 0
      if (any(keep)) {
        counts <- alloc[i, keep]
        clusters <- which(keep)
        expanded_vaf <- c(expanded_vaf, rep(unique_vaf[i], sum(counts)))
        expanded_cluster <- c(expanded_cluster, rep(clusters, counts))
        expanded_freq <- c(expanded_freq, rep(counts, counts))
      }
    }
    final_df <- data.frame(vaf = expanded_vaf, cluster = expanded_cluster, freq = expanded_freq)
  }

  missing_vaf <- !(vaf %in% final_df$vaf)
  if (any(missing_vaf)) {
    row_clusters <- max.col(probs[missing_vaf, , drop = FALSE], ties.method = "first")
    result <- data.frame(vaf = vaf[missing_vaf], cluster = row_clusters, freq = 1)
    final_df <- rbind(final_df, result)
  }

  final_df
}

generate_bootstrap_samples <- function(original_data, n_samples, num_decimal) {
  replicate(
    n = n_samples,
    expr = mean(sample(original_data, size = length(original_data), replace = TRUE))
  )
}

.cliffs_delta_abs <- function(x, y) {
  if (length(x) == 0L || length(y) == 0L) {
    return(NA_real_)
  }
  ## Bit-identical replacement for outer(x, y, "-")-based counting.
  ## Sort y once, then count y < x_i / y <= x_i via findInterval.
  ## The per-x_i counts n_lt/n_gt are integers, and the final value depends only
  ## on the INTEGER totals sum(n_lt) and sum(n_gt). When x has many repeats
  ## (e.g. the ~75-valued round(rbinom(100k)/depth) in peak_test), we evaluate
  ## findInterval on the unique x values and weight by their counts: the totals
  ## sum_u cnt_u*(n_lt_u - n_gt_u) equal the per-element sums exactly (integer
  ## arithmetic), so the result is bit-identical while collapsing 100k -> ~75.
  nx <- length(x); ny <- length(y)
  y_sorted <- sort(y)
  ux <- unique(x)
  if (length(ux) < nx) {
    cnt  <- tabulate(match(x, ux), nbins = length(ux))
    n_lt <- findInterval(ux, y_sorted, left.open = TRUE)
    n_le <- findInterval(ux, y_sorted)
    n_gt <- ny - n_le
    num  <- sum(as.numeric(cnt) * (n_lt - n_gt))      # exact for these integers
    abs(num / (nx * ny))
  } else {
    n_lt <- findInterval(x, y_sorted, left.open = TRUE)  # # y_j  <  x_i
    n_le <- findInterval(x, y_sorted)                    # # y_j <= x_i
    n_gt <- ny - n_le                                    # # y_j  >  x_i
    abs((sum(n_lt) - sum(n_gt)) / (nx * ny))
  }
}

# Count-based two-sample Wilcoxon rank-sum asymptotic p-value (two.sided),
# BIT-IDENTICAL to wilcox_test_p_asym / stats::wilcox.test on the asymptotic
# branch. The statistic, tie-corrected variance, continuity correction and
# 2*min(pnorm) two-sided p are algebraically the standard formulas; computing
# the average ranks from per-value counts (instead of sorting the full x) gives
# the identical numbers when `x` is heavily tied -- e.g. the 100k simulated VAFs
# in peak_test collapse to ~75 unique values. Verified identical() over random
# cases in /tmp/proto_counts.R. Used by the injected lean peak_test.
.wilcox_p_counts <- function(x, y) {
  nx <- length(x); ny <- length(y); n <- nx + ny
  if (nx == 0L || ny == 0L) return(NA_real_)
  uall <- sort(unique(c(x, y)))
  cx <- tabulate(match(x, uall), nbins = length(uall))
  cy <- tabulate(match(y, uall), nbins = length(uall))
  m  <- cx + cy                              # total count at each distinct value
  cum_before <- cumsum(m) - m                # # observations ranked before this group
  avg_rank   <- cum_before + (m + 1) / 2     # average rank within each tied group
  sum_rx     <- sum(cx * avg_rank)           # sum of ranks of x
  statistic  <- sum_rx - nx * (nx + 1) / 2
  tie_sum    <- sum(m^3 - m)
  sigma_sq   <- (nx * ny / 12) * ((nx + ny + 1) - tie_sum / (n * (n - 1)))
  if (!(sigma_sq > 0)) return(NA_real_)
  z    <- statistic - nx * ny / 2
  corr <- if (z > 0) 0.5 else if (z < 0) -0.5 else 0
  z    <- (z - corr) / sqrt(sigma_sq)
  p    <- 2 * min(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE))
  min(p, 1)
}

# Precomputed-counts variants of the two count-based statistics: the heavily
# tied side `x` is supplied as already-unique sorted values `xv` with integer
# counts `xc` (sum(xc) == length(x)). Used by the lean peak_test, where the
# simulated set's value-counts are obtained in ONE pass via tabulate() on the
# raw integer rbinom draws -- avoiding unique()/match() over the 100k doubles
# (twice, once per statistic). Results are bit-identical to .wilcox_p_counts /
# .cliffs_delta_abs called on the expanded vector.
.wilcox_p_counts_pre <- function(xv, xc, y) {
  nx <- sum(xc); ny <- length(y); n <- nx + ny
  if (nx == 0 || ny == 0L) return(NA_real_)
  uall <- sort(unique(c(xv, y)))
  cx <- numeric(length(uall)); cx[match(xv, uall)] <- xc
  cy <- tabulate(match(y, uall), nbins = length(uall))
  m  <- cx + cy
  cum_before <- cumsum(m) - m
  avg_rank   <- cum_before + (m + 1) / 2
  sum_rx     <- sum(cx * avg_rank)
  statistic  <- sum_rx - nx * (nx + 1) / 2
  tie_sum    <- sum(m^3 - m)
  sigma_sq   <- (nx * ny / 12) * ((nx + ny + 1) - tie_sum / (n * (n - 1)))
  if (!(sigma_sq > 0)) return(NA_real_)
  z    <- statistic - nx * ny / 2
  corr <- if (z > 0) 0.5 else if (z < 0) -0.5 else 0
  z    <- (z - corr) / sqrt(sigma_sq)
  min(2 * min(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE)), 1)
}

.cliffs_delta_abs_pre <- function(xv, xc, y) {
  nx <- sum(xc); ny <- length(y)
  if (nx == 0 || ny == 0L) return(NA_real_)
  ys   <- sort(y)
  n_lt <- findInterval(xv, ys, left.open = TRUE)
  n_le <- findInterval(xv, ys)
  n_gt <- ny - n_le
  abs(sum(as.numeric(xc) * (n_lt - n_gt)) / (nx * ny))
}

# Drop-in BIT-IDENTICAL replacement for sum(likelihoodExplore::likbeta(x, s1, s2,
# log=TRUE)). likbeta with length-K shape vectors returns value[k] = sum_i
# dbeta(x[i], s1[k], s2[k], log=TRUE) via plyr::llply over i + Reduce("+");
# inside likbeta `density` is bound to stats::dbeta. So
# sum(stats::dbeta(x, s1[k], s2[k], log=TRUE)) gives the same FP terms in the
# same left-fold order as likbeta's value[k], and indexing s[unique_idx]
# reconstructs vec[k] for every k without changing any FP value or order.
# Dedup is what gives the speedup (TEATIME's p_assign typically has <= 50
# unique values out of K~5000). Verified identical() in
# test/microbench_likbeta.R; ~85x faster on n=5693, K=5693, 50 unique pairs.
.likbeta_sum <- function(x, shape1, shape2) {
  K <- length(shape1)
  if (K == 1L) {
    return(sum(stats::dbeta(x, shape1, shape2, log = TRUE)))
  }
  key <- paste0(shape1, "_", shape2)
  unique_key <- unique(key)
  unique_idx <- match(key, unique_key)
  s <- vapply(seq_along(unique_key), function(j) {
    k_first <- match(unique_key[j], key)
    sum(stats::dbeta(x, shape1[k_first], shape2[k_first], log = TRUE))
  }, numeric(1))
  sum(s[unique_idx])
}

log_likelihood_mixture <- function(data, p_vec, depth) {
  a <- depth * p_vec
  b <- depth - a
  probs <- dbeta_matrix(data, a, b)
  df <- data.frame(probs, vaf = data, check.names = FALSE)
  names(df)[seq_along(a)] <- paste0("prob.", seq_along(a))

  if (length(p_vec) > 1L) {
    assigned <- beta_reassign(df)
    p_assign <- p_vec[assigned$cluster]
    vaf_assign <- assigned$vaf
  } else {
    p_assign <- rep(p_vec, length(data))
    vaf_assign <- data
  }

  # Fast mode uses the dedup-based .likbeta_sum (bit-identical, ~85x faster).
  # Default mode keeps the original likelihoodExplore::likbeta call so the
  # default reference path is unchanged.
  if (isTRUE(getOption("teatime.fast_version", FALSE))) {
    .likbeta_sum(vaf_assign, p_assign * depth, depth - p_assign)
  } else {
    sum(
      likelihoodExplore::likbeta(
        x = vaf_assign,
        shape1 = p_assign * depth,
        shape2 = depth - p_assign
      )
    )
  }
}

compute_AIC <- function(log_likelihood, num_params) {
  -2 * log_likelihood + 2 * num_params
}

compute_BIC <- function(log_likelihood, num_params, sample_size) {
  -2 * log_likelihood + num_params * log(sample_size)
}

# peak_simulate draws a length-`n_sim` reference sample of VAFs from a binomial
# read model: each read drawn from Binomial(depth, vaf), divided by depth and
# rounded to `num_decimal` digits. The resulting empirical distribution is the
# null against which peak_test compares observed cluster VAFs.
peak_simulate <- function(vaf, depth, num_decimal, n_sim = 100000L) {
  round(stats::rbinom(n_sim, depth, vaf) / depth, num_decimal)
}

peak_test_with_simulated <- function(
  df,
  Min_Sample_size,
  k,
  num_decimal,
  simulated,
  temp_check = FALSE
) {
  temp_k <- if (temp_check) 0 else k
  if (nrow(df) < temp_k) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }
  if (k < Min_Sample_size) {
    largest_right_values <- generate_bootstrap_samples(df[seq_len(k), "vaf"], Min_Sample_size, num_decimal)
  } else {
    largest_right_values <- df[seq_len(k), "vaf"]
  }
  largest_right_values <- as.vector(stats::na.omit(largest_right_values))

  # Two-sided Wilcoxon rank-sum p-value between the simulated null and the
  # observed cluster VAFs (asymptotic with continuity correction). Implemented
  # in C++ via the package's wilcox_test_p_asym for tight-loop performance.
  p_val <- wilcox_test_p_asym(as.numeric(simulated), as.numeric(largest_right_values))
  c(
    mean(simulated),
    mean(largest_right_values),
    .cliffs_delta_abs(simulated, largest_right_values),
    p_val
  )
}

# peak_test compares the largest k observed VAFs in `df` against a binomial-
# read null centered at `vaf`. Returns c(mean(simulated), mean(observed),
# absolute Cliff's delta, two-sided Wilcoxon p-value). The early-return guard
# below preserves the empty-input behaviour without drawing from the RNG, so
# nothing further upstream depends on this call's draws when k is too small.
peak_test <- function(
  df,
  Min_Sample_size,
  vaf,
  k,
  depth,
  num_decimal,
  temp_check = FALSE,
  n_sim = NULL,
  tol = NULL,
  batch_size = NULL
) {
  temp_k <- if (temp_check) 0 else k
  if (nrow(df) < temp_k) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }
  n_simulations <- 100000
  simulated <- round(stats::rbinom(n_simulations, depth, vaf) / depth, num_decimal)

  if (k < Min_Sample_size) {
    largest_right_values <- generate_bootstrap_samples(df[seq_len(k), "vaf"], Min_Sample_size, num_decimal)
  } else {
    largest_right_values <- df[seq_len(k), "vaf"]
  }
  largest_right_values <- as.vector(stats::na.omit(largest_right_values))

  # In fast mode the asymptotic Wilcoxon p-value is computed by the C++
  # implementation (bit-identical to stats::wilcox.test on the asymptotic
  # branch) for tight-loop throughput; default mode keeps the R routine.
  p_val <- if (isTRUE(getOption("teatime.fast_version", FALSE))) {
    wilcox_test_p_asym(as.numeric(simulated), as.numeric(largest_right_values))
  } else {
    stats::wilcox.test(simulated, largest_right_values)$p.value
  }
  c(
    mean(simulated),
    mean(largest_right_values),
    .cliffs_delta_abs(simulated, largest_right_values),
    p_val
  )
}
