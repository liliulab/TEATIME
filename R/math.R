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
# stats::dbeta(x, a[j], b[j]). In fast mode the K columns are computed by a
# single vectorised stats::dbeta call (recycled arguments); default mode keeps
# the per-column sapply form. The two paths produce bit-identical output.
dbeta_matrix <- function(x, a, b) {
  if (isTRUE(getOption("teatime.fast_version", FALSE))) {
    matrix(stats::dbeta(rep(x, length(a)),
                        rep(a, each = length(x)),
                        rep(b, each = length(x))),
           nrow = length(x))
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

  # Fast mode expands the (vaf x cluster) allocation matrix into a tidy
  # (vaf, cluster, freq) table in one vectorised pass (row-major order to
  # preserve downstream ordering); default mode walks the matrix row by row.
  if (isTRUE(getOption("teatime.fast_version", FALSE))) {
    nz <- which(alloc > 0, arr.ind = TRUE)
    if (nrow(nz) > 0L) {
      nz <- nz[order(nz[, "row"], nz[, "col"]), , drop = FALSE]
      counts <- alloc[nz]
      final_df <- data.frame(
        vaf     = rep(unique_vaf[nz[, "row"]], counts),
        cluster = rep(nz[, "col"], counts),
        freq    = rep(counts, counts)
      )
    } else {
      final_df <- data.frame(vaf = numeric(0), cluster = integer(0), freq = numeric(0))
    }
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
  ## Sort y once, then for each x_i count y < x_i and y <= x_i via
  ## findInterval. O((n+m) log m) instead of O(n*m).
  nx <- length(x); ny <- length(y)
  y_sorted <- sort(y)
  n_lt <- findInterval(x, y_sorted, left.open = TRUE)  # # y_j  <  x_i
  n_le <- findInterval(x, y_sorted)                    # # y_j <= x_i
  n_gt <- ny - n_le                                    # # y_j  >  x_i
  abs((sum(n_lt) - sum(n_gt)) / (nx * ny))
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

  sum(
    likelihoodExplore::likbeta(
      x = vaf_assign,
      shape1 = p_assign * depth,
      shape2 = depth - p_assign
    )
  )
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
