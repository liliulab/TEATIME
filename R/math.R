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

  final_df <- data.frame(
    vaf = expanded_vaf,
    cluster = expanded_cluster,
    freq = expanded_freq
  )

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
  cmp <- outer(x, y, FUN = "-")
  abs((sum(cmp > 0) - sum(cmp < 0)) / (length(x) * length(y)))
}

log_likelihood_mixture <- function(data, p_vec, depth) {
  a <- depth * p_vec
  b <- depth - a
  probs <- sapply(seq_along(a), function(i) stats::dbeta(data, a[i], b[i]))
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

peak_test <- function(
  df,
  Min_Sample_size,
  vaf,
  k,
  depth,
  num_decimal,
  temp_check = FALSE,
  n_sim = 10000,
  tol = 1e-4,
  batch_size = 2000
) {
  temp_k <- if (temp_check) 0 else k
  if (nrow(df) < temp_k) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }

  if (tol <= 0) {
    simulated <- round(stats::rbinom(n_sim, depth, vaf) / depth, num_decimal)
  } else {
    simulated <- numeric(0)
    running_mean <- NA_real_
    while (length(simulated) < n_sim) {
      draw_n <- min(batch_size, n_sim - length(simulated))
      new_vals <- round(stats::rbinom(draw_n, depth, vaf) / depth, num_decimal)
      simulated <- c(simulated, new_vals)
      new_mean <- mean(simulated)
      if (!is.na(running_mean) && abs(new_mean - running_mean) < tol && length(simulated) >= batch_size * 2) {
        break
      }
      running_mean <- new_mean
    }
  }

  if (k < Min_Sample_size) {
    largest_right_values <- generate_bootstrap_samples(df[seq_len(k), "vaf"], Min_Sample_size, num_decimal)
  } else {
    largest_right_values <- df[seq_len(k), "vaf"]
  }
  largest_right_values <- as.vector(stats::na.omit(largest_right_values))

  w.right.result <- stats::wilcox.test(simulated, largest_right_values)
  c(
    mean(simulated),
    mean(largest_right_values),
    .cliffs_delta_abs(simulated, largest_right_values),
    w.right.result$p.value
  )
}
