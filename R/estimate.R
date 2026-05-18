.vaf_prob_df <- function(vaf_set, p_vec, depth) {
  probs <- sapply(seq_along(p_vec), function(i) stats::dbeta(vaf_set, depth * p_vec[i], depth - depth * p_vec[i]))
  df <- data.frame(probs, vaf = vaf_set, check.names = FALSE)
  names(df)[seq_along(p_vec)] <- paste0("prob.", seq_along(p_vec))
  df
}

compare_real_simu_peak <- function(
  peakdata,
  p,
  vaf_set,
  ctx,
  num_decimal,
  right_df = NULL,
  right_most_vaf = NA_real_,
  mu_small = FALSE,
  bac = FALSE,
  n_wilcox = 50,
  tol = NULL
) {
  depth <- ctx$depth
  row_indices <- seq_len(nrow(peakdata))

  compute_row <- function(i) {
    row <- peakdata[i, , drop = FALSE]
    cell_div_value <- row[["cell.div"]]
    mu_value <- row[["mu_est"]]

    if (is.na(cell_div_value)) {
      largest_right_values <- right_df[seq_len(mu_value), "vaf"]
      p_values <- replicate(n_wilcox, {
        simulated <- round(stats::rbinom(round(mu_value), depth, right_most_vaf) / depth, num_decimal)
        stats::wilcox.test(largest_right_values, simulated)$p.value
      })
      ll <- log_likelihood_mixture(largest_right_values, right_most_vaf, depth)
      return(c(mean(p_values), ll, compute_AIC(ll, 1), compute_BIC(ll, 1, length(largest_right_values))))
    }

    if (!bac) {
      i_values <- if (mu_small) seq_len(20) else seq_len(cell_div_value)
      result_vector <- c(0.5, p / 2, vaf_at_div(i_values, p, ctx))
      df <- beta_reassign(.vaf_prob_df(vaf_set, result_vector, depth))
      data <- df[df$cluster >= 2, "vaf"]
    } else {
      data <- vaf_set
    }

    vaf_list <- c(p / 2, vaf_at_div(seq_len(cell_div_value), p, ctx))
    p_values <- replicate(n_wilcox, {
      sim_list <- sapply(vaf_list, function(vaf) round(stats::rbinom(round(mu_value), depth, vaf) / depth, num_decimal))
      stats::wilcox.test(data, as.vector(t(sim_list)))$p.value
    })
    ll <- log_likelihood_mixture(data, vaf_list, depth)
    c(mean(p_values), ll, compute_AIC(ll, length(vaf_list)), compute_BIC(ll, length(vaf_list), length(data)))
  }

  t(vapply(row_indices, compute_row, numeric(4)))
}

calculate_mu <- function(data, ctx) {
  possible_error <- tryCatch({
    strucchange::breakpoints(cumsum ~ x, data = data, h = 3 / nrow(data))
  }, error = function(e) e)

  if (!inherits(possible_error, "error")) {
    mu_turn <- possible_error
    bf <- strucchange::breakfactor(mu_turn)
    aic_values <- c()
    bic_values <- c()
    slopes <- c()
    for (i in unique(bf)) {
      segment_data <- subset(data, bf == i)
      fm_segment <- stats::lm(cumsum ~ x, data = segment_data)
      coeffs <- stats::coef(fm_segment)
      slope <- coeffs[-1] * (-1) * ctx$beta * log(2)
      slopes <- c(slopes, slope)
      aic_values <- c(aic_values, stats::AIC(fm_segment))
      bic_values <- c(bic_values, stats::BIC(fm_segment))
    }
    return((slopes[which.min(aic_values)] + slopes[which.min(bic_values)]) / 2)
  }

  fm <- stats::lm(cumsum ~ x, data = data)
  summary(fm)$coefficients[2, 1] * (-1) * ctx$beta * log(2)
}

get_slope <- function(df, p, result_vector, ctx, num_decimal = 3) {
  g1 <- df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(count = sum(freq), .groups = "drop")
  g1 <- g1[order(g1$cluster), , drop = FALSE]

  g1 <- g1 %>%
    dplyr::arrange(cluster) %>%
    dplyr::mutate(vaf = result_vector[cluster])

  g1_sub <- g1 %>%
    dplyr::mutate(vaf = format(vaf, nsmall = num_decimal)) %>%
    dplyr::arrange(dplyr::desc(vaf), dplyr::desc(count)) %>%
    dplyr::mutate(cumsum = cumsum(count)) %>%
    dplyr::group_by(vaf) %>%
    dplyr::slice_max(cumsum) %>%
    dplyr::ungroup()

  g1_sub$vaf <- as.numeric(g1_sub$vaf)
  g1_sub$lnf <- 2 * g1_sub$vaf - p
  g1 <- g1_sub[g1_sub$lnf > 0, , drop = FALSE]
  g1$x <- log(g1$lnf)
  g1 <- g1[order(g1$cluster), , drop = FALSE]
  if (nrow(g1) >= 6) {
    g1 <- g1[2:floor(nrow(g1) * 0.5), , drop = FALSE]
  }
  calculate_mu(g1[, c("cumsum", "x"), drop = FALSE], ctx)
}

slope_simu <- function(cell_div, mu, p, ctx, num_decimal = 3) {
  result_vector <- vaf_at_div(seq_len(cell_div), p, ctx)
  a <- ctx$depth * result_vector
  b <- ctx$depth - a
  mu_est_list <- numeric(3)
  for (try_idx in 1:3) {
    all_simulated_vafs_list <- sapply(result_vector, function(vaf) {
      round(stats::rbinom(round(mu), ctx$depth, vaf) / ctx$depth, num_decimal)
    })
    all_simulated_vafs <- as.vector(t(all_simulated_vafs_list))
    probs <- sapply(seq_along(a), function(i) stats::dbeta(all_simulated_vafs, a[i], b[i]))
    df <- data.frame(prob = probs, vaf = all_simulated_vafs)
    df <- beta_reassign(df)
    df <- df[!duplicated(df), , drop = FALSE]
    mu_est_list[try_idx] <- get_slope(df, p, result_vector, ctx, num_decimal)
  }
  mu_est_list
}

slope_method <- function(p, vaf_set, ctx, p_thre, start_div = 1, end_div = NA, num_decimal = 3, ...) {
  min_mu <- 3
  cell_div <- if (is.na(end_div)) start_div:round(length(vaf_set) / min_mu) else start_div:end_div
  mu_list <- round(length(vaf_set) / cell_div)
  mu_data <- data.frame(mu_list = mu_list, cell.div = cell_div)
  unique_mu <- unique(mu_list)
  selected_div <- vapply(unique_mu, function(x) {
    candidates <- mu_data$cell.div[mu_data$mu_list == x]
    if (length(candidates) > 1) sample(candidates, 1) else candidates
  }, numeric(1))
  unique_df <- data.frame(mu_list = unique_mu, cell.div = selected_div)
  div_list <- unique_df$cell.div[
    unique_df$cell.div <= stats::quantile(unique_df$cell.div, 0.75) &
      unique_df$cell.div >= stats::quantile(unique_df$cell.div, 0.25)
  ]
  if (is.na(end_div)) {
    vaf_set <- vaf_set[vaf_set > p / 2]
  }

  collect_rows(function(cell_div_value) {
    mu_est <- length(vaf_set) / cell_div_value
    if (mu_est <= 3) {
      return(NULL)
    }
    tryCatch({
      result_vector <- vaf_at_div(seq_len(cell_div_value), p, ctx)
      df <- beta_reassign(.vaf_prob_df(vaf_set, result_vector, ctx$depth))
      df <- df[!duplicated(df), , drop = FALSE]
      mulist <- c()
      plist <- c()
      mean_list <- NA_real_
      for (try_idx in 1:3) {
        mu_from_real <- get_slope(df, p, result_vector, ctx, num_decimal)
        mu_from_simu <- slope_simu(cell_div_value, mu_est, p, ctx, num_decimal)
        mean_list <- mean(mu_from_simu)
        sd_list <- max(stats::sd(mu_from_simu), 1.5)
        z_score <- (mu_from_real - mean_list) / sd_list
        mulist <- c(mulist, mu_from_real)
        plist <- c(plist, abs(z_score))
      }
      best <- which.min(plist)
      data.frame(
        mu = mu_est,
        mu_real = mulist[best],
        mu_simu = mean_list,
        cell.div = cell_div_value,
        z_score = plist[best],
        p = p
      )
    }, error = function(e) NULL)
  }, as.list(div_list))
}

peak_method <- function(p, vaf_set, ctx, p_thre, celldivlist = NULL, num_decimal = 3, ...) {
  dots <- list(...)
  n_sim_peak <- dots$n_sim_peak %||% 100000
  tol_peak <- dots$tol_peak %||% 0
  n_wilcox <- dots$n_wilcox %||% 50
  collect_data_check <- dots$collect_data_check
  collect <- is.null(collect_data_check)
  min_sample_size <- 6
  right_most_vaf <- vaf_at_div(1, p, ctx)
  right_df <- data.frame(vaf = vaf_set)
  right_df$abs_diff <- abs(right_df$vaf - right_most_vaf)
  right_df$prob <- stats::pbeta(vaf_set, ctx$depth * right_most_vaf, ctx$depth - ctx$depth * right_most_vaf)
  right_df$score <- right_df$prob / max(right_df$prob) - right_df$abs_diff / max(right_df$abs_diff)
  right_df <- right_df[order(-right_df$score), , drop = FALSE]

  if (!collect) {
    collect.data <- collect_data_check
  } else if (is.null(celldivlist)) {
    collect.data <- collect_rows(function(mu) {
      right_result <- peak_test(right_df, min_sample_size, right_most_vaf, mu, ctx$depth, num_decimal, n_sim = n_sim_peak, tol = tol_peak)
      data.frame(
        cell.div = NA,
        mu_est = mu,
        simu.left.mean = NA,
        simu.right.mean = right_result[1],
        left.mean = NA,
        right.mean = right_result[2],
        left.cd = NA,
        right.cd = right_result[3],
        left.p = NA,
        right.p = right_result[4]
      )
    }, as.list(3:length(vaf_set)))
  } else {
    collect.data <- collect_rows(function(cell_div) {
      result_vector <- c(0.5, p / 2, vaf_at_div(seq_len(cell_div), p, ctx))
      df <- beta_reassign(.vaf_prob_df(vaf_set, result_vector, ctx$depth))
      mu1 <- max(nrow(df[df$cluster >= 3, , drop = FALSE]) / cell_div, 1)
      left_df <- df[df$cluster == 2 + cell_div, , drop = FALSE]
      right_sub <- df[df$cluster == 3, , drop = FALSE]
      if (nrow(right_sub) < nrow(left_df)) {
        right_sub <- df[df$cluster == 2, , drop = FALSE]
      }
      if (!(nrow(right_sub) > mu1 && mu1 > 3)) {
        return(NULL)
      }
      left_most_vaf <- vaf_at_div(cell_div, p, ctx)
      left_df$abs_diff <- abs(left_df$vaf - left_most_vaf)
      left_df$prob <- stats::pbeta(left_df$vaf, ctx$depth * left_most_vaf, ctx$depth - ctx$depth * left_most_vaf)
      left_df$score <- left_df$prob / max(left_df$prob) - left_df$abs_diff / max(left_df$abs_diff)
      left_df <- left_df[order(-left_df$score), , drop = FALSE]

      right_sub$abs_diff <- abs(right_sub$vaf - right_most_vaf)
      right_sub$prob <- stats::pbeta(right_sub$vaf, ctx$depth * right_most_vaf, ctx$depth - ctx$depth * right_most_vaf)
      right_sub$score <- right_sub$prob / max(right_sub$prob) - right_sub$abs_diff / max(right_sub$abs_diff)
      right_sub <- right_sub[order(-right_sub$score), , drop = FALSE]

      k <- round(mu1)
      left_result <- peak_test(left_df, min_sample_size, left_most_vaf, k, ctx$depth, num_decimal, n_sim = n_sim_peak, tol = tol_peak)
      right_result <- peak_test(right_sub, min_sample_size, right_most_vaf, k, ctx$depth, num_decimal, n_sim = n_sim_peak, tol = tol_peak)
      data.frame(
        cell.div = cell_div,
        mu_est = mu1,
        simu.left.mean = left_result[1],
        simu.right.mean = right_result[1],
        left.mean = left_result[2],
        right.mean = right_result[2],
        left.cd = left_result[3],
        right.cd = right_result[3],
        left.p = left_result[4],
        right.p = right_result[4]
      )
    }, as.list(celldivlist))
  }

  temp_keep <- collect.data
  if (all(is.na(temp_keep$left.p))) {
    if (!all(is.na(temp_keep$right.p))) {
      collect.data <- temp_keep[temp_keep$right.p > p_thre, , drop = FALSE]
    }
    if (nrow(collect.data) == 0) {
      collect.data <- temp_keep
    }
    collect.data$right_rank <- rank(collect.data$right.cd)
    top_n <- floor((if (collect) (if (nrow(collect.data) > 4) 0.25 else if (nrow(collect.data) > 2) 0.5 else 1) else (if (nrow(collect.data) > 2) 0.5 else 1)) * nrow(collect.data))
    top_pick <- utils::head(collect.data[order(collect.data$right_rank), , drop = FALSE], top_n)
  } else {
    collect.data <- temp_keep[temp_keep$left.p > p_thre | temp_keep$right.p > p_thre, , drop = FALSE]
    collect.data <- collect.data[stats::complete.cases(collect.data[, 3]), , drop = FALSE]
    collect.data$left.cd <- ifelse(is.na(collect.data$left.cd), 10, collect.data$left.cd)
    top_n <- floor((if (nrow(collect.data) > 2) 0.5 else 1) * nrow(collect.data))
    if (all(collect.data$left.cd == 10)) {
      collect.data$right_rank <- rank(collect.data$right.cd)
      top_pick <- utils::head(collect.data[order(collect.data$right_rank), , drop = FALSE], top_n)
    } else {
      collect.data$left_rank <- rank(collect.data$left.cd)
      collect.data$right_rank <- rank(collect.data$right.cd)
      top_left <- utils::head(collect.data[order(collect.data$left_rank), , drop = FALSE], top_n)
      top_right <- utils::head(collect.data[order(collect.data$right_rank), , drop = FALSE], top_n)
      top_pick <- unique(rbind(top_left, top_right))
    }
  }

  top_pick_p <- as.data.frame(compare_real_simu_peak(top_pick, p, vaf_set, ctx, num_decimal, right_df, right_most_vaf, n_wilcox = n_wilcox))
  top_pick$r1 <- ifelse(top_pick_p$V1 < 0.05, 0, top_pick_p$V1)
  top_pick$loglike <- top_pick_p$V2
  top_pick$aic <- top_pick_p$V3
  top_pick$bic <- top_pick_p$V4
  top_pick$top_pick_score <- top_pick$r1 + 1 / top_pick$right_rank

  if (is.null(celldivlist)) {
    top_pick <- top_pick[order(-top_pick$top_pick_score), , drop = FALSE]
  } else {
    if (nrow(top_pick) > 1) {
      top_pick <- top_pick[top_pick$top_pick_score > stats::median(top_pick$top_pick_score), , drop = FALSE]
    }
    top_pick <- top_pick[order(top_pick$right_rank), , drop = FALSE]
  }

  utils::head(top_pick, min(10, nrow(top_pick)))
}

simu_slope <- function(vaf_set, p, ctx) {
  G1.sub <- data.frame(vaf = vaf_set, stringsAsFactors = FALSE)
  G1.sub$count <- 1
  G1.sub$vaf <- format(G1.sub$vaf, nsmall = 3)
  G1.sub <- G1.sub[order(G1.sub$vaf, G1.sub$count, decreasing = TRUE), , drop = FALSE]
  G1.sub$cumsum <- cumsum(G1.sub$count)
  G1.sub <- G1.sub[!duplicated(G1.sub$vaf, fromLast = TRUE), , drop = FALSE]
  G1.sub$vaf <- as.numeric(G1.sub$vaf)
  G1.sub$lnf <- 2 * G1.sub$vaf - p
  G1 <- G1.sub[G1.sub$lnf > 0, , drop = FALSE]
  G1$x <- log(G1$lnf)
  num_rows_to_remove <- round(nrow(G1) * 0.05)
  rows_to_remove <- c(1:num_rows_to_remove, (nrow(G1) - num_rows_to_remove + 1):nrow(G1))
  if (nrow(G1[-rows_to_remove, , drop = FALSE]) > 5) {
    G1 <- G1[-rows_to_remove, , drop = FALSE]
  }
  calculate_mu(G1[, c("cumsum", "x"), drop = FALSE], ctx)
}

mu_find_small <- function(vaf_set, p, ctx) {
  slope1 <- simu_slope(vaf_set, p, ctx)
  slope2 <- simu_slope(vaf_set, p + 0.01, ctx)
  slope3 <- simu_slope(vaf_set, p - 0.01, ctx)
  if (slope1 == slope2 || slope1 == slope3 || slope2 == slope3) {
    noise <- stats::runif(3, min = -1e-5, max = 1e-5)
    slope1 <- slope1 + noise[1]
    slope2 <- slope2 + noise[2]
    slope3 <- slope3 + noise[3]
  }
  mu_collect <- c(slope1, slope2, slope3)
  means_to_test <- 1:3
  test_results <- sapply(means_to_test, function(x) stats::t.test(mu_collect, mu = x)$p.value)
  bac_p <- 1 - p
  sorted_vector <- sort(vaf_set)
  test_results2 <- sapply(means_to_test, function(x) {
    i <- seq_len(round(length(vaf_set) / x))
    diff1 <- sorted_vector[(1 + x):length(sorted_vector)] - sorted_vector[1:(length(sorted_vector) - x)]
    diff2 <- (bac_p / 2 * exp(log(2) * ctx$beta * i)) - (bac_p / 2 * exp(log(2) * ctx$beta * (i + 1)))
    stats::wilcox.test(diff2, diff1)$p.value
  })
  data.frame(mu_est = means_to_test, p_value = test_results, wx_p = test_results2)
}

mu_estimation_small <- function(data, p, ctx, p_thre, num_decimal = 3, ...) {
  dots <- list(...)
  n_wilcox <- dots$n_wilcox %||% 50
  vaf_set <- ctx$main_cluster_vaf
  result_vector <- c(0.5, p / 2, vaf_at_div(seq_len(20), p, ctx))
  df <- beta_reassign(.vaf_prob_df(vaf_set, result_vector, ctx$depth))
  result_df <- mu_find_small(df[df$cluster >= 2, "vaf"], p, ctx)
  mu_small_selection <- result_df[result_df$p_value > 0.05 & result_df$wx_p > 0.05, , drop = FALSE]
  data$mu_est <- data$mu_real
  top_pick <- as.data.frame(compare_real_simu_peak(data, p, vaf_set, ctx, num_decimal, n_wilcox = n_wilcox))
  data$loglike <- top_pick$V2
  data$aic <- top_pick$V3
  data$bic <- top_pick$V4

  pick.cell.div <- data$cell.div
  pick.mu <- data$mu_real
  z_score <- data$z_score
  pick.log <- data$loglike
  pick.bic <- data$bic
  pick.aic <- data$aic

  if (nrow(mu_small_selection) > 0) {
    num_df <- length(df[df$cluster >= 2, "vaf"])
    mu_small_selection$cell.div <- round(num_df / mu_small_selection$mu_est)
    small_mu_pick <- as.data.frame(compare_real_simu_peak(mu_small_selection, p, vaf_set, ctx, num_decimal, mu_small = TRUE, n_wilcox = n_wilcox))
    mu_small_selection$r1 <- ifelse(small_mu_pick$V1 < 0.05, 0, small_mu_pick$V1)
    mu_small_selection$loglike <- small_mu_pick$V2
    mu_small_selection$aic <- small_mu_pick$V3
    mu_small_selection$bic <- small_mu_pick$V4
    mu_small_selection <- mu_small_selection[which(mu_small_selection$bic == min(mu_small_selection$bic)), , drop = FALSE]

    if (mu_small_selection$bic < min(data$bic) * 0.5) {
      pick.cell.div <- round(length(df[df$cluster >= 2, "vaf"]) / mu_small_selection$mu_est)
      pick.mu <- mu_small_selection$mu_est
      pick.bic <- mu_small_selection$bic
      pick.aic <- mu_small_selection$aic
      pick.log <- mu_small_selection$loglike
      z_score <- rep(1, length(mu_small_selection$aic))
    }
  }

  data.frame(
    cell.div = pick.cell.div,
    mu = pick.mu,
    loglike = pick.log,
    bic = pick.bic,
    aic = pick.aic,
    z_score = z_score,
    p = p
  )
}

calculate_border_vaf_fit <- function(ctx) {
  probs <- stats::pbeta(ctx$main_cluster_vaf, ctx$depth * 0.5, ctx$depth - ctx$depth * 0.5)
  combined <- data.frame(main_vaf = ctx$main_cluster_vaf, probs = probs)
  if (nrow(combined[combined$probs <= 0.05, , drop = FALSE]) > 0) {
    max(combined[combined$probs <= 0.05, "main_vaf"])
  } else {
    min(combined$main_vaf)
  }
}

calculate_left_right_most_vaf_fit <- function(border_vaf, ctx, num_decimal = 3) {
  main_vaf <- ctx$main_cluster_vaf
  depth <- ctx$depth
  beta <- ctx$beta
  prob_df <- function(x, s1, s2) {
    probs <- sapply(seq_along(s1), function(j) stats::dbeta(x, s1[j], s2[j]))
    d <- data.frame(probs, vaf = x, check.names = FALSE)
    names(d)[seq_along(s1)] <- paste0("prob.", seq_along(s1))
    d
  }
  m <- RBesT::automixfit(main_vaf, type = "beta", Nc = 2:10, thresh = 0, k = 6, Ninit = min(50, round(length(main_vaf) / 5)), Niter.max = 10000)
  a <- m["a", ]
  b <- m["b", ]
  mean_ab <- a / (a + b)
  approx_vaf <- min(mean_ab)
  approx_vaf_index <- which.min(mean_ab)
  left_most_vaf <- round(approx_vaf, num_decimal)
  inita <- c(0.5 * depth, a[approx_vaf_index])
  initb <- c(0.5 * depth, b[approx_vaf_index])
  df <- beta_reassign(prob_df(main_vaf, inita, initb))
  main_vaf_update <- df[df$cluster < 2, "vaf"]
  right_most_vaf <- approx_vaf
  right_save <- approx_vaf
  iteration <- 1
  while (iteration <= 10000) {
    m <- RBesT::automixfit(main_vaf_update, type = "beta", Nc = 1:10, thresh = 0, k = 6, Ninit = min(50, round(length(main_vaf_update) / 5)), Niter.max = 10000)
    updatea <- m["a", ]
    updateb <- m["b", ]
    mean_ab <- updatea / (updatea + updateb)
    approx_vaf <- min(mean_ab)
    if (approx_vaf >= 0.5 || length(mean_ab) == 1) break
    approx_vaf_index <- which.min(mean_ab)
    currenta <- c(0.5 * depth, updatea[approx_vaf_index])
    currentb <- c(0.5 * depth, updateb[approx_vaf_index])
    df <- beta_reassign(prob_df(main_vaf_update, currenta, currentb))
    main_vaf_update <- df[df$cluster < 2, "vaf"]
    if (approx_vaf < 0.5) {
      right_most_vaf <- approx_vaf
      if (approx_vaf > right_save) right_save <- approx_vaf
    }
    iteration <- iteration + 1
  }
  right_th_vaf <- left_most_vaf + (0.5 - left_most_vaf) / exp(log(2) * beta)
  right_most_vaf <- round(max(right_most_vaf, border_vaf, right_th_vaf, right_save), num_decimal)
  list(
    p = (right_most_vaf * 2 * exp(log(2) * beta) - 1) / (exp(log(2) * beta) - 1),
    start.div = if (right_most_vaf == left_most_vaf) 1 else 2
  )
}

estimate_mu <- function(p, vaf_set, ctx, p_thre, start_div, num_decimal = 3, ...) {
  slope_candidates <- get_mu_method("slope")(p, vaf_set, ctx, p_thre, start_div = start_div, num_decimal = num_decimal, ...)
  slope_candidates <- slope_candidates[slope_candidates$mu_real > 3, , drop = FALSE]
  filtered_slope <- if (nrow(slope_candidates[abs(slope_candidates$z_score) <= 1.96, , drop = FALSE]) > 0) {
    slope_candidates[abs(slope_candidates$z_score) <= 1.96, , drop = FALSE]
  } else {
    slope_candidates
  }

  peak_candidates <- get_mu_method("peak")(p, vaf_set, ctx, p_thre, celldivlist = filtered_slope$cell.div, num_decimal = num_decimal, ...)
  range1 <- ifelse(nrow(peak_candidates) == 1, peak_candidates$mu_est - 1, min(peak_candidates$mu_est))
  range2 <- ifelse(nrow(peak_candidates) == 1, peak_candidates$mu_est + 1, max(peak_candidates$mu_est))
  overlap <- if (range2 > stats::median(filtered_slope$mu)) {
    filtered_slope[filtered_slope$mu_real >= range1 & filtered_slope$mu_real <= range2, , drop = FALSE]
  } else {
    filtered_slope[filtered_slope$mu_real >= range1, , drop = FALSE]
  }
  if (nrow(overlap) == 0) {
    overlap <- filtered_slope
  }
  overlap <- stats::na.omit(overlap)

  all_p_data <- mu_estimation_small(overlap, p, ctx, p_thre, num_decimal, ...)
  if (nrow(all_p_data) >= 3) {
    mean_filtered <- stats::median(all_p_data$mu)
    sd_filtered <- stats::mad(all_p_data$mu)
    if (sd_filtered > 0) {
      z1 <- abs((all_p_data$mu - mean_filtered) / sd_filtered)
      all_p_data <- all_p_data[z1 < 2, , drop = FALSE]
    }
  }
  all_p_data$z_score1 <- 1
  all_p_data$z_score <- abs(all_p_data$z_score)
  all_p_data$bicrank <- rank(all_p_data$bic)
  all_p_data$score <- 1 / all_p_data$bicrank + (1 - all_p_data$z_score) * 0.25
  all_p_data$lowerbound1 <- max(peak_candidates$mu_est)
  all_p_data$lowerbound2 <- min(peak_candidates$mu_est)
  all_p_data[order(-all_p_data$score), , drop = FALSE]
}

filter_mu_estimate <- function(data) {
  if (nrow(data[data$mu > 3, , drop = FALSE]) > 1) {
    data <- data[data$mu > 3, , drop = FALSE]
  }
  if (nrow(data) > 3) {
    mean_filtered <- stats::median(data$cell.div)
    sd_filtered <- stats::mad(data$cell.div)
    data$z_score <- abs((data$cell.div - mean_filtered) / sd_filtered)
    data <- data[data$z_score < 2, , drop = FALSE]
  }
  data
}

find_p_process <- function(second.vaf, mean.a.b, upper_clonal_vaf, clonal.vaf.left, ctx) {
  m <- RBesT::automixfit(
    second.vaf,
    type = "beta",
    Nc = 1:10,
    thresh = 0,
    k = 6,
    Ninit = min(50, round(length(second.vaf) / 3)),
    Niter.max = 10000
  )
  a1 <- m["a", ]
  b1 <- m["b", ]
  mean.a.b1 <- a1 / (a1 + b1)
  second_clone_left_vaf <- min(mean.a.b1)
  clear <- ifelse(length(mean.a.b) > 1, FALSE, TRUE)

  left_limit <- (min(second.vaf) + second_clone_left_vaf) / 2
  if (length(mean.a.b) > 1) {
    right_limit <- max(
      (((upper_clonal_vaf + max(clonal.vaf.left)) / 2 * 2 * exp(log(2) * ctx$beta) - 1) /
        (exp(log(2) * ctx$beta) - 1) / 2),
      (max(second.vaf) + second_clone_left_vaf) / 2
    )
  } else {
    right_limit <- ((max(second.vaf) + upper_clonal_vaf) / 2 * 2 * exp(log(2) * ctx$beta) - 1) /
      (exp(log(2) * ctx$beta) - 1) / 2
    if (right_limit > 0.45 || right_limit < left_limit) {
      right_limit <- (max(second.vaf) + second_clone_left_vaf) / 2
    }
  }
  if (length(second.vaf[second.vaf > right_limit]) <= 6) {
    sorted_vaf <- sort(second.vaf, decreasing = TRUE)
    right_limit <- sorted_vaf[max(round(0.25 * length(sorted_vaf)), 6)]
  }
  my_seq <- c(left_limit, (left_limit + right_limit) / 2, right_limit)
  my_seq <- my_seq[my_seq > 0.5 / (1 + exp(log(2) * ctx$beta))]
  list(my_seq = my_seq, clear = clear)
}

iterate_p_optimize <- function(clear, give.vaf, upper_clonal_vaf, clonal.vaf.left, second.vaf, num_decimal, p_thre, ctx) {
  min_sample_size <- 6
  p <- give.vaf * 2
  reliable <- 1

  suppose_right_vaf <- vaf_at_div(1, p, ctx)
  first_vaf_list <- c()
  if (!clear) {
    a <- if (suppose_right_vaf > upper_clonal_vaf) {
      c(ctx$depth * 0.5, ctx$depth * suppose_right_vaf)
    } else {
      c(ctx$depth * upper_clonal_vaf, ctx$depth * suppose_right_vaf)
    }
    b <- ctx$depth - a
    probs <- sapply(seq_along(a), function(i) stats::dbeta(clonal.vaf.left, a[i], b[i]))
    df <- data.frame(probs, vaf = clonal.vaf.left, check.names = FALSE)
    names(df)[seq_along(a)] <- paste0("prob.", seq_along(a))
    df <- beta_reassign(df)
    first_vaf_list <- df[df$cluster == 2, "vaf"]
  }

  vaf_set <- c(first_vaf_list, second.vaf[second.vaf > give.vaf])
  result_vector <- vaf_at_div(seq_len(20), p, ctx)
  df <- beta_reassign(.vaf_prob_df(vaf_set, result_vector, ctx$depth))
  df <- df[!duplicated(df), , drop = FALSE]

  df_count_freq <- df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(count = sum(freq), .groups = "drop")
  df_count_freq <- as.data.frame(df_count_freq)
  if (min(df_count_freq$cluster) > 1) {
    df_count_freq <- rbind(data.frame(cluster = 1, count = 0), df_count_freq)
    df_count_freq <- df_count_freq[order(df_count_freq$cluster), , drop = FALSE]
  }

  max_mu <- df_count_freq[1, ]$count + df_count_freq[2, ]$count
  if (nrow(df_count_freq) > 6 & df_count_freq[2, "count"] / df_count_freq[3, "count"] <= 4) {
    df_count_freq <- df_count_freq[2:3, , drop = FALSE]
  } else {
    df_count_freq <- df_count_freq[2, , drop = FALSE]
  }
  min_mu <- max(3, min(df_count_freq$count / 4))
  total_count_temp <- length(vaf_set)
  start.div <- round(total_count_temp / max_mu)
  end.div <- round(total_count_temp / min_mu)

  collect.data <- slope_method(p, vaf_set, ctx, p_thre, start_div = start.div, end_div = end.div, num_decimal = num_decimal)
  collect.data <- collect.data[collect.data$mu_real > 3, , drop = FALSE]

  if (nrow(collect.data[abs(collect.data$z_score) <= 1.96, , drop = FALSE]) > 0) {
    collect.data <- collect.data[abs(collect.data$z_score) <= 1.96, , drop = FALSE]
    first_div_vaf_data <- df[df$cluster == as.numeric(df_count_freq$cluster[1]), , drop = FALSE]
    sec.div.vaf <- vaf_at_div(2, p, ctx)
    first_div_vaf_data$abs_diff <- abs(first_div_vaf_data$vaf - sec.div.vaf)
    first_div_vaf_data$prob <- stats::pbeta(first_div_vaf_data$vaf, ctx$depth * sec.div.vaf, ctx$depth - ctx$depth * sec.div.vaf)
    first_div_vaf_data$score <- first_div_vaf_data$prob / max(first_div_vaf_data$prob) - first_div_vaf_data$abs_diff / max(first_div_vaf_data$abs_diff)
    right_df <- first_div_vaf_data[order(-first_div_vaf_data$score), , drop = FALSE]
    celldivlist <- collect.data$cell.div

    collect.data1 <- collect_rows(function(cell.div) {
      mu_suppose <- length(vaf_set) / cell.div
      k <- round(mu_suppose)
      right_result <- peak_test(right_df, min_sample_size, sec.div.vaf, k, ctx$depth, num_decimal)
      data.frame(
        cell.div = cell.div,
        mu_est = mu_suppose,
        simu.right.mean = right_result[1],
        right.mean = right_result[2],
        right.cd = right_result[3],
        right.p = right_result[4]
      )
    }, as.list(celldivlist))

    collect.data2 <- peak_method(p, vaf_set, ctx, p_thre, celldivlist = celldivlist, num_decimal = num_decimal, collect_data_check = collect.data1)
    range1 <- ifelse(nrow(collect.data2) == 1, collect.data2$mu_est - 1, min(collect.data2$mu_est))
    range2 <- ifelse(nrow(collect.data2) == 1, collect.data2$mu_est + 1, max(collect.data2$mu_est))
    over12 <- collect.data[collect.data$mu_real >= range1 & collect.data$mu_real <= range2, , drop = FALSE]
    overlap.pick <- if (nrow(over12) > 0) over12 else NULL

    if (nrow(df_count_freq) > 1) {
      max_mu <- df_count_freq[df_count_freq$cluster == df_count_freq[2, ]$cluster, ]$count
      min_mu <- max(1, df_count_freq[df_count_freq$cluster == 3, ]$count / 2)
      start.div <- round(total_count_temp / max_mu)
      end.div <- round(total_count_temp / min_mu)
      first_div_vaf_data <- df[df$cluster == 3, , drop = FALSE]
      sec.div.vaf <- vaf_at_div(3, p, ctx)
      first_div_vaf_data$abs_diff <- abs(first_div_vaf_data$vaf - sec.div.vaf)
      first_div_vaf_data$prob <- stats::pbeta(first_div_vaf_data$vaf, ctx$depth * sec.div.vaf, ctx$depth - ctx$depth * sec.div.vaf)
      first_div_vaf_data$score <- first_div_vaf_data$prob / max(first_div_vaf_data$prob) - first_div_vaf_data$abs_diff / max(first_div_vaf_data$abs_diff)
      right_df <- first_div_vaf_data[order(-first_div_vaf_data$score), , drop = FALSE]

      collect.data1 <- collect_rows(function(cell.div) {
        mu_suppose <- length(vaf_set) / cell.div
        k <- round(mu_suppose)
        right_result <- peak_test(right_df, min_sample_size, sec.div.vaf, k, ctx$depth, num_decimal)
        data.frame(
          cell.div = cell.div,
          mu_est = mu_suppose,
          simu.right.mean = right_result[1],
          right.mean = right_result[2],
          right.cd = right_result[3],
          right.p = right_result[4]
        )
      }, as.list(celldivlist))

      collect.data3 <- peak_method(p, vaf_set, ctx, p_thre, celldivlist = celldivlist, num_decimal = num_decimal, collect_data_check = collect.data1)
      range1 <- ifelse(nrow(collect.data3) == 1, collect.data3$mu_est - 1, min(collect.data3$mu_est))
      range2 <- ifelse(nrow(collect.data3) == 1, collect.data3$mu_est + 1, max(collect.data3$mu_est))
      over23 <- collect.data[collect.data$mu_real >= range1 & collect.data$mu_real <= range2, , drop = FALSE]
      overlap_values <- intersect(over12$mu_real, over23$mu_real)
      if (nrow(over23) > 0 || nrow(over12) > 0) {
        if (length(overlap_values) > 0) {
          overlap.pick <- over23[over23$mu_real %in% overlap_values, , drop = FALSE]
        } else {
          overlap.pick <- rbind(over12, over23)
        }
      }
    }

    if (is.null(overlap.pick)) {
      overlap.pick <- collect.data
    }
    overlap.pick <- stats::na.omit(overlap.pick)
  } else {
    overlap.pick <- collect.data
    reliable <- 0
  }

  collect.data4 <- mu_find_small(vaf_set, p, ctx)
  mu_small_selection <- collect.data4[collect.data4$p_value > 0.05 & collect.data4$wx_p > 0.05, , drop = FALSE]

  pick.cell.div <- overlap.pick$cell.div
  pick.mu <- overlap.pick$mu_real
  z_score <- overlap.pick$z_score

  overlap.pick$mu_est <- overlap.pick$mu_real
  top_pick <- as.data.frame(compare_real_simu_peak(overlap.pick, p, vaf_set, ctx, num_decimal))
  overlap.pick$loglike <- top_pick$V2
  overlap.pick$aic <- top_pick$V3
  overlap.pick$bic <- top_pick$V4
  pick.log <- overlap.pick$loglike
  pick.bic <- overlap.pick$bic
  pick.aic <- overlap.pick$aic

  if (nrow(mu_small_selection) > 0) {
    num_df <- length(vaf_set)
    mu_small_selection$cell.div <- round(num_df / mu_small_selection$mu_est)
    small_mu_pick <- as.data.frame(compare_real_simu_peak(mu_small_selection, p, vaf_set, ctx, num_decimal, bac = TRUE))
    mu_small_selection$r1 <- ifelse(small_mu_pick$V1 < 0.05, 0, small_mu_pick$V1)
    mu_small_selection$loglike <- small_mu_pick$V2
    mu_small_selection$aic <- small_mu_pick$V3
    mu_small_selection$bic <- small_mu_pick$V4
    mu_small_selection <- mu_small_selection[which(mu_small_selection$bic == min(mu_small_selection$bic)), , drop = FALSE]

    if (mu_small_selection$bic < min(overlap.pick$bic) * 0.5) {
      pick.mu <- mu_small_selection$mu_est
      pick.bic <- mu_small_selection$bic
      pick.aic <- mu_small_selection$aic
      pick.log <- mu_small_selection$loglike
      z_score <- rep(1, length(mu_small_selection$aic))
    }
  }

  data.frame(
    cell.div = pick.cell.div,
    mu = pick.mu,
    loglike = pick.log,
    bic = pick.bic,
    aic = pick.aic,
    p = p,
    z_score = z_score,
    reliable = reliable
  )
}

run_fit <- function(ctx, p_thre) {
  num_decimal <- nchar(as.character(ctx$depth))
  border_vaf <- calculate_border_vaf_fit(ctx)
  tea.result <- calculate_left_right_most_vaf_fit(border_vaf, ctx, num_decimal)
  all.p.data <- estimate_mu(
    p = tea.result$p,
    vaf_set = ctx$main_cluster_vaf,
    ctx = ctx,
    p_thre = p_thre,
    start_div = tea.result$start.div,
    num_decimal = num_decimal,
    n_wilcox = 50,
    n_sim_peak = 100000,
    tol_peak = 0
  )
  list(
    select = data.frame(mu = max(all.p.data$mu), up = max(all.p.data$lowerbound1), p1 = all.p.data$p[which.max(all.p.data$mu)], name = ctx$id),
    all = all.p.data
  )
}

run_bac <- function(ctx, p_thre) {
  num_decimal <- nchar(as.character(ctx$depth))
  second.vaf <- ctx$second_cluster_vaf$bac
  main.vaf <- ctx$main_cluster_vaf
  m <- RBesT::automixfit(second.vaf, type = "beta", Nc = 1:10, thresh = 0, k = 6, Ninit = min(50, round(length(second.vaf) / 3)), Niter.max = 10000)
  a <- m["a", ]
  b <- m["b", ]
  mean.a.b <- a / (a + b)
  left_most_vaf <- min(mean.a.b)
  right_most.vaf <- max(mean.a.b)
  right_vaf_index <- which.max(mean.a.b)
  keep.right.vaf <- right_most.vaf
  probs <- sapply(seq_along(a), function(i) stats::dbeta(second.vaf, a[i], b[i]))
  df <- data.frame(probs, vaf = second.vaf, check.names = FALSE)
  names(df)[seq_along(a)] <- paste0("prob.", seq_along(a))
  if (length(mean.a.b) > 1) {
    df <- beta_reassign(df)
    second.update.vaf <- df[df$cluster == right_vaf_index, "vaf"]
  } else {
    second.update.vaf <- df$vaf
  }

  while (length(mean.a.b) > 1) {
    possible_error <- tryCatch({
      RBesT::automixfit(second.update.vaf, type = "beta", Nc = 1:10, thresh = 0, k = 6, Ninit = min(50, round(length(second.update.vaf) / 3)), Niter.max = 10000)
    }, error = function(e) e)
    if (inherits(possible_error, "error")) break
    m <- possible_error
    a <- m["a", ]
    b <- m["b", ]
    mean.a.b <- a / (a + b)
    left_most_vaf <- min(mean.a.b)
    right_most.vaf <- max(mean.a.b)
    if (length(mean.a.b) < 2) break
    right_vaf_index <- which.max(mean.a.b)
    probs <- sapply(seq_along(a), function(i) stats::dbeta(second.update.vaf, a[i], b[i]))
    df <- data.frame(probs, vaf = second.update.vaf, check.names = FALSE)
    names(df)[seq_along(a)] <- paste0("prob.", seq_along(a))
    df <- beta_reassign(df)
    second.update.vaf <- df[df$cluster == right_vaf_index, "vaf"]
    if (length(second.update.vaf) == 0 || length(unique(second.update.vaf)) == 1) break
  }

  right_most.vaf <- max(right_most.vaf, keep.right.vaf)
  p <- (right_most.vaf * 2 * exp(log(2) * ctx$beta) - 1) / (exp(log(2) * ctx$beta) - 1)
  if (p < 0) stop("bac: non-physical p < 0")

  righta <- c(0.5, right_most.vaf)
  df <- beta_reassign(.vaf_prob_df(main.vaf, righta, ctx$depth))
  vaf_set <- c(df[df$cluster == 2, "vaf"], second.update.vaf)
  pick_mu_cell.div <- peak_method(p, vaf_set, ctx, p_thre, celldivlist = NULL, num_decimal = num_decimal, n_sim_peak = 100000, n_wilcox = 50, tol_peak = 0)

  vaf_set <- c(df[df$cluster == 2, "vaf"], second.vaf)
  result_df <- mu_find_small(vaf_set, p, ctx)
  mu_small_selection <- result_df[result_df$p_value > 0.05 & result_df$wx_p > 0.05, , drop = FALSE]
  pick.mu <- pick_mu_cell.div$mu_est
  pick.bic <- pick_mu_cell.div$bic
  pick.aic <- pick_mu_cell.div$aic
  pick.log <- pick_mu_cell.div$loglike

  if (nrow(mu_small_selection) > 0) {
    num_df <- length(vaf_set)
    mu_small_selection$cell.div <- round(num_df / mu_small_selection$mu_est)
    small_mu_pick <- as.data.frame(compare_real_simu_peak(mu_small_selection, p, vaf_set, ctx, num_decimal, bac = TRUE))
    mu_small_selection$r1 <- ifelse(small_mu_pick$V1 < 0.05, 0, small_mu_pick$V1)
    mu_small_selection$loglike <- small_mu_pick$V2
    mu_small_selection$aic <- small_mu_pick$V3
    mu_small_selection$bic <- small_mu_pick$V4
    mu_small_selection <- mu_small_selection[which(mu_small_selection$bic == min(mu_small_selection$bic)), , drop = FALSE]
    if (mu_small_selection$bic < min(pick_mu_cell.div$bic) * 0.5) {
      pick.mu <- mu_small_selection$mu_est
      pick.bic <- mu_small_selection$bic
      pick.aic <- mu_small_selection$aic
      pick.log <- mu_small_selection$loglike
    }
  }

  bac_pick <- data.frame(cell.div = pick.mu, mu = pick.mu, loglike = pick.log, bic = pick.bic, aic = pick.aic, p = p)
  list(
    select = data.frame(mu = max(bac_pick$mu), up = length(ctx$second_cluster_vaf$bac) / 3, p1 = bac_pick$p[which.max(bac_pick$mu)], name = ctx$id),
    all = bac_pick
  )
}

run_normal <- function(ctx, p_thre) {
  num_decimal <- nchar(as.character(ctx$depth))
  vafdata.summary <- ctx$vafdata_summary
  vafdata.summary.filter <- ctx$vafdata_summary_filter
  vafdata <- ctx$vafdata
  second.vaf <- ctx$second_cluster_vaf$normal
  color <- vafdata.summary$colors[
    vafdata.summary$max %in% vafdata.summary.filter[1, "max"] |
      vafdata.summary$min %in% vafdata.summary.filter[1, "min"]
  ]
  clonal.vaf <- vafdata[vafdata$colors %in% color, "vaf.1"]
  m <- RBesT::automixfit(clonal.vaf, type = "beta", Nc = 1:10, thresh = 0, k = 6, Ninit = min(50, round(length(clonal.vaf) / 3)), Niter.max = 10000)
  a <- m["a", ]
  b <- m["b", ]
  mean.a.b <- a / (a + b)
  upper_clonal_vaf <- min(mean.a.b)
  approx_vaf_index <- which.min(mean.a.b)
  probs <- sapply(seq_along(a), function(i) stats::dbeta(clonal.vaf, a[i], b[i]))
  df <- data.frame(probs, vaf = clonal.vaf, check.names = FALSE)
  names(df)[seq_along(a)] <- paste0("prob.", seq_along(a))
  df <- beta_reassign(df)
  clonal.vaf.left <- df[df$cluster == approx_vaf_index, "vaf"]
  second.try <- nrow(vafdata.summary.filter) > 1
  if (!second.try) {
    temp_vaf <- (upper_clonal_vaf + 0.5) / 2
    second.vaf <- clonal.vaf[clonal.vaf < temp_vaf]
  }
  seq_data <- find_p_process(second.vaf, mean.a.b, upper_clonal_vaf, clonal.vaf.left, ctx)
  my_seq <- seq_data$my_seq
  clear <- seq_data$clear
  if (length(my_seq) <= 1 && second.try) {
    temp_vaf <- (upper_clonal_vaf + 0.5) / 2
    second.vaf <- clonal.vaf[clonal.vaf < temp_vaf]
    seq_data <- find_p_process(second.vaf, mean.a.b, upper_clonal_vaf, clonal.vaf.left, ctx)
    my_seq <- seq_data$my_seq
    clear <- seq_data$clear
  }

  all.p.data.final <- data.frame(
    cell.div = NA, mu = NA, loglike = NA, bic = NA, aic = NA, p = NA,
    z_score = NA, reliable = NA, z_score1 = NA, bicrank = NA, score = 0, clear = NA
  )
  if (length(my_seq) > 0) {
    all.p.data <- collect_rows(function(give.vaf) {
      res <- tryCatch(
        .ok(iterate_p_optimize(clear, give.vaf, upper_clonal_vaf, clonal.vaf.left, second.vaf, num_decimal, p_thre, ctx)),
        error = .err
      )
      if (res$status == "ok") res$value else NULL
    }, as.list(my_seq))

    if (nrow(all.p.data) > 0) {
      all.p.data$mu <- round(all.p.data$mu, num_decimal)
      if (nrow(all.p.data[abs(all.p.data$z_score) <= 1.96, , drop = FALSE]) >= 2) {
        all.p.data <- all.p.data[abs(all.p.data$z_score) <= 1.96, , drop = FALSE]
      }
      all.p.data$z_score1 <- 1
      all.p.data$z_score <- abs(all.p.data$z_score)
      all.p.data$bicrank <- rank(all.p.data$bic)
      all.p.data$score <- 1 / all.p.data$bicrank + (1 - all.p.data$z_score) * 0.25
      all.p.data$clear <- clear
      all.p.data.final <- all.p.data
    }
  }

  list(
    select = data.frame(
      mu = ifelse(nrow(all.p.data.final) > 0, all.p.data.final$mu[which.min(all.p.data.final$bic)], NA),
      up = NA,
      p1 = ifelse(nrow(all.p.data.final) > 0, all.p.data.final$p[which.min(all.p.data.final$bic)], NA),
      name = ctx$id
    ),
    all = all.p.data.final
  )
}

run_estimates <- function(ctx, p_thre = 0.01) {
  estimator_names <- list_estimators()
  preferred_order <- c("fit", "bac", "normal")
  estimator_names <- c(preferred_order[preferred_order %in% estimator_names], setdiff(estimator_names, preferred_order))
  sample_name <- ctx$id

  # Two modes. DEFAULT (ctx$fast_version FALSE): SEQUENTIAL in v1
  # Run.para.estimate.maincluster order (estimator-outer fit,bac,normal /
  # try-inner 1:3) in one continuous RNG stream -> seeded run faithful to
  # v1. fast_version TRUE: dispatch the 9 independent units with mclapply
  # (each fork its own RNG) -> ~2-3x faster, results in v1's distribution
  # but NOT bit-identical (seeded approximation).
  jobs <- expand.grid(try_idx = 1:3, estimator = estimator_names,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  n_jobs <- nrow(jobs)
  run_one <- function(k) {
    estimator_fn <- get_estimator(jobs$estimator[k])
    res <- tryCatch(.ok(estimator_fn(ctx = ctx, p_thre = p_thre)), error = .err)
    if (res$status == "ok") res$value else NULL
  }
  if (isTRUE(ctx$fast_version) && .Platform$OS.type == "unix" && n_jobs > 1L) {
    mc <- getOption("teatime.mc.cores", max(1L, parallel::detectCores(logical = TRUE) - 1L))
    mc <- suppressWarnings(as.integer(mc))
    if (is.na(mc) || mc < 1L) mc <- 1L
    job_out <- parallel::mclapply(seq_len(n_jobs), run_one,
                                  mc.cores = min(n_jobs, mc), mc.set.seed = TRUE)
  } else {
    job_out <- vector("list", n_jobs)
    for (k in seq_len(n_jobs)) job_out[[k]] <- run_one(k)
  }

  # Per-estimator aggregation: identical to the original sequential worker
  # body, fed the 3 try results in try order (1,2,3).
  aggregate_est <- function(estimator_name) {
    all_data <- NULL
    mu_list <- c()
    up_list <- c()
    p_list <- c()
    for (try_idx in 1:3) {
      k <- which(jobs$estimator == estimator_name & jobs$try_idx == try_idx)
      case_result <- job_out[[k]]
      if (is.null(case_result)) next
      case_all <- case_result$all
      if (!is.null(case_all) && nrow(case_all) > 0) {
        if (identical(estimator_name, "normal")) {
          case_all <- case_all[order(-case_all$score), , drop = FALSE]
        } else {
          case_all <- case_all[order(case_all$bic), , drop = FALSE]
        }
        case_all$try <- try_idx
        case_all$name <- sample_name
        all_data <- if (is.null(all_data)) case_all else rbind(all_data, case_all)
        if (identical(estimator_name, "fit")) {
          filtered <- filter_mu_estimate(case_all)
          mu_list <- c(mu_list, filtered$mu)
          up_list <- c(up_list, filtered$lowerbound1)
          p_list <- c(p_list, filtered$p)
        } else if (identical(estimator_name, "bac")) {
          filtered <- filter_mu_estimate(case_all)
          mu_list <- c(mu_list, filtered$mu)
          p_list <- c(p_list, filtered$p)
        } else {
          mu_list <- c(mu_list, case_all$mu)
          p_list <- c(p_list, case_all$p)
        }
      }
    }
    select_df <- data.frame(
      mu = ifelse(length(mu_list) > 0, max(mu_list), NA),
      up = ifelse(identical(estimator_name, "bac"), length(ctx$second_cluster_vaf$bac) / 3, ifelse(length(up_list) > 0, max(up_list), NA)),
      p1 = ifelse(length(mu_list) > 0, p_list[which.max(mu_list)], NA),
      name = sample_name
    )
    if (!is.null(all_data) && nrow(all_data) > 1) {
      all_data <- all_data[!is.na(all_data$mu), , drop = FALSE]
    }
    list(select = select_df, all = all_data)
  }

  results <- lapply(estimator_names, aggregate_est)
  stats::setNames(results, estimator_names)
}

keep_longest_consecutive_rows <- function(df) {
  vec <- df$start
  d <- c(NA, diff(vec))
  breaks <- which(d > 1)
  start <- c(1, breaks)
  end <- c(breaks - 1, length(vec))
  lengths <- end - start + 1
  longest_seq <- which.max(lengths)
  df[start[longest_seq]:end[longest_seq], , drop = FALSE]
}

pick_s <- function(df) {
  df$bicrank <- rank(df$bic)
  df$score <- df$vaf
  if (max(df$start) > 1) {
    df <- keep_longest_consecutive_rows(df)
  }

  pickdata <- df[which(df$bic == min(df$bic)), , drop = FALSE]
  if (max(df$start) >= 1) {
    data.pick.sub <- df
    data.pick.sub$truescore <- 0:(nrow(data.pick.sub) - 1)
    valid_rows <- which(data.pick.sub$truescore <= data.pick.sub$start & data.pick.sub$start > 0)
    if (length(valid_rows) > 0) {
      data.pick.sub <- data.pick.sub[valid_rows, , drop = FALSE]
      pickdata <- df[which(df$bic == min(data.pick.sub$bic)), , drop = FALSE]
    }
  }
  pickdata
}

generate_sequence <- function(n, s, set_max = FALSE, set_max_num = 0) {
  cell.list <- c()
  max_n <- n + 2
  if (set_max) {
    max_n <- set_max_num
  }
  current_val <- 0
  for (i in 1:max_n) {
    current_val <- if (i %% n == 0) (1 + s) * (i / n) else current_val + 1
    cell.list[i] <- current_val
  }
  cell.list
}

s_dataframe_update <- function(cell.list, p, vaf.t1, vaf_set, min.s.detect, ctx, evaluate = FALSE, simulation = FALSE) {
  result_vector <- sapply(cell.list, function(i) p / (2 * exp(log(2) * ctx$beta * i)))
  mix_check <- sapply(result_vector, function(i) ifelse(log((1 - p) / 2 / i) / (log(2) * ctx$beta) > 0, 1, 0))
  cluster_idx <- seq_along(mix_check)
  named_mix_vec <- stats::setNames(mix_check, cluster_idx)

  inita <- ctx$depth * c(vaf.t1, result_vector)
  initb <- ctx$depth - inita
  probs <- sapply(seq_along(inita), function(i) stats::dbeta(vaf_set, inita[i], initb[i]))
  df <- data.frame(prob = probs, vaf = vaf_set)
  df <- beta_reassign(df)

  if (evaluate) {
    p_vec <- c(vaf.t1, result_vector)
    df <- dplyr::mutate(df, p = p_vec[cluster])
  }
  df <- df[df$cluster > 1, , drop = FALSE]

  sample_or_not <- function(data, cluster_id, fraction, mix_vec) {
    if (mix_vec[cluster_id - 1] == 1) {
      dplyr::sample_frac(data, size = fraction)
    } else {
      data
    }
  }

  sampled_df <- df %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_modify(~ sample_or_not(.x, .y$cluster, p, named_mix_vec)) %>%
    dplyr::ungroup()

  if (simulation) {
    cluster.result <- sampled_df %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        count = dplyr::n(),
        mean_vaf = mean(vaf),
        min_vaf = min(vaf),
        .groups = "drop"
      )
    return(cluster.result)
  }

  if (evaluate) {
    sampled_df <- sampled_df %>%
      dplyr::mutate(
        log_likelihood = likelihoodExplore::likbeta(
          x = vaf,
          shape1 = p * ctx$depth,
          shape2 = ctx$depth - p
        )
      )
    cluster.result <- sampled_df %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        count = dplyr::n(),
        mean_vaf = mean(vaf),
        .groups = "drop"
      )
    cluster.result$realvaf <- result_vector[cluster.result$cluster - 1]
    return(list(sampled_df = sampled_df, cluster.result = cluster.result))
  }

  sampled_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      count = dplyr::n(),
      mean_vaf = mean(vaf),
      min_vaf = min(vaf),
      .groups = "drop"
    ) %>%
    dplyr::filter(cluster > 1) %>%
    dplyr::mutate(
      new_s = log(p / (mean_vaf * 2)) / (log(2) * ctx$beta) - 1,
      alt_s = log(p / (min_vaf * 2)) / (log(2) * ctx$beta) - 1,
      new_s = ifelse(new_s < 0, alt_s, new_s)
    ) %>%
    dplyr::filter(new_s <= min.s.detect) %>%
    dplyr::mutate(cluster = cluster - 1)
}

s_update_process <- function(give_n, cluster.result, vaf_set, min.s.detect, vaf.t1, p, ctx, p_thre = 1e-6) {
  s <- cluster.result$new_s[which(cluster.result$cluster == give_n)]
  i <- 1
  while (i <= 100) {
    max_n <- give_n + 2
    if (give_n > 1) {
      cell.list <- generate_sequence(give_n, s)
    } else {
      cell.list <- (1:max_n) * (1 + s)
    }
    cluster.result <- s_dataframe_update(cell.list, p, vaf.t1, vaf_set, min.s.detect, ctx)
    update_s <- cluster.result$new_s[which(cluster.result$cluster == give_n)]
    s <- update_s
    if (abs(update_s - s) < p_thre) break
    i <- i + 1
  }
  update_s
}

evaluate_all_s <- function(svalue.list, vaf.t1, p, ctx) {
  result <- ctx$vafdata
  end.vaf <- p / (2 * exp(log(2) * ctx$beta * (1 + max(svalue.list))))
  vaf_set <- result$vaf.1[which(result$vaf.1 < vaf.t1 & result$vaf.1 > end.vaf)]
  bic.list <- c()
  vaf.list <- c()
  start.list <- c()
  past_diff <- data.frame()
  s.total.num <- length(svalue.list)

  for (i in seq_along(svalue.list)) {
    s <- if (svalue.list[i] > 0) svalue.list[i] else 0
    give_n <- i
    if (give_n > 1) {
      cell.list <- generate_sequence(give_n, s, set_max = TRUE, set_max_num = s.total.num)
    } else {
      cell.list <- (1:s.total.num) * (1 + s)
    }
    result_list <- s_dataframe_update(cell.list = cell.list, p = p, vaf.t1 = vaf.t1, vaf_set = vaf_set, min.s.detect = NA, ctx = ctx, evaluate = TRUE, simulation = FALSE)
    cluster.result <- result_list$cluster.result
    cluster.result$try <- i
    sampled_df <- result_list$sampled_df
    total_log_likelihood <- sum(sampled_df$log_likelihood)
    bic <- compute_BIC(total_log_likelihood, 1, length(vaf_set))
    bic.list <- c(bic.list, bic)
    vaf.list <- c(vaf.list, sum(abs(cluster.result$realvaf - cluster.result$mean_vaf)) / nrow(cluster.result))

    if (i > 1) {
      down <- 0
      for (iter in 1:(i - 1)) {
        past_diff_pick <- past_diff[past_diff$try == iter, , drop = FALSE]
        diff_past <- past_diff_pick[1:iter, , drop = FALSE]$realvaf - past_diff_pick[1:iter, , drop = FALSE]$mean_vaf
        diff_cluster <- cluster.result[1:iter, , drop = FALSE]$realvaf - cluster.result[1:iter, , drop = FALSE]$mean_vaf
        binary_result <- ifelse(diff_past - diff_cluster >= 0, 1, 0)
        down <- down + sum(binary_result)
      }
      start.list <- c(start.list, down)
      past_diff <- rbind(past_diff, cluster.result)
    } else {
      past_diff <- rbind(past_diff, cluster.result)
      start.list <- c(start.list, 0)
    }
  }

  data.frame(
    s = svalue.list,
    bic = bic.list,
    vaf = vaf.list,
    start = start.list
  )
}

get_second_peak_ratio <- function(simu_vaf, vaf_list, cell_list, ctx, pick.ratio = 1) {
  inita <- ctx$depth * vaf_list
  initb <- ctx$depth - inita
  probs <- sapply(seq_along(inita), function(i) stats::dbeta(simu_vaf, inita[i], initb[i]))
  df <- data.frame(prob = probs, vaf = simu_vaf)
  df <- beta_reassign(df)
  df <- df[df$cluster > 1, , drop = FALSE]
  cluster.result <- df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      count = dplyr::n(),
      mean_vaf = mean(vaf),
      min_vaf = min(vaf),
      .groups = "drop"
    )
  as.numeric(cluster.result[cluster.result$cluster == 2, "count", drop = TRUE]) / (cell_list[2] * 100 * pick.ratio)
}

simulate_peak <- function(vaf_list, cell_list, vaf_set, ctx) {
  num_decimal <- nchar(as.character(ctx$depth))
  iterations <- 10
  ratio <- 0
  lowerratio <- 0
  vaf.div <- seq_along(vaf_list)
  for (i in 1:iterations) {
    all_simulated_vafs_list <- sapply(vaf.div, function(idx) {
      vaf <- vaf_list[idx]
      size <- cell_list[idx] * 100
      round(stats::rbinom(size, ctx$depth, vaf) / ctx$depth, num_decimal)
    })
    all_simulated_vafs <- unlist(all_simulated_vafs_list)
    all_simulated_vafs <- all_simulated_vafs[all_simulated_vafs > min(vaf_set)]
    ratio <- ratio + get_second_peak_ratio(all_simulated_vafs, vaf_list, cell_list, ctx, pick.ratio = 1)

    sampled_vafs_list <- lapply(vaf.div, function(idx) {
      vaf_vec <- all_simulated_vafs_list[idx]
      sample_size <- ceiling(length(vaf_vec) / cell_list[idx])
      sample(vaf_vec, size = sample_size)
    })
    sampled_vafs_vector <- unlist(sampled_vafs_list)
    sampled_vafs_vector <- sampled_vafs_vector[sampled_vafs_vector > min(vaf_set)]
    lowerratio <- lowerratio + get_second_peak_ratio(sampled_vafs_vector, vaf_list, cell_list, ctx, pick.ratio = 1 / cell_list[2])
  }
  as.numeric(c(ratio / iterations, lowerratio / iterations))
}

simulate_ratio_peak <- function(p, s, ctx) {
  result <- ctx$vafdata
  vaf.t1 <- p / 2
  give_n <- floor(1 + s)
  if (give_n > 1) {
    cell.list <- generate_sequence(give_n, s, set_max = TRUE, set_max_num = 3)
  } else {
    cell.list <- (1:3) * (1 + s)
  }
  vaf_set <- result$vaf.1[which(result$vaf.1 < vaf.t1)]
  result_vector <- sapply(cell.list, function(i) p / (2 * exp(log(2) * ctx$beta * i)))
  cluster.result <- s_dataframe_update(cell.list = cell.list, p = p, vaf.t1 = vaf.t1, vaf_set = vaf_set, min.s.detect = NA, ctx = ctx, evaluate = FALSE, simulation = TRUE)
  vaf_list <- c(vaf.t1, result_vector)
  cell_list <- c(1, exp(log(2) * ctx$beta * cell.list))
  simu_ratio <- simulate_peak(vaf_list, cell_list, vaf_set, ctx)
  mu.range <- cluster.result[cluster.result$cluster == 2, ]$count / simu_ratio / exp(log(2) * ctx$beta * cell.list[1])
  mu.range
}

get_s <- function(p, ctx, p_thre = 1e-6) {
  result <- ctx$vafdata
  vaf.t1 <- p / 2
  min.s.detect <- log(p / (stats::quantile(result$vaf.1, 0.01) * 2)) / (log(2) * ctx$beta) - 1
  n <- max(2, floor(1 + min.s.detect))
  cell.list <- 1:n
  vaf_set <- result$vaf.1[which(result$vaf.1 < vaf.t1 & result$vaf.1 > stats::quantile(result$vaf.1, 0.05))]
  if (length(vaf_set) == 0) {
    vaf_set <- result$vaf.1[which(result$vaf.1 < vaf.t1)]
  }

  estimate.s.data <- data.frame(s = NA, bic = NA, vaf = NA, start = NA)
  if (length(vaf_set) > 0) {
    cluster.result <- s_dataframe_update(cell.list, p, vaf.t1, vaf_set, min.s.detect, ctx)
    check_div <- if (nrow(cluster.result) == 1) min(cluster.result$cluster) else 1:nrow(cluster.result)
    svalue.list <- sapply(check_div, function(give_n) {
      s_update_process(give_n, cluster.result, vaf_set, min.s.detect, vaf.t1, p, ctx, p_thre)
    })
    if (any(svalue.list > 0)) {
      estimate.s.data <- evaluate_all_s(svalue.list, vaf.t1, p, ctx)
      estimate.s.data <- estimate.s.data[estimate.s.data$s > 0, , drop = FALSE]
    }
  }
  estimate.s.data
}

find_s_from_predict <- function(predict.result, ctx) {
  p.list <- unique(predict.result$p)
  s.results <- lapply(p.list, function(p) {
    df <- get_s(p, ctx, p_thre = 1e-6)
    df$p_value <- p
    data <- pick_s(df)
    murange <- simulate_ratio_peak(p, data$s, ctx)
    data$minmu <- min(murange)
    data$maxmu <- max(murange)
    data
  })
  fit.data <- data.frame()
  combined_df <- do.call(rbind, s.results)
  for (try_idx in seq_along(p.list)) {
    p.pick <- p.list[try_idx]
    fit.check <- predict.result[which(predict.result$p == p.pick), , drop = FALSE]
    if (nrow(fit.check) > 0) {
      minmu <- combined_df[which(combined_df$p_value == p.pick), , drop = FALSE]$minmu
      maxmu <- combined_df[which(combined_df$p_value == p.pick), , drop = FALSE]$maxmu
      fit.check.sub <- fit.check[which(fit.check$mu >= minmu & fit.check$mu <= maxmu), , drop = FALSE]
      fit.check.sub <- stats::na.omit(fit.check.sub)
      if (nrow(fit.check.sub) == 0) {
        minmu <- minmu - 1.96 * 1.5
        fit.check.sub <- fit.check[fit.check$mu >= minmu & fit.check$mu <= maxmu, , drop = FALSE]
        fit.check.sub <- stats::na.omit(fit.check.sub)
      }
      if (nrow(fit.check.sub) > 0) {
        best_fit <- fit.check.sub[which(fit.check.sub$bic == min(fit.check.sub$bic)), , drop = FALSE]
        fit.data <- rbind(fit.data, best_fit)
      }
    }
  }

  if (nrow(fit.data) > 0) {
    min_bic <- min(fit.data$bic)
    min_bic_indices <- which(fit.data$bic == min_bic)
    selected_index <- sample(min_bic_indices, 1)
    selected_row <- fit.data[selected_index, , drop = FALSE]
    selected_row$s <- combined_df[which(combined_df$p_value == selected_row$p), , drop = FALSE]$s
    return(selected_row)
  }
  combined_df
}

run_fitness <- function(estimates, ctx) {
  fit.all <- estimates$fit$all
  inter.all <- estimates$normal$all
  fit.select.one <- estimates$fit$select
  inter.select.one <- estimates$normal$select
  bac.select.one <- estimates$bac$select

  fitmu <- NA
  fitmu_candidate <- 1
  fitcell <- NA
  fits <- NA
  fitp <- NA
  if (!is.null(fit.all) && nrow(fit.all[!is.na(fit.all$cell.div), , drop = FALSE]) > 0) {
    fit.all <- fit.all[!is.na(fit.all$cell.div), , drop = FALSE]
    fit.data <- find_s_from_predict(fit.all, ctx)
    if (!("mu" %in% colnames(fit.data))) {
      mu.fit.pick <- fit.select.one$up
      fit.select.row <- fit.all[fit.all$lowerbound1 == mu.fit.pick, , drop = FALSE]
      fit.select.row <- fit.select.row[sample(nrow(fit.select.row), 1), , drop = FALSE]
      p <- fit.select.row$p
      mu.s.update <- fit.data[fit.data$p_value == p, , drop = FALSE]
      mu.select2 <- mu.s.update[sample(nrow(mu.s.update), 1), , drop = FALSE]$minmu
      df <- get_s(p, ctx, p_thre = 1e-6)
      df$p_value <- p
      data <- pick_s(df)
      fitmu <- mu.fit.pick
      fitmu_candidate <- mu.select2
      fitcell <- fit.select.row$cell.div
      fits <- data$s
      fitp <- p
    } else {
      fitmu <- fit.data$mu
      fitmu_candidate <- 1
      fitcell <- fit.data$cell.div
      fits <- fit.data$s
      fitp <- fit.data$p
    }
  }

  intermu <- NA
  intermu_candidate <- NA
  intertrust <- NA
  intercell <- NA
  inters <- NA
  interp <- NA
  if (!is.null(inter.all) && nrow(inter.all[!is.na(inter.all$cell.div), , drop = FALSE]) > 0) {
    inter.all <- inter.all[!is.na(inter.all$cell.div), , drop = FALSE]
    inter.data <- find_s_from_predict(inter.all, ctx)
    if (!("mu" %in% colnames(inter.data))) {
      mu.inter.pick <- inter.select.one$mu
      inter.select.row <- inter.all[inter.all$mu == mu.inter.pick, , drop = FALSE]
      inter.select.row <- inter.select.row[sample(nrow(inter.select.row), 1), , drop = FALSE]
      p <- inter.select.row$p
      mu.s.update <- inter.data[inter.data$p_value == p, , drop = FALSE]
      mu.select2 <- mu.s.update[sample(nrow(mu.s.update), 1), , drop = FALSE]$minmu
      df <- get_s(p, ctx, p_thre = 1e-6)
      df$p_value <- p
      data <- pick_s(df)
      intermu <- mu.select2
      intermu_candidate <- mu.inter.pick
      intertrust <- (2 * abs(mu.select2 - mu.inter.pick)) / (mu.inter.pick + mu.select2)
      intercell <- inter.select.row$cell.div
      inters <- data$s
      interp <- p
    } else {
      intermu <- inter.data$mu
      intermu_candidate <- 1
      intertrust <- 0
      intercell <- inter.data$cell.div
      inters <- inter.data$s
      interp <- inter.data$p
    }
  }

  data.frame(
    name = ctx$id,
    fitmu = fitmu,
    fitmu_candidate = fitmu_candidate,
    fitcell = fitcell,
    fits = fits,
    fitp = fitp,
    intermu = intermu,
    intermu_candidate = intermu_candidate,
    intertrust = intertrust,
    intercell = intercell,
    inters = inters,
    interp = interp,
    backp = bac.select.one$p1,
    stringsAsFactors = FALSE
  )
}

