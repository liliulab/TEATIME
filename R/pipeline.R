extract_subp <- function(result) {
  result$count <- 1
  meanvaf <- aggregate(result$vaf.1, list(result$colors), mean)
  maxvaf <- aggregate(result$vaf.1, list(result$colors), max)
  minvaf <- aggregate(result$vaf.1, list(result$colors), min)
  sumvaf <- aggregate(result$count, list(result$colors), sum)
  mago.result <- data.frame(max = maxvaf$x, min = minvaf$x, vaf = meanvaf$x, sum = sumvaf$x)

  max_vaf_index <- which.max(mago.result$vaf)
  mago.result <- mago.result[-max_vaf_index, ]
  min_min_index <- which.min(mago.result$min)
  mago.result <- mago.result[-min_min_index, ]
  min_min_index <- which.max(mago.result$min)
  mago.result <- mago.result[min_min_index, ]
  filtered_result <- mago.result[mago.result$sum >= 0, ]
  if (nrow(filtered_result) > 0) {
    total_sum <- sum(filtered_result$sum)
    weighted_mean_vaf <- sum(filtered_result$vaf * filtered_result$sum) / total_sum
    weighted_mean_vaf * 2
  } else {
    0
  }
}

adapt_magos <- function(input) {
  purity <- min(input$purity, 1)
  vafdata <- input$result
  vafdata$vaf.1 <- vafdata$vaf.1 * (2 - purity) / (2 * vafdata$vaf.1 * (1 - purity) + purity)
  vafdata[, c("vaf.1", "depth.1", "colors")]
}

adapt_vcf <- function(input) {
  input[, ncol(input)] <- as.numeric(input[, ncol(input)])
  input <- input[input[, ncol(input)] == 2, , drop = FALSE]
  input <- input[, -ncol(input), drop = FALSE]
  mag <- mag.single.run(input, fold = TRUE)
  adapt_magos(list(purity = mag$purity, result = mag$results))
}

adapt_raw <- function(input) {
  required <- c("vaf.1", "depth.1", "colors")
  if (!all(required %in% names(input))) {
    stop("raw input must contain columns: vaf.1, depth.1, colors")
  }
  input[, required]
}

prepare_data <- function(
  input,
  beta,
  depth = NA,
  input_format = "vcf",
  growth_model = "exponential",
  verbose = FALSE,
  output_folder = "./",
  output_prefix = "TEATIME",
  id = "T01",
  write_final = TRUE,
  score_method = "bic",
  extra = list(),
  seed = NA,
  save_magos = FALSE,
  fast_version = FALSE
) {
  # `save_magos` is only valid for input_format == "vcf" (the only mode
  # where TEATIME runs MAGOS internally). Warn loudly if a user mistakenly
  # passes save_magos with magos / raw input -- those modes have no
  # MAGOS object to save.
  if (!isFALSE(save_magos) && !identical(input_format, "vcf")) {
    warning(sprintf("`save_magos` is only honored when input_format = 'vcf' (got '%s'); ignoring.",
                    input_format))
    save_magos <- FALSE
  }

  # When save_magos is set in vcf mode, intercept the adapter so we can write
  # the intermediate MAGOS clustering to disk. TRUE -> default path
  # <output_folder>/<output_prefix>_MAGOS.rds.  A string -> custom path.
  if (identical(input_format, "vcf") && !isFALSE(save_magos)) {
    input2 <- input
    input2[, ncol(input2)] <- as.numeric(input2[, ncol(input2)])
    input2 <- input2[input2[, ncol(input2)] == 2, , drop = FALSE]
    input2 <- input2[, -ncol(input2), drop = FALSE]
    mag <- mag.single.run(input2, fold = TRUE)
    out_path <- if (is.character(save_magos)) save_magos
                else file.path(output_folder, paste0(output_prefix, "_MAGOS.rds"))
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(mag, out_path)
    cat(sprintf("[TEATIME] saved MAGOS clustering -> %s\n", out_path)); flush.console()
    standardized <- adapt_magos(list(purity = mag$purity, result = mag$results))
  } else {
    adapter <- get_adapter(input_format)
    standardized <- adapter(input)
  }

  if (is.na(depth)) {
    depth <- round(mean(standardized$depth.1))
  }

  magosp <- extract_subp(standardized)

  vafdata.summary <- standardized %>%
    dplyr::group_by(colors) %>%
    dplyr::summarise(
      max = max(vaf.1),
      min = min(vaf.1),
      vaf = mean(vaf.1),
      sum = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(count = 1)

  vafdata.summary <- vafdata.summary[order(vafdata.summary$max, decreasing = TRUE), , drop = FALSE]
  vafdata.summary.filter <- vafdata.summary[which(vafdata.summary$max > 0.25 | vafdata.summary$min > 0.25), , drop = FALSE]
  if (nrow(vafdata.summary.filter) > 1) {
    for (i in 2:nrow(vafdata.summary.filter)) {
      if (vafdata.summary.filter$sum[i] <= 6) {
        vafdata.summary.filter$min[i - 1] <- vafdata.summary.filter$min[i]
        vafdata.summary.filter$count[i] <- 0
      }
    }
  }
  vafdata.summary.filter <- vafdata.summary.filter[vafdata.summary.filter$count > 0, , drop = FALSE]
  vafdata.summary.filter <- vafdata.summary.filter[order(vafdata.summary.filter$max, decreasing = TRUE), , drop = FALSE]

  max_color <- vafdata.summary$colors[which.max(vafdata.summary$max)]
  main.cluster.vaf <- standardized$vaf.1[which(standardized$colors == max_color)]
  second.cluster.vaf <- list()

  if (nrow(vafdata.summary.filter) > 1) {
    second_highest_value <- sort(vafdata.summary.filter$max, decreasing = TRUE)[2]
    second_max_color <- vafdata.summary.filter$colors[vafdata.summary.filter$max %in% second_highest_value]
    second.cluster.vaf$fit <- standardized$vaf.1[which(standardized$colors %in% second_max_color)]
    second.cluster.vaf$bac <- standardized$vaf.1[which(standardized$colors %in% second_max_color)]
    second.cluster.vaf$normal <- standardized$vaf.1[which(standardized$colors %in% second_max_color)]
  } else {
    second.cluster.vaf$fit <- NULL
    vafdata.summary.filter.bac <- rbind(vafdata.summary.filter, vafdata.summary[2, , drop = FALSE])
    second_highest_value <- sort(vafdata.summary.filter.bac$max, decreasing = TRUE)[2]
    second_max_color <- vafdata.summary.filter.bac$colors[vafdata.summary.filter.bac$max %in% second_highest_value]
    second.cluster.vaf$bac <- standardized$vaf.1[which(standardized$colors %in% second_max_color)]

    second_highest_value <- sort(vafdata.summary$max, decreasing = TRUE)[2]
    second_max_color <- vafdata.summary$colors[vafdata.summary$max %in% second_highest_value]
    second.cluster.vaf$normal <- standardized$vaf.1[which(standardized$colors %in% second_max_color)]
  }

  ctx <- list(
    depth = depth,
    beta = beta,
    vafdata = standardized,
    main_cluster_vaf = main.cluster.vaf,
    second_cluster_vaf = second.cluster.vaf,
    vafdata_summary = vafdata.summary,
    vafdata_summary_filter = vafdata.summary.filter,
    magosp = magosp,
    id = id,
    output_folder = output_folder,
    output_prefix = output_prefix,
    write_final = write_final,
    verbose = verbose,
    growth_model = growth_model,
    score_method = score_method,
    extra = extra,
    seed = seed,
    fast_version = fast_version
  )

  required <- c(
    "depth", "beta", "vafdata", "main_cluster_vaf", "second_cluster_vaf",
    "vafdata_summary", "vafdata_summary_filter", "magosp", "id",
    "output_folder", "output_prefix", "write_final", "verbose",
    "growth_model", "score_method", "extra", "seed", "fast_version"
  )
  missing <- setdiff(required, names(ctx))
  if (length(missing) > 0) {
    stop("prepare_data missing required fields: ", paste(missing, collapse = ", "))
  }

  ctx
}

run_rbest <- function(ctx) {
  main.vaf <- ctx$main_cluster_vaf
  possible_error <- tryCatch({
    RBesT::automixfit(
      main.vaf,
      type = "beta",
      Nc = 1:10,
      thresh = 0,
      k = 6,
      Ninit = min(50, round(length(main.vaf) / 5)),
      Niter.max = 10000
    )
  }, error = function(e) e)

  if (!inherits(possible_error, "error")) {
    m <- possible_error
    a <- m["a", ]
    b <- m["b", ]
    mean.a.b <- a / (a + b)
    minvaf <- min(mean.a.b)
    maxvaf <- max(mean.a.b)
    mean.a.b.sub <- mean.a.b[abs(mean.a.b - 0.5) > 0.005]
    len_2 <- length(mean.a.b)
    len <- length(mean.a.b.sub) + 1
  } else {
    minvaf <- mean(main.vaf)
    maxvaf <- mean(main.vaf)
    len_2 <- 1
    len <- 1
  }

  adjustment <- 0.5 - mean(main.vaf)
  main.vaf.adj <- main.vaf + adjustment
  possible_error <- tryCatch({
    RBesT::automixfit(
      main.vaf.adj,
      type = "beta",
      Nc = 1:10,
      thresh = 0,
      k = 6,
      Ninit = min(50, round(length(main.vaf.adj) / 5)),
      Niter.max = 10000
    )
  }, error = function(e) e)

  if (!inherits(possible_error, "error")) {
    m <- possible_error
    a <- m["a", ]
    b <- m["b", ]
    mean.a.b <- a / (a + b)
    minvaf_adj <- min(mean.a.b)
    maxvaf_adj <- max(mean.a.b)
    len_adj <- length(mean.a.b)
    len_adj_2 <- length(mean.a.b)
  } else {
    minvaf_adj <- mean(main.vaf.adj)
    maxvaf_adj <- mean(main.vaf.adj)
    len_adj <- 1
    len_adj_2 <- 1
  }

  data.frame(
    samplename = ctx$id,
    minvaf = minvaf,
    maxvaf = maxvaf,
    len = len,
    len_2 = len_2,
    minvaf_adj = minvaf_adj,
    maxvaf_adj = maxvaf_adj,
    len_adj = len_adj,
    len_adj_2 = len_adj_2,
    stringsAsFactors = FALSE
  )
}

adjust_mu <- function(mu, mu_candidate, times) {
  ifelse(
    mu_candidate != 1 & (mu / mu_candidate > times | mu / mu_candidate < 1 / times),
    0.5 * mu_candidate + 0.5 * mu,
    mu
  )
}

determine_mupick <- function(fitdiff, interdiff, fitdiff2, interdiff2, fitmu, intermu) {
  dplyr::case_when(
    is.na(fitdiff) & !is.na(interdiff) ~ intermu,
    !is.na(fitdiff) & is.na(interdiff) ~ fitmu,
    is.na(fitdiff) & is.na(interdiff) ~ NA_real_,
    fitdiff2 < interdiff2 ~ fitmu,
    TRUE ~ intermu
  )
}

adjust_p <- function(data, magosp, beta, cut = 0.3) {
  data$goodp <- magosp
  data <- data[!(is.na(data$interp) & is.na(data$fitp)), , drop = FALSE]

  data$pickp.alt <- ifelse(
    is.na(data$intermu),
    data$fitp,
    ifelse(
      data$mupick == data$intermu,
      dplyr::coalesce(data$interp, data$fitp),
      dplyr::coalesce(data$fitp, data$interp)
    )
  )

  data$pickp.close <- apply(data, 1, function(row) {
    candidates <- as.numeric(c(row["fitp"], row["interp"], row["backp"]))
    candidates <- candidates[!is.na(candidates)]
    if (length(candidates) == 0) {
      return(NA_real_)
    }
    freq1 <- as.numeric(row["goodp"])
    candidates[which.min(abs(candidates - freq1))]
  })

  data <- data[!is.na(data$pickp.close), , drop = FALSE]
  data$pickp.close <- as.numeric(data$pickp.close)

  data$mupick.new <- ifelse(
    is.na(data$fitp),
    data$intermu,
    ifelse(
      abs(data$pickp.close - data$fitp) < 0.01,
      dplyr::coalesce(data$fitmu, data$intermu),
      dplyr::coalesce(data$intermu, data$fitmu)
    )
  )

  data$mupick.choose <- ifelse(round(abs(data$pickp.alt - data$goodp), 2) >= cut, data$mupick.new, data$mupick)
  data$mupick <- data$mupick.choose
  data$pickp <- ifelse(data$mupick.choose == data$mupick.new, data$pickp.close, data$pickp.alt)
  data$pickp <- ifelse(round(abs(data$goodp - data$pickp), 2) >= cut, data$goodp / 2 + data$pickp / 2, data$pickp)
  data$picks <- ifelse(
    is.na(data$intermu),
    data$fits,
    ifelse(data$mupick == data$intermu, data$inters, data$fits)
  )
  data$pickt1 <- ifelse(
    is.na(data$intermu),
    data$fit_len_diff / data$mupick,
    ifelse(data$mupick == data$intermu, data$inter_len_diff / data$mupick, data$fit_len_diff / data$mupick)
  )
  data$picktend <- ((log(data$pickp) - log(1 - data$pickp)) / (log(2) * beta) + (1 + data$picks) * data$pickt1) / data$picks
  data$picktend <- ifelse(data$picktend > 0, data$picktend, NA)

  valid_idx <- which(data$picktend >= data$pickt1 & data$picktend >= 1 & data$picktend <= 1000)
  invalid_idx <- setdiff(seq_len(nrow(data)), valid_idx)
  data[invalid_idx, c("picks", "mupick", "pickp", "pickt1", "picktend")] <- NA
  data
}

maincluster_candidate_update <- function(ctx) {
  main.vaf <- ctx$main_cluster_vaf
  temp.vaf <- (0.5 + min(main.vaf)) / 2
  iter <- 0
  repeat {
    old.vaf <- if (iter == 0) temp.vaf else new.vaf
    df <- beta_reassign(.vaf_prob_df(main.vaf, c(0.5, old.vaf), ctx$depth))
    df <- df[df$cluster > 1, , drop = FALSE]
    new.vaf <- mean(df$vaf)
    diff <- abs(new.vaf - old.vaf)
    if (diff <= 1 / ctx$depth || iter > 100) {
      break
    }
    iter <- iter + 1
  }
  new.vaf
}

final_process <- function(data.rearrange, rbest_data, ctx) {
  final.data <- data.frame(
    name = ctx$id,
    mu = NA,
    s = NA,
    emergence_time = NA,
    tau = NA,
    p = NA,
    stringsAsFactors = FALSE
  )

  if (identical(data.rearrange$label, "inter")) {
    data <- data.rearrange
    times <- 2
    data <- data %>%
      dplyr::mutate(
        intermu = adjust_mu(intermu, intermu_candidate, times),
        fitmu = adjust_mu(fitmu, fitmu_candidate, times)
      ) %>%
      dplyr::mutate(
        mupick_low_depth = determine_mupick(fitdiff, interdiff, fitdiff2, interdiff2, fitmu, intermu)
      )
    # adjust_p reads `data$mupick` (the original/"low-depth" pick) at line 312:
    #   data$mupick.choose <- ifelse(cond, data$mupick.new, data$mupick)
    # When `data$mupick` is NULL and any element of `cond` is FALSE, R's ifelse
    # errors "replacement has length zero". The rename mupick -> mupick_low_depth
    # was incomplete; restore the column adjust_p expects.
    data$mupick <- data$mupick_low_depth

    data <- adjust_p(data, ctx$magosp, ctx$beta)
    if (nrow(data) > 0) {
      final.data <- data.frame(
        name = ctx$id,
        mu = data$mupick,
        s = data$picks,
        emergence_time = data$pickt1,
        tau = data$picktend / data$pickt1,
        p = data$pickp,
        stringsAsFactors = FALSE
      )
    }
  }

  if (ctx$write_final) {
    header_comment <- paste0(
      "## name: sample ID | mu: mutation rate | ",
      "s: selection coefficient | ",
      "emergence_time: emergence time of the subclone | ",
      "tau: subclone expansion score | ",
      "p: subclonal fraction\n"
    )
    out_path <- file.path(ctx$output_folder, paste0(ctx$output_prefix, ".final.txt"))
    writeLines(header_comment, con = out_path)
    suppressWarnings(write.table(
      final.data,
      file = out_path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      append = TRUE
    ))
  }
  final.data
}

post_process <- function(fitness_result, rbest_result, ctx) {
  data.rearrange <- fitness_result
  second.vaf <- ctx$second_cluster_vaf$normal
  main.vaf <- ctx$main_cluster_vaf
  vaf.all <- c(main.vaf, second.vaf)
  # robust pick: when the candidate vector is empty (no Beta-mix component below
  # close_05_vaf AND empty insert.vaf), assignments like
  # `mean.a.b[which.min(abs(mean.a.b - X))]` evaluate to numeric(0) and crash on
  # data.frame replacement. Return NA in that case (faithful to v1's intent --
  # v1 has the same vulnerable pattern and just escaped the empty case by RNG luck).
  pick_closest <- function(x, target) if (length(x) == 0) NA_real_ else x[which.min(abs(x - target))]

  m <- RBesT::automixfit(vaf.all, type = "beta", Nc = 2:10, thresh = 0, k = 6, Niter.max = 10000)
  a <- m["a", ]
  b <- m["b", ]
  mean.a.b.new <- a / (a + b)
  close_05_vaf <- mean.a.b.new[which.min(abs(mean.a.b.new - 0.5))]
  mean_a_b_filtered <- mean.a.b.new[which(mean.a.b.new > min(main.vaf) & mean.a.b.new < close_05_vaf)]

  if (length(mean_a_b_filtered) == 0) {
    insert.vaf <- maincluster_candidate_update(ctx)
    mean.a.b <- c(mean.a.b.new[which(mean.a.b.new < close_05_vaf)], insert.vaf)
    data.rearrange$clonallen <- 1
  } else {
    mean.a.b <- mean.a.b.new[which(mean.a.b.new < close_05_vaf)]
    data.rearrange$clonallen <- 2
  }

  data.rearrange$minvaf <- if (length(mean.a.b) == 0) NA_real_ else min(mean.a.b)
  data.rearrange$maxvaf <- if (length(mean.a.b) == 0) NA_real_ else max(mean.a.b)
  data.rearrange$closevaf <- pick_closest(mean.a.b, 0.5)

  fitp <- data.rearrange$fitp
  if (!is.na(fitp)) {
    fit1 <- vaf_at_div(1, fitp, ctx)
    data.rearrange$closefitvaf <- pick_closest(mean.a.b, fit1)
    data.rearrange$fitdiff <- abs(data.rearrange$closefitvaf - fit1)
    fit2 <- fitp / 2
    data.rearrange$closefitvaf2 <- pick_closest(mean.a.b, fit2)
    data.rearrange$fitdiff2 <- abs(data.rearrange$closefitvaf2 - fit2)
    df <- beta_reassign(.vaf_prob_df(main.vaf, c(0.5, fit1, fit2), ctx$depth))
    data.rearrange$fit_len_diff <- data.rearrange$fitmu * data.rearrange$fitcell
    data.rearrange$fit_len_ratio <- (data.rearrange$fitmu * data.rearrange$fitcell) / nrow(df[df$cluster > 1, , drop = FALSE])
  } else {
    data.rearrange$closefitvaf <- NA
    data.rearrange$fitdiff <- NA
    data.rearrange$closefitvaf2 <- NA
    data.rearrange$fitdiff2 <- NA
    data.rearrange$fit_len_diff <- NA
    data.rearrange$fit_len_ratio <- NA
  }

  interp <- data.rearrange$interp
  if (!is.na(interp)) {
    inter1 <- vaf_at_div(1, interp, ctx)
    data.rearrange$closeintervaf <- pick_closest(mean.a.b, inter1)
    data.rearrange$interdiff <- abs(data.rearrange$closeintervaf - inter1)
    inter2 <- interp / 2
    data.rearrange$closeintervaf2 <- pick_closest(mean.a.b, inter2)
    data.rearrange$interdiff2 <- abs(data.rearrange$closeintervaf2 - inter2)
    df <- beta_reassign(.vaf_prob_df(main.vaf, c(0.5, inter1, inter2), ctx$depth))
    intermu <- data.rearrange$intermu
    data.rearrange$inter_len_diff <- intermu * data.rearrange$intercell
    data.rearrange$inter_len_ratio <- (intermu * data.rearrange$intercell) / nrow(df[df$cluster > 1, , drop = FALSE])
  } else {
    data.rearrange$closeintervaf <- NA
    data.rearrange$interdiff <- NA
    data.rearrange$closeintervaf2 <- NA
    data.rearrange$interdiff2 <- NA
    data.rearrange$inter_len_diff <- NA
    data.rearrange$inter_len_ratio <- NA
  }

  backp <- data.rearrange$backp
  if (!is.na(backp)) {
    bac1 <- vaf_at_div(1, backp, ctx)
    data.rearrange$closebacvaf <- pick_closest(mean.a.b, bac1)
    data.rearrange$bacdiff <- abs(data.rearrange$closebacvaf - bac1)
    bac2 <- (1 - backp) / (2 * exp(log(2) * ctx$beta * 1))
    data.rearrange$closebacvaf2 <- pick_closest(mean.a.b, bac2)
    data.rearrange$bacdiff2 <- abs(data.rearrange$closebacvaf2 - bac2)
  } else {
    data.rearrange$closebacvaf <- NA
    data.rearrange$bacdiff <- NA
    data.rearrange$closebacvaf2 <- NA
    data.rearrange$bacdiff2 <- NA
  }

  inter_sample <- rbest_result$samplename[which(rbest_result$len_adj > 1)]
  missample <- rbest_result$samplename[which(rbest_result$len_adj == 1)]

  fit_sum <- sum(data.rearrange[c("fitdiff2", "fitdiff2")], na.rm = TRUE)
  if (any(is.na(data.rearrange[c("fitdiff2", "fitdiff2")]))) {
    fit_sum <- Inf
  }
  inter_sum <- sum(data.rearrange[c("interdiff2", "interdiff2")], na.rm = TRUE)
  if (any(is.na(data.rearrange[c("interdiff2", "interdiff2")]))) {
    inter_sum <- Inf
  }
  bac_sum <- sum(data.rearrange[c("bacdiff2", "bacdiff2")], na.rm = TRUE)
  if (any(is.na(data.rearrange[c("bacdiff2", "bacdiff2")]))) {
    bac_sum <- Inf
  }

  data.rearrange$label <- NA_character_
  if (data.rearrange$name %in% inter_sample) {
    data.rearrange$label <- "inter"
  } else if (data.rearrange$name %in% missample) {
    if (inter_sum > bac_sum) {
      data.rearrange$label <- "bac"
    } else {
      data.rearrange$label <- "inter"
    }
  }

  final_process(data.rearrange, rbest_result, ctx)
}

#' Run the TEATIME pipeline
#'
#' Main entry point for TEATIME. Accepts pre-processed MAGOS clustering output
#' or raw VCF-like data, estimates evolutionary parameters, and optionally
#' writes a tab-separated result file.
#'
#' @param input Input data. Either a list with elements `purity` and `result`
#'   (MAGOS object) when `input_format = "magos"`, a `data.frame` with columns
#'   `REF` and `ALT` (when `input_format = "vcf"`), or a `data.frame` already
#'   containing `vaf.1`, `depth.1`, and `colors` (when `input_format = "raw"`).
#' @param beta Numeric. Survival rate (default `0.9`). Used as the `r`
#'   parameter of the growth model.
#' @param depth Integer or `NA`. Mean sequencing depth. If `NA` (default),
#'   computed from the data.
#' @param p_thre Numeric. P-value threshold for breakpoint and fitness tests
#'   (default `0.01`).
#' @param input_format Character. One of `"vcf"` (default), `"magos"`, or
#'   `"raw"`.
#' @param growth_model Character. Name of a registered growth model (default
#'   `"exponential"`). Register custom models with [register_growth_model()].
#' @param verbose Logical. If `TRUE`, prints step-level timing messages
#'   (default `FALSE`).
#' @param debug Logical. If `TRUE`, wraps every pipeline step in error
#'   handling and prints step name, key intermediate values, and elapsed time.
#'   On failure the step name, error message, and available context are printed
#'   before stopping, so the problem can be located immediately (default
#'   `FALSE`).
#' @param output_folder Character. Directory for output files (default `"./"`).
#' @param output_prefix Character. Prefix for output file names (default
#'   `"TEATIME"`).
#' @param id Character. Sample identifier written to the result table (default
#'   `"T01"`).
#' @param write_final Logical. If `TRUE` (default), writes
#'   `<output_prefix>.final.txt` to `output_folder`.
#' @param seed Integer or `NA`. Random seed for reproducibility (default `123`).
#'   Set to `NA` to disable seeding and run all estimators three independent
#'   times for stochastic robustness.
#' @param extra Named list of additional parameters passed to custom growth
#'   models via `ctx$extra` (default `list()`). Ignored by the built-in
#'   exponential model.
#' @param fast_version Logical (default `FALSE`). `FALSE` is the reference
#'   pipeline. `TRUE` swaps the hot inner kernels for C++ implementations
#'   (Wilcoxon p-value, dbeta matrix, beta-reassignment) and uses vectorised
#'   rbinom in `slope_simu`, giving roughly a 5-7x per-sample speedup while
#'   keeping the same single-stream RNG order as default. Output is identical
#'   to default mode on the deterministic majority of inputs; small drift on
#'   borderline samples can occur because the C++ kernels round at the last
#'   floating-point bit.
#'
#' @return A one-row `data.frame` with columns:
#'   \describe{
#'     \item{name}{Sample ID (`id` argument).}
#'     \item{mu}{Mutation rate.}
#'     \item{s}{Selection coefficient.}
#'     \item{emergence_time}{Emergence time of the subclone.}
#'     \item{tau}{Subclone expansion score.}
#'     \item{p}{Subclonal fraction.}
#'   }
#'   All numeric columns are `NA` when the pipeline cannot produce a reliable
#'   estimate.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using a pre-computed MAGOS object
#' input <- list(purity = magos.114$purity, result = magos.114$results)
#' result <- TEATIME.run(input, beta = 0.9, depth = 150, seed = 42)
#' }
TEATIME.run <- function(
  input,
  beta = 0.9,
  depth = NA,
  p_thre = 0.01,
  input_format = "vcf",
  growth_model = "exponential",
  verbose = FALSE,
  output_folder = "./",
  output_prefix = "TEATIME",
  id = "T01",
  write_final = TRUE,
  seed = 123,
  extra = list(),
  debug = FALSE,
  save_magos = FALSE,
  fast_version = FALSE
) {
  # Default mode (fast_version = FALSE) dispatches to the embedded af3e64d v1
  # source so the output is bit-identical to the v1 release. Fast mode keeps
  # the v2 modular pipeline (registry / Rcpp / mclapply). The fork is taken
  # only for input_format == "vcf", which is v1's only accepted input shape.
  if (!isTRUE(fast_version) && identical(input_format, "vcf")) {
    return(.run_v1_default(
      input = input, beta = beta, depth = depth, p_thre = p_thre,
      output_folder = output_folder, output_prefix = output_prefix,
      id = id, write_final = write_final, seed = seed, debug = debug
    ))
  }

  # Single switch read by the inner speedup paths (Rcpp Wilcoxon, vectorised
  # beta_reassign, get_slope cache in slope_method). When `fast_version = TRUE`
  # these use their accelerated implementations; otherwise every inner function
  # behaves exactly as in the reference release.
  prev_fast_opt <- getOption("teatime.fast_version", FALSE)
  options(teatime.fast_version = isTRUE(fast_version))
  on.exit(options(teatime.fast_version = prev_fast_opt), add = TRUE)

  if (!is.na(seed)) {
    set.seed(seed)
    # v1 (af3e64d TEATIME.r:3183-3186): always advance the seed once before
    # the prepare-data step, regardless of input format. Matches v1 RNG path.
    seed <- seed + 1
    set.seed(seed)
  }

  # When debug = FALSE, silence the noisy "max(empty)" / "ties in p-value"
  # warnings that come from estimate.R's many max() / wilcox calls -- they
  # are harmless to the final result. Wrap the WHOLE pipeline body so every
  # step is silenced. Set debug = TRUE to see them again.
  .step <- function(name, expr, info_fn = NULL) {
    if (!debug) return(suppressWarnings(expr))
    t0  <- proc.time()[["elapsed"]]
    cat(sprintf("[DEBUG] %-20s ... ", name))
    result <- tryCatch(withCallingHandlers(expr,
      error = function(e) {
        cs <- sys.calls()
        cat(sprintf("FAILED (%.1fs)\n", proc.time()[["elapsed"]] - t0))
        cat(sprintf("[DEBUG]   error  : %s\n", conditionMessage(e)))
        cat("[DEBUG]   trace  :\n")
        for (i in seq_along(cs)) {
          s <- paste(deparse(cs[[i]]), collapse = " ")
          if (nchar(s) > 200) s <- paste0(substr(s, 1, 200), "...")
          cat(sprintf("[DEBUG]    %2d: %s\n", i, s))
        }
      }),
      error = function(e) {
        if (!is.null(info_fn)) {
          info <- tryCatch(info_fn(), error = function(e2) NULL)
          if (!is.null(info)) cat(sprintf("[DEBUG]   context: %s\n", info))
        }
        stop(e)
      })
    elapsed <- proc.time()[["elapsed"]] - t0
    info_str <- if (!is.null(info_fn)) tryCatch(info_fn(result), error = function(e) "") else ""
    cat(sprintf("OK (%.1fs)%s\n", elapsed, if (nzchar(info_str)) paste0("  |  ", info_str) else ""))
    result
  }

  # All-NA fallback writer: guarantees a result (and a .final.txt when
  # write_final = TRUE) even if some pipeline step errors. A failed sample
  # is reported as NA rather than crashing the caller / leaving no file.
  .na_result <- function() {
    fd <- data.frame(name = id, mu = NA, s = NA,
                     emergence_time = NA, tau = NA, p = NA,
                     stringsAsFactors = FALSE)
    if (isTRUE(write_final)) {
      hdr <- paste0(
        "## name: sample ID | mu: mutation rate | s: selection coefficient | ",
        "emergence_time: emergence time of the subclone | ",
        "tau: subclone expansion score | p: subclonal fraction\n")
      out_path <- file.path(output_folder, paste0(output_prefix, ".final.txt"))
      tryCatch({
        dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
        writeLines(hdr, con = out_path)
        suppressWarnings(write.table(fd, file = out_path, sep = "\t",
                                     row.names = FALSE, quote = FALSE,
                                     append = TRUE))
      }, error = function(e2) NULL)
    }
    fd
  }

  tryCatch({
    ctx <- .step("prepare_data", prepare_data(
      input = input, beta = beta, depth = depth,
      input_format = input_format, growth_model = growth_model,
      verbose = verbose, output_folder = output_folder,
      output_prefix = output_prefix, id = id,
      write_final = write_final, seed = seed, extra = extra,
      save_magos = save_magos, fast_version = fast_version
    ), info_fn = function(r = NULL) {
      if (is.null(r)) return(sprintf("input_format=%s  beta=%s  depth=?", input_format, beta))
      sprintf("depth=%d  n_mut=%d  n_clusters=%d  main_vaf=%.3f  magosp=%.3f",
              r$depth, length(r$vafdata$vaf.1),
              length(unique(r$vafdata$colors)),
              mean(r$main_cluster_vaf), r$magosp)
    })

    rbest <- .step("run_rbest", run_rbest(ctx),
      info_fn = function(r = NULL) if (!is.null(r)) sprintf("label=%s", r$label %||% "?") else NULL)

    estimates <- .step("run_estimates", run_estimates(ctx, p_thre = p_thre),
      info_fn = function(r = NULL) {
        if (is.null(r)) return(NULL)
        fit_n   <- if (!is.null(r$fit$all))    nrow(r$fit$all)    else 0L
        bac_n   <- if (!is.null(r$bac$all))    nrow(r$bac$all)    else 0L
        normal_n <- if (!is.null(r$normal$all)) nrow(r$normal$all) else 0L
        sprintf("fit=%d rows  bac=%d rows  normal=%d rows", fit_n, bac_n, normal_n)
      })

    fitness <- .step("run_fitness", run_fitness(estimates, ctx),
      info_fn = function(r = NULL) {
        if (is.null(r)) return(NULL)
        sprintf("fitmu=%s  intermu=%s  fitp=%s  backp=%s",
                round(r$fitmu,    3), round(r$intermu,  3),
                round(r$fitp,     4), round(r$backp,    4))
      })

    .step("post_process", post_process(fitness, rbest, ctx),
      info_fn = function(r = NULL) {
        if (is.null(r) || !is.data.frame(r)) return(NULL)
        sprintf("mu=%s  s=%s  emerge=%s  tau=%s  p=%s",
                r$mu, round(r$s, 3), r$emergence_time,
                round(r$tau, 3), round(r$p, 4))
      })
  }, error = function(e) {
    if (debug) {
      cat(sprintf("[DEBUG] pipeline failed: %s -- returning NA result\n",
                  conditionMessage(e)))
    }
    .na_result()
  })
}
