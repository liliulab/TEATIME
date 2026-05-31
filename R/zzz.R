.teatime_ref_env <- NULL
.teatime_prod_env <- NULL
.teatime_prod_fast_env <- NULL

# Cache the VCF reference environment on first use.
.load_reference <- function() {
  if (!is.null(.teatime_ref_env)) return(.teatime_ref_env)
  env <- new.env(parent = asNamespace("TEATIME"))
  magos_path   <- system.file("extdata/reference_source/MAGOS.R",   package = "TEATIME")
  teatime_path <- system.file("extdata/reference_source/TEATIME.r", package = "TEATIME")
  if (!nzchar(magos_path) || !nzchar(teatime_path)) {
    stop("TEATIME: reference source missing from inst/extdata/reference_source/",
         call. = FALSE)
  }
  suppressMessages(sys.source(magos_path, envir = env, keep.source = FALSE))
  src <- readLines(teatime_path)
  src <- gsub('^source\\(file\\.path\\("R", "MAGOS\\.R"\\)\\)', '# (already loaded)', src)
  tmp <- tempfile(fileext = ".r"); on.exit(unlink(tmp), add = TRUE)
  writeLines(src, tmp)
  suppressMessages(sys.source(tmp, envir = env, keep.source = FALSE))
  utils::assignInMyNamespace(".teatime_ref_env", env)
  env
}

# Default-mode dispatcher: run the reference implementation into a scratch
# tempdir so all intermediate per-step text files land there and are cleaned up
# on exit. Only the harvested result is written at the user's output_folder
# when write_final = TRUE.
.run_default_dispatch <- function(input, beta, depth, p_thre,
                                  output_folder, output_prefix,
                                  id, write_final, seed, debug) {
  ref <- .load_reference()
  scratch <- tempfile("teatime_scratch_"); dir.create(scratch, recursive = TRUE)
  on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)
  run_seed <- if (is.na(seed)) NA else as.integer(seed)
  call_args <- list(
    input.file    = input,
    beta          = beta,
    depth         = depth,
    p_thre        = p_thre,
    magos_object  = TRUE,
    output.folder = scratch,
    output.prefix = output_prefix,
    id            = id,
    steps         = 0:5,
    write_final   = TRUE,
    debug_mode    = isTRUE(debug),
    purity_set    = 0,
    seed          = run_seed
  )
  if (isTRUE(debug)) {
    do.call(ref$TEATIME.run, call_args)
  } else {
    suppressWarnings(suppressMessages(
      utils::capture.output(do.call(ref$TEATIME.run, call_args))
    ))
  }

  scratch_final <- file.path(scratch, paste0(output_prefix, ".final.txt"))
  raw_df <- if (file.exists(scratch_final)) {
    utils::read.table(scratch_final, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
  } else {
    data.frame(name = id, mu = NA_real_, s = NA_real_,
               t1 = NA_real_, tend = NA_real_, p = NA_real_,
               stringsAsFactors = FALSE)
  }

  out_df <- .normalize_output_schema(raw_df)

  if (isTRUE(write_final)) {
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    }
    final_path <- file.path(output_folder, paste0(output_prefix, ".final.txt"))
    hdr <- paste0(
      "## name: sample ID | mu: mutation rate | s: selection coefficient | ",
      "emergence_time: emergence time of the subclone | ",
      "tau: subclone expansion score | p: subclonal fraction\n")
    tryCatch({
      writeLines(hdr, con = final_path)
      suppressWarnings(utils::write.table(
        out_df, file = final_path, sep = "\t",
        row.names = FALSE, quote = FALSE, append = TRUE))
    }, error = function(e2) NULL)
  }
  out_df
}

# Build (and cache) the production estimation environment. `fast = TRUE` injects
# the accelerated kernels. Used by the MAGOS dispatch.
.load_production <- function(fast = FALSE) {
  if (!fast && !is.null(.teatime_prod_env))      return(.teatime_prod_env)
  if ( fast && !is.null(.teatime_prod_fast_env)) return(.teatime_prod_fast_env)

  env <- new.env(parent = asNamespace("TEATIME"))

  # The production functions use unqualified names from several CRAN packages;
  # copy those exports into the dispatch env so they resolve here.
  for (pkg in c("RBesT", "dplyr", "tidyr", "parallel",
                "likelihoodExplore", "magrittr", "Matrix", "igraph")) {
    ns <- asNamespace(pkg)
    for (nm in getNamespaceExports(pkg)) {
      if (!exists(nm, envir = env, inherits = FALSE)) {
        assign(nm, get(nm, envir = ns), envir = env)
      }
    }
  }

  # `cliffDelta` shim over `.cliffs_delta_abs` (avoids an rcompanion dependency).
  # Returns signif(|Cliff's delta|, 3); grouping follows the factor-level order of
  # `group`, matching rcompanion (A = x[g == levels(g)[1]]).
  env$cliffDelta <- function(formula, data, ...) {
    vars   <- all.vars(formula)
    val_v  <- as.numeric(data[[vars[1]]])
    grp_v  <- data[[vars[2]]]
    levs   <- levels(factor(grp_v))
    a      <- val_v[grp_v == levs[1]]
    b      <- val_v[grp_v == levs[2]]
    signif(.cliffs_delta_abs(a, b), 3)
  }

  # Copy the production functions (`.prod_defs`) into the dispatch env, re-pointing
  # each closure to `env` so their names resolve here and fast-kernel overrides can
  # shadow them.
  for (nm in ls(.prod_defs, all.names = TRUE)) {
    obj <- get(nm, envir = .prod_defs)
    if (is.function(obj)) environment(obj) <- env
    assign(nm, obj, envir = env)
  }

  # Fast variant: override the hot primitives in `env` with accelerated
  # equivalents. Each override produces identical numbers and preserves the RNG
  # draw order, so fast mode matches default mode under the same seed.
  if (isTRUE(fast)) {
    .inject_fast_kernels(env)
    utils::assignInMyNamespace(".teatime_prod_fast_env", env)
  } else {
    utils::assignInMyNamespace(".teatime_prod_env", env)
  }
  env
}

# Override the bundle's hot primitives in-place with accelerated, bit-identical
# implementations. Each override is gated/structured to match its faithful
# counterpart exactly (verified in tests/). Deterministic kernels (no RNG) are
# free to swap; RNG-batched function overrides preserve draw order + count.
.inject_fast_kernels <- function(env) {
  # (1) beta_reassign -> the package's math.R beta_reassign, which uses the C++
  #     beta_reassign_core under the teatime.fast_version option (set by the
  #     fast dispatcher during the run). Bit-identical incl. row order (verified).
  pkg_beta_reassign <- beta_reassign  # resolves to math.R version in the namespace
  env$beta_reassign <- function(df) pkg_beta_reassign(df)

  # (1b) log_likelihood_mixture -> math.R version (dbeta_matrix C++ + C++
  #      beta_reassign + .likbeta_sum dedup). Bit-identical (rel diff 0) to the
  #      bundle's dplyr+likbeta version but ~1000x faster on large clusters,
  #      where the bundle's per-row likelihoodExplore::likbeta over tens of
  #      thousands of rows is THE dominant cost (~25s/call on 20K muts).
  pkg_llm <- log_likelihood_mixture
  env$log_likelihood_mixture <- function(data, p_vec, depth) pkg_llm(data, p_vec, depth)

  # (1c) dbeta -> deduplicating wrapper. The bundle computes dbeta(vaf_set, a, b)
  #      over tens of thousands of VAFs that take only ~depth distinct values
  #      (e.g. 410 unique of 20009). Deduplicating the x argument is bit-identical
  #      (same values, expanded by match()) and collapses the dbeta work by the
  #      redundancy factor. Only dedup for scalar shapes + large x; otherwise
  #      delegate to stats::dbeta unchanged.
  env$dbeta <- function(x, shape1, shape2, ncp, log = FALSE) {
    if (!missing(ncp)) return(stats::dbeta(x, shape1, shape2, ncp, log = log))
    if (length(shape1) != 1L || length(shape2) != 1L || length(x) < 200L) {
      return(stats::dbeta(x, shape1, shape2, log = log))
    }
    ux <- unique(x)
    stats::dbeta(ux, shape1, shape2, log = log)[match(x, ux)]
  }

  # (2) wilcox.test -> asymptotic C++ p-value, but ONLY on the asymptotic
  #     branch (n>=50 or ties). stats::wilcox.test takes the EXACT branch when
  #     both samples are < 50 with no ties; there the C++ value differs, so we
  #     fall back. Returns a minimal list exposing $p.value (the only field the
  #     bundle reads from wilcox.test results; cliffDelta's $statistic use is
  #     covered by the separate cliffDelta shim).
  env$wilcox.test <- function(x, y = NULL, ...) {
    if (is.null(y)) return(stats::wilcox.test(x, ...))
    use_asym <- length(x) >= 50L || length(y) >= 50L ||
      anyDuplicated(c(x, y)) > 0L
    if (use_asym) {
      list(p.value = wilcox_test_p_asym(as.numeric(x), as.numeric(y)))
    } else {
      stats::wilcox.test(x, y, ...)
    }
  }

  # (3) ks.test -> C++ statistic + R's own asymptotic p-value, but ONLY on the
  #     asymptotic branch. stats::ks.test sets exact <- (n.x*n.y < 10000) when
  #     exact is unspecified; below that threshold it uses psmirnov_exact, so we
  #     defer to stats::ks.test there (small, and RNG-neutral). At/above 10000
  #     it goes asymptotic -> psmirnov_asymp, which depends only on the D
  #     statistic (ties drop out). ks_stat_2s_cpp computes D bit-identically and
  #     stats::psmirnov(exact=FALSE) gives the same scalar p-value, so the whole
  #     result is bit-identical to stats::ks.test on this branch. This is the
  #     dominant cost inside simulate_from_estimation's 50-try loop on large
  #     clusters (the order/sort/cumsum on ~10^4-element vectors).
  env$ks.test <- function(x, y, ...) {
    nx <- length(x); ny <- length(y)
    if (nx < 1L || ny < 1L || nx * ny < 10000) return(stats::ks.test(x, y, ...))
    D <- ks_stat_2s_cpp(as.numeric(x), as.numeric(y))
    pval <- stats::psmirnov(D, sizes = c(nx, ny), alternative = "two.sided",
                            exact = FALSE, lower.tail = FALSE)
    list(statistic = c(D = D), p.value = min(1, max(0, pval)))
  }

  # (4) simulate_from_estimation -> batch the per-try rbinom into a single draw.
  #     The bundle's loop does, for try in 1:ntry:
  #         sapply(vaf_list, \(v) round(rbinom(round(mu), depth, v)/depth, nd))
  #     i.e. per try it draws round(mu) binomials at each vaf in order, giving a
  #     (mu_int x K) matrix read row-major via as.vector(t(.)). Because rbinom
  #     consumes the stream element-by-element with each element's own (size,
  #     prob) -- size = depth is constant here -- one rbinom over the tiled prob
  #     vector rep(rep(vaf_list, each=mu_int), times=ntry) draws the IDENTICAL
  #     values in the IDENTICAL order as the ntry separate sapply calls. Each
  #     try's segment is reshaped with the same matrix(.,nrow=mu_int);
  #     as.vector(t(.)) so the simulated vector fed to ks/wilcox is bit-identical.
  #     The per-try ks.test/wilcox.test (now the C++ kernels above) and the final
  #     log_likelihood_mixture / AIC / BIC are unchanged. This removes the only
  #     RNG-consuming for-loop in the hot path without altering the draw stream.
  sfe <- function(compare_set, cell.div, mu, p, depth, beta = 0.8, num_decimal) {
    vaf_list <- c(p / 2)
    for (i in 1:cell.div) {
      cell.num <- exp(log(2) * beta * i)
      vaf_list <- c(vaf_list, p / 2 + (1 - p) / (2 * cell.num))
    }
    ntry    <- getOption("teatime.sim_ntry", 50L)
    mu_int  <- round(mu)
    K       <- length(vaf_list)
    per_try <- mu_int * K
    prob_vec <- rep(rep(vaf_list, each = mu_int), times = ntry)
    draws    <- round(rbinom(per_try * ntry, depth, prob_vec) / depth, num_decimal)
    p_v_1_values <- numeric(ntry)
    p_v_2_values <- numeric(ntry)
    # Batch the 50-try KS: compare_set is fixed across tries, so sort it once and
    # compute every try's D in one C++ pass (ks_d_batch_cpp). The KS statistic
    # depends only on each segment's multiset, so the per-try reshape is
    # irrelevant to D -> we feed the raw `draws`. Gated to the asymptotic regime
    # (nx*ny >= 10000, matching ks.test's exact<-(nx*ny<10000)); below that the
    # per-try stats::ks.test exact branch is kept so small (MMRF) samples stay
    # bit-identical. The p-values come from R's own psmirnov (vectorised).
    nx <- length(compare_set)
    use_batch <- (nx >= 1L && per_try >= 1L && nx * per_try >= 10000)
    if (use_batch) {
      # Both KS and Wilcoxon are multiset statistics (the per-try reshape is
      # irrelevant), and compare_set is fixed across tries -> compute all ntry
      # at once with one sort of compare_set. nx*per_try >= 10000 guarantees the
      # asymptotic branch for both (ks.test exact<-(nx*ny<10000); a sample with
      # max(nx,ny) >= 100 is >= 50 so wilcox is asymptotic too), so this is
      # bit-identical to the per-try injected ks.test / wilcox.test path.
      D_vec <- ks_d_batch_cpp(compare_set, draws, ntry)
      ks_p_all <- pmin(1, pmax(0, stats::psmirnov(
        D_vec, sizes = c(nx, per_try), alternative = "two.sided",
        exact = FALSE, lower.tail = FALSE)))
      wx_p_all <- wilcox_p_batch_cpp(compare_set, draws, ntry)
    }
    for (try in 1:ntry) {
      if (use_batch) {
        p_v_1_values[try] <- ks_p_all[try]
        p_v_2_values[try] <- wx_p_all[try]
      } else {
        seg <- draws[((try - 1L) * per_try + 1L):(try * per_try)]
        all_simulated_vafs <- as.vector(t(matrix(seg, nrow = mu_int, ncol = K)))
        p_v_1_values[try] <- ks.test(compare_set, all_simulated_vafs)$p.value
        p_v_2_values[try] <- wilcox.test(compare_set, all_simulated_vafs)$p.value
      }
    }
    mean_p_v_1 <- mean(p_v_1_values)
    mean_p_v_2 <- mean(p_v_2_values)
    log_likelihood <- log_likelihood_mixture(compare_set, vaf_list, depth)
    num_params  <- length(vaf_list)
    sample_size <- length(compare_set)
    AIC_value <- compute_AIC(log_likelihood, num_params)
    BIC_value <- compute_BIC(log_likelihood, num_params, sample_size)
    c(mean_p_v_1, mean_p_v_2, log_likelihood, AIC_value, BIC_value)
  }
  environment(sfe) <- env   # resolve ks.test/wilcox.test/llm/compute_* in env
  env$simulate_from_estimation <- sfe

  # (5) peak_test: count-based Wilcoxon p + Cliff's delta. The simulated set
  #     round(rbinom(n, depth, vaf)/depth) has only ~depth distinct values, so both
  #     statistics are computed from value-counts (.wilcox_p_counts /
  #     .cliffs_delta_abs) instead of the full vector. The rbinom draw and
  #     generate_bootstrap_samples are kept verbatim, preserving the RNG stream.
  #     The approximation branch (teatime.approx) replaces the draw with its exact
  #     dbinom distribution.
  pt <- function(df, Min_Sample_size, vaf, k, depth, num_decimal) {
    if (nrow(df) < k) return(c(NA, NA, NA, NA))
    n_simulations <- getOption("teatime.peak_nsim", 100000)
    vals_int <- round((0:depth) / depth, num_decimal)
    if (isTRUE(getOption("teatime.approx", FALSE))) {
      # APPROXIMATION MODE (mathematical substitution, NOT fewer simulations):
      # round(rbinom(depth, vaf)/depth) has the exact PMF dbinom(0:depth, depth, vaf).
      # Use n_simulations * PMF as the *expected* value-counts -> this is the
      # n_simulations -> Inf limit of the Monte-Carlo (noise-free, deterministic,
      # consumes NO RNG and skips generating 100k draws). Verified to match the
      # 100k-sim mean/Cliff/Wilcoxon-p to ~1e-3 (PoC), since 100k sampling noise
      # is exactly what the analytic form removes.
      probs    <- dbinom(0:depth, depth, vaf)
      cnt_int  <- round(n_simulations * probs)
      mean_sim <- sum(vals_int * probs)
    } else {
      raw <- rbinom(n_simulations, depth, vaf)                 # integers 0..depth (RNG-identical)
      cnt_int  <- tabulate(raw + 1L, nbins = depth + 1L)       # value-counts in one pass
      mean_sim <- mean(round(raw / depth, num_decimal))        # bit-identical mean()
    }
    if (k < Min_Sample_size) {
      largest_right_values <- generate_bootstrap_samples(df[1:k, ]$vaf, Min_Sample_size, num_decimal)
    } else {
      largest_right_values <- df[1:k, ]$vaf
    }
    largest_right_values <- as.vector(na.omit(largest_right_values))
    # aggregate counts by rounded value (rowsum -> sorted unique values + summed
    # counts); shared by both statistics. Equals unique(values)+counts.
    keep <- cnt_int > 0
    ag <- rowsum(as.numeric(cnt_int[keep]), group = vals_int[keep], reorder = TRUE)
    xv <- as.numeric(rownames(ag)); xc <- ag[, 1]
    w_p <- .wilcox_p_counts_pre(xv, xc, largest_right_values)
    cd  <- signif(.cliffs_delta_abs_pre(xv, xc, largest_right_values), 3)
    c(mean_sim, mean(largest_right_values), abs(cd), w_p)
  }
  environment(pt) <- env   # resolve generate_bootstrap_samples (env) + helpers (ns)
  env$peak_test <- pt

  invisible(env)
}

# Production dispatcher (MAGOS input). Runs the estimation pipeline in a scratch
# tempdir (per-sample intermediates land there and are cleaned up), computes Rbest
# per sample, and applies the magosp/adjust_p correction so the returned schema
# matches TEATIME.run's contract.
.run_production_dispatch <- function(input, beta, depth, p_thre,
                                     output_folder, output_prefix,
                                     id, write_final, seed, debug,
                                     cached_intermediates = NULL,
                                     fast = FALSE) {
  prod <- .load_production(fast = fast)

  # In fast mode, the injected kernels (beta_reassign C++, dbeta_matrix C++,
  # wilcox C++, log_likelihood_mixture, dbeta dedup — the math.R fast gates)
  # read getOption("teatime.fast_version"). Set it for the duration of this
  # run; restore on exit.
  #
  # NOTE: fast mode does NOT cap the mu-candidate sweep and does NOT reduce
  # peak_nsim / sim_ntry. Monte-Carlo size (100000), repeat count (50) and the
  # full candidate list stay faithful, so fast is *bit-identical* to a fresh
  # seeded production run — speed comes only from the accelerated deterministic
  # kernels, not from doing less work. teatime.max_candidates / teatime.peak_nsim
  # / teatime.sim_ntry remain available for manual opt-in tolerance if a caller
  # explicitly wants to trade accuracy for speed.
  if (isTRUE(fast)) {
    prev_fast <- getOption("teatime.fast_version", FALSE)
    options(teatime.fast_version = TRUE)
    on.exit(options(teatime.fast_version = prev_fast), add = TRUE)
  }

  scratch <- tempfile("teatime_prod_"); dir.create(scratch, recursive = TRUE)
  on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)

  # Layout the on-disk structure the reference scripts expect:
  #   <scratch>/<id>/All_<id>.rds
  #   <scratch>/Rbest.classify.rds         (Rbest classification table)
  sample_dir <- file.path(scratch, id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  m_obj <- list(purity = input$purity, result = input$result)
  saveRDS(m_obj, file.path(sample_dir, paste0("All_", id, ".rds")))

  savepath <- paste0(scratch, "/")
  dpcode   <- "mean"
  depth_int <- as.integer(round(depth))

  # Per-step wall-time accumulator. Emitted when options(teatime.timing=TRUE)
  # and attached to the result as attr(.,"step_times") so the slowest stage
  # is visible for targeting speedups. `timed()` also quiets the wrapped call
  # (capture.output returns the captured text, so we assign the value inside).
  .timings <- list()
  timed <- function(name, expr_call) {
    t0 <- proc.time()[["elapsed"]]
    res <- NULL
    if (isTRUE(debug)) {
      res <- force(expr_call)
    } else {
      suppressWarnings(suppressMessages(
        invisible(utils::capture.output({ res <- force(expr_call) }))
      ))
    }
    dt <- proc.time()[["elapsed"]] - t0
    prev <- if (is.null(.timings[[name]])) 0 else .timings[[name]]
    .timings[[name]] <<- prev + dt
    if (isTRUE(getOption("teatime.timing", FALSE))) {
      cat(sprintf("[TEATIME stage] %-20s %8.2fs\n", name, dt)); flush.console()
    }
    res
  }
  silent <- function(expr_call) timed("misc", expr_call)

  # Seed BEFORE any RNG (Rbest_classify -> automixfit is stochastic).
  # Approximation mode is deterministic: the analytic kernels remove the heavy
  # Monte-Carlo, and we pin the residual stochastic *choices* (BIC tie-breaks,
  # automixfit init, sample_frac) with a fixed seed so the result is repeatable
  # run-to-run with no seed needed. An explicit user seed still wins.
  if (!is.na(seed)) {
    set.seed(seed)
  } else if (isTRUE(getOption("teatime.approx", FALSE))) {
    set.seed(0L)
  }

  # Rbest classification is ALWAYS computed per-sample from this sample's own
  # MAGOS clustering (the bundle's Rbest_classify reads All_<id>.rds and writes
  # Rbest.classify.rds). The pipeline is fully self-contained: each sample gets
  # its own Rbest, never a precomputed cross-sample table.
  timed("Rbest_classify", prod$Rbest_classify(id, savepath))

  # Optional shortcut: if the caller points `cached_intermediates` at a
  # directory holding a previously-written `all.sample.all.guess.depth*.rds`
  # for this sample, skip the upstream re-estimation and run only Stage 5+6 on
  # that cached intermediate. When absent, the full pipeline runs from scratch.
  cached_file <- if (!is.null(cached_intermediates)) {
    file.path(cached_intermediates, paste0("all.sample.all.guess.depth", dpcode, ".rds"))
  } else NA_character_

  if (!is.na(cached_file) && file.exists(cached_file)) {
    file.copy(cached_file, file.path(sample_dir, basename(cached_file)),
              overwrite = TRUE)
    cached_start <- file.path(cached_intermediates,
                              paste0("all.sample.all.guess.start.choose.depth", dpcode, ".rds"))
    if (file.exists(cached_start)) {
      file.copy(cached_start, file.path(sample_dir, basename(cached_start)),
                overwrite = TRUE)
    }
  } else {
    timed("fit.all.run",     prod$fit.all.run(id, depth_int, beta, savepath, dpcode))
    timed("bac.all.run",     prod$bac.all.run(id, depth_int, beta, savepath, dpcode))
    timed("normal.all.run",  prod$normal.all.run(id, depth_int, beta, savepath, dpcode))
    timed("All.guess.update", prod$All.guess.update(id, savepath, depth_int, beta, dpcode))
  }
  dr <- timed("Post_process", prod$Post_process(id, savepath, depth_int, beta, dpcode))
  if (!is.data.frame(dr) || nrow(dr) == 0 || is.na(dr$label[1]) || dr$label[1] != "inter") {
    out_na <- .na_v257_row(id)
    .emit_timings(.timings, id, getOption("teatime.timing", FALSE))
    attr(out_na, "step_times") <- .timings
    return(out_na)
  }
  timed("Inter_post_process", prod$Inter_post_process(id, dr, savepath, depth_int, beta, dpcode))

  out_rds <- file.path(sample_dir, paste0("final.result.depth", dpcode, ".rds"))
  if (!file.exists(out_rds)) {
    out_na <- .na_v257_row(id)
    .emit_timings(.timings, id, getOption("teatime.timing", FALSE))
    attr(out_na, "step_times") <- .timings
    return(out_na)
  }
  r <- readRDS(out_rds)

  # Stage-C magosp/adjust_p correction.
  r <- timed("stage_c", {
    magosp <- .stageC_magosp(input$result, input$purity)
    .stageC_adjust_p(r, magosp, beta = beta)
  })

  # A NA picktend (no valid emergence timeline) means no usable result: return the
  # all-NA row and write no final output.
  if (is.na(r$picktend[1])) {
    out_na <- .na_v257_row(id)
    .emit_timings(.timings, id, getOption("teatime.timing", FALSE))
    attr(out_na, "step_times") <- .timings
    return(out_na)
  }

  out_df <- data.frame(
    name           = id,
    mu             = r$mupick[1],
    s              = r$picks[1],
    emergence_time = r$pickt1[1],
    tau            = r$picktend[1] / r$pickt1[1],
    p              = r$pickp[1],
    stringsAsFactors = FALSE
  )
  .emit_timings(.timings, id, getOption("teatime.timing", FALSE))
  attr(out_df, "step_times") <- .timings

  if (isTRUE(write_final)) {
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    }
    final_path <- file.path(output_folder, paste0(output_prefix, ".final.txt"))
    hdr <- paste0(
      "## name: sample ID | mu: mutation rate | s: selection coefficient | ",
      "emergence_time: emergence time of the subclone | ",
      "tau: subclone expansion score | p: subclonal fraction\n")
    tryCatch({
      writeLines(hdr, con = final_path)
      suppressWarnings(utils::write.table(
        out_df, file = final_path, sep = "\t",
        row.names = FALSE, quote = FALSE, append = TRUE))
    }, error = function(e) NULL)
  }
  out_df
}

# Fast-mode dispatcher: same control flow + RNG order as the production
# dispatch, but with the injected accelerated kernels (fast = TRUE). Verified
# bit-identical to a fresh seeded production run on the same MAGOS input.
.run_fast_dispatch <- function(input, beta, depth, p_thre,
                               output_folder, output_prefix,
                               id, write_final, seed, debug,
                               cached_intermediates = NULL) {
  .run_production_dispatch(
    input = input, beta = beta, depth = depth, p_thre = p_thre,
    output_folder = output_folder, output_prefix = output_prefix,
    id = id, write_final = write_final, seed = seed, debug = debug,
    cached_intermediates = cached_intermediates,
    fast = TRUE
  )
}

# Pretty-print a per-step wall-time breakdown (descending), with the slowest
# step highlighted. Called when options(teatime.timing=TRUE).
.emit_timings <- function(timings, id, on) {
  if (!isTRUE(on) || length(timings) == 0) return(invisible(NULL))
  tt <- sort(unlist(timings), decreasing = TRUE)
  total <- sum(tt)
  cat(sprintf("\n[TEATIME timing] %s  (total %.1fs)\n", id, total))
  for (nm in names(tt)) {
    cat(sprintf("  %-20s %8.2fs  %5.1f%%\n", nm, tt[[nm]],
                if (total > 0) 100 * tt[[nm]] / total else 0))
  }
  invisible(NULL)
}

.na_v257_row <- function(id) {
  data.frame(
    name = id, mu = NA_real_, s = NA_real_,
    emergence_time = NA_real_, tau = NA_real_, p = NA_real_,
    stringsAsFactors = FALSE
  )
}

# Subclonal fraction from the MAGOS purity-adjusted clustering: drop the highest-VAF
# and lowest-min clusters, then take 2 * the size-weighted mean VAF of the rest.
.stageC_magosp <- function(result, purity) {
  purity <- min(purity, 1)
  result$vaf.1 <- result$vaf.1 * (2 - purity) / (2 * result$vaf.1 * (1 - purity) + purity)
  result$count <- 1
  meanvaf <- stats::aggregate(result$vaf.1, list(result$colors), mean)
  maxvaf  <- stats::aggregate(result$vaf.1, list(result$colors), max)
  minvaf  <- stats::aggregate(result$vaf.1, list(result$colors), min)
  sumvaf  <- stats::aggregate(result$count,  list(result$colors), sum)
  mago    <- data.frame(max = maxvaf$x, min = minvaf$x,
                        vaf = meanvaf$x, sum = sumvaf$x)
  mago <- mago[-which.max(mago$vaf), , drop = FALSE]
  mago <- mago[-which.min(mago$min), , drop = FALSE]
  filt <- mago[mago$sum >= 0, , drop = FALSE]
  if (nrow(filt) > 0) sum(filt$vaf * filt$sum) / sum(filt$sum) * 2 else 0
}

# Final pick adjustment: choose the mu/p candidate closest to the MAGOS-derived
# subclonal fraction, then derive picks/pickt1/picktend.
.stageC_adjust_p <- function(row, magosp, beta = 0.9, cut = 0.3) {
  row$goodp <- magosp
  if (is.na(row$interp) && is.na(row$fitp)) return(row)
  row$pickp.alt <- if (is.na(row$intermu)) row$fitp
    else if (row$mupick == row$intermu) dplyr::coalesce(row$interp, row$fitp)
    else dplyr::coalesce(row$fitp, row$interp)
  cand <- c(row$fitp, row$interp, row$backp)
  cand <- cand[!is.na(cand)]
  if (length(cand) == 0) return(row)
  row$pickp.close <- cand[which.min(abs(cand - row$goodp))]
  row$mupick.new <- if (is.na(row$fitp)) row$intermu
    else if (abs(row$pickp.close - row$fitp) < 0.01)
      dplyr::coalesce(row$fitmu, row$intermu)
    else dplyr::coalesce(row$intermu, row$fitmu)
  row$mupick.choose <- ifelse(round(abs(row$pickp.alt - row$goodp), 2) >= cut,
                              row$mupick.new, row$mupick)
  row$mupick <- row$mupick.choose
  row$pickp  <- ifelse(row$mupick.choose == row$mupick.new,
                       row$pickp.close, row$pickp.alt)
  row$pickp  <- ifelse(round(abs(row$goodp - row$pickp), 2) >= cut,
                       row$goodp / 2 + row$pickp / 2, row$pickp)
  row$picks  <- if (is.na(row$intermu)) row$fits
                else if (row$mupick == row$intermu) row$inters else row$fits
  row$pickt1 <- if (is.na(row$intermu)) row$fit_len_diff / row$mupick
                else if (row$mupick == row$intermu) row$inter_len_diff / row$mupick
                else row$fit_len_diff / row$mupick
  row$picktend <- ((log(row$pickp) - log(1 - row$pickp)) / (log(2) * beta)
                   + (1 + row$picks) * row$pickt1) / row$picks
  # picktend must be a physical timeline: > 0, >= pickt1, and within [1, 1000].
  # Otherwise there is no valid subclone end-time, so set it to NA.
  pt <- row$picktend
  if (is.na(pt) || pt <= 0 || pt < 1 || pt > 1000 ||
      (!is.na(row$pickt1) && pt < row$pickt1)) {
    row$picktend <- NA_real_
  }
  row
}

# Translate (name, mu, s, t1, tend, p) to (name, mu, s, emergence_time, tau, p).
#   emergence_time = t1
#   tau            = tend / t1
.normalize_output_schema <- function(df) {
  if (!all(c("t1", "tend") %in% names(df))) return(df)
  data.frame(
    name = df$name,
    mu = df$mu,
    s = df$s,
    emergence_time = df$t1,
    tau = ifelse(is.na(df$t1) | df$t1 == 0, NA_real_, df$tend / df$t1),
    p = df$p,
    stringsAsFactors = FALSE
  )
}

.onLoad <- function(libname, pkgname) {
  ## Growth models
  .register_builtin_registries()

  ## Mu methods
  register_mu_method("slope", slope_method)
  register_mu_method("peak",  peak_method)

  ## Estimators
  register_estimator("fit",    run_fit)
  register_estimator("bac",    run_bac)
  register_estimator("normal", run_normal)

  ## Input adapters
  register_adapter("magos", adapt_magos)
  register_adapter("vcf",   adapt_vcf)
  register_adapter("raw",   adapt_raw)

  invisible(NULL)
}
