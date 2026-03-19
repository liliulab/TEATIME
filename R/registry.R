`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.teatime_registry_env <- new.env(parent = emptyenv())
.teatime_registry_env$growth <- new.env(parent = emptyenv())
.teatime_registry_env$estimator <- new.env(parent = emptyenv())
.teatime_registry_env$mu_method <- new.env(parent = emptyenv())
.teatime_registry_env$adapter <- new.env(parent = emptyenv())

.validate_fn_signature <- function(fn, expected_args, label) {
  if (!is.function(fn)) {
    stop(sprintf("%s must be a function.", label), call. = FALSE)
  }

  fn_args <- names(formals(fn))
  if (length(fn_args) < length(expected_args) ||
      !identical(fn_args[seq_along(expected_args)], expected_args)) {
    stop(
      sprintf(
        "%s must start with arguments: %s",
        label,
        paste(expected_args, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.register_entry <- function(registry_env, name, fn, expected_args, label) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  .validate_fn_signature(fn, expected_args, label)
  assign(name, fn, envir = registry_env)
  invisible(fn)
}

.get_entry <- function(registry_env, name, label) {
  if (!exists(name, envir = registry_env, inherits = FALSE)) {
    stop(sprintf("Unknown %s: %s", label, name), call. = FALSE)
  }
  get(name, envir = registry_env, inherits = FALSE)
}

#' Register a custom growth model
#'
#' Adds a growth model function to the TEATIME registry so it can be referenced
#' by name in [TEATIME.run()] via the `growth_model` argument.
#'
#' The function must accept at least `i` (cell divisions, integer), `p` (clonal
#' proportion, numeric), and `beta` (survival rate, numeric) as named arguments,
#' and must return a single numeric VAF value. If the function also declares a
#' `ctx` argument, the full pipeline context list is passed automatically —
#' use this to access any additional model parameters stored in `ctx$extra`.
#'
#' @param name A single character string naming the model (e.g. `"exponential"`).
#' @param fn A function accepting at least `i`, `p`, `beta` (named arguments,
#'   any order). Optionally declare a `ctx` argument to receive the full
#'   pipeline context for models that require extra parameters.
#'   Must return a single numeric VAF value.
#'
#' @return `fn` invisibly.
#' @export
#'
#' @examples
#' # Simple model — only needs i, p, beta
#' register_growth_model("exponential", function(i, p, beta) {
#'   p / 2 + (1 - p) / (2 * exp(log(2) * beta * i))
#' })
#'
#' # Model with extra parameters — access them via ctx$extra
#' register_growth_model("logistic", function(i, p, beta, ctx) {
#'   K <- ctx$extra$carrying_capacity %||% 1
#'   p / 2 + (1 - p) / (2 * (1 + (K / p - 1) * exp(-beta * i)))
#' })
register_growth_model <- function(name, fn) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  if (!is.function(fn)) stop("growth model must be a function.", call. = FALSE)
  assign(name, fn, envir = .teatime_registry_env$growth)
  invisible(fn)
}

#' Register a custom mutation rate estimator
#'
#' Adds an estimator function to the TEATIME registry. Estimators are the
#' internal workers that produce candidate (mu, s, t1, p) tuples for a given
#' clonal proportion.
#'
#' @param name A single character string naming the estimator.
#' @param fn A function whose first two arguments are `ctx` (the pipeline
#'   context list) and `p_thre` (the p-value threshold). Must return a named
#'   list with at least `all` (a data frame of candidates) and `select`
#'   (a one-row summary).
#'
#' Built-in estimators (simplified logic):
#' ```r
#' # "fit" — scans clonal fractions and estimates mu + s for the fit case
#' register_estimator("fit", function(ctx, p_thre) {
#'   results <- lapply(seq(0.3, 0.9, by = 0.1), function(p) {
#'     mu_df <- get_mu_method("slope")(p, ctx$main_cluster_vaf, ctx, p_thre)
#'     data.frame(p = p, mu = mu_df$mu, s = NA_real_, cell.div = NA_real_)
#'   })
#'   all <- do.call(rbind, results)
#'   list(all = all, select = all[which.max(all$mu), , drop = FALSE])
#' })
#'
#' # "bac" — background case: uses second cluster VAF to anchor the estimate
#' # "normal" — normal case: uses all subclonal clusters combined
#' ```
#'
#' Example — add a new estimator that uses a fixed known clonal fraction:
#' ```r
#' register_estimator("fixed_p", function(ctx, p_thre) {
#'   p     <- ctx$extra$known_p
#'   mu_df <- get_mu_method("peak")(p, ctx$main_cluster_vaf, ctx, p_thre)
#'   row   <- data.frame(p = p, mu = mu_df$mu, s = NA_real_, cell.div = NA_real_)
#'   list(all = row, select = row)
#' })
#' ```
#'
#' @return `fn` invisibly.
#' @export
register_estimator <- function(name, fn) {
  .register_entry(
    registry_env = .teatime_registry_env$estimator,
    name = name,
    fn = fn,
    expected_args = c("ctx", "p_thre"),
    label = "estimator"
  )
}

#' Register a custom mu estimation method
#'
#' Adds a mu-estimation method to the TEATIME registry. Mu methods estimate
#' the mutation rate for a given clonal proportion and VAF set. The function
#' receives `(p, vaf_set, ctx, p_thre)` and must return a `data.frame` with
#' at least columns `mu`, `p`, and `score`.
#'
#' Built-in methods (simplified logic):
#' ```r
#' # "slope" — estimates mu from the slope of the sorted VAF distribution
#' register_mu_method("slope", function(p, vaf_set, ctx, p_thre) {
#'   sorted <- sort(vaf_set)
#'   x      <- seq_along(sorted) / length(sorted)
#'   slope  <- coef(lm(sorted ~ x))[2]
#'   mu_est <- abs(slope) * ctx$depth * 2
#'   data.frame(mu = mu_est, p = p, score = abs(slope))
#' })
#'
#' # "peak" — estimates mu from the position of the dominant VAF peak
#' register_mu_method("peak", function(p, vaf_set, ctx, p_thre) {
#'   peak_vaf <- as.numeric(names(which.max(table(round(vaf_set, 2)))))
#'   mu_est   <- ctx$depth * abs(peak_vaf - p / 2) * 4
#'   data.frame(mu = mu_est, p = p, score = length(vaf_set))
#' })
#' ```
#'
#' Example — add a new method based on median VAF deviation:
#' ```r
#' register_mu_method("median_ratio", function(p, vaf_set, ctx, p_thre) {
#'   mu_est <- max(1, ctx$depth * abs(median(vaf_set) - p / 2) * 4)
#'   data.frame(mu = mu_est, p = p, score = 1)
#' })
#' ```
#'
#' @param name A single character string naming the method.
#' @param fn A function whose first four arguments are `p`, `vaf_set`, `ctx`,
#'   and `p_thre`.
#'
#' @return `fn` invisibly.
#' @export
register_mu_method <- function(name, fn) {
  .register_entry(
    registry_env = .teatime_registry_env$mu_method,
    name = name,
    fn = fn,
    expected_args = c("p", "vaf_set", "ctx", "p_thre"),
    label = "mu method"
  )
}

#' Register a custom input format adapter
#'
#' Adds an adapter function to the TEATIME registry. Adapters normalise raw
#' input data into the standardised `data.frame` with columns `vaf.1`,
#' `depth.1`, and `colors` that the pipeline requires.
#'
#' Built-in adapters (simplified logic):
#' ```r
#' # "magos" — corrects VAF for tumour purity from a MAGOS result object
#' register_adapter("magos", function(input) {
#'   purity        <- min(input$purity, 1)
#'   vafdata       <- input$result
#'   vafdata$vaf.1 <- vafdata$vaf.1 * (2 - purity) /
#'                    (2 * vafdata$vaf.1 * (1 - purity) + purity)
#'   vafdata[, c("vaf.1", "depth.1", "colors")]
#' })
#'
#' # "vcf"  — runs MAGOS internally on a REF/ALT/CN data.frame
#' # "raw"  — passes through a data.frame already in vaf.1/depth.1/colors form
#' ```
#'
#' Example — add a new adapter for a custom TSV with different column names:
#' ```r
#' register_adapter("custom_tsv", function(input) {
#'   data.frame(
#'     vaf.1   = input$frequency,
#'     depth.1 = input$coverage,
#'     colors  = input$cluster
#'   )
#' })
#'
#' my_data <- read.table("clusters.tsv", header = TRUE, sep = "\t")
#' TEATIME.run(my_data, input_format = "custom_tsv", ...)
#' ```
#'
#' @param name A single character string naming the format (e.g. `"magos"`).
#' @param fn A function whose first argument is `input` (the raw input object).
#'   Must return a `data.frame` with columns `vaf.1`, `depth.1`, and `colors`.
#'
#' @return `fn` invisibly.
#' @export
register_adapter <- function(name, fn) {
  .register_entry(
    registry_env = .teatime_registry_env$adapter,
    name = name,
    fn = fn,
    expected_args = c("input"),
    label = "adapter"
  )
}

get_growth_model <- function(name) {
  .get_entry(.teatime_registry_env$growth, name, "growth model")
}

get_estimator <- function(name) {
  .get_entry(.teatime_registry_env$estimator, name, "estimator")
}

get_mu_method <- function(name) {
  .get_entry(.teatime_registry_env$mu_method, name, "mu method")
}

get_adapter <- function(name) {
  .get_entry(.teatime_registry_env$adapter, name, "adapter")
}

list_growth_models <- function() {
  sort(ls(envir = .teatime_registry_env$growth, all.names = TRUE))
}

list_estimators <- function() {
  sort(ls(envir = .teatime_registry_env$estimator, all.names = TRUE))
}

list_mu_methods <- function() {
  sort(ls(envir = .teatime_registry_env$mu_method, all.names = TRUE))
}

list_adapters <- function() {
  sort(ls(envir = .teatime_registry_env$adapter, all.names = TRUE))
}

vaf_at_div <- function(i, p, ctx) {
  beta <- ctx$beta
  if (is.null(beta) || !is.numeric(beta) || length(beta) != 1L || is.na(beta)) {
    stop("ctx$beta must be a single numeric value.", call. = FALSE)
  }

  model_fn <- get_growth_model(ctx$growth_model %||% "exponential")
  fn_args   <- names(formals(model_fn))

  # Always supply i, p, beta; also supply ctx if the function declares it
  call_args <- list(i = i, p = p, beta = beta)
  if ("ctx" %in% fn_args) call_args$ctx <- ctx

  do.call(model_fn, call_args)
}

.register_builtin_growth_models <- function() {
  register_growth_model(
    "exponential",
    function(i, p, beta) {
      p / 2 + (1 - p) / (2 * exp(log(2) * beta * i))
    }
  )
}

.register_builtin_registries <- function() {
  .register_builtin_growth_models()
}

