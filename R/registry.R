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
#' Adds a mu-estimation method to the TEATIME registry. Mu methods compute
#' candidate mutation-rate values for a given clonal proportion and VAF set.
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

