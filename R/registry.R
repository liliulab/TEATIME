`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.teatime_registry_env <- new.env(parent = emptyenv())
.teatime_registry_env$growth    <- new.env(parent = emptyenv())
.teatime_registry_env$estimator <- new.env(parent = emptyenv())
.teatime_registry_env$mu_method <- new.env(parent = emptyenv())
.teatime_registry_env$adapter   <- new.env(parent = emptyenv())

.validate_fn_signature <- function(fn, expected_args, label) {
  if (!is.function(fn)) {
    stop(sprintf("%s must be a function.", label), call. = FALSE)
  }
  fn_args <- names(formals(fn))
  if (length(fn_args) < length(expected_args) ||
      !identical(fn_args[seq_along(expected_args)], expected_args)) {
    stop(
      sprintf("%s must start with arguments: %s",
              label, paste(expected_args, collapse = ", ")),
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
#' The function must accept `i` (cell divisions), `p` (subclonal fraction), and
#' `beta` (survival rate) as named arguments and return a single numeric VAF
#' value. Declare an optional `ctx` argument to receive the full pipeline
#' context — use `ctx$extra` to pass in additional model parameters via the
#' `extra` argument of [TEATIME.run()].
#'
#' @param name A single character string naming the model (e.g. `"exponential"`).
#' @param fn A function accepting at least `i`, `p`, `beta` as named arguments,
#'   and optionally `ctx`. Must return a single numeric VAF value.
#'
#' @return `fn` invisibly.
#' @export
#'
#' @examples
#' # Built-in exponential model (registered automatically on load)
#' register_growth_model("exponential", function(i, p, beta) {
#'   p / 2 + (1 - p) / (2 * exp(log(2) * beta * i))
#' })
#'
#' # Logistic model — growth slows as the tumour approaches carrying capacity K
#' # Extra parameters are passed in via extra = list(...) in TEATIME.run()
#' register_growth_model("logistic", function(i, p, beta, ctx) {
#'   K <- ctx$extra$carrying_capacity %||% 1
#'   p / 2 + (1 - p) / (2 * (1 + (K / p - 1) * exp(-beta * i)))
#' })
#'
#' \dontrun{
#' TEATIME.run(
#'   input        = input.file,
#'   input_format = "magos",
#'   growth_model = "logistic",
#'   extra        = list(carrying_capacity = 2)
#' )
#' }
register_growth_model <- function(name, fn) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  if (!is.function(fn)) stop("growth model must be a function.", call. = FALSE)
  assign(name, fn, envir = .teatime_registry_env$growth)
  invisible(fn)
}

# ── Internal-only registries ─────────────────────────────────────────────────

register_estimator <- function(name, fn) {
  .register_entry(.teatime_registry_env$estimator, name, fn,
                  c("ctx", "p_thre"), "estimator")
}

register_mu_method <- function(name, fn) {
  .register_entry(.teatime_registry_env$mu_method, name, fn,
                  c("p", "vaf_set", "ctx", "p_thre"), "mu method")
}

register_adapter <- function(name, fn) {
  .register_entry(.teatime_registry_env$adapter, name, fn,
                  c("input"), "adapter")
}

get_growth_model <- function(name) .get_entry(.teatime_registry_env$growth,    name, "growth model")
get_estimator    <- function(name) .get_entry(.teatime_registry_env$estimator, name, "estimator")
get_mu_method    <- function(name) .get_entry(.teatime_registry_env$mu_method, name, "mu method")
get_adapter      <- function(name) .get_entry(.teatime_registry_env$adapter,   name, "adapter")

list_growth_models <- function() sort(ls(envir = .teatime_registry_env$growth,    all.names = TRUE))
list_estimators    <- function() sort(ls(envir = .teatime_registry_env$estimator, all.names = TRUE))
list_mu_methods    <- function() sort(ls(envir = .teatime_registry_env$mu_method, all.names = TRUE))
list_adapters      <- function() sort(ls(envir = .teatime_registry_env$adapter,   all.names = TRUE))

vaf_at_div <- function(i, p, ctx) {
  beta <- ctx$beta
  if (is.null(beta) || !is.numeric(beta) || length(beta) != 1L || is.na(beta)) {
    stop("ctx$beta must be a single numeric value.", call. = FALSE)
  }
  model_fn  <- get_growth_model(ctx$growth_model %||% "exponential")
  call_args <- list(i = i, p = p, beta = beta)
  if ("ctx" %in% names(formals(model_fn))) call_args$ctx <- ctx
  do.call(model_fn, call_args)
}

.register_builtin_growth_models <- function() {
  register_growth_model("exponential", function(i, p, beta) {
    p / 2 + (1 - p) / (2 * exp(log(2) * beta * i))
  })
}

.register_builtin_registries <- function() {
  .register_builtin_growth_models()
}
