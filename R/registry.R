`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.teatime_registry_env <- new.env(parent = emptyenv())
.teatime_registry_env$growth    <- new.env(parent = emptyenv())
.teatime_registry_env$estimator <- new.env(parent = emptyenv())
.teatime_registry_env$mu_method <- new.env(parent = emptyenv())
.teatime_registry_env$adapter   <- new.env(parent = emptyenv())
.teatime_registry_env$ensemble  <- new.env(parent = emptyenv())

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

#' Register an ensemble strategy
#'
#' Ensembles wrap a per-try `check_fn(ctx, p_thre)` that returns a candidate
#' data frame and aggregate across `n_tries` repeated calls. Used by
#' `run_fit` / `run_bac` / `run_normal` so the multi-try aggregation policy
#' can be swapped without touching estimator code.
#'
#' @param name A character string.
#' @param fn A function with signature
#'   `function(check_fn, ctx, p_thre, n_tries)`.
#' @return `fn` invisibly.
#' @export
register_ensemble <- function(name, fn) {
  .register_entry(.teatime_registry_env$ensemble, name, fn,
                  c("check_fn", "ctx", "p_thre"), "ensemble")
}

get_growth_model <- function(name) .get_entry(.teatime_registry_env$growth,    name, "growth model")
get_estimator    <- function(name) .get_entry(.teatime_registry_env$estimator, name, "estimator")
get_mu_method    <- function(name) .get_entry(.teatime_registry_env$mu_method, name, "mu method")
get_adapter      <- function(name) .get_entry(.teatime_registry_env$adapter,   name, "adapter")
get_ensemble     <- function(name) .get_entry(.teatime_registry_env$ensemble,  name, "ensemble")

list_growth_models <- function() sort(ls(envir = .teatime_registry_env$growth,    all.names = TRUE))
list_estimators    <- function() sort(ls(envir = .teatime_registry_env$estimator, all.names = TRUE))
list_mu_methods    <- function() sort(ls(envir = .teatime_registry_env$mu_method, all.names = TRUE))
list_adapters      <- function() sort(ls(envir = .teatime_registry_env$adapter,   all.names = TRUE))
list_ensembles     <- function() sort(ls(envir = .teatime_registry_env$ensemble,  all.names = TRUE))

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

# Default ensemble. Runs `check_fn(ctx, p_thre)` `n_tries` times and aggregates
# the per-try candidate frames. Stochastic re-fits inside `check_fn` (e.g.
# `RBesT::automixfit`) give different candidates each pass; the ensemble picks:
#   mu1 = mu of the highest-scored row across the maxima of each try (main)
#   mu2 = max over a per-try mix of medians and maxima (broader pool)
#   mu3 = the same broader pool's row whose try had the lowest BIC
#   p1  = p paired with the mu1 pick
#   p2  = p of the row closest to the per-try median
# `check_fn` must pre-sort so its row[1] is the per-try best candidate.
.register_builtin_ensembles <- function() {
  register_ensemble("score_3try", function(check_fn, ctx, p_thre, n_tries = 3) {
    mulist  <- numeric(0)
    mulist_pool <- numeric(0)
    biclist <- numeric(0)
    plist   <- numeric(0)
    plist1  <- numeric(0)
    data.all <- data.frame()

    for (try_i in seq_len(n_tries)) {
      check <- tryCatch(check_fn(ctx, p_thre), error = function(e) NULL)
      if (is.null(check) || !is.data.frame(check) || nrow(check) == 0) next

      check$try <- try_i
      data.all  <- rbind(data.all, check)

      if (!is.na(check$mu[1])) {
        mulist_pool <- c(mulist_pool, stats::median(check$mu, na.rm = TRUE))
        mulist_pool <- c(mulist_pool, max(check$mu, na.rm = TRUE))
        max_mu_row  <- which(check$mu == max(check$mu, na.rm = TRUE))[1]
        if (!is.na(max_mu_row)) biclist <- c(biclist, check$bic[max_mu_row])

        top    <- check[1, , drop = FALSE]
        mulist <- c(mulist, top$mu)
        plist  <- c(plist,  top$p)

        med_ref <- stats::median(check$mu, na.rm = TRUE)
        nearest <- which.min(abs(check$mu - med_ref))
        plist1  <- c(plist1, check$p[nearest])
      }
    }

    if (length(mulist) > 0) {
      idx_main <- which.max(mulist)
      mu1      <- mulist[idx_main]
      p1       <- plist[idx_main]

      idx_med  <- which.max(mulist_pool)
      mu2      <- mulist_pool[idx_med]
      p2       <- plist1[min(length(plist1), idx_med)]

      idx_bic  <- if (length(biclist) > 0) which.min(biclist) else integer(0)
      mu3      <- if (length(idx_bic)) mulist_pool[idx_bic] else NA_real_
    } else {
      mu1 <- NA_real_; mu2 <- NA_real_; mu3 <- NA_real_
      p1  <- NA_real_; p2  <- NA_real_
    }

    list(
      select = data.frame(
        mu  = mu1,
        up  = NA_real_,
        p1  = p1,
        name = ctx$id,
        stringsAsFactors = FALSE
      ),
      candidate = data.frame(
        mu1 = mu1, mu2 = mu2, mu3 = mu3,
        p1  = p1,  p2  = p2,
        stringsAsFactors = FALSE
      ),
      all = data.all
    )
  })
}

.register_builtin_registries <- function() {
  .register_builtin_growth_models()
  .register_builtin_ensembles()
}
