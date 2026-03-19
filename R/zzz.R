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
