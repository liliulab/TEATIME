.teatime_ref_env <- NULL

# Lazily source the bundled reference implementation into a private env on first
# use. Default-mode TEATIME.run dispatches into this env; fast mode keeps the
# in-package modular pipeline.
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
