.teatime_v1_env <- NULL

# Source the embedded af3e64d v1 reference (TEATIME.r + MAGOS.R) into a private
# env on first request. Default-mode TEATIME.run dispatches into this env so the
# default-mode output is bit-identical to the v1 release. The fast-mode pipeline
# remains the v2 modular implementation.
.load_v1_ref <- function() {
  if (!is.null(.teatime_v1_env)) return(.teatime_v1_env)
  env <- new.env(parent = asNamespace("TEATIME"))
  magos_path <- system.file("extdata/v1_reference/MAGOS.R",  package = "TEATIME")
  teatime_path <- system.file("extdata/v1_reference/TEATIME.r", package = "TEATIME")
  if (!nzchar(magos_path) || !nzchar(teatime_path)) {
    stop("TEATIME: v1 reference source missing from inst/extdata/v1_reference/",
         call. = FALSE)
  }
  suppressMessages(sys.source(magos_path, envir = env, keep.source = FALSE))
  src <- readLines(teatime_path)
  # The bundled v1 TEATIME.r begins with source(file.path("R","MAGOS.R")); we've
  # already loaded MAGOS.R into the same env above, so neutralize that line.
  src <- gsub('^source\\(file\\.path\\("R", "MAGOS\\.R"\\)\\)', '# (v1 MAGOS already loaded)', src)
  tmp <- tempfile(fileext = ".r"); on.exit(unlink(tmp), add = TRUE)
  writeLines(src, tmp)
  suppressMessages(sys.source(tmp, envir = env, keep.source = FALSE))
  utils::assignInMyNamespace(".teatime_v1_env", env)
  env
}

# Default-mode dispatcher: translate v2's TEATIME.run args to the v1 reference
# signature, invoke v1's TEATIME.run, and read its .final.txt back as the
# result. v1 writes name/mu/s/t1/tend/p; v2 default mode returns those same
# columns so downstream consumers see exactly v1's output.
.run_v1_default <- function(input, beta, depth, p_thre,
                            output_folder, output_prefix,
                            id, write_final, seed, debug) {
  v1 <- .load_v1_ref()
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  }
  v1_seed <- if (is.na(seed)) NA else as.integer(seed)
  call_args <- list(
    input.file    = input,
    beta          = beta,
    depth         = depth,
    p_thre        = p_thre,
    magos_object  = TRUE,
    output.folder = output_folder,
    output.prefix = output_prefix,
    id            = id,
    steps         = 0:5,
    write_final   = isTRUE(write_final),
    debug_mode    = isTRUE(debug),
    purity_set    = 0,
    seed          = v1_seed
  )
  if (isTRUE(debug)) {
    do.call(v1$TEATIME.run, call_args)
  } else {
    suppressWarnings(suppressMessages(
      utils::capture.output(do.call(v1$TEATIME.run, call_args))
    ))
  }
  final_path <- file.path(output_folder, paste0(output_prefix, ".final.txt"))
  if (file.exists(final_path)) {
    v1_df <- utils::read.table(final_path, header = TRUE, sep = "\t",
                               stringsAsFactors = FALSE)
    return(.v1_to_v2_schema(v1_df, final_path, write_final = isTRUE(write_final)))
  }
  data.frame(name = id, mu = NA_real_, s = NA_real_,
             emergence_time = NA_real_, tau = NA_real_, p = NA_real_,
             stringsAsFactors = FALSE)
}

# Translate v1's (name, mu, s, t1, tend, p) to v2's
# (name, mu, s, emergence_time, tau, p). Mapping (from pipeline.R post_process):
#   emergence_time = t1
#   tau            = tend / t1
# If write_final is TRUE, also rewrite .final.txt so downstream readers see the
# v2 schema regardless of which mode produced the file.
.v1_to_v2_schema <- function(v1_df, final_path, write_final) {
  if (!all(c("t1", "tend") %in% names(v1_df))) return(v1_df)
  v2_df <- data.frame(
    name = v1_df$name,
    mu = v1_df$mu,
    s = v1_df$s,
    emergence_time = v1_df$t1,
    tau = ifelse(is.na(v1_df$t1) | v1_df$t1 == 0,
                 NA_real_, v1_df$tend / v1_df$t1),
    p = v1_df$p,
    stringsAsFactors = FALSE
  )
  if (isTRUE(write_final)) {
    hdr <- paste0(
      "## name: sample ID | mu: mutation rate | s: selection coefficient | ",
      "emergence_time: emergence time of the subclone | ",
      "tau: subclone expansion score | p: subclonal fraction\n")
    tryCatch({
      writeLines(hdr, con = final_path)
      suppressWarnings(utils::write.table(
        v2_df, file = final_path, sep = "\t",
        row.names = FALSE, quote = FALSE, append = TRUE))
    }, error = function(e2) NULL)
  }
  v2_df
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
