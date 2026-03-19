#' Example data bundled with TEATIME
#'
#' Three files are provided in `inst/extdata` for testing and validation.
#' Access them with [system.file()]:
#'
#' \describe{
#'   \item{`exampledata.rds`}{Raw somatic mutation data (REF/ALT counts) used
#'     as the primary test input.}
#'   \item{`MAGOS.rds`}{Pre-computed MAGOS clustering result derived from
#'     `exampledata.rds`. Use for reproducible validation runs without
#'     re-running MAGOS.}
#'   \item{`TEATIME.final.txt`}{Reference output from the original TEATIME
#'     pipeline on the same data. Used to verify v2 results match within
#'     tolerance.}
#' }
#'
#' @examples
#' magos <- readRDS(system.file("extdata", "MAGOS.rds", package = "TEATIME"))
#' input <- list(purity = magos$purity, result = magos$results)
#' result <- TEATIME.run(input, beta = 0.9, depth = round(mean(magos$results$depth.1)),
#'                       seed = 123, write_final = FALSE)
#' print(result)
#'
#' @name example-data
NULL
