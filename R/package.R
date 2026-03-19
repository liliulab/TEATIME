#' TEATIME: Estimating EvolutionAry Events Through Single-Timepoint Sequencing
#'
#' TEATIME derives evolutionary parameters — mutation rate, selection
#' coefficient, and subclone emergence timing — from cancer sequencing data
#' by analysing the composition of variant allele frequency (VAF) clusters.
#'
#' @section Main function:
#' The single public entry point is [TEATIME.run()].
#'
#' @section Extending TEATIME:
#' Custom growth models, estimators, mu methods, and input adapters can be
#' registered at runtime via the registry API:
#' - [register_growth_model()]
#' - [register_estimator()]
#' - [register_mu_method()]
#' - [register_adapter()]
#'
#' @keywords internal
"_PACKAGE"

## Suppress R CMD check NOTE for dplyr/data.table non-standard evaluation columns
utils::globalVariables(c(
  ## dplyr / pipeline column references
  "colors", "vaf.1", "cluster", "freq", "vaf", "count",
  "mean_vaf", "min_vaf", "new_s", "alt_s",
  ## post-process column references
  "fitmu", "fitmu_candidate", "intermu", "intermu_candidate",
  "fitdiff", "fitdiff2", "interdiff", "interdiff2"
))

#' @import igraph
#' @import Matrix
#' @importFrom dplyr group_by summarise_all summarize pull select contains funs
#' @importFrom stats aggregate dbeta pbeta quantile rbeta sd var
#' @importFrom utils combn flush.console write.table
NULL
