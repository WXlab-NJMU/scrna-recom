#' @include utils.R
NULL

#' Project Data Integration
#'
#' @param csv The header is "project,path", the path must included matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir The ouput directory for rds and plot
#' @param used Methods for integration: SeuratCCA, SCTransform, Harmony, Liger, default is SeuratCCA
#' @return Return a combined object
#' 
#' @concept integration
#' @export Integration
#'
Integration <- function(csv, outdir, used = "SeuratCCA") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE) 
  }
  obj <- 1
  class(obj) <- used
  print(used)
  UseMethod(generic = 'Integration', object = obj)
}

