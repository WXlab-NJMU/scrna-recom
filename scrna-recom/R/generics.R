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
#'
Integration <- function(csv, outdir, used = "SeuratCCA") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE) 
  }
  UseMethod(generic = 'Integration', 
            object = structure(1, class = used))
}


#' Cell Type Annotation
#'
#' @param input Seurat input object
#' @param outdir The ouput directory for rds and plot
#' @param used Methods for cell type annotation: scCATCH, singleR
#' @return Return a combined object
#' 
#' @concept cell type annotation
#' @export AnnotateCellType
#'
AnnotateCellType <- function(input, outdir, used = "scCATCH") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE) 
  }
  UseMethod(generic = 'AnnotateCellType', 
            object = structure(1, class = used))
}


#' Trajectories
#'
#' @param object The header is "project,path", the path must included matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir The ouput directory for rds and plot
#' @param used Methods for  trajectories: Monocole3, scVelo
#' @return Return a combined object
#' 
#' @concept trajectory
#' @export Trajectories
#'
Trajectories <- function(csv, outdir, used = "Monocole3") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE) 
  }
  UseMethod(generic = 'AnnotateCellType', 
            object = structure(1, class = used))
}


