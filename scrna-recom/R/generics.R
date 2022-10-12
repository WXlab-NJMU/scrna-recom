#' @include utils.R
NULL


#' @export
setGeneric("basic_analysis", function(object) {
  standardGeneric("basic_analysis")
})

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
Integration <- function (csv, outdir, project, used = "SeuratCCA", dim) {
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


#' Trajectory Analysis
#'
#' @param object The header is "project,path", the path must included matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir The ouput directory for rds and plot
#' @param used Methods for  trajectories: Monocole3, scVelo
#' @return Return a combined object
#'
#' @concept trajectory
#' @export Trajectory
#'
Trajectory <- function(csv, outdir, used = "Monocole3") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  UseMethod(generic = 'Trajectory',
            object = structure(1, class = used))
}

#' Metabolic Analysis
#'
#' @param object The header is "project,path", the path must included matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir The ouput directory for rds and plot
#' @param used Methods for anlysis: scMetabolism, scFEA
#' @return Return a combined object
#'
#' @concept Metabolic Analysis
#' @export MetabolicAnalysis
#'
MetabolicAnalysis <- function(csv, outdir, used = "scMetabolism") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  UseMethod(generic = 'MetabolicAnalysis',
            object = structure(1, class = used))
}


#' Cell Communication
#'
#' @param infile Input file path
#' @param outdir The ouput directory for rds and plot
#' @param used Methods for anlysis: CellChat
#' @return Return a combined object
#'
#' @concept Cell Communication
#' @export CellCommunication
#'
CellCommunication <- function(infile, outdir, used = "CellChat") {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  UseMethod(generic = 'CellCommunication',
            object = structure(1, class = used))
}
