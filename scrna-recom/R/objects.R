check_input <- function(object) {
  errors <- character()
  # check indir and included files
  if (dir.exists(object@indir)) {
    matrix_file <- file.path(object@indir, "matrix.mtx")
    if (! file.exists(matrix_file)) {
      errors <- c(errors, paste0(matrix_file, " not existed!!!")) 
    }
    barcode_file <- file.path(object@indir,"barcodes.tsv")
    if (! file.exists(barcode_file)) {
      errors <- c(errors, paste0(barcode_file, " not existed!!! ")) 
    }
  } else {
    errors <- c(errors, paste0(object@indir," not existed!!!"))
  }
  # check project name
  if (is.na(object@project)){
    errors <- c(errors, "project name not defined!!!") 
  }
  if (length(errors) == 0) TRUE else errors
}
#' Input class for basic analysis
#'
#' @slot indir cellranger count matrix directory
#' @slot outdir output directory
#' @slot steps list of steps(default: 1,2,3,4,5,8): 
#'   1-QualityControl, 2-Normalization, 3-FeatureSelection, 
#'   4-Scaling, 5-dimensionalReduction, 6-RemoveDoublet,
#'   7-DifferentialExpression, 8-AutoCellType
#' @slot maxrna qc thresold: max rna counts
#' @slot minrna qc thresold: min rna counts
#' @slot mincell qc thresold: min cell counts
#' @slot maxmt qc thresold: max mt percentage, default: 5
#' @slot project project name
#' @export
#' 
setClass(Class = "Input",
         representation(indir = "character",  outdir = "character", steps = "numeric",
                        maxrna = "numeric", minrna = "numeric", mincell = "numeric",
                        maxmt = "numeric", project = "character",
                        sctransform = "logical", pc = "numeric"),
         prototype(indir = NA_character_,  outdir = NA_character_, steps = c(1,2,3,4,5,8),
                   maxrna = 2500, minrna = 200, mincell = 3, maxmt = 5, project = NA_character_,
                   sctransform = FALSE, pc = NA_integer_),
         #contains = "",
         validity = check_input)
