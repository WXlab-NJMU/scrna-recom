#' TryCatch Using WithCallingHandlers
#' 
TryCatchWithCallingHandlers <- function(expression){
  withCallingHandlers(
    expression,
    warning = function(w){ 
      message("warning:\n", w)
      lobstr::cst()
      invokeRestart("muffleWarning")},
    error = function(e){ 
      lobstr::cst()
      message("error:\n", e) },
    finally = { message("...") })
}

#' Merge Data to List of Seurat Object
#'
#' @description
#' Merge data included in the csvfile to one list
#'
#' @param csv The header must be "project,path"
#' @return Returns List consisted of multiple Seurat objects
#'
#' @concept utility
#' @importFrom Seurat Read10X CreateSeuratObject
#' @export
#'
MergeFileData <- function(csv) {
  inputs <- read.csv(csv)
  data.list <- apply(inputs, 1, FUN = function(item) {
    data <- Seurat::Read10X(data.dir = item[["path"]]) %>%
      Seurat::CreateSeuratObject(project = item[["project"]], min.cells = 3, min.features = 200)
    data[["dataset"]] <- item[["project"]]
    return(data)
    })
  return(data.list)
}

#' Read files to Seurat Object
#'
#' @description
#' Merge data included in the csvfile to one list
#'
#' @param csv The header must be "project,path"
#' @return Returns List consisted of multiple Seurat objects
#'
#' @concept utility
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom SeuratDisk Connect Convert LoadH5Seurat
#' @importFrom tools file_ext
#' @export
#'
ReadtoSeuratObject <- function(infile) {
  message("Reading: ", infile)
  if (file_test('-f', infile)) {
    extension <- tools::file_ext(infile)
    # Convert from SingleCellExperiment
    if (extension == "rds") {
      input <- readRDS(file = infile)
      final <- Seurat::as.Seurat(input, counts = "counts", data = "logcounts")
    } else if (extension == "loom") {
      input <- SeuratDisk::Connect(filename = infile, mode = "r")
      final <- Seurat::as.Seurat(input)
    } else if (extension == "h5ad") {
      SeuratDisk::Convert(infile, dest = "h5seurat", overwrite = TRUE)
      converted.file = sub(".h5ad",'.h5seurat', infile)
      final <- SeuratDisk::LoadH5Seurat(converted.file)
    } else {
      warning("Input is not supported !!!")
    }
  } else if (file_test('-d', infile)) {
    input <- Seurat::Read10X(infile, strip.suffix=TRUE)
    final <- Seurat::CreateSeuratObject(input)
  }
  final
}

#' Convert Seurat Object to Other Format
#'
#' @description
#' convert a seurat object to other format object
#'
#' @param csv The header must be "project,path"
#' @return Returns List consisted of multiple Seurat objects
#'
#' @concept utility
#' @importFrom Seurat as.SingleCellExperiment 
#' @importFrom SeuratDisk SaveH5Seurat Convert as.loom
#' @export
#' 
OutFromSeuratToOthers <- function(obj, format, outfile = NULL) {
  if (format == "SingleCellExperiment") {
    final <- Seurat::as.SingleCellExperiment(obj)
  } else if (format == "loom") {
    message("always remember to close loom files when done: `final$close_all()` ")
    final <- SeuratDisk::as.loom(obj, filename = outfile)
    # always remember to close loom files when done: `final$close_all()`
  } else if (format == "AnnData") {
    outfile <- paste0(outfile, '.h5seurat')
    # convert to h5seurat
    SeuratDisk::SaveH5Seurat(obj, filename = outfile, overwrite=TRUE)
    # convert to h5ad
    SeuratDisk::Convert(outfile, dest = "h5ad", overwrite=TRUE)
    final <- sub('.h5seurat', ".h5ad", outfile)
  }
  final
}

