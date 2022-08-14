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