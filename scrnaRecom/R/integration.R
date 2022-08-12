#' @include preprocess.R
NULL

#' @param x desc
#' @rdname SeuratIntegration
#' @export
#' @method SeuratIntegration default
#' @concept integration
#'
SeuratIntegration.default <- function(
  anchors,
  vars = NULL,
  ...
) {
  inputs = read.csv(csv)
  projects <- lapply(inputs, FUN = function(item) {
                       data <- Read10X(data.dir = item["path"])
                       x <- CreateSeuratObject(counts = data,
                                               project = item["project"],
                                               min.cells = 3,
                                               min.features = 200)
                       x <- NormalizeData(x)
                       x <- FindVariableFeatures(x, selection.method = "vst",
                                                nfeatures = 2000)
                       return(x)
  })
  features <- SelectIntegrationFeatures(object.list = projects)
  anchors <- FindIntegrationAnchors(object.list = projects,
                                    anchor.features = features)
  combined <- IntegrateData(anchorset = anchors)
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  return(combined)
}

SeuratIntegration.SCTransform <- function(csv){
  inputs = read.csv(csv)
  projects <- lapply(inputs, function(item) {
                       data <- Read10X(data.dir = item["path"])
                       x <- CreateSeuratObject(counts = data,
                                               project = item["project"],
                                               min.cells = 3,
                                               min.features = 200)
                       x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
                       x <- SCTransform(x, vars.to.regress = "percent.mt", method = "glmGamPoi")
                       return(x)}
  )
  features <- SelectIntegrationFeatures(object.list = projects, nfeatures = 3000)
  pojects <- PrepSCTIntegration(object.list = projects, anchor.features = features)
  pojects <- lapply(projects, FUN = RunPCA, features = features)
  anchors <- FindIntegrationAnchors(object.list = pojects, normalization.method = "SCT",
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
  combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
}


