#' @include utils.R
NULL

#' Seurat Integration CCA
#'
#' @description
#' merge multiple projects data into one object
#' 
#' @seealso 
#' https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
#' 
#' @param csv The header is "project,path", the path must included matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir
#'
#' @rdname integration
#' @order 1
#' @import Seurat
#' @export
#' @concept integration
#' @method SeuratIntegration default
#' @concept integration
#' 
SeuratIntegration.default <- function(
  csv,
  outdir,
  ...
) {
  projects <- MergeFileData(csv) %>% 
    lapply(FUN = function(x) {
      NormalizeData(x) %>% 
        FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
  features <- SelectIntegrationFeatures(object.list = projects)
  combined.data <- projects %>% 
    FindIntegrationAnchors(anchor.features = features) %>% 
    IntegrateData() %>% 
    ScaleData() %>% 
    RunPCA(npcs = 30, verbose = FALSE) %>% 
    RunUMAP(reduction = "pca", dims = 1:30)
  dir.create(outdir, recursive = TRUE)
  saveRDS(combined.data, file.path(outdir, "integration.seurat-cca.rds"))
  plot(file.path(outdir, "integration.seurat-cca.pdf"))
  p1 <- DimPlot(combined.data, reduction = "pca", group.by = "dataset")
  print(p1)
  p2 <- DimPlot(combined.data, reduction = "umap", group.by = "dataset") 
  print(p2)
  dev.off()
  return(combined.data)
}

#' Seurat SCTransform Integration 
#'
#' @description
#' merge multiple projects data into one object
#' 
#' @seealso 
#' https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
#'
#' @param csv The header is "project,path", the path must included matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir
#'
#' @rdname integration
#' @order 2
#' @import Seurat
#' @export
#' @concept integration
#' @method Integration SCTransform
#' @concept integration
#' 
Integration.SCTransform <- function(
  csv,
  outdir,
  ...  
){
  projects <- MergeFileData(csv) %>% lapply(SCTransform)
  features <- Seurat::SelectIntegrationFeatures(object.list = projects, nfeatures = 5000)
  pojects <- Seurat::PrepSCTIntegration(object.list = projects, anchor.features = features)
  anchors <- Seurat::FindIntegrationAnchors(object.list = pojects, normalization.method = "SCT", anchor.features = features)
  combined.data <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT") %>%
    Seurat::RunPCA(verbose = FALSE)  %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:30)
  dir.create(outdir, recursive = TRUE)
  saveRDS(combined.data, file.path(outdir, "integration.sctransform.rds"))
  plot(file.path(outdir, "integration.sctransform.pdf"))
  p1 <- Seurat::DimPlot(combined.data, reduction = "pca", group.by = "dataset")
  print(p1)
  p2 <- Seurat::DimPlot(combined.data, reduction = "umap", group.by = "dataset") 
  print(p2)
  dev.off()
  return(combined.data)
}

#' Seurat Integration Using Harmony 
#'
#' @description
#' merge multiple projects data into one object
#' 
#' @source 
#' https://portals.broadinstitute.org/harmony/articles/quickstart.html
#' 
#' @param csv The header is "project,path", the path must included matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir
#'
#' @rdname integration
#' @order 3
#' @import Seurat
#' @import harmony
#' @export
#' @concept integration
#' @method Integration Harmony
#' @concept integration
#' 
Integration.Harmony <- function(csv, outdir){
  projects <- MergeFileData(csv)
  combined.data <- merge(projects[[1]], tail(projects, length(projects)-1)) %>%
    Seurat::NormalizeData()  %>% 
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    Seurat::ScaleData() %>% 
    Seurat::RunPCA(npcs = 20) %>%
    harmony::RunHarmony(merged, "dataset")
  dir.create(outdir, recursive = TRUE)
  saveRDS(combined.data, file.path(outdir, "integration.harmony.rds"))
  plot(file.path(outdir, "integration.harmony.pdf"))
  p1 <- Seurat::DimPlot(object = combined.data, reduction = "pca", group.by = "dataset")
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", group.by = "dataset")
  print(p1 + p2)
  p3 <- Seurat::DimPlot(object = combined.data, reduction = "harmony", group.by = "dataset")
  p4 <- Seurat::VlnPlot(object = combined.data, features = "harmony_1", group.by = "dataset")
  print(p3 + p4)
  dev.off()
  return(combined.data)
}


#' Seurat Integration Using Harmony 
#'
#' @description
#' merge multiple projects data into one object
#' 
#' @source 
#' https://github.com/welch-lab/liger
#' https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html
#' 
#' 
#' @param csv The header is "project,path", the path must included matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir
#'
#' @rdname integration
#' @order 3
#' @import rliger
#' @import dplyr
#' @import Seurat
#' @import SeuratWrappers
#' @export
#' @concept integration
#' @method Integration Liger
#' @concept integration
#' 
Integration.Liger <- function(csv, outdir){
  projects <- MergeFileData(csv)
  combined.data <- merge(projects[[1]], tail(projects, length(projects)-1)) %>%
    Seurat::NormalizeData()  %>% 
    Seurat::FindVariableFeatures() %>% 
    Seurat::ScaleData(split.by = "dataset", do.center = FALSE) %>% 
    SeuratWrappers::RunOptimizeALS(k = 20, lambda = 5, split.by = "dataset")  %>% 
    SeuratWrappers::RunQuantileNorm(split.by = "dataset") %>%
    Seurat::FindNeighbors(reduction = "iNMF", dims = 1:20) %>%
    Seurat::FindClusters(resolution = 0.55)
  combined.data <- Seurat::RunUMAP(dims = 1:ncol(combined.data[["iNMF"]]), reduction = "iNMF")
  dir.create(outdir, recursive = TRUE)
  saveRDS(combined.data, file.path(outdir, "integration.liger.rds"))
  plot(file.path(outdir, "integration.liger.pdf"))
  p1 <- Seurat::DimPlot(combined.data, group.by = c("dataset", "ident"), ncol = 3)
  print(p1)
  dev.off()
  return(combined.data)
}
