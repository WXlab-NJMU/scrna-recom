#' @include utils.R
#' @import dplyr Seurat SeuratWrappers
#' @importFrom grDevices dev.off pdf
#' @importFrom utils read.csv tail
NULL

#' @section Seurat CCA:
#' * source code: <https://github.com/satijalab/seurat>
#' * quick start: <https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1>
#' @md
#'
#' @import Seurat
#' @export
#' @rdname Integration
#' @method Integration SeuratCCA
#' @concept integration
#'
Integration.SeuratCCA <- function(object, outdir, project, used, dim){
  #all.genes <- rownames(object)
  object <- Seurat::NormalizeData(object) %>% 
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 4000)
  projects <- Seurat::SplitObject(object, split.by = "orig.ident")
  anchors <- Seurat::SelectIntegrationFeatures(object.list = projects)
  combined.data <- projects %>%
    Seurat::FindIntegrationAnchors(anchor.features = anchors) %>%
    Seurat::IntegrateData() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = dim) %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:dim) %>%
    Seurat::FindNeighbors() %>%
    Seurat::FindClusters()
  prefix <- file.path(outdir, sprintf("%s.debatch.seurat-cca.dim=%d", project, dim))
  saveRDS(combined.data, paste0(prefix, ".rds"))
  pdf(paste0(prefix, ".pdf"))
  p1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("orig.ident"))
  print(p1)
  p1_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("seurat_clusters"))
  print(p1_1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", group.by = c("orig.ident", "seurat_clusters"))
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p2 + p3)
  p7 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"))
  print(p7)
  p7_1 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"))
  print(p7_1)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"))
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p8 + p9)
  dev.off()
  return(combined.data)
}

#' @section Seurat Large Data:
#' For very large datasets, Seruat provide two options that can improve efficiency and runtimes:
#' * Reciprocal PCA (RPCA): faster and more conservative
#' * Reference-based integration
#' - quick start: <https://satijalab.org/seurat/articles/integration_large_datasets.html>
#' - RPCA: <https://satijalab.org/seurat/articles/integration_rpca.html>
#' @md
#'
#' @param reference Sample IDs to used, dim as integration
#' @import Seurat
#' @export
#' @rdname Integration
#' @method Integration SeuratRPCA
#' @concept integration
#'
Integration.SeuratRPCA <- function(object, outdir, project, used, dim, reference = c(1,2)){
  #projects <- MergeFileData(csv, outdir, mincell, minrna, maxrna, maxmt)
  #projects <- lapply(projects, FUN = function(x) {
  #  Seurat::NormalizeData(x) %>%
  #    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  #})
  object <- object %>% Seurat::NormalizeData() %>% 
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 4000)
  projects <- Seurat::SplitObject(object, split.by = "orig.ident")
  anchors <- Seurat::SelectIntegrationFeatures(object.list = projects)
  # difference with SeruatCCA, run PCA first then
  projects <- lapply(projects, FUN = function(x) {
    x <- ScaleData(x, features = anchors)
    x <- RunPCA(x, npcs = 50, features = anchors)
  })
  anchors <- FindIntegrationAnchors(projects, dims = 1:50,
                                    reference = reference, reduction = "rpca")
  combined.data <- IntegrateData(anchorset = anchors, dims = 1:50) %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = dim) %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:dim) %>%
    Seurat::FindNeighbors() %>%
    Seurat::FindClusters()
  DefaultAssay(combined.data) <- "integrated"
  prefix <- file.path(outdir, sprintf("%s.debatch.seurat-rpca.dim=%d", project, dim))
  saveRDS(combined.data, paste0(prefix, ".rds"))
  pdf(paste0(prefix, ".pdf"))
  p1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("orig.ident"))
  print(p1)
  p1_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("seurat_clusters"))
  print(p1_1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", group.by = c("orig.ident", "seurat_clusters"))
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p2 + p3)
  p7 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"))
  print(p7)
  p7_1 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"))
  print(p7_1)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"))
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p8 + p9)
  dev.off()
  return(combined.data)
}


#' @section SCTransform Integration:
#' * source code: <https://github.com/satijalab/seurat>
#' * quick start: <https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1>
#'
#' @import Seurat
#' @export
#' @method Integration SCTransform
#' @rdname Integration
#'
Integration.SCTransform <- function(object, outdir, project, used, dim){
  object <- Seurat::SCTransform(object)
  projects <- Seurat::SplitObject(object, split.by = "orig.ident")
  features <- Seurat::SelectIntegrationFeatures(object.list = projects, nfeatures = 4000)
  pojects <- Seurat::PrepSCTIntegration(object.list = projects, anchor.features = features)
  anchors <- Seurat::FindIntegrationAnchors(object.list = pojects, normalization.method = "SCT", anchor.features = features)
  combined.data <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT") %>%
    Seurat::RunPCA(npcs = dim)  %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:dim) %>%
    Seurat::FindNeighbors() %>%
    Seurat::FindClusters()
  prefix <- file.path(outdir, sprintf("%s.debatch.seurat-sctransform.dim=%d", project, dim))
  saveRDS(combined.data, paste0(prefix, ".rds"))
  pdf(paste0(prefix, ".pdf"))
  p1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("orig.ident"))
  print(p1)
  p1_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("seurat_clusters"))
  print(p1_1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", group.by = c("orig.ident", "seurat_clusters"))
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p2 + p3)
  p7 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"))
  print(p7)
  p7_1 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"))
  print(p7_1)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"))
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p8 + p9)
  dev.off()
  return(combined.data)
}


#' @section Integration Using Harmony:
#' * source code: <https://github.com/immunogenomics/harmony>
#' * quick start: <https://portals.broadinstitute.org/harmony/articles/quickstart.html>
#' @md
#'
#' @import Seurat
#' @importFrom harmony RunHarmony
#' @export
#' @rdname Integration
#'
Integration.Harmony <- function(object, outdir, project, used, dim){
  all.genes <- rownames(object)
  combined.data <- object %>%
    Seurat::NormalizeData()  %>%
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = dim) %>%
    Seurat::FindNeighbors() %>%
    Seurat::FindClusters()
  combined.data <- harmony::RunHarmony(combined.data, group.by.vars = c("orig.ident"))
  combined.data <- Seurat::RunUMAP(combined.data,
                                   dims = 1:ncol(combined.data[["harmony"]]), reduction = "harmony")
  combined.data <- Seurat::FindNeighbors(combined.data, reduction = "harmony") %>% Seurat::FindClusters()
  prefix <- file.path(outdir, sprintf("%s.debatch.harmony.dim=%d", project, dim))
  saveRDS(combined.data, paste0(prefix, ".rds"))
  pdf(paste0(prefix, ".pdf"))
  p1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("orig.ident"))
  print(p1)
  p1_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("seurat_clusters"))
  print(p1_1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", group.by = c("orig.ident", "seurat_clusters"))
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p2 + p3)
  p4 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "harmony", group.by = c("orig.ident"))
  print(p4)
  p4_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "harmony", group.by = c("seurat_clusters"))
  print(p4_1)
  p5 <- Seurat::VlnPlot(object = combined.data, features = "harmony_1", group.by = c("orig.ident", "seurat_clusters"))
  p6 <- Seurat::VlnPlot(object = combined.data, features = "harmony_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p5 + p6)
  p7 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"))
  print(p7)
  p7_1 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"))
  print(p7_1)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"))
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p8 + p9)
  dev.off()
  return(combined.data)
}


#' @section Integration Using Liger:
#' * source code: <https://github.com/welch-lab/liger>
#' * quick start: <https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html>
#'
#' @import Seurat
#' @importFrom SeuratWrappers RunOptimizeALS RunQuantileNorm
#' @export
#' @method Integration Liger
#' @rdname Integration
#'
Integration.Liger <- function(object, outdir, project, used, dim){
  #projects <- MergeFileData(csv, outdir, mincell, minrna, maxrna, maxmt)
  #combined.data <- merge(projects[[1]], tail(projects, length(projects)-1))
  all.genes <- rownames(object)
  combined.data <- object
  combined.data <- combined.data %>%
    Seurat::NormalizeData()  %>%
    Seurat::FindVariableFeatures(nfeatures = 4000) %>%
    Seurat::ScaleData(split.by = "orig.ident", do.center = FALSE) %>%
    Seurat::RunPCA(npcs = dim) %>%
    SeuratWrappers::RunOptimizeALS(k = dim, split.by = "orig.ident")  %>%
    SeuratWrappers::RunQuantileNorm(split.by = "orig.ident")
  # run pca to display, liger only need scaled data
  combined.data <- Seurat::RunUMAP(combined.data,
                                   dims = 1:ncol(combined.data[["iNMF"]]), reduction = "iNMF")
  combined.data <- Seurat::FindNeighbors(combined.data, reduction = "iNMF") %>% Seurat::FindClusters()
  prefix <- file.path(outdir, sprintf("%s.debatch.liger.dim=%d", project, dim))
  saveRDS(combined.data, paste0(prefix, ".rds"))
  pdf(paste0(prefix, ".pdf"))
  p1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("orig.ident"))
  print(p1)
  p1_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("seurat_clusters"))
  print(p1_1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", group.by = c("orig.ident", "seurat_clusters"))
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p2 + p3)
  p4 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "iNMF", group.by = c("orig.ident"))
  print(p4)
  p4_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "iNMF", group.by = c("seurat_clusters"))
  print(p4_1)
  p5 <- Seurat::VlnPlot(object = combined.data, features = "iNMF_1", group.by = c("orig.ident", "seurat_clusters"))
  p6 <- Seurat::VlnPlot(object = combined.data, features = "iNMF_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p5 + p6)
  p7 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"))
  print(p7)
  p7_1 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"))
  print(p7_1)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"))
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"))
  print(p8 + p9)
  dev.off()
  return(combined.data)
}
