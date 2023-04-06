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
Integration.SeuratCCA <- function(object, outdir, project, used, dims,
                                  nfeatures = 2000, k = 20, resolution = 0.5,
                                  plot.features = c("nFeature_RNA", "percent.mt", "percent.rb")){
  prefix <- file.path(outdir, sprintf("%s.debatch.seurat-cca", project))
  pdf(paste0(prefix, ".pdf"))
  # do cca
  object <- Seurat::NormalizeData(object) %>%
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures)
  projects <- Seurat::SplitObject(object, split.by = "orig.ident")
  anchors <- Seurat::SelectIntegrationFeatures(object.list = projects)
  combined.data <- projects %>%
    Seurat::FindIntegrationAnchors(anchor.features = anchors) %>%
    Seurat::IntegrateData() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = 50)
  # determine the optimal dims
  p <- Seurat::ElbowPlot(combined.data, ndims = 50, reduction = "pca")
  data <- tibble(dims=head(p$data$dims,-1),
                stdev=head(p$data$stdev,-1),
                slope= - diff(p$data$stdev)/diff(p$data$dims))
  print(data)
  write.csv(data, paste0(prefix, ".elbow.csv"), quote = F)
  if (dims == "auto"){
    opt_dim <- determineOptimalDims(p$data)
    print(paste0("Optimal dimensional: ", opt_dim))
  } else {
    opt_dim <- as.integer(dims)
  }
  p <- p + ggplot2::geom_vline(xintercept = opt_dim, color = "red") +
    ggplot2::geom_text(x=c(opt_dim + 2), y=c(2), label=paste0("dim=",opt_dim))
  print(p)
  combined.data <- combined.data %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:opt_dim, min_dist = 0.1) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:opt_dim, k.param = k) %>%
    Seurat::FindClusters(resolution = resolution)
  p1 <- Seurat::DimPlot(object = combined.data, reduction = "pca", group.by = c("orig.ident", "seurat_clusters"),shuffle = TRUE, raster = T, label = T, repel = T) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "after integration")
  print(p1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2", group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  print(p2 + p3)
  mycolors <- scicolors(length(unique(combined.data@meta.data$orig.ident)))
  p7 <- (Seurat::DimPlot(combined.data, cols = mycolors, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"), raster=T) %>% AddTag()) +
    ggplot2::labs(caption = "after integration")
  print(p7)
  mycolors <- scicolors(length(unique(combined.data@meta.data$seurat_clusters)))
  p7_2 <- (Seurat::DimPlot(combined.data, cols = mycolors, shuffle = TRUE, reduction = "umap", split.by = c("orig.ident"), raster=T) %>% AddTag()) +
    ggplot2::labs(title = "after integration")
  print(p7_2)
  p7_3 <- PlotCellRatio(combined.data, "orig.ident", "seurat_clusters")
  print(p7_3)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"), raster = T)
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"), raster = T)
  print(p8 + p9)
  # plot features
  if (! .hasSlot(combined.data@meta.data, "percent.mt")) {
    combined.data[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^MT-")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.hb")) {
    combined.data[["percent.hb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^HB[AB]")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.rb")){
    combined.data[["percent.rb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^RP[SL]")
  }
  for (feature in plot.features){
    p <- Seurat::FeaturePlot(combined.data, features = feature, reduction = "umap", raster = T) &
      ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  dev.off()
  saveRDS(combined.data, paste0(prefix, ".rds"))
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
Integration.SeuratRPCA <- function(object, outdir, project, used, dims, ref.samples = c(0, 1),
                                   nfeatures = 2000, k = 20, resolution = 0.5,
                                   plot.features = c("nFeature_RNA", "percent.mt", "percent.rb")){
  prefix <- file.path(outdir, sprintf("%s.debatch.seurat-rpca", project))
  pdf(paste0(prefix, ".pdf"))

  # do rpca
  object <- object %>% Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures)
  projects <- Seurat::SplitObject(object, split.by = "orig.ident")
  anchors <- Seurat::SelectIntegrationFeatures(object.list = projects)
  # difference with SeruatCCA, run PCA first then
  projects <- lapply(projects, FUN = function(x) {
    x <- Seurat::ScaleData(x, features = anchors)
    x <- Seurat::RunPCA(x, npcs = 50, features = anchors)
  })
  anchors <- Seurat::FindIntegrationAnchors(projects, dims = 1:50, reduction = "rpca", reference = ref.samples)
  combined.data <- IntegrateData(anchorset = anchors, dims = 1:50) %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = 50)
  # determine the optimal dims
  p <- Seurat::ElbowPlot(combined.data, ndims = 50, reduction = "pca")
  data <- tibble(dims=head(p$data$dims,-1),
                stdev=head(p$data$stdev,-1),
                slope= - diff(p$data$stdev)/diff(p$data$dims))
  print(data)
  write.csv(data, paste0(prefix, ".elbow.csv"), quote = F)
  if (dims == "auto"){
    opt_dim <- determineOptimalDims(p$data)
    print(paste0("Optimal dimensional: ", opt_dim))
  } else {
    opt_dim <- as.integer(dims)
  }
  p <- p + ggplot2::geom_vline(xintercept = opt_dim, color = "red") +
    ggplot2::geom_text(x=c(opt_dim + 2), y=c(2), label=paste0("dim=",opt_dim))
  print(p)
  combined.data <- combined.data %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:opt_dim, min_dist = 0.1) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:opt_dim, k.param = k) %>%
    Seurat::FindClusters(resolution = resolution)
  DefaultAssay(combined.data) <- "integrated"
  p1 <- Seurat::DimPlot(object = combined.data, reduction = "pca", group.by = c("orig.ident", "seurat_clusters"),
                        shuffle = TRUE, label = T, repel = T, raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "after integration")
  print(p1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  print(p2 + p3)
  mycolors <- scicolors(length(unique(combined.data@meta.data$orig.ident)))
  p7 <- (Seurat::DimPlot(combined.data, cols = mycolors, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"), raster=T) %>% AddTag()) +
    ggplot2::labs(caption = "after integration")
  print(p7)
  mycolors <- scicolors(length(unique(combined.data@meta.data$seurat_clusters)))
  p7_2 <- (Seurat::DimPlot(combined.data, cols = mycolors, shuffle = TRUE, reduction = "umap", split.by = c("orig.ident"), raster=T) %>% AddTag()) +
    ggplot2::labs(title = "after integration")
  print(p7_2)
  p7_3 <- PlotCellRatio(combined.data, "orig.ident", "seurat_clusters")
  print(p7_3)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T)
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T)
  print(p8 + p9)
  # plot features
  if (! .hasSlot(combined.data@meta.data, "percent.mt")) {
    combined.data[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^MT-")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.hb")) {
    combined.data[["percent.hb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^HB[AB]")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.rb")){
    combined.data[["percent.rb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^RP[SL]")
  }
  for (feature in plot.features){
    p <- Seurat::FeaturePlot(combined.data, features = feature, reduction = "umap", raster = T) &
      ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  dev.off()
  saveRDS(combined.data, paste0(prefix, ".rds"))
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
Integration.SCTransform <- function(object, outdir, project, used, dims,
                                    nfeatures = 2000,k = 20, resolution = 0.5,
                                    plot.features = c("nFeature_RNA", "percent.mt", "percent.rb")){
  prefix <- file.path(outdir, sprintf("%s.debatch.seurat-sctransform", project, dim))
  pdf(paste0(prefix, ".pdf"))

  object <- Seurat::SCTransform(object)
  projects <- Seurat::SplitObject(object, split.by = "orig.ident")
  features <- Seurat::SelectIntegrationFeatures(object.list = projects, nfeatures = nfeatures)
  pojects <- Seurat::PrepSCTIntegration(object.list = projects, anchor.features = features)
  anchors <- Seurat::FindIntegrationAnchors(object.list = pojects, normalization.method = "SCT", anchor.features = features)
  combined.data <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT") %>%
    Seurat::RunPCA(npcs = 50)
  # determine the optimal dims
  p <- Seurat::ElbowPlot(combined.data, ndims = 50, reduction = "pca")
  data <- tibble(dims=head(p$data$dims,-1),
                stdev=head(p$data$stdev,-1),
                slope= - diff(p$data$stdev)/diff(p$data$dims))
  print(data)
  write.csv(data, paste0(prefix, ".elbow.csv"), quote = F)
  if (dims == "auto"){
    opt_dim <- determineOptimalDims(p$data)
    print(paste0("Optimal dimensional: ", opt_dim))
  } else {
    opt_dim <- as.integer(dims)
  }
  p <- p + ggplot2::geom_vline(xintercept = opt_dim, color = "red") +
    ggplot2::geom_text(x=c(opt_dim + 2), y=c(2), label=paste0("dim=",opt_dim))
  print(p)
  combined.data <- combined.data %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:opt_dim, min_dist = 0.1) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:opt_dim, k.param = k) %>%
    Seurat::FindClusters(resolution = resolution)
  p1 <- Seurat::DimPlot(object = combined.data, reduction = "pca", group.by = c("orig.ident", "seurat_clusters"),
                        shuffle = TRUE, label = T, repel = T, raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "after integration")
  print(p1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  print(p2 + p3)
  mycolors <- scicolors(length(unique(combined.data@meta.data$orig.ident)))
  p7 <- (Seurat::DimPlot(combined.data, cols = mycolors, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"), raster=T) %>% AddTag()) +
    ggplot2::labs(caption = "after integration")
  print(p7)
  mycolors <- scicolors(length(unique(combined.data@meta.data$seurat_clusters)))
  p7_2 <- (Seurat::DimPlot(combined.data, cols = mycolors, shuffle = TRUE, reduction = "umap", split.by = c("orig.ident"), raster=T) %>% AddTag()) +
    ggplot2::labs(title = "after integration")
  print(p7_2)
  p7_3 <- PlotCellRatio(combined.data, "orig.ident", "seurat_clusters")
  print(p7_3)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"), raster = T)
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"), raster = T)
  print(p8 + p9)
  # plot features
  if (! .hasSlot(combined.data@meta.data, "percent.mt")) {
    combined.data[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^MT-")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.hb")) {
    combined.data[["percent.hb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^HB[AB]")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.rb")){
    combined.data[["percent.rb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^RP[SL]")
  }
  for (feature in plot.features){
    p <- Seurat::FeaturePlot(combined.data, features = feature, reduction = "umap") & ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  dev.off()
  saveRDS(combined.data, paste0(prefix, ".rds"))
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
Integration.Harmony <- function(object, outdir, project, used, dims,
                                nfeatures = 2000, k = 20, resolution = 0.5,
                                plot.features = c("nFeature_RNA", "percent.mt", "percent.rb")){
  prefix <- file.path(outdir, sprintf("%s.debatch.harmony", project))
  pdf(paste0(prefix, ".pdf"))

  # before
  combined.data <- object %>%
    Seurat::NormalizeData()  %>%
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures) %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = 50)
  # determine the optimal dims
  p <- Seurat::ElbowPlot(combined.data, ndims = 50, reduction = "pca")
  data <- tibble(dims=head(p$data$dims,-1),
                stdev=head(p$data$stdev,-1),
                slope= - diff(p$data$stdev)/diff(p$data$dims))
  print(data)
  write.csv(data, paste0(prefix, ".elbow.csv"), quote = F)
  if (dims == "auto"){
    opt_dim <- determineOptimalDims(p$data)
    print(paste0("Optimal dimensional: ", opt_dim))
  } else {
    opt_dim <- as.integer(dims)
  }
  p <- p + ggplot2::geom_vline(xintercept = opt_dim, color = "red") +
    ggplot2::geom_text(x=c(opt_dim + 2), y=c(2), label=paste0("dim=",opt_dim))
  print(p)
  combined.data <- combined.data %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:opt_dim, min_dist = 0.1) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:opt_dim, k.param = k) %>%
    Seurat::FindClusters(resolution = resolution)
  p1 <- (Seurat::DimPlot(object = combined.data, reduction = "pca",
                        group.by = c("orig.ident", "seurat_clusters"),
                        shuffle = TRUE, label = T, repel = T, raster = T) %>%
        AddTag()) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "before integration")
  print(p1)
  p1_1 <- Seurat::DimPlot(combined.data, reduction = "umap", ncol = 3,
                          split.by = "orig.ident", raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(caption = "before integration")
  print(p1_1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "before integration")
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "before integration")
  print(p2 + p3)

  # do harmony
  combined.data <- harmony::RunHarmony(combined.data, group.by.vars = c("orig.ident"))
  combined.data <- combined.data %>%
    Seurat::RunUMAP(reduction = "harmony", dims = 1:ncol(combined.data[["harmony"]])) %>%
    Seurat::FindNeighbors(reduction = "harmony", dims = 1:opt_dim, k.param = k) %>%
    Seurat::FindClusters(resolution = resolution)
  p4 <- Seurat::DimPlot(object = combined.data, reduction = "harmony", group.by = c("orig.ident"),
                        shuffle = TRUE, raster = T)  &
    ggplot2::labs(title = "after integration")
  print(p4)
  p5 <- Seurat::VlnPlot(object = combined.data, features = "harmony_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  p6 <- Seurat::VlnPlot(object = combined.data, features = "harmony_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  print(p5 + p6)
  mycolors <- scicolors(length(unique(combined.data@meta.data$orig.ident)))
  p7 <- (Seurat::DimPlot(combined.data, cols = mycolors, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"), raster=T) %>% AddTag()) +
    ggplot2::labs(caption = "after integration")
  print(p7)
  mycolors <- scicolors(length(unique(combined.data@meta.data$seurat_clusters)))
  p7_2 <- (Seurat::DimPlot(combined.data, cols = mycolors, shuffle = TRUE, reduction = "umap", split.by = c("orig.ident"), raster=T) %>% AddTag()) +
    ggplot2::labs(title = "after integration")
  print(p7_2)
  p7_3 <- PlotCellRatio(combined.data, "orig.ident", "seurat_clusters")
  print(p7_3)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T)
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T)
  print(p8 + p9)
  # plot features
  if (! .hasSlot(combined.data@meta.data, "percent.mt")) {
    combined.data[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^MT-")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.hb")) {
    combined.data[["percent.hb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^HB[AB]")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.rb")){
    combined.data[["percent.rb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^RP[SL]")
  }
  for (feature in plot.features){
    p <- (Seurat::FeaturePlot(combined.data, features = feature, reduction = "umap") %>% AddTag()) &
      ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  dev.off()
  saveRDS(combined.data, paste0(prefix, ".rds"))
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
Integration.Liger <- function(object, outdir, project, used, dims,
                              nfeatures = 2000, k = 20, resolution = 0.5,
                              plot.features = c("nFeature_RNA", "percent.mt", "percent.rb")){
  prefix <- file.path(outdir, sprintf("%s.debatch.liger", project))
  pdf(paste0(prefix, ".pdf"))

  combined.data <- object %>%
    Seurat::NormalizeData()  %>%
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures) %>%
    Seurat::ScaleData(split.by = "orig.ident", do.center = FALSE) %>%
    Seurat::RunPCA(npcs = 50)
  # determine the optimal dims
  p <- Seurat::ElbowPlot(combined.data, ndims = 50, reduction = "pca")
  data <- tibble(dims=head(p$data$dims,-1),
                stdev=head(p$data$stdev,-1),
                slope= - diff(p$data$stdev)/diff(p$data$dims))
  print(data)
  write.csv(data, paste0(prefix, ".elbow.csv"), quote = F)
  if (dims == "auto"){
    opt_dim <- determineOptimalDims(p$data)
    print(paste0("Optimal dimensional: ", opt_dim))
  } else {
    opt_dim <- as.integer(dims)
  }
  p <- p + ggplot2::geom_vline(xintercept = opt_dim, color = "red") +
    ggplot2::geom_text(x=c(opt_dim + 2), y=c(2), label=paste0("dim=",opt_dim))
  print(p)
  combined.data <- combined.data %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:opt_dim, min_dist = 0.1) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:opt_dim, k.param = k) %>%
    Seurat::FindClusters(resolution = resolution)
  p1 <- Seurat::DimPlot(object = combined.data, reduction = "pca",
                        group.by = c("orig.ident", "seurat_clusters"),
                        shuffle = TRUE, label = T, repel = T) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "before integration")
  print(p1)
  p1_1 <- Seurat::DimPlot(combined.data, reduction = "umap", ncol = 3,
                          split.by = "orig.ident", raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(caption = "before integration")
  print(p1_1)
  p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "before integration")
  p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "before integration")
  print(p2 + p3)

  # do liger
  combined.data <- combined.data %>%
    SeuratWrappers::RunOptimizeALS(k = opt_dim, split.by = "orig.ident")  %>%
    SeuratWrappers::RunQuantileNorm(split.by = "orig.ident")
  combined.data <- combined.data %>%
    Seurat::RunUMAP(dims = 1:ncol(combined.data[["iNMF"]]), reduction = "iNMF") %>%
    Seurat::FindNeighbors(reduction = "iNMF", k.param = k) %>%
    Seurat::FindClusters(resolution = resolution)
  p4 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "iNMF",
                        group.by = c("orig.ident"), raster = T) &
    ggplot2::labs(title = "after integration")
  print(p4)
  p5 <- Seurat::VlnPlot(object = combined.data, features = "iNMF_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  p6 <- Seurat::VlnPlot(object = combined.data, features = "iNMF_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T) &
    ggplot2::labs(caption = "after integration")
  print(p5 + p6)
  p7 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap",
                        group.by = c("orig.ident"), raster = T) &
    ggplot2::labs(title = "after integration")
  print(p7)
  p7_1 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap",
                          group.by = c("seurat_clusters"), raster = T) &
    ggplot2::labs(title = "after integration")
  print(p7_1)
  p7_2 <- Seurat::DimPlot(combined.data, reduction = "umap", ncol = 3,
                          split.by = "orig.ident", raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(caption = "after integration")
  print(p7_2)
  p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T)
  p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2",
                        group.by = c("orig.ident", "seurat_clusters"), raster = T)
  print(p8 + p9)
  # plot features
  if (! .hasSlot(combined.data@meta.data, "percent.mt")) {
    combined.data[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^MT-")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.hb")) {
    combined.data[["percent.hb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^HB[AB]")
  }
  if (! .hasSlot(combined.data@meta.data, "percent.rb")){
    combined.data[["percent.rb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^RP[SL]")
  }
  for (feature in plot.features){
    p <- Seurat::FeaturePlot(combined.data, features = feature, reduction = "umap") & ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  dev.off()
  saveRDS(combined.data, paste0(prefix, ".rds"))
  return(combined.data)
}
