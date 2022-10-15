library(dplyr)


#' @section Clustering
#'
#' @import Seurat
#' @import ggplot2
#' @export
#' @method Clustering
#' @rdname clustering
#'
clustering <- function (input, outdir, project,
                        nfeatures = 2000, plot.features = c("nFeature_RNA", "percent.mt", "percent.rb"),
                        reduction = "pca", dims = 30, resolution = 0.5){
  # create outdir
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  prefix <- file.path(outdir, sprintf("%s.cluster.fs=%d.ds=%s.dims=%d.reso=%.2f", 
                                      project, nfeatures, reduction, dims, resolution))
  pdf(paste0(prefix, ".pdf"))
  # prepare
  assay <- input@active.assay
  if (! .hasSlot(input@meta.data, "percent.mt")) {
    input[["percent.mt"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "^MT-")  
  }
  if (! .hasSlot(input@meta.data, "percent.hb")) {
    input[["percent.hb"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "^HB[AB]")
  }
  if (! .hasSlot(input@meta.data, "percent.rb")){
    input[["percent.rb"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "^RP[SL]")  
  }
  ## norm
  if (! .hasSlot(input@assays[[assay]], "data")) input <- Seurat::NormalizeData(input)
  ## features
  if (! .hasSlot(input@assays[[assay]], "var.features")) {
    input <- Seurat::FindVariableFeatures(input, selection.method = "vst", nfeatures = nfeatures)
  }
  if (reduction == "pca") {
    ## scale
    if (! .hasSlot(input@assays[[assay]], "scale.data")) input <- Seurat::ScaleData(input)
    # pca
    input <- Seurat::RunPCA(input, npcs = dims)
    p1 <- Seurat::DimPlot(input, reduction = reduction)
    print(p1)
    p2 <- Seurat::VizDimLoadings(input, dims = 1:9, reduction = reduction) & 
      ggplot2::theme(axis.text=ggplot2::element_text(size=5), 
                     axis.title=ggplot2::element_text(size=8,face="bold"))
    print(p2)
    p3 <- Seurat::DimHeatmap(input, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
    print(p3)
    p4 <- Seurat::ElbowPlot(input)
    print(p4)
  } else if (reduction == "harmony"){
    if (! .hasSlot(input@assays[[assay]], "scale.data")) input <- Seurat::ScaleData(input)
    input <- Seurat::RunPCA(input, npcs = dims)
    input <- harmony::RunHarmony(input, group.by.vars = c("orig.ident"))
  } else if (reduction == "iNMF") {
    ## scale
    if (! .hasSlot(input@assays[[assay]], "scale.data")) {
      Seurat::ScaleData(input, split.by = "orig.ident", do.center = FALSE)
    }
    ## reduct
    input <- SeuratWrappers::RunOptimizeALS(input, k = dims, split.by = "orig.ident") %>%
      SeuratWrappers::RunQuantileNorm(split.by = "orig.ident")
  }
  # cluster
  input <- Seurat::FindNeighbors(input, reduction = reduction, dims = 1:dims)
  input <- Seurat::FindClusters(input, resolution = resolution)
  # umap
  input <- Seurat::RunUMAP(input, dims = 1:dims, verbose = F)
  write.table(table(input@meta.data$seurat_clusters), 
              paste0(prefix, ".clusters.tsv"),
              quote = FALSE, row.names = FALSE)
  p5 <- Seurat::DimPlot(input, label.size = 4, repel = T,label = T)
  print(p5)
  # plot features
  for (feature in plot.features){
    p <- Seurat::FeaturePlot(input, features = feature) & ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  dev.off()
  saveRDS(input, paste0(prefix, ".rds"))
}