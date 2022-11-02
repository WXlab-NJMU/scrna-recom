library(dplyr)


#' Clustering on seurat object
#'
#' @import Seurat
#' @import ggplot2
#' @export
#' @param input Input seurat object
#' @param outdir Output folder
#' @param project Project name
#' @param nfeatures Number of variable features to used
#' @param plot.features Features to plot on UMAP, default is c("nFeature_RNA", "percent.mt", "percent.rb")
#' @param reduction Reduction method, current support is pca, harmony, iNMF (liger)
#' @param dim Dimensions to use for clustering
#' @param resolution Resolution to use for clustering
#' @rdname clustering
#'
clustering <- function (input, outdir, project, 
                        nfeatures = 2000, plot.features = c("nFeature_RNA", "percent.mt", "percent.rb"),
                        reduction = "pca", dims = 30, k =20, resolution = 0.8){
  # create outdir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.cluster.fs=%d.ds=%s.dims=%d.reso=%.2f", 
                                      project, nfeatures, reduction, dims, resolution))
  pdf(paste0(prefix, ".pdf"))
  # prepare
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
  assay <- input@active.assay
  if (assay == "RNA") {
    input <- Seurat::NormalizeData(input)  
  }
  ## features
  input <- Seurat::FindVariableFeatures(input, selection.method = "vst", nfeatures = as.numeric(nfeatures))
  top10 <- head(Seurat::VariableFeatures(input), 10)
  p <- Seurat::VariableFeaturePlot(input)
  p <- Seurat::LabelPoints(plot = p, points = top10, 
                           repel = TRUE, xnudge = 0, ynudge = 0) + ggplot2::theme(legend.position="bottom")
  print(p)
  if (reduction == "pca") {
    # pca
    input <- Seurat::ScaleData(input, features= rownames(input))
    input <- Seurat::RunPCA(input, npcs = dims, 
                            features = Seurat::VariableFeatures(object = input))
    p1 <- Seurat::DimPlot(input, reduction = reduction, group.by = c("orig.ident"),
                          shuffle = TRUE, raster = T) 
    print(p1)
    p2 <- Seurat::VizDimLoadings(input, dims = 1:9, reduction = reduction) & 
      ggplot2::theme(axis.text=ggplot2::element_text(size=5), 
                     axis.title=ggplot2::element_text(size=8,face="bold"))
    print(p2)
    p3 <- Seurat::DimHeatmap(input, reduction = reduction, 
                             dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
    print(p3)
    p4 <- Seurat::ElbowPlot(input, reduction = reduction)
    print(p4)
  } else if (reduction == "harmony"){
    input <- Seurat::ScaleData(input)
    if (!("harmony" %in% names(input@reductions))){
      input <- Seurat::RunPCA(input, npcs = dims)
      input <- harmony::RunHarmony(input, group.by.vars = c("orig.ident"))
    }
    p1 <- Seurat::DimPlot(input, reduction = reduction, group.by = c("orig.ident"), 
                          shuffle = TRUE, raster = T) 
    print(p1)
    p2 <- Seurat::VizDimLoadings(input, dims = 1:9, reduction = reduction) & 
      ggplot2::theme(axis.text=ggplot2::element_text(size=5), 
                     axis.title=ggplot2::element_text(size=8,face="bold"))
    print(p2)
    p3 <- Seurat::DimHeatmap(input, reduction = reduction, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
    print(p3)
    p4 <- Seurat::ElbowPlot(input, reduction = reduction)
    print(p4)
  } else if (reduction == "iNMF") {
    ## scale
    input <- Seurat::ScaleData(input, split.by = "orig.ident", do.center = FALSE)
    ## reduct
    if (!("iNMF" %in% names(input@reductions))){
      input <- SeuratWrappers::RunOptimizeALS(input, k = dims, split.by = "orig.ident") %>%
        SeuratWrappers::RunQuantileNorm(split.by = "orig.ident")  
    }
    p1 <- Seurat::DimPlot(input, reduction = reduction, group.by = c("orig.ident"), 
                          shuffle = TRUE, raster = T)     
    print(p1)
  }
  # cluster
  input <- Seurat::FindNeighbors(input, reduction = reduction, k.param = k, dims = 1:dims)
  input <- Seurat::FindClusters(input, resolution = resolution)
  # umap
  input <- Seurat::RunUMAP(input, dims = 1:dims, reduction = reduction)
  write.table(table(input@meta.data$seurat_clusters, input@meta.data$orig.ident), 
              paste0(prefix, ".clusters.tsv"),
              quote = FALSE, row.names = FALSE)
  p5 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"),
                        label.size = 5, repel = T,label = T, raster = T)
  print(p5)
  p6 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"), raster = T)
  print(p6)
  p7 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", split.by = "orig.ident", raster = T) 
  print(p7)
  # plot features
  for (feature in plot.features){
    p <- Seurat::FeaturePlot(input, reduction = "umap", features = feature, raster = T) & ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  # markers
  markers <- Seurat::FindAllMarkers(input, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC) -> top20
  write.csv(top20, paste0(prefix, ".top20_genes.csv"))
  markers %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC) -> top3
  Seurat::DotPlot(input, features = unique(top3$gene[1:20])) & 
    Seurat::NoLegend() & 
    ggplot2::labs(title = "Top3 markers") &
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust=1)) -> p  
  print(p)
  
  dev.off()
  saveRDS(input, paste0(prefix, ".rds"))
  return(input)
}