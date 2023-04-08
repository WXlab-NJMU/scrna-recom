#' @include utils.R
#' @include objects.R
#' @include generics.R
NULL

library(Seurat)
QualityControl <- function (indir, outdir, project,
                            mincell, minrna, maxrna, maxmt,
                            genecol = 2) {
  pbmc.data <- Seurat::Read10X(data.dir = indir,
                               gene.column = genecol, strip.suffix = TRUE)
  pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data,
                             project = project,
                             min.cells = mincell,
                             min.features = minrna)
  pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "(^MT|:MT)-")
  #plot
  pdf(file.path(outdir,"01_qc_before.pdf"))
  plot1 <- Seurat::VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(plot1)
  plot2 <- Seurat::FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot3 <- Seurat::FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot2 + plot3)
  dev.off()
  pbmc <- subset(pbmc,
                 subset = nFeature_RNA > minrna
                 & nFeature_RNA < maxrna
                 & percent.mt < maxmt)
  outrds = file.path(outdir,"01_after_qc.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

Normalization <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  pbmc <- Seurat::NormalizeData(pbmc,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)
  outrds = file.path(outdir, "02_normalization.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

FeatureSelection <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(pbmc), 10)
  pdf(file.path(outdir, "03_top10_features.pdf"))
  plot1 <- Seurat::VariableFeaturePlot(pbmc)
  print(plot1)
  plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot2)
  dev.off()
  outrds = file.path(outdir, "03_feature_selection.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

Scale <- function(input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  all.genes <- rownames(pbmc)
  #pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- Seurat::ScaleData(pbmc) # faster! use 2000 instead of all genes
  outrds = file.path(outdir, "04_scale.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

DimensionReduction <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  # pca
  pbmc <- Seurat::RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pdf(file.path(outdir, "05_pca.pdf"))
  plot1 <- Seurat::VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
  plot2 <- Seurat::DimPlot(pbmc, reduction = "pca")
  plot3 <- Seurat::DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
  print(plot1 + plot2 + plot3)
  dev.off()
  # cluster
  pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
  pbmc <- Seurat::FindClusters(pbmc, resolution = 2)
  # umap
  pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
  pdf(file.path(outdir, "05_umap.pdf"))
  plot4 <- Seurat::DimPlot(pbmc, reduction = "umap")
  print(plot4)
  dev.off()
  # determine dimensionality
  pbmc <- Seurat::JackStraw(pbmc, num.replicate = 100)
  pbmc <- Seurat::ScoreJackStraw(pbmc, dims = 1:20)
  pdf(file.path(outdir, "05_determine_dimensionality.pdf"))
  plot5 <- Seurat::JackStrawPlot(pbmc, dims = 1:15)
  print(plot5)
  dev.off()
  outrds = file.path(outdir,"05_dimension_reduction.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

DifferentExpression <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  pbmc.markers <- Seurat::FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  pbmc.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
  pdf(file.path(outdir,"06_different_expression.heatmap.pdf"))
  pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
  plot1 <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
  print(plot1)
  plot2 <- Seurat::VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
  print(plot2)
  plot3 <- Seurat::FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14",
                                          "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
  print(plot3)
  dev.off()
  outrds = file.path(outdir,"06_different_expression.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

CellType <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  # cell type
  new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
                       "FCGR3A+ Mono", "NK", "DC", "Platelet")
  names(new.cluster.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.cluster.ids)
  pdf(file.path(outdir, "07_celltype.pdf"))
  plot1 <- Seurat::DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  print(plot1)
  dev.off()
  outrds = file.path(outdir, "07_cell_type.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}


#' Basic Analysis Using Seurat
#'
#' @description
#' from cellranger matrix to cluster, including normalization, scaling,
#' dimension reduction, clustering
#'
#' @importFrom dplyr if_else
#' @import Seurat
#' @export
#'
setMethod("basic_analysis",
          signature("Input"),
          function(object) {
            dir.create(object@outdir, showWarnings = FALSE)
            cat("Steps: ", object@steps)
            # qc
            rds <- QualityControl(object@indir, object@outdir,
                                  object@project, object@mincell,
                                  object@minrna, object@maxrna, object@maxmt,
                                  object@genecol
            )
            # norm
            rds <- dplyr::if_else(2 %in% object@steps, Normalization(rds, object@outdir), rds)
            # feature select
            rds <- dplyr::if_else(3 %in% object@steps, FeatureSelection(rds, object@outdir), rds)
            # scale
            rds <- dplyr::if_else(4 %in% object@steps, Scale(rds, object@outdir), rds)
            # dimension reduction
            rds <- dplyr::if_else(5 %in% object@steps, DimensionReduction(rds, object@outdir), rds)
            # remove doublet
            rds <- dplyr::if_else(6 %in% object@steps, RemoveDoublet(rds, object@outdir), rds)
            # different expression
            rds <- dplyr::if_else(7 %in% object@steps, DifferentExpression(rds, object@outdir), rds)
            # cell type
            rds <- dplyr::if_else(8 %in% object@steps, CellType(rds, object@outdir), rds)
          }
)
