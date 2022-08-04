'scRNA Pipeline Using Seurat.

Usage:
  seurat.R qc <indir> <outdir>
              [--maxRNA=<int>] [--minRNA=<int>] [--mt=<float>]
              [--minCell=<int>] [--project=<char>]

  seurat.R norm <input> <outdir>

  seurat.R feat <input> <outdir>

  seurat.R scale <input> <outdir>

  seurat.R reduct <input> <outdir>

  seurat.R diff <input> <outdir>

  seurat.R type <input> <outdir>

  seurat.R workflow <indir> <outdir> [--steps=<step1,step2,...>]

  seurat.R -h | --help
  seurat.R --version

Options:
  -h --help       Show this screen.
  --version       Show version.
  --maxRNA=<int>  nFeature_RNA maximum [default: 2500]
  --minRNA=<int>  nFeature_RNA minimum [default: 200]
  --minCell=<int> Cell minimum [default: 3]
  --mt=<float>    Percent of maximum mt genes [default: 5].
  --project=<char>
                  Project name [default: pbmc].

  --steps=<step1,step2,...>   Workflow:
                qc,norm,diff,reduct,scale,diff,type [default: all]
' -> doc

library(docopt)
arguments <- docopt(doc, version = '1.0.0')
print(arguments)

library(Seurat)
library(dplyr)

QualityControl <- function (indir, outdir) {
  pbmc.data <- Read10X(data.dir = indir)
  pbmc <- CreateSeuratObject(counts = pbmc.data,
                             project = arguments$project,
                             min.cells = as.numeric(arguments$minCell),
                             min.features = as.numeric(arguments$minRNA))
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  #plot
  pdf(file.path(outdir,"01_qc_before_metrics_violin_plot_all_samples.pdf"))
  plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(plot1)
  dev.off()
  pdf(file.path("01_qc_countRNA_vs_mt_and_feature.pdf"))
  plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot2 + plot3)
  dev.off()
  pbmc <- subset(pbmc,
                 subset = nFeature_RNA > as.numeric(arguments$minRNA)
                 & nFeature_RNA < as.numeric(arguments$maxRNA)
                 & percent.mt < as.numeric(arguments$mt))
  outrds = file.path(outdir,"01_after_qc.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

Normalization <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  pbmc <- NormalizeData(pbmc,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)
  outrds = file.path(outdir, "02_normalization.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

FeatureSelection <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(pbmc), 10)
  pdf(file.path(outdir, "03_top10_features.pdf"))
  plot1 <- VariableFeaturePlot(pbmc)
  print(plot1)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
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
  pbmc <- ScaleData(pbmc) # faster! use 2000 instead of all genes
  outrds = file.path(outdir, "04_scale.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

DimensionReduction <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  # pca
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pdf(file.path(outdir, "05_pca.pdf"))
  plot1 <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
  plot2 <- DimPlot(pbmc, reduction = "pca")
  plot3 <- DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
  print(plot1 + plot2 + plot3)
  dev.off()
  # cluster
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  # umap
  pbmc <- RunUMAP(pbmc, dims = 1:10)
  pdf(file.path(outdir, "05_umap.pdf"))
  plot4 <- DimPlot(pbmc, reduction = "umap")
  print(plot4)
  dev.off()
  # determine dimensionality
  pbmc <- JackStraw(pbmc, num.replicate = 100)
  pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
  pdf(file.path(outdir, "05_determine_dimensionality.pdf"))
  plot5 <- JackStrawPlot(pbmc, dims = 1:15)
  print(plot5)
  dev.off()
  outrds = file.path(outdir,"05_dimension_reduction.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

DifferentExpression <- function (input, outdir) {
  pbmc <- readRDS(normalizePath(input))
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  pbmc.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
  pdf(file.path(outdir,"06_different_expression.heatmap.pdf"))
  pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
  plot1 <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
  print(plot1)
  plot2 <- VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
  print(plot2)
  plot3 <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14",
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
  plot1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  print(plot1)
  dev.off()
  outrds = file.path(outdir, "07_cell_type.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}

steps = list()
if (arguments$workflow) {
  if (arguments$steps == "all"){
    steps = c("qc", "norm", "feat", "scale", "reduct", "diff", "type")
  } else {
    steps = unlist(strsplit(arguments$steps, split=','))
  }
}
print(steps)
dir.create(arguments$outdir, showWarnings = FALSE)
rds <- ''
if (arguments$qc | 'qc' %in% steps) {
  rds <- arguments$indir %>% QualityControl(arguments$outdir)
}
if (arguments$norm | 'norm' %in% steps) {
  rds <- rds %>% Normalization(arguments$outdir)
}
if (arguments$feat | 'feat' %in% steps) {
  rds <- rds %>% FeatureSelection(arguments$outdir)
}
if (arguments$scale | 'scale' %in% steps) {
  rds <- rds %>% Scale(arguments$outdir)
}
if (arguments$reduct | 'reduct' %in% steps) {
  rds <- rds %>% DimensionReduction(arguments$outdir)
}
if (arguments$diff | 'diff' %in% steps) {
  rds <- rds %>% DifferentExpression(arguments$outdir)
}
if (arguments$type | 'type' %in% steps) {
  rds <- rds %>% CellType(arguments$outdir)
}

