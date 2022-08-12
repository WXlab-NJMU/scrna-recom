#' @import
NULL


#' Seurat Standard Prepocessing
#'
#' @description
#' Preprocessing steps: QualityControl, Normalize, FeatureSelection, Scaling
#'
#' @param indir the 10X genomics data folder,
#'    including matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir output results folder
#' @return Returns a processed seurat object rds
#'
#' @import Seurat
#' @export
#' @rdname preprocess
#' @order 1
#' @concept preprocess
#'
SeuratPreprocess <- function(indir, outdir) {
  rds <- indir %>%
    QualityControl(outdir) %>%
    Normalization(outdir) %>%
    FeatureSelection(outdir) %>%
    Scale(outdir)
}

QualityControl <- function (project, indir, outdir) {
  pbmc.data <- Read10X(data.dir = indir)
  pbmc <- CreateObject(counts = pbmc.data,
                             project = project,
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

#' Seurat SCTransform Prepocessing
#'
#' @description
#' SCTransform replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
#' and removes confounding sources of variation, such as mitochondrial percentage
#'
#' @references
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1
#'
#' @param indir the 10X genomics data folder,
#'    including matrix.mtx, genes.tsv, bardcodes.tsv
#' @param outdir
#'
#' @rdname preprocess
#' @import Seurat
#' @import sctransform
#' @export
#' @concept preprocess
#'
SeuratSCTransform <- function(indir, outdir) {
  pbmc_data <- Read10X(data.dir = indir)
  pbmc <- CreateSeuratObject(counts = pbmc_data)
  pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
  # option of method = "glmGamPoi" could improve the speed
  pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  outrds = file.path(outdir, "sctransformed.rds")
  saveRDS(pbmc, outrds)
  return(outrds)
}
