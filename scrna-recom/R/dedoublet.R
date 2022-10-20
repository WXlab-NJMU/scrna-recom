library(dplyr)
#' Remove doublet on seurat object
#'
#' @import Seurat
#' @import DoubletFinder
#' @import ggplot2
#' @export
#' @param input Input seurat object
#' @param outdir Output folder
#' @param project Project name
#' @param nfeatures Number of variable features to used
#' @param dims Dimensions to use for clustering
#' @param cores Cores to compute 
#'
remove_doublet <- function (input, outdir, project, 
                            nfeatures = 2000, dims = 50, cores = 10) {
  # create outdir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.dedoublet.dims=%d", project, dims))

  # check input seurat object 
  ## NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE
  ## only support assay=RNA
  if (! "pca" %in% names(input@reductions)) {
    if (input@active.assay == "RNA") input <- Seurat::NormalizeData(input)
    ## features
    input <- Seurat::FindVariableFeatures(input, selection.method = "vst", nfeatures = nfeatures) 
    input <- Seurat::ScaleData(input) %>%
      Seurat::RunPCA(npcs = dims, features = Seurat::VariableFeatures(object = input)) %>%
      Seurat::FindNeighbors(reduction = "pca", dims = 1:dims) %>%
      Seurat::FindClusters()
  }
  if (! "umap" %in% names(input@reductions))  input <- Seurat::RunUMAP(input, reduction = "pca", dims = 1:dims)
  if (! "tsne" %in% names(input@reductions)) input <- Seurat::RunTSNE(input, reduction="pca", dims=1:dims)
  
  # pk
  sweep.res <- DoubletFinder::paramSweep_v3(input, PCs = 1:dims, sct = FALSE, num.cores = cores)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  # nExp
  annotations <- input@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)      ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(input@meta.data))   ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  # final
  pANN1 <-  paste(0.25, mpK, nExp_poi, sep="_")
  pANN2 <-  paste(0.25, mpK, nExp_poi.adj, sep="_")
  input <- DoubletFinder::doubletFinder_v3(input, PCs = 1:dims, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  input <- DoubletFinder::doubletFinder_v3(input, PCs = 1:dims, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = paste("pANN", 0.25, mpK, nExp_poi, sep="_"), sct = FALSE)
  classify.pANN1 <- paste0("DF.classifications_", pANN1)
  classify.pANN2 <- paste0("DF.classifications_", pANN2)
  pdf(paste0(prefix, ".pdf"))
  p1 <- Seurat::DimPlot(input, reduction = "pca", shuffle = TRUE, raster = T, label.size = 6,
                        group.by = c(classify.pANN1, classify.pANN2)) & 
    ggplot2::labs(caption = "Before doublet removal") &
    ggplot2::theme(plot.title = element_text(size=10), legend.position="bottom")
  print(p1)
  p2 <- Seurat::DimPlot(input, reduction = "umap", shuffle = TRUE, raster = T, label.size = 6,
                        group.by = c(classify.pANN1, classify.pANN2)) & 
    ggplot2::labs(caption = "Before doublet removal") &
    ggplot2::theme(plot.title = element_text(size=10), legend.position="bottom")
  print(p2)
  p3 <- Seurat::DimPlot(input, reduction = "tsne", shuffle = TRUE, raster = T, label.size = 6,
                  group.by = c(classify.pANN1, classify.pANN2)) &
    ggplot2::labs(caption = "Before doublet removal") &
    ggplot2::theme(plot.title = element_text(size=10), legend.position="bottom")
  print(p3)
  #saveRDS(input, paste0(prefix, ".before.rds"))
  # keep only singlet
  classify <- input@meta.data[classify.pANN2]
  singlet <- Seurat::Cells(input)[classify == "Singlet"]
  output <- subset(input, cells = singlet)
  p4 <- Seurat::DimPlot(output, reduction = "pca", shuffle = TRUE, raster = T) & 
    Seurat::NoLegend() &
    ggplot2::labs(title = project, caption = "After doublet removal")
  print(p4)
  p5 <- Seurat::DimPlot(output, reduction = "umap", shuffle = TRUE, raster = T) & 
    Seurat::NoLegend() &
    ggplot2::labs(title = project, caption = "After doublet removal")
  p6 <- Seurat::DimPlot(output, reduction = "tsne", shuffle = TRUE, raster = T) & 
    Seurat::NoLegend() &
    ggplot2::labs(title = project, caption = "After doublet removal") 
  print(p5 + p6)
  dev.off()
  saveRDS(output, paste0(prefix, ".after.rds"))
  return(output)
}


#' Remove doublet on grouped seurat object
#'
#' @import Seurat
#' @import DoubletFinder
#' @import ggplot2
#' @export
#' @param input Input seurat object
#' @param outdir Output folder
#' @param project Project name
#' @param nfeatures Number of variable features to used
#' @param dims Dimensions to use for clustering
#'
group_remove_doublet <- function(input, outdir, project,
                                 nfeatures = 2000, dims = 50, cores = 10) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.dedoublet", project))
  obj.list <- Seurat::SplitObject(input, split.by = "orig.ident")
  samples <- names(obj.list)
  obj.dedoublet <- lapply(seq(obj.list), FUN = function(i){
    sample <- names(obj.list)[i]
    obj <- obj.list[[sample]]
    remove_doublet(obj, outdir, sample, nfeatures = nfeatures, dims = dims, cores = cores)
  })
  output <- merge(obj.dedoublet[[1]], tail(obj.dedoublet, length(obj.dedoublet)-1), add.cell.ids = samples, project = project)
  saveRDS(output, paste0(prefix, ".rds"))
  return(output)
}
