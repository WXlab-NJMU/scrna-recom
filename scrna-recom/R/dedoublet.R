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
#' @param anno.used Used annotation, SingleR or scCATCH, default is NULL
#'
remove_doublet <- function (input, outdir, project, 
                            nfeatures = 2000, dims = 30, anno.used = NULL) {
  # create outdir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.dedoublet.dims=%d", 
                                      project, dims))
  pdf(paste0(prefix, ".pdf"))

  # check input seurat object 
  ## NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE
  if (! "pca" %in% names(input@reductions)) {
    if (input@active.assay == "RNA") input <- Seurat::NormalizeData(input)
    ## features
    input <- Seurat::FindVariableFeatures(input, selection.method = "vst", nfeatures = nfeatures)
    input <- Seurat::ScaleData(input)
    input <- Seurat::RunPCA(input, npcs = dims, 
                            features = Seurat::VariableFeatures(object = input))
  }
  if (! "tsne" %in% names(input@reductions)) input <- Seurat::RunTSNE(input, reduction="pca", dims=1:dims)
  
  # pk
  sweep.res <- DoubletFinder::paramSweep_v3(input, PCs = 1:dims, sct = FALSE)
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
  p1 <- Seurat::DimPlot(input, reduction = "umap", shuffle = TRUE, raster = T, label.size = 6,
                  group.by = c(classify.pANN1, classify.pANN2)) & ggplot2::theme(
                    plot.title = element_text(size=10), legend.position="bottom")
  print(p1)
  p2 <- Seurat::DimPlot(input, reduction = "umap", shuffle = TRUE, raster = T,
                  split.by = "orig.ident", group.by = c(classify.pANN2)) & Seurat::NoLegend()  
  print(p2)
  p3 <- Seurat::DimPlot(input, reduction = "tsne", shuffle = TRUE, raster = T, label.size = 6,
                  group.by = c(classify.pANN1, classify.pANN2)) & ggplot2::theme(
                    plot.title = element_text(size=10), legend.position="bottom")
  print(p3)
  p4 <- Seurat::DimPlot(input, reduction = "tsne", shuffle = TRUE, raster = T,
                  split.by = "orig.ident", group.by = c(classify.pANN2)) & Seurat::NoLegend()
  print(p4)
  dev.off()
  saveRDS(input, paste0(prefix, ".before.rds"))
  # keep only singlet
  classify <- input@meta.data[classify.pANN2]
  singlet <- Seurat::Cells(input)[classify == "Singlet"]
  output <- subset(input, cells = singlet)
  p5 <- Seurat::DimPlot(output, reduction = "umap", shuffle = TRUE, raster = T, 
                        label.size = 6, repel = T, label = T, group.by = "seurat_clusters") & 
    ggplot2::labs(title = "After doublet removal") &
    ggplot2::theme(legend.position="bottom", legend.text = ggplot2::element_text(size=8))
  print(p5)
  if (! is.null(anno.used)){
    p5_1 <- Seurat::DimPlot(output, reduction = "umap", shuffle = TRUE, raster = T, 
                          label.size = 6, repel = T, label = T, group.by = paste0(anno.used, ".cluster_type")) & 
      ggplot2::labs(title = "After doublet removal") &
      ggplot2::theme(legend.position="bottom", legend.text = ggplot2::element_text(size=8))
    print(p5_1)   
  }
  p6 <- Seurat::DimPlot(output, reduction = "tsne", shuffle = TRUE, raster = T, 
                        label.size = 6, repel = T, label = T, group.by = "seurat_clusters") & 
    ggplot2::labs(title = "After doublet removal") &
    ggplot2::theme(legend.position="bottom", legend.text = ggplot2::element_text(size=8))
  print(p6)
  if (! is.null(anno.used)){
    p6_1 <- Seurat::DimPlot(output, reduction = "tsne", shuffle = TRUE, raster = T, 
                            label.size = 6, repel = T, label = T, group.by = paste0(anno.used, ".cluster_type")) & 
      ggplot2::labs(title = "After doublet removal") &
      ggplot2::theme(legend.position="bottom", legend.text = ggplot2::element_text(size=8))
    print(p6_1)   
  }
  saveRDS(output, paste0(prefix, ".after.rds"))
  return(output)
}