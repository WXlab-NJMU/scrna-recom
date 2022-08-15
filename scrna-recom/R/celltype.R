#' @section scCATCH:
#' Cluster-based Toolkit for Cellular Heterogeneity (scCATCH) 
#' from cluster potential marker genes identification to cluster annotation 
#' based on evidence-based score by matching the potential marker genes 
#' with known cell markers in tissue-specific cell taxonomy reference database (CellMatch).
#' 
#' Support species is 'Human' or 'Mouse'
#' * source code: <https://github.com/ZJUFanLab/scCATCH>
#' * quick start: <https://raw.githack.com/ZJUFanLab/scCATCH/master/vignettes/tutorial.html>
#' @md
#'
#' @import scCATCH
#' @export
#' @rdname AnnotateCellType
#' @method AnnotateCellType scCATCH
#' @concept cell type annotation
#' 
#' @param species Value: Human or 'Mouse'
#' @param tissue Tissue name
#' @param strict Use the most strict condition to identify marker genes, default is FALSE
#' @param geneinfo NCBI gene info for a given species, see [scCATCH::demo_geneinfo()]
#' @param cellmatch Known markers of human and mouse, other species see [scCATCH::demo_marker()]
#' @param cell_min_pct Include the gene detected in at least this many cells in each cluster, default: 0.25
#' @param logfc	Include the gene with at least this fold change of average gene expression compared to every other clusters, default is 0.25
#' @param pvalue Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters, default is 0.05

AnnotateCellType.scCATCH <- function(input, outdir, 
                                     species = "Human",
                                     tissue = NULL,
                                     strict = FALSE,
                                     geneinfo = scCATCH::geneinfo,
                                     cellmatch = scCATCH::cellmatch,
                                     cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05) {
  # test
  input <- pbmc
  geneinfo = scCATCH::geneinfo
  cellmatch = scCATCH::cellmatch
  species = "Human"
  tissue = "Blood"
  strict = FALSE
  # revise gene symbol for scRNA data
  expr.norm.data <- input[['RNA']]@data
  expr.norm.data <- scCATCH::rev_gene(data = expr.norm.data, species = species, 
                             data_type = "data", geneinfo = geneinfo)
  obj <- createscCATCH(data = input[['RNA']]@data, cluster = as.character(Seurat::Idents(input)))
  # set the condition to identify marker genes, 1 for strict, 2 for loose
  method = ifelse(strict == TRUE, "1", "2")
  obj <- findmarkergene(object = obj, species = species, tissue = tissue, 
                        marker = cellmatch, use_method = method,  
                        cell_min_pct = cell_min_pct, logfc = logfc, pvalue = pvalue)
  obj <- findcelltype(object = obj)
  # write output
  saveRDS(out, file.path(outdir, "annotated_cell_type.scCATCH.rds"))
  write.table(obj@celltype, file.path(outdir, "annotated_cell_type.scCATCH.tsv"))
  #png(file.path(outdir,"annotated_cell_type.scCATCH.pdf"))
  #dev.off()
}

#' @section SingleR:
#' * source code: <https://github.com/dviraran/SingleR>
#' * quick start: <https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html>
#' ## Two mode:
#' * use built-in reference: HumanPrimaryCellAtlasData、BlueprintEncodeData、MouseRNAseqData、
#'   DatabaseImmuneCellExpressionData、NovershternHematopoieticData、MonacoImmuneData，
#'   see details in `browseVignettes("celldex")` 
#' * use pre-labelled dataset to annotate: (todo) 
#'   See `browseVignettes("SingleR::trainSingleR")`
#' ## Paramenters
#' * input: 
#' @md
#' 
#' @import SingleR SingleR plotScoreHeatmap plotDeltaDistribution
#' @import scater plotHeatmap
#' @import celldex HumanPrimaryCellAtlasData
#' @export
#' @rdname AnnotateCellType
#' @method AnnotateCellType SingleR
#' @concept cell type annotation
#' 
AnnotateCellType.SingleR <- function(input, used,
                                     reference = NULL) {
  hpca <- celldex::HumanPrimaryCellAtlasData()
  # scale data: input@assays[["integrated"]]@scale.data
  # use norm data
  out <- SingleR::SingleR(test = input@assays[["RNA"]]@data, 
                          ref = hpca, labels = hpca$label.main)
  table(out$labels)
  saveRDS(out, file.path(outdir, "annotated_cell_type.SingleR.rds"))
  png(file.path(outdir,"annotated_cell_type.SingleR.pdf"))
  p1 <- SingleR::plotScoreHeatmap(out)
  print(p1)
  p2 <- SingleR::plotDeltaDistribution(out, ncol = 3)
  print(p2)
  all.markers <- metadata(out)$de.genes
  p3 <- scater::plotHeatmap(out, order_columns_by="labels", 
                            features=unique(unlist(all.markers$beta)))
  print(p3)
  dev.off()
  return(out)
}
