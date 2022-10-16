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
#' @param tissue Tissue name in scCATCH database, default Blood is used
#' @param strict Use the most strict condition to identify marker genes, default is FALSE
#' @param geneinfo NCBI gene info for a given species, see [scCATCH::demo_geneinfo()]
#' @param cellmatch Known markers of human and mouse, other species see [scCATCH::demo_marker()]
#' @param cell_min_pct Include the gene detected in at least this many cells in each cluster, default: 0.25
#' @param logfc	Include the gene with at least this fold change of average gene expression compared to every other clusters, default is 0.25
#' @param pvalue Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters, default is 0.05

AnnotateCellType.scCATCH <- function(input, outdir, project, used,
                                     species = "Human",
                                     tissue = NULL,
                                     strict = FALSE,
                                     geneinfo = scCATCH::geneinfo,
                                     cellmatch = scCATCH::cellmatch,
                                     cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05) {
  #geneinfo = scCATCH::geneinfo
  #cellmatch = scCATCH::cellmatch
  # revise gene symbol for scRNA data
  expr.norm.data <- input[['RNA']]@data
  expr.norm.data <- scCATCH::rev_gene(data = expr.norm.data, species = species,
                             data_type = "data", geneinfo = geneinfo)
  obj <- scCATCH::createscCATCH(data = input[['RNA']]@data,
                                cluster = as.character(Seurat::Idents(input)))
  # set the condition to identify marker genes, 1 for strict, 2 for loose
  method = ifelse(strict == TRUE, "1", "2")
  obj <- scCATCH::findmarkergene(object = obj, species = species, tissue = tissue,
                        marker = cellmatch, use_method = method,
                        cell_min_pct = cell_min_pct, logfc = logfc, pvalue = pvalue)
  obj <- scCATCH::findcelltype(object = obj)
  # write output
  clusters <- obj@celltype
  clusters <- sapply(input@active.ident,
                     function(x) ifelse(x %in% clusters$cluster,
                                        clusters[clusters$cluster == x,]$cell_type,
                                        "unknown"))
  meta.cluster <- data.frame(clusters, row.names = Seurat::Cells(input))
  input <- Seurat::AddMetaData(input, metadata = meta.cluster,
                               col.name = "scCATCH.cluster_type")

  prefix <- file.path(outdir, sprintf("%s.celltype.scCATCH.tissue=%s",
                                      project, sub(" ","_",tissue)))
  saveRDS(input, paste0(prefix, ".rds"))
  write.table(obj@celltype, paste0(prefix, ".detail.csv"), sep = ",")
  pdf(paste0(prefix, ".pdf"))
  p5 <- Seurat::DimPlot(input, shuffle = TRUE, label = T, repel = T,
                        reduction = "umap", group.by = c("scCATCH.cluster_type"))
  p5 <- p5 + ggplot2::theme(legend.position="bottom",
                            legend.text = ggplot2::element_text(size=8),
                            legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p5)
  p5_1 <- Seurat::DimPlot(input, shuffle = TRUE, label = T, repel = T, label.size = 3,
                          reduction = "umap", 
                          group.by = c("scCATCH.cluster_type"), split.by = "orig.ident")
  p5_1 <- p5_1 + ggplot2::theme(legend.position="bottom",
                                legend.text = ggplot2::element_text(size=8),
                                legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p5_1)
  p7 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"))
  print(p7)
  dev.off()
}

#' @section SingleR:
#' * source code: <https://github.com/dviraran/SingleR>
#' * quick start: <https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html>
#' ## Two mode:
#' * use built-in reference: HumanPrimaryCellAtlasData、BlueprintEncodeData、MouseRNAseqData、
#'   DatabaseImmuneCellExpressionData、NovershternHematopoieticData、MonacoImmuneData
#'   see details in `browseVignettes("celldex")` or https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_General-purpose_references
#' * use pre-labelled dataset to annotate: (todo)
#'   See `browseVignettes("SingleR::trainSingleR")`
#' ## Paramenters
#' * input:
#' @md
#' @param reference Reference dataset in celldex,
#'   HumanPrimaryCellAtlasData、BlueprintEncodeData、MouseRNAseqData、
#'   DatabaseImmuneCellExpressionData、NovershternHematopoieticData、MonacoImmuneData
#' @param level Label levels: main (broad), fine (fine-grained), ont (standard in Cell Ontology)
#' @importFrom  SingleR SingleR plotScoreHeatmap plotDeltaDistribution
#' @importFrom  scater plotHeatmap
#' @importFrom  celldex HumanPrimaryCellAtlasData
#' @export
#' @rdname AnnotateCellType
#' @method AnnotateCellType SingleR
#' @concept cell type annotation
#'
AnnotateCellType.SingleR <- function(input, outdir, project, used,
                                     reference = NULL, level = NULL) {
  if (reference == "HumanPrimaryCellAtlasData") {
    refdata <- celldex::HumanPrimaryCellAtlasData() #general, 158 entries
  }else if (reference == "BlueprintEncodeData"){
    refdata <- celldex::BlueprintEncodeData() # pure stroma and immune cells, 43 entries
  }else if (reference == "MouseRNAseqData"){
    refdata <- celldex::MouseRNAseqData() # 28 entries
  }else if (reference == "ImmGenData"){
    refdata <- celldex::ImmGenData() # exhaustive coverage, need to remove certain samples, 356 entries
  }else if (reference == "DatabaseImmuneCellExpressionData"){
    refdata <- celldex::DatabaseImmuneCellExpressionData() # CD4+ T cell subsets, 15 entires
  } else if (reference == "NovershternHematopoieticData"){
    refdata <- celldex::NovershternHematopoieticData() # greatest resolution for myeloid and progenitor cells, 38 entries
  } else if (reference == "MonacoImmuneData"){
    refdata <- celldex::MonacoImmuneData() #best covers all of the bases for a typical PBMC sample, 29 entries
  }
  if (level == "main") {
    labels <- refdata$label.main
  } else if (level == "fine") {
    labels <- refdata$label.fine
  } else if (level == "ont") {
    labels <- refdata$label.ont
  }
  # scale data: input@assays[["integrated"]]@scale.data
  # use norm data: input@assays[["integrated"]]@data
  out.cluster <- SingleR::SingleR(test = input@assays[["RNA"]]@data, clusters = input@active.ident,
                          ref = refdata, labels = labels)
  out.barcode <- SingleR::SingleR(test = input@assays[["RNA"]]@data,
                                  ref = refdata, labels = labels)
  table(out.cluster$labels)
  table(out.barcode$labels)

  input <- Seurat::AddMetaData(input, metadata = factor(out.barcode$pruned.labels), col.name = "SingleR.cell_type")
  meta.cluster <- sapply(input@active.ident, function(x) out.cluster$pruned.labels[as.integer(x)+1])
  input <- Seurat::AddMetaData(input, metadata = meta.cluster, col.name = "SingleR.cluster_type")

  prefix <- file.path(outdir, sprintf("%s.celltype.SingleR.dataset=%s", project, reference))
  saveRDS(input, paste0(prefix, ".rds"))
  write.table(out.cluster, paste0(prefix, ".clusters.detail.csv"), sep = ",")
  write.table(out.barcode, paste0(prefix, ".barcodes.detail.csv"), sep = ",")
  pdf(paste0(prefix, ".pdf"))
  p4 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", group.by = c("SingleR.cell_type"))
  p4 <- p4 + ggplot2::theme(legend.position="bottom",
                            legend.text = ggplot2::element_text(size=8),
                            legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p4)
  p4_1 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", 
                        group.by = c("SingleR.cell_type"), split.by = "orig.ident")
  p4_1 <- p4_1 + ggplot2::theme(legend.position="bottom",
                                legend.text = ggplot2::element_text(size=8),
                                legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p4_1)
  p5 <- Seurat::DimPlot(input, shuffle = TRUE, label = T, repel = T, label.size = 3,
                        reduction = "umap", group.by = c("SingleR.cluster_type"))
  p5 <- p5 + ggplot2::theme(legend.position="bottom",
                            legend.text = ggplot2::element_text(size=8),
                            legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p5)
  p5_1 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", 
                          group.by = c("SingleR.cluster_type"), split.by = "orig.ident")
  p5_1 <- p5_1 + ggplot2::theme(legend.position="bottom",
                                legend.text = ggplot2::element_text(size=8),
                                legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p5_1)
  p7 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"))
  p7 <- p7 + ggplot2::theme(legend.position="bottom",
                            legend.text = ggplot2::element_text(size=8),
                            legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p7)

  #p1 <- SingleR::plotScoreHeatmap(out.cluster)
  #print(p1)
  #p2 <- SingleR::plotDeltaDistribution(out.cluster, ncol = 3)
  #print(p2)
  #all.markers <- metadata(out.cluster)$de.genes
  #p3 <- scater::plotHeatmap(out.cluster, order_columns_by="labels",
  #                          features=unique(unlist(all.markers$beta)))
  #print(p3)
  dev.off()
}
