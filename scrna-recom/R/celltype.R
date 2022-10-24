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
                                     cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05,
                                     plot.features = c("nFeature_RNA", "percent.mt", "percent.rb")) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.celltype.scCATCH.tissue=%s", project, sub(" ","_",tissue)))
  # Human or Mouse
  species <- stringr::str_to_title(species)
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
  write.table(clusters, paste0(prefix, ".clusters.detail.tsv"))
  meta.clusters <- sapply(input@active.ident,
                          function(x) ifelse(x %in% clusters$cluster,
                                            clusters[clusters$cluster == x,]$cell_type,
                                            "unknown"))
  meta.cluster <- data.frame(meta.clusters, row.names = Seurat::Cells(input))
  input <- Seurat::AddMetaData(input, metadata = meta.cluster,
                               col.name = "scCATCH.cluster_type")
  write.table(obj@celltype, paste0(prefix, ".detail.csv"), sep = ",")
  pdf(paste0(prefix, ".pdf"))
  p5 <- Seurat::DimPlot(input, shuffle = TRUE, raster = T, label = T, repel = T, label.size = 4,
                        reduction = "umap", group.by = c("scCATCH.cluster_type"))
  p5 <- p5 + ggplot2::theme(legend.position="bottom",
                            legend.text = ggplot2::element_text(size=8),
                            legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p5)
  p5_1 <- Seurat::DimPlot(input, shuffle = TRUE, raster = T,
                          reduction = "umap",
                          group.by = c("scCATCH.cluster_type"), split.by = "orig.ident")
  p5_1 <- p5_1 + ggplot2::theme(legend.position="bottom",
                                legend.text = ggplot2::element_text(size=8),
                                legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p5_1)
  p7 <- Seurat::DimPlot(input, shuffle = TRUE, raster = T, label = T, repel = T, label.size = 3,
                        reduction = "umap", group.by = c("seurat_clusters", "orig.ident")) & Seurat::NoLegend()
  print(p7)
  # plot 9 features in one
  max.id <- length(plot.features)
  for (i in seq(1, ceiling(max.id/9))){
    start.id <- 9 * (i-1) + 1
    end.id <- ifelse(i*9 > max.id, max.id, i*9)
    features <- plot.features[start.id:end.id]
    p <- Seurat::FeaturePlot(input, raster = T, features = features, ncol = 3,
                             reduction = "umap") & Seurat::NoLegend() & ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  dev.off()
  saveRDS(input, paste0(prefix, ".rds"))
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
                                     reference = NULL, level = NULL, species = "Human",
                                     plot.features = c("nFeature_RNA", "percent.mt", "percent.rb")) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.celltype.SingleR.dataset=%s", project, reference))
  # Human or Mouse
  species <- stringr::str_to_title(species)
  # human
  hpc <- celldex::HumanPrimaryCellAtlasData() #general, 158 entries
  bpe <- celldex::BlueprintEncodeData() # pure stroma and immune cells, 43 entries
  ice <- celldex::DatabaseImmuneCellExpressionData() # CD4+ T cell subsets, 15 entires
  nh <- celldex::NovershternHematopoieticData() # greatest resolution for myeloid and progenitor cells, 38 entries
  mi <- celldex::MonacoImmuneData() #best covers all of a typical PBMC sample, 29 entries
  # mouse
  mr <- celldex::MouseRNAseqData() # 28 entries
  mig <- elldex::ImmGenData() # mouse, exhaustive coverage, 356 entries
  
  if (reference == "HumanPrimaryCellAtlasData") {
    refdata <- hpc
  } else if (reference == "BlueprintEncodeData"){
    refdata <- bpe
  } else if (reference == "DatabaseImmuneCellExpressionData"){
    refdata <- ice
  } else if (reference == "NovershternHematopoieticData"){
    refdata <- nh
  } else if (reference == "MonacoImmuneData"){
    refdata <- mi
  } else if (reference == "MouseRNAseqData"){
    refdata <- mr
  } else if (reference == "ImmGenData"){
    refdata <- mig
  } else if (reference == "combined"){
    if (species == "Human"){
      refdata <- list(HPC=hpc, BPE=bpe, ICE=ice, NH=nh, MI=mi)
    } else if (species == "Mouse"){
      refdata <- list(MR=mr, MIG=mig)
    }
  }
  
  if (level == "main") {
    labels <- if (reference == "combined") lapply(refdata, function(x) x$label.main) else refdata$label.main
  } else if (level == "fine") {
    labels <- if (reference == "combined") lapply(refdata, function(x) x$label.fine) else refdata$label.fine
  } else if (level == "ont") {
    labels <- if (reference == "combined") lapply(refdata, function(x) x$label.ont) else refdata$label.ont
  }

  # scale data: input@assays[["integrated"]]@scale.data
  # use norm data: input@assays[["integrated"]]@data
  out.cluster <- SingleR::SingleR(test = input@assays[["RNA"]]@data, clusters = input@active.ident,
                          ref = refdata, labels = labels)
  out.barcode <- SingleR::SingleR(test = input@assays[["RNA"]]@data,
                                  ref = refdata, labels = labels)
  table(out.cluster$labels)
  table(out.barcode$labels)

  input <- Seurat::AddMetaData(input, metadata = out.barcode$pruned.labels, col.name = "SingleR.cell_type")
  meta.cluster <- sapply(input@active.ident, function(x) out.cluster$pruned.labels[as.integer(x)])
  input <- Seurat::AddMetaData(input, metadata = meta.cluster, col.name = "SingleR.cluster_type")
  saveRDS(input, paste0(prefix, ".rds"))
  write.table(out.cluster, paste0(prefix, ".clusters.detail.tsv"), row.names = F)
  write.table(out.barcode, paste0(prefix, ".barcodes.detail.tsv"), row.names = F)
  pdf(paste0(prefix, ".pdf"))
  p4 <- Seurat::DimPlot(input, shuffle = TRUE, raster = T,
                        reduction = "umap", group.by = c("SingleR.cell_type"))
  p4 <- p4 + ggplot2::theme(legend.position="bottom",
                            legend.text = ggplot2::element_text(size=8),
                            legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p4)
  p4_1 <- Seurat::DimPlot(input, shuffle = TRUE, raster = T, reduction = "umap",
                        group.by = c("SingleR.cell_type"), split.by = "orig.ident")
  p4_1 <- p4_1 + ggplot2::theme(legend.position="bottom",
                                legend.text = ggplot2::element_text(size=8),
                                legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p4_1)
  p5 <- Seurat::DimPlot(input, shuffle = TRUE, raster = T, label = T, repel = T, label.size = 3,
                        reduction = "umap", group.by = c("SingleR.cluster_type"))
  p5 <- p5 + ggplot2::theme(legend.position="bottom",
                            legend.text = ggplot2::element_text(size=8),
                            legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p5)
  p5_1 <- Seurat::DimPlot(input, shuffle = TRUE, raster = T, reduction = "umap",
                          group.by = c("SingleR.cluster_type"), split.by = "orig.ident")
  p5_1 <- p5_1 + ggplot2::theme(legend.position="bottom",
                                legend.text = ggplot2::element_text(size=8),
                                legend.key.size = ggplot2::unit(0.1, 'cm'))
  print(p5_1)
  p7 <- Seurat::DimPlot(input, shuffle = TRUE, raster = T, label = T, repel = T, label.size = 3,
                        reduction = "umap", group.by = c("seurat_clusters", "orig.ident")) & Seurat::NoLegend()
  print(p7)
  # plot 9 features in one
  max.id <- length(plot.features)
  for (i in seq(1, ceiling(max.id/9))){
    start.id <- 9 * (i-1) + 1
    end.id <- ifelse(i*9 > max.id, max.id, i*9)
    features <- plot.features[start.id:end.id]
    p <- Seurat::FeaturePlot(input, features = features, ncol = 3, raster = T,
                             reduction = "umap") & Seurat::NoLegend() & ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  #for (feature in plot.features){
  #  p <- Seurat::FeaturePlot(input, features = feature, reduction = "umap") & ggplot2::theme(plot.title = ggplot2::element_text(size=10))
  #  print(p)
  #}


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

#' @section CellMarker:
#' * homepage: <http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.html>
#' @md
#' @param tissue tissue name in cellmarker database
#' @import dplyr
#' @import ggplot2
#' @import Seurat
#' @export
#' @rdname AnnotateCellType
#' @method AnnotateCellType CellMarker
#' @concept cell type annotation
#'
AnnotateCellType.CellMarker <- function(input, outdir, project, used, 
                                        species = "Human",  tissue = NULL) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.celltype.CellMarker.tissue=%s", project, tissue))
  pdf(paste0(prefix, ".pdf"))
  # Human or Mouse
  species <- stringr::str_to_title(species)
  markers <- cellmarkers %>% filter(species == species) %>% filter(tissue_class == tissue) 
    
  cells <- markers$cell_name
  for (cell in cells){
    genes <- unlist(strsplit(markers[markers["cell_name"] == cell,]$marker_genes, split="\\|"))
    genes <- intersect(genes, rownames(input))
    cat(cell, ":", genes, "\n")
    if (length(genes) > 1){
      input[[cell]] <- sqrt(colMeans(as.matrix(input@assays$RNA@data)[genes,]))  
    } else if(length(genes) == 1){
      input[[cell]] <- as.matrix(input@assays$RNA@data)[genes,] 
    } else {
      cells <- cells[cells != cell]
    }
  }
  p1 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"),
                        label.size = 5, repel = T,label = T, raster = T)
  print(p1)
  p2 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", split.by = "orig.ident", raster = T) 
  print(p2)
  # plot 9 different cells in one
  max.id <- length(cells)
  for (i in seq(1, ceiling(max.id/9))){
    start.id <- 9 * (i-1) + 1
    end.id <- ifelse(i*9 > max.id, max.id, i*9)
    features <- cells[start.id:end.id]
    p <- Seurat::FeaturePlot(input, features = features, ncol = 3, raster = T,reduction = "umap") & 
      Seurat::NoLegend() & 
      ggplot2::theme(plot.title = ggplot2::element_text(size=10)) 
    print(p)
  }
  dev.off()
  #saveRDS(input, paste0(prefix, ".rds"))
  return(input)
}

#' @section SelfMarker:
#' * homepage: <http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.html>
#' @md
#' @param marker.file csv file, including cell_name,marker_genes
#' @import dplyr
#' @import ggplot2
#' @import Seurat
#' @export
#' @rdname AnnotateCellType
#' @method AnnotateCellType SelfMarker
#' @concept cell type annotation
#'
AnnotateCellType.SelfMarker <- function(input, outdir, project, used, 
                                        marker.file = NULL) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.celltype.SelfMarker", project))
  pdf(paste0(prefix, ".pdf"))
  selfmarkers <- readxl::read_excel(marker.file)
  cells <- selfmarkers$cell_name
  for (cell in cells){
    genes <- unlist(strsplit(selfmarkers[selfmarkers["cell_name"] == cell,]$marker_genes, split=","))
    genes <- intersect(genes, rownames(input))
    cat(cell, ":", genes, "\n")
    if (length(genes) > 1){
      input[[cell]] <- sqrt(colMeans(as.matrix(input@assays$RNA@data)[genes,]))  
    } else if(length(genes) == 1){
      input[[cell]] <- as.matrix(input@assays$RNA@data)[genes,] 
    } else {
      cells <- cells[cells != cell]
    }
  }
  p1 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"),
                        label.size = 5, repel = T,label = T, raster = T)
  print(p1)
  p2 <- Seurat::DimPlot(input, shuffle = TRUE, reduction = "umap", split.by = "orig.ident", raster = T) 
  print(p2)
  # plot 9 different cells in one
  max.id <- length(cells)
  for (i in seq(1, ceiling(max.id/9))){
    start.id <- 9 * (i-1) + 1
    end.id <- ifelse(i*9 > max.id, max.id, i*9)
    features <- cells[start.id:end.id]
    p <- Seurat::FeaturePlot(input, features = features, ncol = 3, raster = T,reduction = "umap") & 
      Seurat::NoLegend() & 
      ggplot2::theme(plot.title = ggplot2::element_text(size=10)) 
    print(p)
  }
  dev.off()
  #saveRDS(input, paste0(prefix, ".rds"))
  return(input)
}
