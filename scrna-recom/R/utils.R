library(tidyverse)
#' TryCatch Using WithCallingHandlers
#'
TryCatchWithCallingHandlers <- function(expression){
  withCallingHandlers(
    expression,
    warning = function(w){
      message("warning:\n", w)
      lobstr::cst()
      invokeRestart("muffleWarning")},
    error = function(e){
      lobstr::cst()
      message("error:\n", e) },
    finally = { message("...") })
}

#' Merge Data to List of Seurat Object
#'
#' @description
#' Merge data included in the csvfile to one list
#' @import Seurat
#' @param csv The header must be "project,path"
#' @return Returns List consisted of multiple Seurat objects
#'
#' @concept utility
#' @importFrom Seurat Read10X CreateSeuratObject
#' @export
#'
MergeFileData <- function(csv, outdir,
                          mincell, minrna, maxrna, maxmt,
                          genecol = 2) {
  inputs <- read.csv(csv)
  data.list <- apply(inputs, 1, FUN = function(item) {
    indir  <- item[["path"]]
    sample  <- item[["sample"]]
    data <- Seurat::Read10X(data.dir = indir, gene.column = genecol, strip.suffix = TRUE) %>%
      Seurat::CreateSeuratObject(project = sample)
    data[["dataset"]] <- sample
  })
  return(data.list)
}

#' Read files to Seurat Object
#'
#' @description
#' Merge data included in the csvfile to one list
#'
#' @param csv The header must be "project,path"
#' @return Returns List consisted of multiple Seurat objects
#'
#' @concept utility
#' @importFrom Seurat Read10X CreateSeuratObject
#' @importFrom SeuratDisk Connect Convert LoadH5Seurat
#' @importFrom tools file_ext
#' @export
#'
ReadtoSeuratObject <- function(infile, genecol = 2) {
  message("Reading: ", infile)
  if (file_test('-f', infile)) {
    extension <- tools::file_ext(infile)
    # Convert from SingleCellExperiment
    if (extension == "rds") {
      input <- readRDS(file = infile)
      final <- Seurat::as.Seurat(input, counts = "counts", data = "logcounts")
    } else if (extension == "loom") {
      input <- SeuratDisk::Connect(filename = infile, mode = "r")
      final <- Seurat::as.Seurat(input)
    } else if (extension == "h5ad") {
      SeuratDisk::Convert(infile, dest = "h5seurat", overwrite = TRUE)
      converted.file = sub(".h5ad",'.h5seurat', infile)
      final <- SeuratDisk::LoadH5Seurat(converted.file)
    } else {
      warning("Input is not supported !!!")
    }
  } else if (file_test('-d', infile)) {
    input <- Seurat::Read10X(infile, gene.column = genecol, strip.suffix=TRUE)
    final <- Seurat::CreateSeuratObject(input)
  }
  final
}


#' Read H5ad to Seurat Object
#'
#' @description
#' Merge data included in the csvfile to one list
#'
#' @param infile input is h5ad file
#' @return Returns Seurat objects
#'
#' @concept utility
#' @import Seurat
#' @import SeuratDisk
#' @importFrom reticulate import
#' @export
#'
ConvertH5adToSeurat <- function(infile, outfile){
  scanpy <- reticulate::import("scanpy")
  pandas <- reticulate::import("pandas")

  adata <- scanpy$read(infile)
  meta <- adata$obs
  gene <- adata$var

  adata2 <- t(as.matrix(adata$X))
  genes <- rownames(gene)
  barcodes <- rownames(meta)
  rownames(adata2) <- genes # gene 21606
  colnames(adata2) <- barcodes # barcodes 12262

  seu <- CreateSeuratObject(adata2)
  seu <- AddMetaData(seu, meta)
  seu <- NormalizeData(seu)%>% ScaleData()
  seu$RNA@var.features <- genes[adata$var[["highly_variable"]]]


  pca.embeddings <- adata$obsm[["X_pca"]]
  rownames(pca.embeddings) <- barcodes
  pca.loadings = adata$varm[["PCs"]]
  rownames(pca.loadings) <- genes
  pca.stdev <- as.vector(adata$uns[["pca"]]$variance)

  pca.dr <- CreateDimReducObject(
    embeddings = pca.embeddings,
    loadings = pca.loadings,
    stdev = pca.stdev,
    key = "PC",
    assay = "RNA"
  )
  seu@reductions[["pca"]] <- pca.dr

  umap.embeddings <- adata$obsm[["X_umap"]]
  rownames(umap.embeddings) <- barcodes

  umap.dr <- CreateDimReducObject(
    embeddings = umap.embeddings,
    key = "UMAP",
    assay = "RNA"
  )
  seu@reductions[["umap"]] <- umap.dr
  saveRDS(seu, outfile)
  seu
}

#' Convert Seurat Object to Other Format
#'
#' @description
#' convert a seurat object to other format object
#'
#' @param csv The header must be "project,path"
#' @return Returns List consisted of multiple Seurat objects
#'
#' @concept utility
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SeuratDisk SaveH5Seurat Convert as.loom
#' @export
#'
OutFromSeuratToOthers <- function(obj, format, outfile = NULL) {
  if (format == "SingleCellExperiment") {
    final <- Seurat::as.SingleCellExperiment(obj)
  } else if (format == "loom") {
    message("always remember to close loom files when done: `final$close_all()` ")
    final <- SeuratDisk::as.loom(obj, filename = outfile)
    # always remember to close loom files when done: `final$close_all()`
  } else if (format == "AnnData") {
    outfile <- paste0(outfile, '.h5seurat')
    # convert to h5seurat
    SeuratDisk::SaveH5Seurat(obj, filename = outfile, overwrite=TRUE)
    # convert to h5ad
    SeuratDisk::Convert(outfile, dest = "h5ad", overwrite=TRUE)
    final <- sub('.h5seurat', ".h5ad", outfile)
  }
  final
}

#' Science Colors
#'
#' @description
#' Science Color Platte
#' @importFrom ggsci pal_npg pal_aaas
#' @import dplyr
#' @import purrr
#' @importFrom purrr flatten map
#'
#' @param x length
#' @return Returns vector of fixed length colors
#'
#' @concept utility
scicolors <- function(x, journal="nature") {
  nround  <- ifelse(x %% 8 == 0, x %/% 8, (x %/% 8) + 1)
  if (journal=="science"){
    #mycolors  <- map(rounds, ~ggsci::pal_aaas("default", alpha = .)(8)) %>% flatten()
    mycolors  <- rep(ggsci::pal_aaas("default", alpha = 0.8)(10)[seq(8)], nround)
  } else if  (journal=="nature"){
    mycolors  <- rep(ggsci::pal_npg("nrc", alpha = 0.8)(10)[seq(8)], nround)
  }
  mycolors  <- unlist(mycolors)[1:x]
  return(mycolors)
}

#' Add Umap Tag
#'
#' @description
#' Add umap tag on the topbottom
#' @import tidyr
#' @import Seurat
#' @import ggplot2
#' @import stringr
#' @import dplyr
#' @import tidyr
#'
#' @param p seurat plot
#' @param aspect.ratio default is 1
#' @param strip.font.size default is 12
#' @param legend.position default is right, set none if not needed
#' @param legend.font.size default is 8
#' @param line.length default is 3
#' @param line.color default is black
#' @param line.lab.size default is 3
#' @return Returns plot object
#'
#' @concept utility
#' @export
AddTag <- function(p, aspect.ratio = 1, strip.font.size = 8,
                       legend.position = "right", legend.font.size = 5,
                       line.length=3, line.color = "black", line.lab.size =2) {

  data <- p$data
  x.min <- min(data[1]) - 2
  y.min <- min(data[2]) - 2
  x.mid <- x.min + 2*line.length/3
  y.mid <- y.min + 2*line.length/3
  line.tib <- tibble(x=c(x.min, x.min, x.min, x.min + line.length),
                     y=c(y.min, y.min + line.length, y.min, y.min),
                     group=c(1,1,2,2))
  label.tib <- tibble(x=c(x.mid, x.min - 0.5),
                      y=c(y.min - 0.5, y.mid),
                      angle = c(0, 90),
                      label=str_replace_all(names(data)[1:2], '_', ''))

  p <- p + theme_void() +
           theme( plot.title = element_blank(),
                  strip.text = element_text(colour = "black", face = "bold", size = strip.font.size),
                  legend.position = legend.position, legend.title = element_blank(),
                  legend.text = element_text(size=legend.font.size),
                  legend.key.size=ggplot2::unit(5,'mm'),
                  aspect.ratio = aspect.ratio)
  p <- p + geom_line(data = line.tib, aes(x=x,y=y,group=group),
                     arrow = arrow(length = ggplot2::unit(0.05,"inches"),
                                   ends = "last",type = "closed")) +
           geom_text(data = label.tib,
                     aes(x = x,y = y,angle = angle,label = label),
                     fontface = "italic",color = line.color, size = line.lab.size)
  return(p)
}

#' Plot Cell Type Ratio
#'
#' @description
#' bar plot of cell ratios
#'
#' @param p seurat plot
#' @param input seurat object
#' @return Returns plot object
#'
#' @import ggalluvial
#' @import ggplot2
#' @importFrom ggsci scale_fill_npg
#' @import Seurat
#' @import dplyr
#'
#' @concept utility
#' @export
PlotCellRatio <- function(input, sample = "orig.ident", celltype = "seurat_clusters"){
  ratio.info <- input@meta.data %>% group_by(.data[[sample]],.data[[celltype]]) %>%
    summarise(num = n()) %>% mutate(prop = round(num/sum(num), digits = 4)) %>%
    mutate(id = paste0(.data[[celltype]],": ", prop*100))
  ratio.info[[celltype]] <- factor(ratio.info[[celltype]], levels = levels(input))
  p <- ggplot2::ggplot(ratio.info,
                       ggplot2::aes_string(x = sample,y = "prop",fill = celltype,
                                           alluvium = celltype, stratum = celltype,
                                           label = "id")) +
    ggalluvial::geom_flow(width = .5, curve_type = "quintic", alpha = 0.6) + geom_stratum(width = .5) +
    ggplot2::geom_text(stat = ggalluvial::StatStratum, size = 4, min.y = 0.01) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::scale_fill_manual(values=scicolors(length(unique(ratio.info[[celltype]])))) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = ggplot2::rel(1.2),color = 'black'),
                 axis.title = ggplot2::element_text(size = ggplot2::rel(1.5),color = 'black'),
                 legend.text = ggplot2::element_text(size = ggplot2::rel(1.2),color = 'black'),
                 legend.title = ggplot2::element_text(size = ggplot2::rel(1.5),color = 'black')) +
    ggplot2::xlab('') + ggplot2::ylab('Cell Ratio')
    p
}

#' Plot Cell Marker Expression
#'
#' @description
#' bar plot of cell ratios
#'
#' @param input seurat object
#' @param genes seurat markergenes vector
#' @return Returns plot object
#'
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @export
PlotCellMarkerExpression <- function(input, genes, kw = "celltype"){
  blues <- c("#BDCEEE","#1552BE")
  reds <- c("#FCC1C9","#F80B39")
  p <- DotPlot(input,group.by = 'seurat_clusters',features = genes, assay = "RNA",
               cols = blues)
  ann <- input@meta.data[,c('seurat_clusters',kw)] %>% distinct()
  colnames(ann) <- c("id","celltype")
  p$data <- left_join(p$data,ann)

  p <- p + facet_grid(facets = celltype~., scales = "free_y",
                      space = "free_y", switch = "x") +
        theme(panel.spacing = ggplot2::unit(x = 0.5,  units = "lines"),
              strip.background =element_blank(),
              strip.text.y =element_text(angle = 0)) +
        annotate(geom = 'segment', y = Inf,yend = -Inf,x = Inf, xend = Inf,size=1) +
        rotate_x_text()
  return(p)
}

#' Determine Optimize Dimensionality
#'
#' @description
#' determine optiomal dimensionality
#'
#' @param input elbow plot data
#' @return Returns optimal dims
#'
#' @import tidyr
#' @import dplyr
#' @export
determineOptimalDims <- function(data){
  data2 <- tibble(dims=head(data$dims,-1),
                stdev=head(data$stdev,-1),
                slope= - diff(data$stdev)/diff(data$dims))
  dims <- data2 %>% filter(slope > 0.05) %>%
    arrange(desc(dims)) %>% .$dims
  optimal_dim <- dims[1] + 1
}


#' Determine PCS by stdev
#'
#' @description
#' determine optimal pcs
#'
#' @param obj.seu seurat object
#' @return Returns optimal dims
#'
#' @import tidyr
#' @import dplyr
#' @export
DeterminePCS <- function(obj.seu){
  pct <- obj.seu[["pca"]]@stdev / sum( obj.seu[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)

  co1 <- which(cumu > 90 & pct < 5)[1] # first pc with cumu > 90
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
              decreasing = T)[1] + 1 # last diff < 0.01
  pcs <- min(co1, co2)
  pcs

  plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
  p <- ggplot2::ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
         geom_text() +
         geom_vline(xintercept = 90, color = "red") +
         geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
         theme_bw()
  print(p)
  return(pcs)
}
