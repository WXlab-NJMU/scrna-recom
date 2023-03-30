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
MergeFileData <- function(csv, outdir, mincell, minrna, maxrna, maxmt) {
  inputs <- read.csv(csv)
  data.list <- apply(inputs, 1, FUN = function(item) {
    indir  <- item[["path"]]
    sample  <- item[["sample"]]
    data <- Seurat::Read10X(data.dir = indir) %>% Seurat::CreateSeuratObject(project = sample)
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
ReadtoSeuratObject <- function(infile) {
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
    input <- Seurat::Read10X(infile, strip.suffix=TRUE)
    final <- Seurat::CreateSeuratObject(input)
  }
  final
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
  alphas  <- ifelse(x %% 8 == 0, x %/% 8, (x %/% 8) + 1)
  rounds <- rev(seq(0.3, 0.3*alphas, 0.2))
  if (journal=="science"){
    mycolors  <- map(rounds, ~ggsci::pal_aaas("default", alpha = .)(8)) %>% flatten()
  } else if  (journal=="nature"){
    mycolors  <- map(rounds, ~ggsci::pal_npg("nrc", alpha = .)(8))  %>% flatten()
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
addTag <- function(p, aspect.ratio = 1, strip.font.size = 8,
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
plotCellRatio <- function(input, sample = "orig.ident", celltype = "seurat_clusters"){
  ratio.info <- input@meta.data %>% group_by(.data[[sample]],.data[[celltype]]) %>%
    summarise(num = n()) %>% mutate(rel_num = num/sum(num))
  p <- ggplot2::ggplot(ratio.info,ggplot2::aes_string(x = sample,y = "rel_num")) +
    ggplot2::geom_col(ggplot2::aes_string(fill = celltype), width =  0.5) +
    ggalluvial::geom_flow(ggplot2::aes_string(stratum = celltype, alluvium = celltype, fill = celltype),
                          width = .5, curve_type = "quintic", alpha = 0.4) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggsci::scale_fill_npg(alpha = 0.7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = ggplot2::rel(1.2),color = 'black'),
                   axis.title = ggplot2::element_text(size = ggplot2::rel(1.5),color = 'black'),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(1.2),color = 'black'),
                   legend.title = ggplot2::element_text(size = ggplot2::rel(1.5),color = 'black')) +
    ggplot2::xlab('') + ggplot2::ylab('Cell percent ratio')

  print(p)
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
#' @importFrom ggsci pal_npg
plotCellMarkerExpression <- function(input, genes, kw = "celltype"){
  p <- DotPlot(input,group.by = 'seurat_clusters',features = genes, assay = "RNA",
               cols = ggsci::pal_npg("nrc", alpha = 0.7)(3))
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
determineOptimalDims <- function(data){
  data2 <- tibble(dims=head(data$dims,-1),
                stdev=head(data$stdev,-1),
                slope= - diff(data$stdev)/diff(data$dims))
  optimal_dim <- data2 %>% filter(slope < 0.02) %>% head(3) %>% tail(1) %>% .$dims
  if (length(optimal_dim) == 0){
    optimal_dim <- 20
  } else{
    optimal_dim <- optimal_dim + 1
  }
  return(optimal_dim)
}
