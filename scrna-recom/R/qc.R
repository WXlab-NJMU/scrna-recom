#' @include utils.R
#' @include objects.R
#' @include generics.R
NULL

library(dplyr)

#' Quality Control
#'
#' @description
#' from cellranger matrix to cluster, including normalization, scaling,
#' dimension reduction, clustering
#'
#' @import Seurat
#' @rdname qc
#' @export
#'
qc <- function (indir, outdir, project, 
                max.counts, min.counts, max.genes, min.genes,
                min.cells, max.mt, max.hb) {
  pbmc.data <- Seurat::Read10X(data.dir = indir)
  pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data,
                                     project = project,
                                     min.cells = min.cells)
  total.cells = ncol(x = pbmc)
  total.genes = nrow(x = pbmc)
  print(project)
  cat("total cell counts: ", total.cells, "\n")
  cat("total gene counts: ", total.genes, "\n")
  pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc[["percent.hb"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^HB[AB]")
  # filt counts
  pbmc <- subset(pbmc, subset = nCount_RNA > min.counts & nCount_RNA < max.counts)
  count.filt = total.cells - ncol(x = pbmc)
  count.range = sprintf("%i-%i", min.counts, max.counts)
  sprintf("filt cells using counts(%s): %i", count.range, count.filt)
  # filt genes
  pbmc <- subset(pbmc, subset = nFeature_RNA > min.genes & nFeature_RNA < max.genes)
  gene.filt = total.genes - nrow(x = pbmc)
  gene.range = sprintf("%i-%i", min.genes, max.genes)
  sprintf("filt genes using genes(%s): %i", gene.range, gene.filt)
  # filt mt
  pbmc <- subset(pbmc, subset = percent.mt < max.mt)
  mt.filt = total.cells - count.filt - ncol(x = pbmc)
  sprintf("filt cells using mt<%i: %i", max.mt, mt.filt)
  # filt hb
  print(subset(pbmc, subset = percent.hb > 0))
  pbmc <- subset(pbmc, subset = percent.hb < max.hb)
  hb.filt = total.cells - count.filt - mt.filt - ncol(x = pbmc)
  sprintf("filt cells using hb<%i: %i", max.hb, hb.filt)
  # final stat
  cat("final cell counts: ", ncol(pbmc), "\n")
  cat("final gene counts: ", nrow(pbmc), "\n")
  detail <- data.frame(
    item = c("sample", "total.cells", "total.genes", 
             "counts.thresold", "counts.filtcells", 
             "genes.thresold", "filtgenes",
             "mt.thresold", "mt.filtcells",
             "hb.thresold", "hb.filtcells",
             "final.cells", "final.genes"),
    info = c(project, total.cells, total.genes, 
             count.range, count.filt, 
             gene.range, gene.filt, 
             max.mt, mt.filt, 
             max.hb, hb.filt,
             ncol(pbmc), nrow(pbmc))
  )
  #detail <- t(detail)
  write.csv(detail, file.path(outdir, paste0(project,".qc.stat.csv")), 
            row.names = FALSE, col.names = FALSE, quote = FALSE )
  return(pbmc) 
}

#' group quality control
#' 
#' @import Seurat
#' @import dplyr
#' @rdname qc
#' @export
#'
group_qc <- function (csv, outdir, project, 
                      max.counts, min.counts, max.genes, min.genes,
                      min.cells, max.mt, max.hb) {
  # create outdir
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  # read raw data
  inputs <- read.csv(csv)
  samples <- inputs[["sample"]]
  sample.data <- apply(inputs, 1, FUN = function(item) {
    indir  <- item[["path"]]
    sample  <- item[["sample"]]
    data <- Seurat::Read10X(data.dir = indir) %>% Seurat::CreateSeuratObject(project = sample)
    data[["percent.mt"]] <- Seurat::PercentageFeatureSet(data, pattern = "^MT-")
    data[["percent.hb"]] <- Seurat::PercentageFeatureSet(data, pattern = "^HB[AB]")
    data
  })
  raw.data <- merge(sample.data[[1]], tail(sample.data, length(sample.data)-1), 
                    add.cell.ids = samples, project = project)
  
  # plot before qc
  pdf(file.path(outdir, paste0(project,".qc_before.pdf")))
  plot1 <- Seurat::VlnPlot(raw.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 2)
  #plot1 <- Seurat::VlnPlot(raw.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(plot1)
  plot2 <- Seurat::FeatureScatter(raw.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle=TRUE)
  print(plot2)
  plot3 <- Seurat::FeatureScatter(raw.data, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle=TRUE)
  print(plot3)
  plot4 <- Seurat::FeatureScatter(raw.data, feature1 = "nCount_RNA", feature2 = "percent.hb", shuffle=TRUE)
  print(plot4)
  dev.off()
  # qc 
  sample.qcdata <- apply(inputs, 1, FUN = function(item) {
    indir  <- item[["path"]]
    sample  <- item[["sample"]]
    #print(item[["qc"]])
    if (! is.na(item[["qc"]])) {
      text <- unlist(strsplit(item[["qc"]], split = '&'))
      thresolds <- setNames(
        lapply(text, function(x) unlist(strsplit(x, '='))[2]),
        lapply(text, function(x) unlist(strsplit(x, '='))[1])
      )
      print(thresolds)
      max.counts <- ifelse("max.counts" %in% names(thresolds), thresolds$max.counts, max.counts)
      min.counts <- ifelse("min.counts" %in% names(thresolds), thresolds$min.counts, min.counts)
      max.genes <- ifelse("max.genes" %in% names(thresolds), thresolds$max.genes, max.genes)
      min.genes <- ifelse("min.genes" %in% names(thresolds), thresolds$min.genes, min.genes)
      min.cells <- ifelse("min.cells" %in% names(thresolds), thresolds$min.cells, min.cells)
      max.mt <- ifelse("max.mt" %in% names(thresolds), thresolds$max.mt, max.mt)
      max.hb <- ifelse("max.hb" %in% names(thresolds), thresolds$max.hb, max.hb)
    }
    print(sample)
    print(c("max.counts", "min.counts", "max.genes", "min.genes", "min.cells", "max.mt", "max.hb"))
    print(c(max.counts, min.counts, max.genes, min.genes, min.cells, max.mt, max.hb))
    qc(indir, outdir, sample, 
       as.numeric(max.counts), as.numeric(min.counts), 
       as.numeric(max.genes), as.numeric(min.genes), 
       as.numeric(min.cells), as.numeric(max.mt), as.numeric(max.hb))
  })
  qc.data <- merge(sample.qcdata[[1]], tail(sample.qcdata, length(sample.qcdata)-1), 
                    add.cell.ids = samples, project = project)
  # plot after qc
  pdf(file.path(outdir, paste0(project,".qc_after.pdf")))
  plot1 <- Seurat::VlnPlot(qc.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), ncol = 2)
  print(plot1)
  plot2 <- Seurat::FeatureScatter(qc.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle=TRUE)
  print(plot2)
  plot3 <- Seurat::FeatureScatter(qc.data, feature1 = "nCount_RNA", feature2 = "percent.hb", shuffle=TRUE)
  print(plot3)
  plot4 <- Seurat::FeatureScatter(qc.data, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle=TRUE)
  print(plot4)
  dev.off()
  # merge stat
  stat <- do.call(cbind, lapply(samples, function(x) read.csv(file.path(outdir, paste0(x,".qc.stat.csv")))))
  cols <- append(1, seq(1, length(samples))*2 )
  write.csv(stat[cols], file.path(outdir, paste0(project,".qc.stat.csv")), 
            row.names = FALSE, quote = FALSE )
  for (sample in samples){ 
    file.remove(file.path(outdir, paste0(sample,".qc.stat.csv")))
  }
}




