#' @include utils.R
NULL

library(dplyr)
library(ggplot2)

#' Remove doublet on seurat object
#'
#' @import Seurat
#' @import DoubletFinder
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom purrr detect_index
#' @importFrom ggsci scale_fill_npg
#' @export
#' @param input Input seurat object
#' @param outdir Output folder
#' @param project Project name
#' @param nfeatures Number of variable features to used
#' @param dims Dimensions to use for clustering
#' @param cores Cores to compute
#'
remove_doublet <- function (input, outdir, project, dims,
                            nfeatures = 2000, resolution = 1,
                            cores = 10) {
  # create outdir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.dedoublet", project))
  pdf(paste0(prefix, ".pdf"))

  # check input seurat object
  ## NormalizeData, FindVariableGenes, ScaleData, RunPCA, and RunTSNE
  ## only support assay=RNA
  #if (! "pca" %in% names(input@reductions)) {
  #  if (input@active.assay == "RNA") input <- Seurat::NormalizeData(input)
  #  ## features
  #  input <- Seurat::FindVariableFeatures(input, selection.method = "vst", nfeatures = nfeatures)
  #  input <- Seurat::ScaleData(input) %>%
  #    Seurat::RunPCA(npcs = 50, features = Seurat::VariableFeatures(object = input))
  #}
  input <- Seurat::NormalizeData(input) %>%
    Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures) %>%
    Seurat::ScaleData()
  input <- Seurat::RunPCA(input, npcs = 50,
                          features = Seurat::VariableFeatures(object = input))
  # determine the optimal dims
  p <- Seurat::ElbowPlot(input, ndims = 50, reduction = "pca")
  data <- tibble(dims=head(p$data$dims,-1),
                stdev=head(p$data$stdev,-1),
                slope= - diff(p$data$stdev)/diff(p$data$dims))
  print(data)
  write.csv(data, paste0(prefix, ".elbow.csv"), quote = F)
  if (dims == "auto"){
    opt_dim <- determineOptimalDims(p$data)
    DeterminePCS(input)
    print(paste0("Optimal dimensional: ", opt_dim))
  } else {
    opt_dim <- as.integer(dims)
  }
  p <- p + ggplot2::geom_vline(xintercept = opt_dim, color = "red") +
    ggplot2::geom_text(x=c(opt_dim + 2), y=c(2), label=paste0("dim=",opt_dim))
  print(p)
  saveRDS(input, paste0(prefix, ".before.rds"))
  input <- input %>%
      Seurat::FindNeighbors(reduction = "pca", dims = 1:opt_dim) %>%
      Seurat::FindClusters(resolution = resolution)
  input <- Seurat::RunUMAP(input, reduction = "pca", dims = 1:opt_dim, min_dist = 0.1)
  input <- Seurat::RunTSNE(input, reduction="pca", dims=1:opt_dim, check_duplicates = FALSE)
  #if (! "umap" %in% names(input@reductions))  input <- Seurat::RunUMAP(input, reduction = "pca", dims = 1:dims)
  #if (! "tsne" %in% names(input@reductions)) input <- Seurat::RunTSNE(input, reduction="pca", dims=1:dims)

  # pk
  sweep.res <- DoubletFinder::paramSweep_v3(input, PCs = 1:opt_dim, sct = FALSE, num.cores = cores)
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
  input <- DoubletFinder::doubletFinder_v3(input, PCs = 1:opt_dim, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  input <- DoubletFinder::doubletFinder_v3(input, PCs = 1:opt_dim, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = paste("pANN", 0.25, mpK, nExp_poi, sep="_"), sct = FALSE)
  classify.pANN1 <- paste0("DF.classifications_", pANN1)
  classify.pANN2 <- paste0("DF.classifications_", pANN2)
  p1 <- Seurat::DimPlot(input, reduction = "umap", shuffle = TRUE, raster = T, label.size = 6,
                        group.by = c(classify.pANN1, classify.pANN2)) &
    ggplot2::labs(caption = "Before doublet removal") &
    ggplot2::theme(plot.title = ggplot2::element_text(size=10), legend.position="bottom")
  print(p1)
  p2 <- Seurat::DimPlot(input, reduction = "tsne", shuffle = TRUE, raster = T, label.size = 6,
                  group.by = c(classify.pANN1, classify.pANN2)) &
    ggplot2::labs(caption = "Before doublet removal") &
    ggplot2::theme(plot.title = ggplot2::element_text(size=10), legend.position="bottom")
  print(p2)
  #saveRDS(input, paste0(prefix, ".before.rds"))
  # keep only singlet
  count <- t(table(input@meta.data[c(classify.pANN2, "seurat_clusters")]))
  cluster <- as.numeric(attributes(count)$dimnames$seurat_clusters)
  freq <- cbind(cluster, count, prop.table(count) * 100)
  colnames(freq) <- c("cluster","doublet", "singlet", "percent.doublet","percent.singlet")
  count.doublets <- sum(freq[,2])
  percent.doublets <- sum(freq[,2])/(sum(freq[,2]) + sum(freq[,3]))
  p3 <- ggplot(as.data.frame(freq)) +
    geom_bar(aes(x=factor(cluster), y=percent.doublet),stat="identity") +
    ggsci::scale_fill_npg(alpha = 0.7) +
    labs(x="cluster",
         title = "Doublets distribution among clusters",
         subtitle = sprintf("Total doublets: %d(%.2f%%)", count.doublets, percent.doublets * 100))
  print(p3)
  plot.tib <- tibble(FetchData(input, "seurat_clusters"), FetchData(input, classify.pANN2))
  names(plot.tib) <- c("Cluster", "Status")
  plot.tib <- plot.tib %>% group_by(Cluster, Status) %>% summarise(n=n()) %>%
    pivot_wider(names_from = Status, values_from = n) %>% replace(is.na(.), 0)%>%
    mutate(Total = Singlet + Doublet)  %>%
    mutate(Singlet=Singlet/Total)  %>%
    mutate(Doublet=1 - Singlet) %>% select(Cluster, Singlet, Doublet) %>%
    pivot_longer(c(Singlet, Doublet),names_to = "Status", values_to = "Percent") %>%
    rowwise()  %>% mutate(LabelPos = if_else(Status == "Singlet", Percent/2, 1 - Percent/2))
  p4 <- ggplot(data=plot.tib, aes(fill=Status, y = Percent, x = Cluster)) +
    geom_bar(position="fill", stat="identity") +
    ggsci::scale_fill_npg(alpha = 0.7) +
    geom_text(data = plot.tib %>% filter(Status == "Doublet"),
              aes(label = sprintf("%.0f",Percent*100),y=LabelPos),size = 3) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Doublets among clusters",
         subtitle = sprintf("Total doublets: %d(%.2f%%)", count.doublets, percent.doublets * 100))
  print(p4)
  freq <- rbind(freq, Total=c(max(freq[,1]) + 1,
                              sum(freq[,2]),
                              sum(freq[,3]),
                              percent.doublets, 1-percent.doublets))
  freq <- round(freq,4)
  write.csv(freq, paste0(prefix, ".stat.csv"), quote = FALSE)
  classify <- input@meta.data[classify.pANN2]
  singlet <- Seurat::Cells(input)[classify == "Singlet"]
  output <- subset(input, cells = singlet)
  p5 <- Seurat::DimPlot(output, label = TRUE, reduction = "umap", shuffle = TRUE, raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(title = project, caption = "After doublet removal")
  p6 <- Seurat::DimPlot(output, label = TRUE, reduction = "tsne", shuffle = TRUE, raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(title = project, caption = "After doublet removal")
  print(p5)
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
group_remove_doublet <- function(input, outdir, project, dims,
                                 nfeatures = 2000, resolution = 1,
                                 cores = 10) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.dedoublet", project))
  obj.list <- Seurat::SplitObject(input, split.by = "orig.ident")
  samples <- names(obj.list)
  obj.dedoublet <- lapply(seq(obj.list), FUN = function(i){
    sample <- names(obj.list)[i]
    obj <- obj.list[[sample]]
    remove_doublet(obj, outdir, sample,
                   nfeatures = nfeatures, dims = dims, resolution = resolution,
                   cores = cores)
  })
  output <- merge(obj.dedoublet[[1]], tail(obj.dedoublet, length(obj.dedoublet)-1), add.cell.ids = samples, project = project)
  saveRDS(output, paste0(prefix, ".rds"))
  stat <- do.call("rbind",
                  lapply(samples, function (sample){
                    file <- file.path(outdir, sprintf("%s.dedoublet.stat.csv", sample))
                    total <- tail(read.csv(file), 1)
                    total$X <- NULL
                    total
                  }))
  row.names(stat) <- samples
  write.csv(stat, paste0(prefix, ".stat.csv"), quote = F)
  return(output)
}
