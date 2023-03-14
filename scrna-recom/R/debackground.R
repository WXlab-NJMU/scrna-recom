#' @include cluster.R
#' @import Seurat dplyr ggplot2 SoupX DropletUtils
NULL


#' Remove background RNA using SoupX
#'
#' @description
#' tutorial: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html
#' output is the count matrix
#'
#' @param raw cellranger output raw_feature_bc_matrix folder
#' @param filtered cellranger output filtered_feature_bc_matrix folder
#' @param project project name
#' @param outdir output folder
#'
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @import SoupX
#' @import DropletUtils
#' @rdname background-removal
#' @export
#'
background.removal <- function(raw, filtered, outdir, project) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.bkremoval.SoupX", project))
  pdf(paste0(prefix, ".pdf"))
  tod <- Seurat::Read10X(raw)
  toc <- Seurat::Read10X(filtered)
  sc = SoupX::SoupChannel(tod, toc)

  # add cluster and reduction
  folder.clustering <- file.path(outdir, "clustering")
  project.before <- paste0(project, ".before")
  Seurat::Read10X(filtered) %>% Seurat::CreateSeuratObject(project = project.before) %>%
    clustering(folder.clustering, project.before) -> obj.before
  Seurat::FindAllMarkers(obj.before, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
    group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC) -> top3.before
  meta <- obj.before@meta.data
  sc = SoupX::setClusters(sc, setNames(meta$seurat_cluster, rownames(meta)))
  sc = SoupX::setDR(sc, as.data.frame(obj.before@reductions$umap[[]]), reductName = "UMAP")
  # perform SoupX
  sc = SoupX::autoEstCont(sc)
  out = SoupX::adjustCounts(sc)
  DropletUtils:::write10xCounts(file.path(outdir, "soupx_filtered_matrix"), out)
  ## difference
  #cntSoggy = rowSums(sc$toc > 0)
  #cntStrained = rowSums(out > 0)
  #mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 20)
  #write.csv(mostZeroed, paste0(prefix, ".mostzeroed.csv", quote=F))

  # Recluster using Seurat
  project.after <- paste0(project, ".after")
  Seurat::CreateSeuratObject(out, project = project.after) %>%
    clustering(folder.clustering, project.after) -> obj.after
  saveRDS(obj.after, paste0(prefix, ".rds"))
  Seurat::FindAllMarkers(obj.after, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
    group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC) -> top3.after
  # plot
  top3genes.before <- unique(top3.before$gene)
  top3genes.after <- unique(top3.after$gene)
  Seurat::DimPlot(obj.before, reduction = "umap", group.by = c("seurat_clusters"),
                  label =T, shuffle = TRUE, raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "before SoupX") -> p1
  Seurat::DimPlot(obj.after, reduction = "umap", group.by = c("seurat_clusters"),
                  label =T, shuffle = TRUE, raster = T) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "after SoupX") -> p2
  print(p1+p2)
  Seurat::DotPlot(obj.before, features = top3genes.before[1:18]) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "before SoupX", caption = "Top3 markers") &
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6, angle = 45, hjust=1),
                   axis.title.x = ggplot2::element_blank()) -> p1
  Seurat::DotPlot(obj.after, features = top3genes.after[1:18]) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "after SoupX", caption = "Top3 markers") &
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6, angle = 45, hjust=1),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) -> p2
  print(p1+p2)
  Seurat::DotPlot(obj.after, features = top3genes.before[1:18]) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "after SoupX", caption = "Top3 before") &
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6, angle = 45, hjust=1),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) -> p1_2
  print(p1 + p1_2)
  Seurat::DotPlot(obj.before, features = top3genes.after[1:18]) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "before SoupX", caption = "Top3 after") &
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6, angle = 45, hjust=1),
                   axis.title.x = ggplot2::element_blank()) -> p2_1
  print(p2_1 + p2)
  dev.off()
  return(obj.after)
}
