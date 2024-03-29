#' @include utils.R
NULL

#' Clustering on seurat object
#'
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @importFrom reticulate source_python
#' @importFrom SeuratDisk Convert SaveH5Seurat
#' @importFrom ggsci pal_npg
#' @export
#' @param input Input seurat object
#' @param outdir Output folder
#' @param project Project name
#' @param nfeatures Number of variable features to used
#' @param plot.features Features to plot on UMAP, default is c("nFeature_RNA", "percent.mt", "percent.rb")
#' @param reduction Reduction method, current support is pca, harmony, iNMF (liger)
#' @param dim Dimensions to use for clustering
#' @param resolution Resolution to use for clustering
#' @param do.sct whether to use sctransform instead of normdata and scale data
#' @rdname clustering
#'
clustering <- function (input, outdir, project, dims,
                        nfeatures = 2000, plot.features = c("nFeature_RNA", "percent.mt", "percent.rb"),
                        reduction = "pca", k =20, resolution = 2, do.sct = FALSE){
  # create outdir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.cluster.ds=%s.reso=%.2f",
                                      project, reduction, resolution))
  pdf(paste0(prefix, ".pdf"), 16,9)
  # prepare
  if (! .hasSlot(input@meta.data, "percent.mt")) {
    input[["percent.mt"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "(^MT|:MT)-")
  }
  if (! .hasSlot(input@meta.data, "percent.hb")) {
    input[["percent.hb"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "(^HB|:HB)[AB]")
  }
  if (! .hasSlot(input@meta.data, "percent.rb")){
    input[["percent.rb"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "(^RP|:RP)[SL]")
  }
  sample.cols <- scicolors(length(unique(input@meta.data$orig.ident)))
  assay <- input@active.assay
  # norm, scale, find variables
  if (do.sct) {
    if (assay == "integrated"){
        sample.ls <- lapply(X = SplitObject(input, split.by = "orig.ident"), FUN = SCTransform)
        features <- SelectIntegrationFeatures(object.list = sample.ls, nfeatures = 3000)
        sample.ls <- PrepSCTIntegration(object.list = sample.ls, anchor.features = features)
        anchors <- FindIntegrationAnchors(object.list = sample.ls,
                                          normalization.method = "SCT",
                                          anchor.features = features)
        input <- IntegrateData(anchorset = anchors, normalization.method = "SCT") %>%
            Seurat::RunPCA(npcs = 50)

    } else {
        input <- SCTransform(input, assay = assay, vst.flavor = "v2", verbose = FALSE, method = "glmGamPoi")  %>%
            Seurat::RunPCA(npcs = 50)
    }
  } else {
    ## norm
    if (assay == "RNA") input <- Seurat::NormalizeData(input)
    ## features
    input <- Seurat::FindVariableFeatures(input, selection.method = "vst", nfeatures = as.numeric(nfeatures))
    top10 <- head(Seurat::VariableFeatures(input), 10)
    p <- Seurat::VariableFeaturePlot(input)
    p <- Seurat::LabelPoints(plot = p, points = top10,
                             repel = TRUE, xnudge = 0, ynudge = 0) +
      ggplot2::theme(legend.position="bottom")
    print(p)
    if (reduction == "pca") {
      # pca
      input <- Seurat::ScaleData(input, features= rownames(input))
      input <- Seurat::RunPCA(input, npcs = 50,
                              features = Seurat::VariableFeatures(object = input))
      p1 <- Seurat::DimPlot(input, cols = sample.cols, reduction = reduction, group.by = c("orig.ident"),
                            shuffle = TRUE, raster = T)
      print(p1)
      p2 <- Seurat::VizDimLoadings(input, dims = 1:9, reduction = reduction) &
        ggplot2::theme(axis.text=ggplot2::element_text(size=5),
                       axis.title=ggplot2::element_text(size=8,face="bold"))
      print(p2)
      p3 <- Seurat::DimHeatmap(input, reduction = reduction,
                               dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
      print(p3)
    } else if (reduction == "harmony"){
      input <- Seurat::ScaleData(input)
      if (!("harmony" %in% names(input@reductions))){
        input <- Seurat::RunPCA(input, npcs = 50)
        input <- harmony::RunHarmony(input, group.by.vars = c("orig.ident"))
      }
      p1 <- Seurat::DimPlot(input, cols = sample.cols, reduction = reduction, group.by = c("orig.ident"),
                            shuffle = TRUE, raster = T)
      print(p1)
      p2 <- Seurat::VizDimLoadings(input, dims = 1:9, reduction = reduction) &
        ggplot2::theme(axis.text=ggplot2::element_text(size=5),
                       axis.title=ggplot2::element_text(size=8,face="bold"))
      print(p2)
      p3 <- Seurat::DimHeatmap(input, reduction = reduction, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
      print(p3)
      p4 <- Seurat::ElbowPlot(input, ndims = dims, reduction = reduction)
      print(p4)
    } else if (reduction == "iNMF") {
      ## scale
      input <- Seurat::ScaleData(input, split.by = "orig.ident", do.center = FALSE)
      ## reduct
      if (!("iNMF" %in% names(input@reductions))){
        input <- SeuratWrappers::RunOptimizeALS(input, k = 20, split.by = "orig.ident") %>%
          SeuratWrappers::RunQuantileNorm(split.by = "orig.ident")
      }
      p1 <- Seurat::DimPlot(input, cols = sample.cols, reduction = reduction, group.by = c("orig.ident"),
                            shuffle = TRUE, raster = T)
      print(p1)
    }
  }

  # determine the optimal dims
  p <- Seurat::ElbowPlot(input, ndims = 50, reduction = "pca")
  data <- tibble(dims=head(p$data$dims,-1),
                stdev=head(p$data$stdev,-1),
                slope= - diff(p$data$stdev)/diff(p$data$dims))
  write.csv(data, paste0(prefix, ".elbow.csv"), quote = F)
  opt_dim <- determineOptimalDims(p$data)
  DeterminePCS(input)
  if (dims == "auto"){
    print(paste0("Optimal dimensional: ", opt_dim))
  } else {
    opt_dim <- as.integer(dims)
  }
  p <- p + ggplot2::geom_vline(xintercept = opt_dim, color = "red") +
    ggplot2::geom_text(x=c(opt_dim + 2), y=c(2), label=paste0("dim=",opt_dim))
  print(p)
  input <- input %>%
      Seurat::FindNeighbors(reduction = "pca", dims = 1:opt_dim) %>%
      Seurat::FindClusters(resolution = resolution)
  input <- Seurat::RunUMAP(input, reduction = "pca", dims = 1:opt_dim, min_dist = 0.1)
  cluster_counts <- table(input@meta.data$seurat_clusters)
  cluster_props <- prop.table(cluster_counts)
  cluster_stats <- t(bind_rows(cluster_counts, cluster_props))
  colnames(cluster_stats) <- c("count", "percentage")
  cluster_stats[,"percentage"] <- round(cluster_stats[,"percentage"]*100, digits = 2)
  write.table(cluster_stats,
              paste0(prefix, ".clusters.tsv"),
              quote = FALSE, row.names = TRUE)
# cluster region
  Idents(input) <- 'seurat_clusters'
  chunks <- length(levels(Seurat::Idents(input)))
  if (chunks > 15) { chunks <- 15 }
  p <- map(split(levels(Seurat::Idents(input)), 1:chunks),
           function(x) {Seurat::DimPlot(input, label = T, label.size = 2,
                                        cells.highlight = Seurat::CellsByIdentities(input, idents = x)) +
                        Seurat::NoAxes(legend.position='top', legend.text = element_text(size = 8)) })
  p2 <- patchwork::wrap_plots(plots = p, ncol=5)
  print(p2)

  cluster.cols <- scicolors(length(unique(input@meta.data$seurat_clusters)))
  p5 <- Seurat::DimPlot(input, cols = cluster.cols, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"),
                        label.size = 5, repel = T,label = T, raster = T) %>% AddTag()
  print(p5)
  p6 <- Seurat::DimPlot(input, cols = sample.cols, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"), raster = T) %>% AddTag()
  print(p6)
  p7 <- Seurat::DimPlot(input, cols = cluster.cols, shuffle = TRUE, reduction = "umap", split.by = "orig.ident", raster = T) %>% AddTag()
  print(p7)
  p8 <- PlotCellRatio(input, sample = "orig.ident", celltype = "seurat_clusters")
  print(p8)
  # cell cycle
  if (length(intersect(rownames(input), Seurat::cc.genes$s.genes)) > 0) {
      input <- Seurat::CellCycleScoring(input, assay = "RNA", s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes, set.ident = FALSE)
      p9 <- Seurat::VlnPlot(input, assay = "RNA",features = c("S.Score", "G2M.Score"), group.by = "seurat_clusters", ncol = 1) &
        ggplot2::theme(axis.title.x=element_blank())
      print(p9)
      p10 <- Seurat::VlnPlot(input, assay = "RNA",features = c("PCNA", "TOP2A",  "MKI67"), group.by = "seurat_clusters", ncol = 1 ) +
        ggplot2::labs(caption="Cell Cycle Markers") & ggplot2::theme(axis.title.x=element_blank())
      print(p10)
      p11 <- Seurat::DotPlot(input, assay = "RNA", features = intersect(rownames(input), Seurat::cc.genes$s.genes), group.by = "seurat_clusters") +
        Seurat::RotatedAxis() + Seurat::NoLegend() +
        ggplot2::ggtitle("Phase S")  + ggplot2::xlab("Genes") + ggplot2::ylab("Clusters")
      print(p11)
      p12 <- Seurat::DotPlot(input, assay = "RNA", features = intersect(rownames(input), Seurat::cc.genes$g2m.genes), group.by = "seurat_clusters") +
        Seurat::RotatedAxis() + Seurat::NoLegend() +
        ggplot2::ggtitle("Phase G2M") + ggplot2::xlab("Genes") + ggplot2::ylab("Clusters")
      print(p12)
  }

  # plot features
  for (feature in plot.features){
    p <- Seurat::FeaturePlot(input, reduction = "umap",
                             features = feature, raster = T) &
         Seurat::NoAxes() &
         ggplot2::theme(plot.title = ggplot2::element_text(size=10))
    print(p)
  }
  # markers
  if (do.sct) {
    input <- Seurat::PrepSCTFindMarkers(input)
    markers <- Seurat::FindAllMarkers(input, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  } else {
    markers <- Seurat::FindAllMarkers(input, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  }
  markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC) -> top50
  write.csv(top50, paste0(prefix, ".top50_genes.csv"))
  markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC) -> top3
  #Seurat::DotPlot(input, cols = ggsci::pal_npg("nrc")(1),
  Seurat::DotPlot(input, assay = "RNA",
                  features = unique(top3$gene[1:30])) &
    ggplot2::labs(title = "Top3 markers") &
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust=1)) -> p
  print(p)
  p <- Seurat::DoHeatmap(input,  assay = "RNA", slot = "data", group.by = "seurat_clusters",
                         features = unique(top3$gene[1:50])) + Seurat::NoLegend() +
       ggplot2::scale_fill_gradientn(colors = c("blue", "white", "red"))
  print(p)

  dev.off()
  saveRDS(input, paste0(prefix, ".rds"))
  # other figures using scanpy
  input.h5seurat = paste0(prefix, ".h5Seurat")
  input.h5ad = paste0(prefix, ".h5ad")
  SeuratDisk::SaveH5Seurat(input, filename = input.h5seurat, overwrite=TRUE)
  SeuratDisk::Convert(input.h5seurat, dest = "h5ad", overwrite=TRUE)

  reticulate::source_python(system.file("extdata", "scanpy-marker-plot.py", package = "scrnaRecom"))
  plot_clusters(input.h5ad, file.path(outdir, paste0(project, "-scanpy-plots")))
  return(input)
}


#' Rename cluster and plot cluster markers
#'
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @importFrom reticulate source_python
#' @importFrom SeuratDisk Convert SaveH5Seurat
#' @importFrom ggsci pal_npg
#' @export
#' @param input Input seurat object
#' @param outdir Output folder
#' @param project Project name
#' @param groupfile Group file, "cluster,cellType"
#' @param markerfile Marker file, "cellType,gene"
#' @param key Cell type name, default is "celltype_fine"
#' @rdname clustering
#'
renameClusterPlotMarkers <- function (input, outdir, project,
                                      groupfile, markerfile,
                                      key = "celltype_fine"){
  # create outdir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  prefix <- file.path(outdir, sprintf("%s.cluster.%s", project, key))
  pdf(paste0(prefix, ".pdf"), 16,9)
  # read group
  groups.df <- read.csv(groupfile)
  groups <- groups.df$cellType
  names(groups) <-  groups.df$cluster
  # read marker
  markers.df <- read.csv(markerfile)
  #input <- Seurat::RenameIdents(input, groups)
  #input[[key]] = Idents(input)
  input <- AddMetaData(input, col.name = key,
                       metadata = sapply(input@meta.data$seurat_clusters, function(x) groups[[x]]))
  write.table(input@meta.data, paste0(prefix, ".clusters.tsv"),
              quote = FALSE, row.names = FALSE)
  saveRDS(input, paste0(prefix, ".rds"))
  seu.ls <- Seurat::SplitObject(input, split.by = key)
  for (celltype in names(seu.ls)){
    seu <- seu.ls[[celltype]]
    saveRDS(seu, file.path(outdir, paste0(celltype, ".rds")))
  }
  # change cluster order
  Idents(input) <- key
  order <-  markers.df %>% dplyr::distinct(cellType) %>% .$cellType
  message("cluster order: ", paste(order,  delimiter =""))
  Idents(input) <- factor(Idents(input), levels = order)
  message("seurat ident levels: ", paste(levels(Idents(input)),  delimiter =""))

  sample.cols <- scicolors(length(unique(input@meta.data$orig.ident)))
  cluster.cols <- scicolors(length(unique(groups.df$cellType)))
  p5 <- Seurat::DimPlot(input, cols = cluster.cols, shuffle = TRUE, reduction = "umap",
                        label.size = 5, repel = T,label = T, raster = T) %>% AddTag()
  print(p5)
  p6 <- Seurat::DimPlot(input, cols = sample.cols, group.by = "orig.ident",
                        shuffle = TRUE, reduction = "umap", raster = T) %>% AddTag()
  print(p6)
  p7 <- Seurat::DimPlot(input, cols = cluster.cols, split.by = "orig.ident",
                        shuffle = TRUE, reduction = "umap", raster = T) %>% AddTag()
  print(p7)
  p8 <- PlotCellRatio(input, sample = "orig.ident", celltype = key)
  print(p8)
  p <- map(levels(Seurat::Idents(input)),
           function(x) {Seurat::DimPlot(input, label = T, label.size = 2,
                                        cells.highlight = Seurat::CellsByIdentities(input, idents = x)) +
                        Seurat::NoAxes(legend.position='top', legend.text = element_text(size = 8)) })
  p2 <- patchwork::wrap_plots(plots = p, ncol=4)
  print(p2)

  # marker expression
  markers <- markers.df %>% dplyr::distinct(gene) %>% .$gene
  Seurat::DotPlot(input, assay = "RNA", features = markers) &
    Seurat::NoLegend() &
    ggplot2::labs(title = "Cell Markers") &
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust=1)) -> p
  print(p)
  p <- Seurat::DoHeatmap(input,  assay = "RNA", slot = "data",
                         features = markers) + Seurat::NoLegend() +
       ggplot2::scale_fill_gradientn(colors = c("blue", "white", "red"))
  print(p)
  # cell cycle
  if (length(intersect(rownames(input), Seurat::cc.genes$s.genes)) > 0) {
      p10 <- Seurat::VlnPlot(input, assay = "RNA",features = c("PCNA", "TOP2A",  "MKI67"), ncol = 1 ) +
        ggplot2::labs(caption="Cell Cycle Markers") & ggplot2::theme(axis.title.x=element_blank())
      print(p10)
      p11 <- Seurat::DotPlot(input, assay = "RNA", features = intersect(rownames(input), cc.genes$s.genes)) +
        Seurat::RotatedAxis() + Seurat::NoLegend() +
        ggplot2::ggtitle("Phase S")  + ggplot2::xlab("Genes") + ggplot2::ylab("Clusters")
      print(p11)
      p12 <- Seurat::DotPlot(input, assay = "RNA", features = intersect(rownames(input), cc.genes$g2m.genes)) +
        Seurat::RotatedAxis() + Seurat::NoLegend() +
        ggplot2::ggtitle("Phase G2M") + ggplot2::xlab("Genes") + ggplot2::ylab("Cell Type")
    print(p12)
  }
  # plot features
  plot.features = c("nFeature_RNA", "percent.rb")
  if (! .hasSlot(input@meta.data, "percent.mt")) {
    input[["percent.mt"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "(^MT|:MT)-")
  }
  if (! .hasSlot(input@meta.data, "percent.hb")) {
    input[["percent.hb"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "(^HB|:HB)[AB]")
  }
  if (! .hasSlot(input@meta.data, "percent.rb")){
    input[["percent.rb"]] <- Seurat::PercentageFeatureSet(input, assay = "RNA", pattern = "(^RP|:RP)[SL]")
  }
  p <- Seurat::FeaturePlot(input, reduction = "umap", features = plot.features, raster = T) &
       Seurat::NoAxes() &
       ggplot2::theme(plot.title = ggplot2::element_text(size=10))
  print(p)
  dev.off()

  # other figures using scanpy
  input.h5seurat = paste0(prefix, ".h5Seurat")
  input.h5ad = paste0(prefix, ".h5ad")
  SeuratDisk::SaveH5Seurat(input, filename = input.h5seurat, overwrite=TRUE)
  SeuratDisk::Convert(input.h5seurat, dest = "h5ad", overwrite=TRUE)

  reticulate::source_python(system.file("extdata", "scanpy-marker-plot.py", package = "scrnaRecom"))
  plot_cluster_markers(input.h5ad, markerfile, file.path(outdir, paste0(project, "-scanpy-plots")), cluster = key)
  return(input)
}
