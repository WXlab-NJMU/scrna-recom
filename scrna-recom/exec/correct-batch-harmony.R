library(argparser)
p <- arg_parser("scRNA Batch Correction")
#p <- add_argument(p, "csv", help="input matrix folder", type="character")
p <- add_argument(p, "input", help="input seurat rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--method", help="SeuratCCA, SeuratRPCA, Harmony, Liger",
                  type="character", default = "SeuratCCA")
p <- add_argument(p, "--dim", type="numeric", default=30,
                  help="feature nums, npc in Seurat::RunPCA or k in rliger::optimizeALS, default is 30")
p <- add_argument(p, "--nfeatures", type="numeric", default=2000,
                  help="number of variable features to use for scaledata and pca, default is 2000")
p <- add_argument(p, "--kparam", type="numeric", default=20,
                  help="k.param of knn in FindNeighbor, default is 20")
p <- add_argument(p, "--resolution", type="numeric", default=0.5,
                  help="resolution of cluster, default is 0.5")
p <- add_argument(p, "--plot", nargs='*', 
                  default = c("nFeature_RNA", "percent.mt", "percent.rb"),
                  help="features to plot on umap")
argv <- parse_args(p)
print(argv$plot)

# prepare
object <- readRDS(argv$input)
outdir <- argv$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
project <- argv$project
dim = argv$dim
nfeatures = argv$nfeatures
k = argv$k
resolution = argv$resolution
plot.features = argv$plot

prefix <- file.path(outdir, sprintf("%s.debatch.harmony.dim=%d", project, dim))
pdf(paste0(prefix, ".pdf"))

# before
combined.data <- object %>%
  Seurat::NormalizeData()  %>%
  Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures) %>%
  Seurat::ScaleData() %>%
  Seurat::RunPCA(npcs = dim) %>%
  Seurat::RunUMAP(reduction = "pca", dims = 1:dim) %>%
  Seurat::FindNeighbors(reduction = "pca", dims = 1:dim, k.param = k) %>%
  Seurat::FindClusters(resolution = resolution)
p1 <- Seurat::DimPlot(
  object = combined.data, reduction = "pca", 
  group.by = c("orig.ident", "seurat_clusters"),
  shuffle = TRUE, label = T, repel = T) & Seurat::NoLegend() & ggplot2::labs(title = "before integration")
print(p1)
p1_1 <- Seurat::DimPlot(combined.data, reduction = "umap", ncol = 3, 
                        split.by = "orig.ident") + Seurat::NoLegend() & ggplot2::labs(caption = "before integration")
print(p1_1)
p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", 
                      group.by = c("orig.ident", "seurat_clusters")) & ggplot2::labs(caption = "before integration")
p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2", 
                      group.by = c("orig.ident", "seurat_clusters")) & ggplot2::labs(caption = "before integration")
print(p2 + p3)

# do harmony
combined.data <- harmony::RunHarmony(combined.data, group.by.vars = c("orig.ident")) 
combined.data <- combined.data %>%
  Seurat::RunUMAP(reduction = "harmony", dims = 1:ncol(combined.data[["harmony"]])) %>%
  Seurat::FindNeighbors(reduction = "harmony", dims = 1:dim, k.param = k) %>%
  Seurat::FindClusters(resolution = resolution)
p4 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, 
                      group.by = c("orig.ident"),
                      reduction = "harmony")  & ggplot2::labs(title = "after integration")
print(p4)
p5 <- Seurat::VlnPlot(object = combined.data, features = "harmony_1", 
                      group.by = c("orig.ident", "seurat_clusters")) & ggplot2::labs(caption = "after integration")
p6 <- Seurat::VlnPlot(object = combined.data, features = "harmony_2", 
                      group.by = c("orig.ident", "seurat_clusters")) & ggplot2::labs(caption = "after integration")
print(p5 + p6)
p7 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", 
                      group.by = c("orig.ident")) & ggplot2::labs(title = "after integration")
print(p7)
p7_1 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", 
                        group.by = c("seurat_clusters")) & ggplot2::labs(title = "after integration")
print(p7_1)
p7_2 <- Seurat::DimPlot(combined.data, reduction = "umap", ncol = 3, 
                        split.by = "orig.ident") + Seurat::NoLegend() & ggplot2::labs(caption = "after integration")
print(p7_2)
p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"))
p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"))
print(p8 + p9)
# plot features
if (! .hasSlot(combined.data@meta.data, "percent.mt")) {
  combined.data[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^MT-")  
}
if (! .hasSlot(combined.data@meta.data, "percent.hb")) {
  combined.data[["percent.hb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^HB[AB]")
}
if (! .hasSlot(combined.data@meta.data, "percent.rb")){
  combined.data[["percent.rb"]] <- Seurat::PercentageFeatureSet(combined.data, assay = "RNA", pattern = "^RP[SL]")  
}
for (feature in plot.features){
  p <- Seurat::FeaturePlot(combined.data, features = feature, reduction = "umap") & ggplot2::theme(plot.title = ggplot2::element_text(size=10))
  print(p)
}
dev.off()
saveRDS(combined.data, paste0(prefix, ".rds"))
