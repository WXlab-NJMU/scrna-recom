library(argparser)
library(Seurat)
library(dplyr)
library(harmony)
p <- arg_parser("scRNA Batch Correction")
#p <- add_argument(p, "csv", help="input matrix folder", type="character")
p <- add_argument(p, "input", help="input seurat rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--method", help="SeuratCCA, SeuratRPCA, Harmony, Liger",
                  type="character", default = "SeuratCCA")
argv <- parse_args(p)
#print(argv)

object <- readRDS(argv$input)
outdir <- argv$outdir
project <- argv$project

combined.data <- object %>%
  Seurat::NormalizeData()  %>%
  Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  Seurat::ScaleData() %>%
  Seurat::RunPCA(npcs = 50) %>%
  Seurat::FindNeighbors() %>%
  Seurat::FindClusters()
combined.data <- harmony::RunHarmony(combined.data, group.by.vars = c("seurat_clusters", "orig.ident"))
combined.data <- Seurat::RunUMAP(combined.data,
                                 dims = 1:ncol(combined.data[["harmony"]]), reduction = "harmony")
combined.data <- Seurat::FindNeighbors(combined.data, reduction = "harmony") %>% Seurat::FindClusters()
saveRDS(combined.data, file.path(outdir, paste0(project, ".harmony.rds")))
pdf(file.path(outdir, paste0(project, ".harmony.pdf")))
p1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("orig.ident"))
print(p1)
p1_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "pca", group.by = c("seurat_clusters"))
print(p1_1)
p2 <- Seurat::VlnPlot(object = combined.data, features = "PC_1", group.by = c("orig.ident", "seurat_clusters"))
p3 <- Seurat::VlnPlot(object = combined.data, features = "PC_2", group.by = c("orig.ident", "seurat_clusters"))
print(p2 + p3)
p4 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "harmony", group.by = c("orig.ident"))
print(p4)
p4_1 <- Seurat::DimPlot(object = combined.data, shuffle = TRUE, reduction = "harmony", group.by = c("seurat_clusters"))
print(p4_1)
p5 <- Seurat::VlnPlot(object = combined.data, features = "harmony_1", group.by = c("orig.ident", "seurat_clusters"))
p6 <- Seurat::VlnPlot(object = combined.data, features = "harmony_2", group.by = c("orig.ident", "seurat_clusters"))
print(p5 + p6)
p7 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("orig.ident"))
print(p7)
p7_1 <- Seurat::DimPlot(combined.data, shuffle = TRUE, reduction = "umap", group.by = c("seurat_clusters"))
print(p7_1)
p8 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_1", group.by = c("orig.ident", "seurat_clusters"))
p9 <- Seurat::VlnPlot(object = combined.data, features = "UMAP_2", group.by = c("orig.ident", "seurat_clusters"))
print(p8 + p9)
dev.off()
