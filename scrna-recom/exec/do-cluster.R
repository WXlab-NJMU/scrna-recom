#!/usr/bin/env Rscript
library(argparser)
library(Seurat)
p <- arg_parser("scRNA-seq Dimension Reduction and Clustering using Seurat")
p <- add_argument(p, "infile", help="input SeuratObject rds", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")

p <- add_argument(p, "--nfeatures", help="number of variable features", type="numeric", default = 2000)
p <- add_argument(p, "--method", help="reduction method, pca harmony or iNMF", type="character", default = "pca")
p <- add_argument(p, "--dims", help="dimensions, that is npcs for pca", type="numeric", default = 30)
p <- add_argument(p, "--resolution", help="resolution of clustering", type="numeric", default = 0.8)
p <- add_argument(p, "--kparam", help="k of knn in FindNeighbor", type="numeric", default = 20)
p <- add_argument(p, "--plot", nargs='*', type="character",
                  default = c("nFeature_RNA", "percent.mt", "percent.rb"),
                  help="features to plot on umap")
argv <- parse_args(p)
print(argv$plot)
input <- readRDS(argv$infile)

library(scrnaRecom)
scrnaRecom::clustering(input, argv$outdir, argv$project, 
                       nfeatures = argv$nfeatures, plot.features = argv$plot,
                       reduction = argv$method, dims = argv$dims, resolution = argv$resolution)