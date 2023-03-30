#!/usr/bin/env Rscript
library(argparser)
p <- arg_parser("scRNA-seq Batch Correction using Seurat, Liger or Harmony")
#p <- add_argument(p, "csv", help="input matrix folder", type="character")
p <- add_argument(p, "input", help="input SeuratObject rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--method", help="SeuratCCA, SeuratRPCA, Harmony, Liger",
                  type="character", default = "SeuratCCA")
p <- add_argument(p, "--dims", default="auto",
                  help="feature nums, npc in Seurat::RunPCA or k in rliger::optimizeALS, default is 30")
p <- add_argument(p, "--nfeatures", type="numeric", default=2000,
                  help="number of variable features to use for scaledata and pca, default is 2000")
p <- add_argument(p, "--kparam", type="numeric", default=20,
                  help="k.param of knn in FindNeighbor, default is 20")
p <- add_argument(p, "--resolution", type="numeric", default=0.5,
                  help="resolution of cluster, default is 0.5")
argv <- parse_args(p)

library(scrnaRecom)
object <- readRDS(argv$input)
Integration(object, argv$outdir, argv$project, argv$method, argv$dims,
            nfeatures = argv$nfeatures,
            resolution = argv$resolution, k = argv$kparam,
            )
