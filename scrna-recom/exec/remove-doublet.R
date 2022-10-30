#!/usr/bin/env Rscript
library(argparser)
p <- arg_parser("scRNA-seq Doublet Removal using DoubletFinder")
p <- add_argument(p, "input", help="input SeuratObject rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--dims", type="numeric", default=30,
                  help="number of dimensionality used in reduction")
p <- add_argument(p, "--resolution", type="numeric", default=0.8,
                  help=" use a value above 1 if you want to obtain a larger number of clusters.")
p <- add_argument(p, "--nfeatures", type="numeric", default=2000,
                  help="number of variable features to use for scaledata and pca")
p <- add_argument(p, "--cores", type="numeric", default=10,
                  help="number of cores used in parallel")

argv <- parse_args(p)
print(argv$plot)

library(scrnaRecom)
input <- readRDS(argv$input)
group_remove_doublet(
  input, argv$outdir, argv$project, 
  nfeatures = argv$nfeatures, dims = argv$dims, resolution = argv$resolution, 
  cores = argv$cores
)

