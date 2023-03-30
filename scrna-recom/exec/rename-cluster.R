#!/usr/bin/env Rscript
library(argparser)
library(Seurat)
p <- arg_parser("scRNA-seq Dimension Reduction and Clustering using Seurat")
p <- add_argument(p, "infile", help="input SeuratObject rds", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "groupfile", help="group file", type="character")
p <- add_argument(p, "markerfile", help="marker file", type="character")
p <- add_argument(p, "key", default="celltype_fine", help="celltype title", type="character")
argv <- parse_args(p)
input <- readRDS(argv$infile)

library(scrnaRecom)
scrnaRecom::renameClusterPlotMarkers(input, argv$outdir, argv$project,
                                     argv$groupfile, argv$markerfile,
                                     argv$key)
