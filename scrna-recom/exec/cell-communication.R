#!/usr/bin/env Rscript
library(argparser)
library(tictoc)
p <- arg_parser("scRNA-seq Dimension Reduction and Clustering using Seurat")
p <- add_argument(p, "infile", help="input SeuratObject rds", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")

p <- add_argument(p, "--method", help="current support is CellChat",
                  type="character", default = "CellChat")
p <- add_argument(p, "--species", help="human or mouse",
                  default = "human")
p <- add_argument(p, "--cores", help="cores to do parallel",
                  default = 20)
argv <- parse_args(p)
tic("Runing Interaction Analysis")
input <- readRDS(argv$infile)
library(scrnaRecom)
CellCommunication(input, argv$outdir, argv$project, argv$method,
                  species = argv$species, cores = argv$cores
)
toc()
