#!/usr/bin/env Rscript
library(argparser)
p <- arg_parser("scRNA Basic Data Analysis Using Seurat")
p <- add_argument(p, "indir", help="input matrix folder", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--maxrna", help="nFeature_RNA maximum", type="numeric", default = 2500)
p <- add_argument(p, "--minrna", help="nFeature_RNA minimum", type="numeric", default = 200)
p <- add_argument(p, "--mincell", help="cell minimum", type="numeric", default = 3)
p <- add_argument(p, "--maxmt", help="percent of maximum mt genes", type="numeric", default = 5)
p <- add_argument(p, "--sctransform", help="whether to use sctransform instead", type="logical", default = FALSE)
p <- add_argument(p, "--pc", help="rerun with fixed pc", type="numeric")
p <- add_argument(p, "--steps", nargs='*', type="numeric", default = c(1,2,3,4,5,8),
                  help="1-qc,2-norm, 3-feature select, 4-scale, 5-dimension reduce, 6-remove doublet, 7-diff, 8-type")
argv <- parse_args(p)
#print(argv)

library(scrnaRecom)
pc <- ifelse(is.na(argv$pc), NA_integer_, argv$pc)
input <- new("Input", indir = argv$indir,  outdir = argv$outdir, project = argv$project, 
             maxrna = argv$maxrna, minrna = argv$minrna, mincell = argv$mincell, maxmt = argv$maxmt, 
             steps = argv$steps, sctransform = argv$sctransform, pc = pc)
basic_analysis(input)
print(input)