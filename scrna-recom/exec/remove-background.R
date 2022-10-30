#!/usr/bin/env Rscript
library(argparser)
p <- arg_parser("scRNA-seq Background RNA Removal using SoupX")
p <- add_argument(p, "raw", help="cellranger raw_feature_bc_matrix folder", type="character")
p <- add_argument(p, "filtered", help="cellranger filtered_feature_bc_matrix folder", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
argv <- parse_args(p)
print(argv$plot)

library(scrnaRecom)
background.removal(argv$raw, argv$filtered, argv$outdir, argv$project)
