library(argparser)
p <- arg_parser("scRNA Batch Correction")
#p <- add_argument(p, "csv", help="input matrix folder", type="character")
p <- add_argument(p, "input", help="input seurat rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--method", help="SeuratCCA, SeuratLargeData, Harmony, Liger", 
                  type="character", default = "SeuratCCA")
argv <- parse_args(p)
#print(argv)

library(scrnaRecom)
object <- readRDS(argv$input)
Integration(object, argv$outdir, argv$project, argv$method)
