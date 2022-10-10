library(argparser)
p <- arg_parser("scRNA Batch Correction")
#p <- add_argument(p, "csv", help="input matrix folder", type="character")
p <- add_argument(p, "input", help="input seurat rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "--method", help="SeuratCCA, SeuratLargeData, Harmony, Liger", 
                  type="character", default = "SeuratCCA")
p <- add_argument(p, "--maxrna", help="nFeature_RNA maximum", type="numeric", default = 2500)
p <- add_argument(p, "--minrna", help="nFeature_RNA minimum", type="numeric", default = 200)
p <- add_argument(p, "--mincell", help="cell minimum", type="numeric", default = 3)
p <- add_argument(p, "--maxmt", help="percent of maximum mt genes", type="numeric", default = 5)
argv <- parse_args(p)
#print(argv)

library(scrnaRecom)
object <- readRDS(argv$input)
Integration(object, argv$outdir, argv$method,
            argv$mincell, argv$minrna, argv$maxrna, argv$maxmt)
