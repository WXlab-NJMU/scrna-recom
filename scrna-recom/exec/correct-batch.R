library(argparser)
p <- arg_parser("scRNA Batch Correction")
#p <- add_argument(p, "csv", help="input matrix folder", type="character")
p <- add_argument(p, "input", help="input seurat rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--method", help="SeuratCCA, SeuratRPCA, Harmony, Liger",
                  type="character", default = "SeuratCCA")
p <- add_argument(p, "--dim", type="numeric", default=30,
                  help="feature nums, npc in Seurat::RunPCA or k in rliger::optimizeALS")
p <- add_argument(p, "--nfeatures", type="numeric", default=2000,
                  help="number of variable features to use for scaledata and pca")
argv <- parse_args(p)
#print(argv)

library(scrnaRecom)
object <- readRDS(argv$input)
Integration(object, argv$outdir, argv$project, argv$method, argv$dim, nfeatures = argv$nfeatures)
