library(argparser)
p <- arg_parser("scRNA Batch Correction")
#p <- add_argument(p, "csv", help="input matrix folder", type="character")
p <- add_argument(p, "input", help="input seurat rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--anno", type="character",
                  help="anno used: SingleR or scCATCH, default is null")
p <- add_argument(p, "--dims", type="numeric", default=30,
                  help="feature nums, npcs in Seurat::RunPCA, default is 30")
p <- add_argument(p, "--nfeatures", type="numeric", default=2000,
                  help="number of variable features to use for scaledata and pca, default is 2000")

argv <- parse_args(p)
print(argv$plot)

library(scrnaRecom)
input <- readRDS(argv$input)
remove_doublet(input, argv$outdir, argv$project, 
               nfeatures = argv$nfeatures, dims = argv$dims, anno.used = argv$anno)

