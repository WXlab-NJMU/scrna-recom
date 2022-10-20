library(argparser)
p <- arg_parser("scRNA Doublet Removal")
p <- add_argument(p, "input", help="input seurat rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--dims", type="numeric", default=50,
                  help="npcs in Seurat::RunPCA, default is 50")
p <- add_argument(p, "--nfeatures", type="numeric", default=2000,
                  help="number of variable features to use for scaledata and pca, default is 2000")

argv <- parse_args(p)
print(argv$plot)

library(scrnaRecom)
input <- readRDS(argv$input)
remove_doublet(input, argv$outdir, argv$project, 
               nfeatures = argv$nfeatures, dims = argv$dims)

