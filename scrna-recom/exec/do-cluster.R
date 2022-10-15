library(argparser)
library(Seurat)
p <- arg_parser("scRNA Basic Data Analysis Using Seurat")
p <- add_argument(p, "infile", help="input seurat object", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")

p <- add_argument(p, "--nfeatures", help="number of variable features", type="numeric", default = 2000)
p <- add_argument(p, "--method", help="reduction method", type="character", default = "pca")
p <- add_argument(p, "--dims", help="dimensions", type="numeric", default = 30)
p <- add_argument(p, "--resolution", help="resolution of cluster", type="numeric", default = 0.5)
p <- add_argument(p, "--plot", nargs='*', 
                  default = c("nFeature_RNA", "percent.mt", "percent.rb"),
                  help="features to plot on umap")
argv <- parse_args(p)
#print(argv)
input <- readRDS(argv$infile)

library(scrnaRecom)
scrnaRecom::clustering(input, argv$outdir, argv$project, 
                       nfeatures = argv$nfeatures, plot.features = argv$plot,
                       reduction = argv$method, dims = argv$dims, resolution = argv$resolution)