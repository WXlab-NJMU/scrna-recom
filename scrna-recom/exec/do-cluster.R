library(argparser)
p <- arg_parser("scRNA Basic Data Analysis Using Seurat")
p <- add_argument(p, "indir", help="input matrix folder", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")

p <- add_argument(p, "--nfeatures", help="number of variable features", type="numeric", default = 2000)
p <- add_argument(p, "--reduction", help="reduction method", type="character", default = "pca")
p <- add_argument(p, "--dim", help="dimensions", type="numeric", default = 30)
p <- add_argument(p, "--resolution", help="resolution of cluster", type="numeric", default = 0.5)
p <- add_argument(p, "--plot", nargs='*', type="character", 
                  default = c("nFeature_RNA", "percent.mt", "percent.rb"),
                  help="features to plot on umap")

input <- readRDS(argv$infile)
clustering(input, argv$outdir, argv$project, 
           nfeatures = argv$nfeatures, plot.features = argv$plot,
           reduction = argv$reduction, dims = argv$dims, resolution = argv$resolution)