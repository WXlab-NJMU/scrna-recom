library(argparser)
p <- arg_parser("scRNA-seq quality control using Seurat")
p <- add_argument(p, "csv", help="csv file including sample, path, qc(specific to sample)", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--max.genes", help="nFeature_RNA maximum", type="numeric", default = 5000)
p <- add_argument(p, "--min.genes", help="nFeature_RNA minimum", type="numeric", default = 200)
p <- add_argument(p, "--max.counts", help="nCount_RNA maximum", type="numeric", default = 40000)
p <- add_argument(p, "--min.counts", help="nCount_RNA minimum", type="numeric", default = 500)
p <- add_argument(p, "--min.cells", help="cell minimum", type="numeric", default = 3)
p <- add_argument(p, "--max.mt", help="percent of maximum mt genes", type="numeric", default = 20)
p <- add_argument(p, "--max.hb", help="percent of maximum hb genes", type="numeric", default = 10)
argv <- parse_args(p)
#print(argv)

library(scrnaRecom)
scrnaRecom::group_qc(
  argv$csv, argv$outdir, argv$project, 
  argv$max.counts, argv$min.counts, argv$max.genes, argv$min.genes,
  argv$min.cells, argv$max.mt, argv$max.hb)