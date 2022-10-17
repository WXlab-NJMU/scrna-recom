library(argparser)
p <- arg_parser("scRNA Batch Correction")
#p <- add_argument(p, "csv", help="input matrix folder", type="character")
p <- add_argument(p, "input", help="input seurat rds file", type="character")
p <- add_argument(p, "outdir", help="output result folder", type="character")
p <- add_argument(p, "project", help="project name", type="character")
p <- add_argument(p, "--method", help="SingleR or scCATCH",
                  type="character", default = "SingleR")
p <- add_argument(p, "--reference", type="character",
                  default="MonacoImmuneData",
                  help="set when using SingleR,
                  HumanPrimaryCellAtlasData (general)、BlueprintEncodeData (pure stroma and immune)、MouseRNAseqData ( low-resolution bulk tissues)、
                  DatabaseImmuneCellExpressionData (exhaustive coverage)、DatabaseImmuneCellExpressionData(CD4+ T cell subsets)、
                  NovershternHematopoieticData (greatest resolution for myeloid and progenitor cells)、
                  MonacoImmuneData (best covers all of the bases for a typical PBMC sample)")
p <- add_argument(p, "--level", type="character", default="main",
                  help="set when using SingleR,
                  main (broad), fine (fine-grained), ont (standard in Cell Ontology)")
p <- add_argument(p, "--species", type="character", default="Human",
                  help="set when using scCATCH,
                  Human or Mouse")
p <- add_argument(p, "--tissue", type="character", 
                  help="set when using scCATCH,
                  Tissue name in scCATCH database")
p <- add_argument(p, "--minpct", type="numeric", default = 0.25, 
                  help="set when using scCATCH,
                  Include the gene detected in at least this many cells in each cluster")
p <- add_argument(p, "--logfc", type="numeric", default = 0.25, 
                  help="set when using scCATCH,
                  Include the gene with at least this fold change of average gene expression compared to every other clusters")
p <- add_argument(p, "--pvalue", type="numeric", default = 0.05, 
                  help="set when using scCATCH,
                  Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters")
p <- add_argument(p, "--plot", nargs='*', type="character",
                  default = c("nFeature_RNA", "percent.mt", "percent.rb"),
                  help="features to plot on umap")
argv <- parse_args(p)
print(argv$plot)

library(scrnaRecom)
object <- readRDS(argv$input)
if (argv$method == "SingleR"){
  AnnotateCellType(object, argv$outdir, argv$project, argv$method, 
                   reference = argv$reference, level = argv$level, 
                   plot.features = argv$plot)
} else if (argv$method == "scCATCH"){
  AnnotateCellType(object, argv$outdir, argv$project, argv$method, 
                   species = argv$species, tissue = argv$tissue, 
                   cell_min_pct = argv$minpct, logfc = argv$logfc, pvalue = argv$pvalue,
                   plot.features = argv$plot)
} else {
  print("Method not supported!!!")
  exit()
}