#' @include utils.R
#' 
NULL


#' @section scMetabolism:
#' * source code: <https://github.com/wu-yc/scMetabolism>
#' * quick start: <https://github.com/wu-yc/scMetabolism>
#' @md
#'
#' @param infile Input file path, supported format is h5ad, loom, rds
#' @param outdir Output path
#' @param pathway Interested Pathway, such as "Glycolysis / Gluconeogenesis"
#' @import Seurat
#' @import scMetabolism
#' @export
#' @rdname MetabolicAnalysis
#' @concept Metabolic Analysis
#' 
MetabolicAnalysis.scMetabolism <- function(infile, outdir, pathways = NULL){
  # 1. Quantify single-cell metabolism with Seurat
  ## method supports VISION (default), AUCell, ssgsea, and gsva
  ## ncores is the number of threads of parallel computation
  obj.seurat <- ReadtoSeuratObject(infile)
  countexp.Seurat<-scMetabolism::sc.metabolism.Seurat(obj = obj.seurat, method = "VISION", 
                                        imputation = F, ncores = 2, metabolism.type = "KEGG")
  outfile <- file.path(outdir, "metabolism.scmetabolism_score_matrix.tsv")
  write.table(countexp.Seurat@assays$METABOLISM$score, outfile)
  
  # 2. Visualize
  ## dimplot interested pathway
  if (!is.null(pathways)) {
    for (pathway in pathways){
      scMetabolism::DimPlot.metabolism(
        obj = countexp.Seurat, pathway = pathway, 
        dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)
    }
    scMetabolism::DotPlot.metabolism(obj = countexp.Seurat, pathway = pathways, 
                       phenotype = "ident", norm = "y")
    scMetabolism::BoxPlot.metabolism(obj = countexp.Seurat, pathway = pathways, 
                       phenotype = "ident", ncol = 1)    
  }
}