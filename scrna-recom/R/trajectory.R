#' @include utils.R
NULL

#' @section Monocole3:
#' * source code: <https://github.com/cole-trapnell-lab/monocle3>
#' * quick start: <https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/>
#' * seuratwrapper: <http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html>
#' @md
#' 
#' @param reduction_method Currently "UMAP", "tSNE", "PCA", "LSI", and "Aligned" are supported.
#' @param root.key Key used to find root
#' @param root.value Value of root key
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom monocle3 cluster_cells plot_cells learn_graph principal_graph order_cells
#' @importFrom SeuratWrappers as.cell_data_set
#' @importFrom Seurat as.Seurat
#' @importFrom patchwork wrap_plots
#' @export
#' @rdname Trajectory
#' @method Trajectory Monocole3
#' 
Trajectory.Monocole3 <- function(input, outdir, used, 
                                 reduction = "UMAP", 
                                 root.key = NULL, root.value = NULL) {
  obj.seurat <- ReadtoSeuratObject(infile)
  #obj.seurat <- pbmc
  cds <- SeuratWrappers::as.cell_data_set(obj.seurat)
  cds <- monocle3::cluster_cells(cds, reduction_method = reduction)
  outpdf <- file.path(outdir, "trajectory.monocle3_output.pdf")
  pdf(outpdf)
  p1 <- monocle3::plot_cells(cds, show_trajectory_graph = FALSE)
  p2 <- monocle3::plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
  p3 <- monocle3::plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = TRUE)
  patchwork::wrap_plots(p1, p2, p3)
  
  integrated.sub <- subset(Seurat::as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
  cds <- SeuratWrappers::as.cell_data_set(integrated.sub)
  cds <- monocle3::learn_graph(cds)
  #genes <- c("ceh-36", "dmd-6")
  #monocle3::plot_cells(cds, genes=genes, label_cell_groups=FALSE, show_trajectory_graph=FALSE)
  monocle3::plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
  # get root
  get_root_node <- function(cds, key, value){
    cell_ids <- which(SummarizedExperiment::colData(cds)[, key] == value)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[
      as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
  }
  # plot pseudotime with root
  if (!is.nan(root.key)){
    cds <- monocle3::order_cells(cds, 
      root_pr_nodes = get_root_node(cds, key = root.key, value = root.value))
    # plot root
    monocle3::plot_cells(
      cds, color_cells_by = "pseudotime", 
      label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
    # plot pseudotime in reduced dimension
    #cds <- monocle3::reduce_dimension(cds, reduction_method = reduction)
    integrated.sub <- Seurat::as.Seurat(cds, assay = NULL)
    Seurat::FeaturePlot(integrated.sub, "monocle3_pseudotime")
  }
  dev.off()
  saveRDS(cds, file = file.path(outdir, "trajectory.monocle3_final.rds"))
}

#' @section scVelo:
#' * source code: <https://github.com/dviraran/SingleR>
#' * quick start: <https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html>
#' @md
#' 
#' @export
#' @rdname Trajectory
#' @method Trajectory scVelo
#' 
Trajectory.scVelo <- function(input, outdir, used) {
  print("python command line")  
}

