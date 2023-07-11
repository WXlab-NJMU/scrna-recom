#' @include utils.R
#' @import patchwork
NULL


#' @section CellChat:
#' * source code: <https://github.com/wu-yc/scMetabolism>
#' * quick start: <https://github.com/wu-yc/scMetabolism>
#' @md
#'
#' @param input Input file, current supported format is seurat object
#' @param outdir Output path
#' @param species Species, only support human or mouse
#' @import NMF ggalluvial CellChat dplyr
#' @importFrom Seurat GetAssayData Idents
#' @export
#' @concept Cell to Cell Communication
#' @rdname CellCommunication
#'
CellCommunication.CellChat <- function(input, outdir, project, used,
                                       species = "human", cores = 20){
  if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
  outpdf <- glue::glue("{outdir}/{project}.interaction.cellchat.output.pdf")
  pdf(outpdf)
  data.input <- Seurat::GetAssayData(input, assay = "RNA", slot = "data")
  labels <- Seurat::Idents(input)
  meta <- data.frame(group = labels, row.names = names(labels))
  ## add C in seurat cluster
  ##meta <- data.frame(labels = paste0("C", Seurat::Idents(obj.seurat)),
  #                   row.names = names(Seurat::Idents(obj.seurat)))
  cellchat <- CellChat::createCellChat(object = data.input,
                                       meta = meta, group.by = "group")
  # CellChat database
  ## current human and mouse are available
  if (species == "human"){
    CellChatDB <- CellChat::CellChatDB.human
  }else if (species == "mouse"){
    CellChatDB <- CellChat::CellChatDB.mouse
  } else {
    print("Only support human or mouse!!!")
  }
  CellChat::showDatabaseCategory(CellChatDB)
  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)
  # use a subset of CellChatDB for cell-cell communication analysis
  cellchat@DB <- CellChatDB
  # preprocessing the expression data
  cellchat <- CellChat::subsetData(cellchat)
  #future::plan("multisession", workers = cores) # do parallel
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  saveRDS(cellchat, glue::glue("{outdir}/{project}.interaction.cellchat.init.rds"))

  # Part II: Inference of cell-cell communication network
  message("Part II: Inference of cell-cell communication network")
  ## Compute the communication probability and infer cellular communication network
  ## default is trimean (25%), set (type="truncatedMean", trim = 0.1) to get more results
  cellchat <- CellChat::computeCommunProb(cellchat)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 5)
  df.net <- CellChat::subsetCommunication(cellchat)
  write.table(df.net, glue::glue("{outdir}/{project}.interaction.cellchat.inferred_network_ligands.tsv"))
  df.netp <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
  write.table(df.netp, glue::glue("{outdir}/{project}.interaction.cellchat.inferred_network_pathway.tsv"))
  ## Calculate the aggregated cell-cell communication network
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  print(cellchat@netP$pathways)
  saveRDS(cellchat, glue::glue("{outdir}/{project}.interaction.cellchat.infer.rds"))
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), mar=c(0.2,0.6,0.6,0.5), xpd=TRUE)
  CellChat::netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

  mat <- cellchat@net$weight
  par(mfrow = c(4,3), mar=c(0.8,0.5,0.8,0.5), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    celltype <- rownames(mat)[i]
    message(celltype)
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    CellChat::netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name =  celltype)
  }
  # Part III: Visualization of cell-cell communication network
  message("Part III: Visualization of cell-cell communication network")
  ## plot a specified pathway
  pathways.show <- cellchat@netP$pathways
  message("Enriched pathways: ", paste(pathways.show, delimiter = " "))
  # set recevier to interest
  #vertex.receiver = seq(1,2)
  for (pathway in pathways.show){
    #par(mfrow = c(1,2), xpd=TRUE)
    #p1 <- CellChat::netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver)
    #print(p1)
    par(mfrow = c(2,1), mar=c(0.2,0.5,0.2,0.5))
    CellChat::netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
    CellChat::netVisual_aggregate(cellchat, signaling = pathway, layout = "chord", pt.title = 8)
    #plot base on single object, could not use par
    par(mar=c(0.2,0.5,0.2,0.5))
    CellChat::netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds")
    ### Compute the contribution
    par(mar=c(0.2,0.5,0.2,0.5))
    CellChat::netAnalysis_contribution(cellchat, signaling = pathway,
                                       title = paste0("Contribution of L-R pairs in ", pathway))
    #### show one ligand-receptor pair
    pairLR <- CellChat::extractEnrichedLR(cellchat, signaling = pathway, geneLR.return = FALSE)
    p <- CellChat::netVisual_individual(cellchat, signaling = pathway, pairLR.use = pairLR[,1], 
                                   layout = "circle", nCol = 3, signaling.name = pairLR[,1])
    print(p)
    ### plot expression distribution
    ###  only inferred significant communications, set `enriched.only = FALSE` for more
    p <- CellChat::plotGeneExpression(cellchat, signaling = pathway)
    print(p)
  }

  p <- CellChat::netVisual_bubble(cellchat, remove.isolate = TRUE, angle.x = 45)
  print(p)
  #CellChat::netVisual_bubble(cellchat, signaling = cellchat@netP$pathways, remove.isolate = FALSE)

  #pairLR.use <- CellChat::extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways)
  #CellChat::netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
  #CellChat::netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
  #CellChat::netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)


  # Part IV: Systems analysis of cell-cell communication network
  cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  saveRDS(cellchat, glue::glue("{outdir}/{project}.interaction.cellchat.system.rds"))
  par(mfrow = c(4,2), mar=c(0.2,0.5,0.2,0.5))
  CellChat::netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                              width = 8, height = 4, font.size = 10)
  ## Visualize the dominant senders (sources) and receivers (targets) in a 2D space
  gg1 <- CellChat::netAnalysis_signalingRole_scatter(cellchat)
  gg2 <- CellChat::netAnalysis_signalingRole_scatter(cellchat, signaling = cellchat@netP$pathways)
  print(gg1 + gg2 + patchwork::plot_layout(ncol = 1))
  ## Identify signals contributing most to outgoing or incoming signaling
  par(mfrow = c(2,1), mar=c(0.2,0.5,0.2,0.5))
  ht1 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", 
                                                     font.size = 8, height = 12, width = 10)
  ht2 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                                     font.size = 8, height = 12, width = 10)
  print(ht1)
  print(ht2)
  
  ### Identify and visualize outgoing communication pattern of secreting cells
  library(NMF)
  library(ggalluvial)
  p <- CellChat::selectK(cellchat, pattern = "outgoing")
  print(p)
  data <- p$data
  data[["slope"]] <- c(-diff(data$score),0)
  print(data)
  #### use the count of drop suddenly in last step
  nPatterns <- data$k[
    which(sapply(seq(16), function(x) data$slope[x] > 0 && data$slope[x] > data$slope[x+1]))][1]
  if (nPatterns > 5) nPatterns = 5
  plot(0,type='n',axes=FALSE,ann=FALSE)
  cellchat <- CellChat::identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  p1 <- CellChat::netAnalysis_river(cellchat, pattern = "outgoing")
  p2 <- CellChat::netAnalysis_dot(cellchat, pattern = "outgoing")
  print(p1 + p2 + patchwork::plot_layout(ncol = 1))
  ### Identify and visualize incoming communication pattern of target cells
  p <- CellChat::selectK(cellchat, pattern = "incoming")
  print(p)
  data <- p$data %>% mutate(slope = c(-diff(data$score),0))
  nPatterns <- data$k[
    which(sapply(seq(16), function(x) data$slope[x] > 0 && data$slope[x] > data$slope[x+1]))][1]
  if (nPatterns > 5) nPatterns = 5
  #### Plot pattern
  plot.new()
  cellchat <- CellChat::identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  p1 <- CellChat::netAnalysis_river(cellchat, pattern = "incoming")
  p2 <- CellChat::netAnalysis_dot(cellchat, pattern = "incoming")
  print(p1 + p2 + patchwork::plot_layout(ncol = 1))
  ## Manifold and classification learning analysis of signaling network
  ### Identify signaling groups based on their functional similarity
  cellchat <- CellChat::computeNetSimilarity(cellchat, type = "functional")
  cellchat <- CellChat::netEmbedding(cellchat, type = "functional")
  cellchat <- CellChat::netClustering(cellchat, type = "functional", do.parallel = T, nCores = cores)
  p1 <- CellChat::netVisual_embedding(cellchat, type = "functional", label.size = 3.5, 
                                      title = "functional similarity")
  ### Identify signaling groups based on structure similarity
  cellchat <- CellChat::computeNetSimilarity(cellchat, type = "structural")
  cellchat <- CellChat::netEmbedding(cellchat, type = "structural")
  cellchat <- CellChat::netClustering(cellchat, type = "structural", do.parallel = T, nCores = cores)
  p2 <- CellChat::netVisual_embedding(cellchat, type = "structural", label.size = 3.5, 
                                      title = "structure similarity")
  print(p1 + p2 + patchwork::plot_layout(ncol = 1))

  # Part V: Save the CellChat object
  saveRDS(cellchat, file = glue::glue("{outdir}/{project}.interaction.cellchat_final.rds"))
  dev.off()
}
