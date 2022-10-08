#' @include utils.R
#' @import patchwork
NULL


#' @section CellChat:
#' * source code: <https://github.com/wu-yc/scMetabolism>
#' * quick start: <https://github.com/wu-yc/scMetabolism>
#' @md
#'
#' @param infile Input file path, supported format is h5ad, loom, rds
#' @param outdir Output path
#' @param species Species, only support human or mouse
#' @import NMF ggalluvial CellChat
#' @importFrom Seurat GetAssayData Idents
#' @export
#' @concept Cell to Cell Communication
#' @rdname CellCommunication
#' 
CellCommunication.CellChat <- function(infile, outdir, species = NULL){
  outpdf <- file.path(outdir, "communication.cellchat_output.pdf")
  pdf(outpdf)
  data.input <- Seurat::GetAssayData(obj.seurat, assay = "RNA", slot = "data")
  ## add C in seurat cluster
  meta <- data.frame(labels = paste0("C", Seurat::Idents(obj.seurat)), 
                     row.names = names(Seurat::Idents(obj.seurat))) 
  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "labels")
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
  CellChatDB.use <- CellChat::subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  # preprocessing the expression data
  cellchat <- CellChat::subsetData(cellchat)
  future::plan("multiprocess", workers = 4) # do parallel
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  
  
  # Part II: Inference of cell-cell communication network
  ## Compute the communication probability and infer cellular communication network
  cellchat <- CellChat::computeCommunProb(cellchat)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  df.net <- CellChat::subsetCommunication(cellchat)
  write.table(df.net, "CellChat.inferred_network_ligands.tsv")
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  df.netp <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
  write.table(df.netp, "CellChat.inferred_network_pathway.tsv")
  ## Calculate the aggregated cell-cell communication network
  cellchat <- CellChat::aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  p1 <- CellChat::netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  p2 <- CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  # Part III: Visualization of cell-cell communication network
  ## plot a specified pathway
  pathways.show <- head(cellchat@netP$pathways)
  vertex.receiver = seq(1,4) 
  for (pathway in pathways.show){
    par(mfrow = c(1,2), xpd=TRUE)
    p1 <- CellChat::netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver)
    p2 <- CellChat::netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
    #plot base on single object, could not use par
    CellChat::netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds")
    
    ### Compute the contribution
    CellChat::netAnalysis_contribution(cellchat, signaling = pathway)
    pairLR <- CellChat::extractEnrichedLR(cellchat, signaling = pathway, geneLR.return = FALSE)
    #### show one ligand-receptor pair
    LR.show <- pairLR[1,]
    vertex.receiver = seq(1,4)
    CellChat::netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR.show, layout = "circle")
    CellChat::netVisual_individual(cellchat, signaling = pathway, pairLR.use = LR.show, layout = "chord")
    
    ### plot expression distribution
    ###  only inferred significant communications, set `enriched.only = FALSE` for more
    CellChat::plotGeneExpression(cellchat, signaling = pathway)
  }
  
  CellChat::netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
  CellChat::netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = cellchat@netP$pathways, remove.isolate = FALSE)
  
  pairLR.use <- CellChat::extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways)
  CellChat::netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
  CellChat::netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
  CellChat::netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
  
  
  # Part IV: Systems analysis of cell-cell communication network
  cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  CellChat::netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  ## Visualize the dominant senders (sources) and receivers (targets) in a 2D space
  gg1 <- CellChat::netAnalysis_signalingRole_scatter(cellchat)
  gg2 <- CellChat::netAnalysis_signalingRole_scatter(cellchat, signaling = cellchat@netP$pathways)
  gg1 + gg2
  ## Identify signals contributing most to outgoing or incoming signaling
  ht1 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  ht1 + ht2
  
  ### Identify and visualize outgoing communication pattern of secreting cells
  library(NMF)
  library(ggalluvial)
  CellChat::selectK(cellchat, pattern = "outgoing")
  #### use the count of drop suddenly in last step
  nPatterns = 3
  cellchat <- CellChat::identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  CellChat::netAnalysis_river(cellchat, pattern = "outgoing")
  CellChat::netAnalysis_dot(cellchat, pattern = "outgoing")
  ### Identify and visualize incoming communication pattern of target cells
  CellChat::selectK(cellchat, pattern = "incoming")
  CellChat::netAnalysis_river(cellchat, pattern = "incoming")
  CellChat::netAnalysis_dot(cellchat, pattern = "incoming")
  ## Manifold and classification learning analysis of signaling network
  ### Identify signaling groups based on their functional similarity
  cellchat <- CellChat::computeNetSimilarity(cellchat, type = "functional")
  cellchat <- CellChat::netEmbedding(cellchat, type = "functional")
  cellchat <- CellChat::netClustering(cellchat, type = "functional")
  CellChat::netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
  ### Identify signaling groups based on structure similarity
  cellchat <- CellChat::computeNetSimilarity(cellchat, type = "structural")
  cellchat <- CellChat::netEmbedding(cellchat, type = "structural")
  cellchat <- CellChat::netClustering(cellchat, type = "structural")
  CellChat::netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
  
 
  # Part V: Save the CellChat object
  saveRDS(cellchat, file = file.path(outdir, "communication.cellchat_final.rds"))
  dev.off()
}