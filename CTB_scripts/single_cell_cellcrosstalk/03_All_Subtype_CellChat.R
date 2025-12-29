rm(list = ls())
library(Seurat)
library(qs)
library(harmony)
library(tidyverse)
library(dittoSeq)
library(cowplot)
library(CellChat)
library(paletteer)
library(patchwork)


sce <- qread("data/ctb_scobj_v4.qs")
DimPlot(sce, group.by = "subtype", raster = T)
sce$disease <- factor(sce$disease, levels = c("Normal","CTB_BL","CTB_T"))

NS <- subset(sce,disease=="Normal")
BL <- subset(sce,disease=='CTB_BL')
BT <- subset(sce,disease=='CTB_T')
rm(sce);gc()

KS_cellchat <- function(input_obj,
                       assay= NULL,
                       group.by = NULL,
                       workers,
                       species=c('human','mouse'),
                       CellChatDB.use=NULL,#a character vector, which is a subset of c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact","Non-protein Signaling")
                       PPIuse=F,
                       type="triMean",# c("triMean", "truncatedMean", "thresholdedMean", "median")
                       min.cells = 10
){
  
  
  
  cellchat.obj = createCellChat(input_obj, assay = assay, group.by = group.by)
  
  if(species=='human'){
    
    CellChatDB <- CellChatDB.human
    ppi = PPI.human
  }
  
  if(species =="mouse"){
    
    CellChatDB <- CellChatDB.mouse
    ppi = PPI.mouse
  }
  
  
  if(is.null(CellChatDB.use)){
    
    cellchat.obj@DB <- CellChatDB
    
  }else{
    
    CellChatDB <- subsetDB(CellChatDB, search = CellChatDB.use, key = "annotation")
    cellchat.obj@DB <- CellChatDB
  }
  
  cellchat.obj <- subsetData(cellchat.obj) 
  future::plan("multisession", workers = workers) 
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  if(PPIuse==F){
    
    cellchat.obj <- computeCommunProb(cellchat.obj, type = type)
    
  }else{
    
    cellchat.obj <- projectData(cellchat.obj, ppi)
    cellchat.obj <- computeCommunProb(cellchat.obj, raw.use=F, type = type)
  }
  
  
  cellchat.obj <- filterCommunication(cellchat.obj, min.cells = min.cells)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  
  cellchat.obj <- aggregateNet(cellchat.obj)
  
  return(cellchat.obj)
  
}



BL.cellchat <- KS_cellchat(BL, 
                           assay = 'RNA', 
                           group.by = "subtype",
                           CellChatDB.use="Secreted Signaling",
                           workers=1, 
                           species='human')
saveRDS(BL.cellchat, file = "CellChat/BL.cellchat.rds")


NS.cellchat <- KS_cellchat(NS, 
                           assay = 'RNA', 
                           group.by = "subtype",
                           CellChatDB.use="Secreted Signaling",
                           workers=1, 
                           species='human')
saveRDS(NS.cellchat, file = "CellChat/NS.cellchat.rds")



rm(list = ls())
gc()
setwd("~/OneDrive_4_2025-3-8")
NS.cellchat  <- readRDS("CellChat/NS.cellchat.rds")
BL.cellchat  <- readRDS("CellChat/BL.cellchat.rds")
T.cellchat  <- readRDS("CellChat/T.cellchat.rds")


#merge cellchat obj of different group
object.list <- list(Normal = NS.cellchat, CTB_BL = BL.cellchat, CTB_T = T.cellchat)
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
}
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

object.list <- list(Normal = NS.cellchat, CTB_BL = BL.cellchat)
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
}
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

object.list <- list(CTB_T = T.cellchat, CTB_BL = BL.cellchat)
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
}
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", 
               color.use = c("#00A0E9","#F47252"), 
               stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison",
               color.use = c("#00A0E9","#F47252"), 
               stacked = F, do.stat = TRUE)
gg1 + gg2


gg1 <- compareInteractions(cellchat, 
                           show.legend = F, 
                           group = c(1,2,3), 
                           color.use = c("#009944", "#F47252","#00A0E9"))
gg2 <- compareInteractions(cellchat, 
                           show.legend = F, 
                           group = c(1,2,3), 
                           measure = "weight", 
                           color.use = c("#009944","#F47252","#00A0E9"))
gg1 + gg2
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, comparison = c(1, 2), weight.scale = T)
netVisual_diffInteraction(cellchat, comparison = c(1, 2), weight.scale = T, measure = "weight")



weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, 
                   weight.scale = T, 
                   label.edge= F, 
                   edge.weight.max = weight.max[2], 
                   edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

