


rm(list = ls())
gc()
setwd("/home/data/gz0436/OneDrive_4_2025-3-8/Pyscenic")
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(data.table)
library(ggrepel)
library(philentropy) # install.packages("philentropy")
source("R/regulon_specificity.R")

seu <- qs::qread("~/OneDrive_4_2025-3-8/data/ctb_my_t_fb.qs")
DefaultAssay(seu) <- "AUCell"
# rows are cells and columns are features (regulons)
rasMat <- t(seu[["AUCell"]]@data)
dim(rasMat)

ctMat <- calIndMat(seu$celltype)
dim(ctMat)

rssMat <- calRSSMat(rasMat, ctMat)
dim(rssMat)

## Overall
plot.list <- lapply(levels(seu$celltype), function(xx) {
  PlotRegulonRank(rssMat, xx, topn=5)
})
cowplot::plot_grid(plotlist = plot.list[1:6], ncol = 6)
dev.off()
cowplot::plot_grid(plotlist = plot.list[7:12], ncol = 6)
dev.off()
cowplot::plot_grid(plotlist = plot.list[13:18], ncol = 6)
dev.off()
cowplot::plot_grid(plotlist = plot.list[19:22], ncol = 6)
dev.off()

## Specific cases
PlotRegulonRank(rssMat, "Th17/Treg")
FeaturePlot(seu, reduction = "umap", features = "RORC(38g)")
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "Th17/Treg") +
  DimPlot2(seu, reduction = "umap", regulon = "RORC(38g)")
dev.off()
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "Th17") +
  DimPlot2(seu, reduction = "umap", regulon = "RORC(38g)")
dev.off()
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "Treg") +
  DimPlot2(seu, reduction = "umap", regulon = "FOXP3(41g)")
dev.off()
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "CD4 Naïve") +
  DimPlot2(seu, reduction = "umap", regulon = "LEF1(11g)")
dev.off()
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "CD4 Naïve") +
  DimPlot2(seu, reduction = "umap", regulon = "FOXP1(61g)", threshold=0.15)
dev.off()
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "pDC") +
  DimPlot2(seu, reduction = "umap", regulon = "SPIB(153g)")
dev.off()

DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "SPP1 Mac") +
  DimPlot2(seu, reduction = "umap", regulon = "MTF1(34g)")
dev.off()
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "SPP1 Mac") +
  DimPlot2(seu, reduction = "umap", regulon = "MITF(58g)")
dev.off()
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "TREM2 Mac") +
  DimPlot2(seu, reduction = "umap", regulon = "MITF(58g)")
dev.off()
DimPlot2(seu, reduction = "umap", group.by = "celltype", group.highlight = "TREM2 Mac") +
  DimPlot2(seu, reduction = "umap", regulon = "MITF(58g)")
dev.off()

DefaultAssay(seu) <- "AUCell"
VlnPlot(seu, group.by = "celltype", features = c("MITF(58g)"), split.by = "disease",
        split.plot = TRUE, pt.size = 0, cols = c("red","blue")) + ylab("TF activity")
VlnPlot(seu, group.by = "celltype", features = c("MTF1(34g)"), split.by = "disease",
        split.plot = TRUE, pt.size = 0, cols = c("red","blue")) + ylab("TF activity")


unique(seu$celltype)

gc()
regulons <- clusterProfiler::read.gmt("output/02-ctb_my_t_fb.regulons.gmt")


regulons[grep("IRF7", regulons$term),]
regulons[grep("IRF1", regulons$term),]
regulons[grep("NCX1", regulons$term),]
regulons[grep("MTF1", regulons$term),]
regulons[grep("MTF1", regulons$gene),]   # CEBPB  IRF1  YBX1  MBD2
regulons[grep("CEBPB", regulons$gene),]
regulons[grep("IRF1", regulons$gene),]  # IRF7  STAT1
regulons[grep("YBX1", regulons$gene),]  # MTF1 SPP1
regulons[grep("MBD2", regulons$gene),]
grep('MTF1',regulons[grep("YBX1", regulons$term),]$gene)


regulons[grep("SPP1", regulons$gene),]  # HMGA1
regulons[grep("TREM2", regulons$gene),]  # SPI1  TFEB
regulons[grep("HMGA1", regulons$gene),]   

## TNF(Th1/17) -- > IRF7(SPP1 Mac) -->  IRF1 -->  MTF1 -->  downstream target genes
## TNF(Th1/17) -- > IRF7(SPP1 Mac) -->  IRF1 -->  CEBPB/MBD2 -->  MTF1  -->  downstream target genes

regulons[grep("IL15", regulons$gene),]





