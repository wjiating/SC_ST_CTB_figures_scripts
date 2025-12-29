rm(list = ls())
gc()
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(ggraph)
library(tidygraph)
setwd("/home/data/gz0436/OneDrive_4_2025-3-8/Pyscenic")
source("R/IO.R")

data <- LoadpySCENICOutput(regulon.gmt = "output/02-ctb_my_t_fb.regulons.gmt",
                           adj.mat.file = "output/01-step1_adj.tsv")
head(data)
summary(data$importance)
# data <- subset(data, importance > 1)

hallmarks <- clusterProfiler::read.gmt("resource/h.all.v2022.1.Hs.symbols.gmt")
genes <- subset(hallmarks, term == "HALLMARK_MTORC1_SIGNALING")$gene
source("R/network_plot.R")

regulon.clusters <- readRDS("output/03-2.ctb_my_t_fb.regulon_modules.rds")
table(regulon.clusters$cluster)

sub("\\([0-9]+g\\)", "", subset(regulon.clusters, cluster=="M1")$regulon)
sub("\\([0-9]+g\\)", "", subset(regulon.clusters, cluster=="M8")$regulon)
subset(regulon.clusters, cluster=="M8")$regulon

for(module in paste0("M", 1:8)) {
  cat("\n=== Module", module, "===\n")
  regs <- subset(regulon.clusters, cluster == module)$regulon
  print(regs)
}

data[grep("IRF7",data$TF),]
data[grep("IRF1",data$TF),]
head(data)

RegulonGraphVis(data, tf.show = c("IRF7","IRF1","MTF1"),
                targets.show = c("MT1B","MT2A"),
                layout = "kk",
                prop = 0.5)

RegulonGraphVis(data, tf.show = sub("\\([0-9]+g\\)", "", subset(regulon.clusters, cluster=="M8")$regulon),
                targets.show = genes, 
                # layout = "circle",
                n = 20)


