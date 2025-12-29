rm(list = ls())
library(Seurat)
library(ggrepel)
library(qs)
library(semla)
library(reshape2)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
gc()

name_list = list('CTB1NT','CTB1T','CTB3NT','CTB3T','CTB7NT','CTB7T')

SCC <- qread('./semla_output/SCC.qs')
table(SCC$orig.ident)

# Plot with semla
cols <- viridis::rocket(11, direction = -1)
cols <- RColorBrewer::brewer.pal(11,"Spectral")|>rev()

p1 <- MapFeatures(SCC, 
                  features = c("CD68","CXCL9","SPP1"), 
                  section_number = 5, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[5]]) # & ThemeLegendRight()
p2 <- MapFeatures(SCC, 
                  features = c("CD68","CXCL9","SPP1"), 
                  section_number = 6, 
                  arrange_features = "row",
                  blend = TRUE,
                  # scale_alpha=TRUE,
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[6]]) # & ThemeLegendRight()
p1
p2

p1 <- MapFeatures(SCC, 
                  features = c("LYZ","CTSS","MMP9","CXCL9","CXCL10","CHIT1",'SPP1'), 
                  section_number = 5, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  # blend = TRUE,
                  # max_cutoff = 1.8, 
                  # min_cutoff = 0.1,
                  colors = cols,
                  ncol = 7,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[5]]) # & ThemeLegendRight()
p2 <- MapFeatures(SCC, 
                  features = c("LYZ","CTSS","MMP9","CXCL9","CXCL10","CHIT1",'SPP1'), 
                  section_number = 6, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  colors = cols,
                  ncol = 7,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[6]]) # & ThemeLegendRight()
p1
p2

p1 <- MapFeatures(SCC, 
                  features = c("APOC1","FBP1","NEAT1","IFI27","ANXA1",'MMP12'), 
                  section_number = 5, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  # blend = TRUE,
                  # max_cutoff = 1.8, 
                  # min_cutoff = 0.1,
                  colors = cols,
                  ncol = 6,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[5]]) # & ThemeLegendRight()
p2 <- MapFeatures(SCC, 
                  features = c("APOC1","FBP1","NEAT1","IFI27","ANXA1",'MMP12'), 
                  section_number = 6, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  colors = cols,
                  ncol = 6,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[6]]) # & ThemeLegendRight()
p1
p2

p1 <- MapFeatures(SCC, 
                  features = c("COL1A2","COL3A1","SPARC","IGKC","IGHG4",'IGHG1'), 
                  section_number = 3, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  # blend = TRUE,
                  # max_cutoff = 1.8, 
                  # min_cutoff = 0.1,
                  colors = cols,
                  ncol = 6,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[3]]) # & ThemeLegendRight()
p2 <- MapFeatures(SCC, 
                  features = c("COL1A2","COL3A1","SPARC","IGKC","IGHG4",'IGHG1'), 
                  section_number = 4, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  colors = cols,
                  ncol = 6,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[4]])
dev.off()
p1
p2

p1 <- MapFeatures(SCC, 
                  features = c("PI16","POSTN","ELN","IGHA1","IGLL5"), 
                  section_number = 1, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  # blend = TRUE,
                  # max_cutoff = 1.8, 
                  # min_cutoff = 0.1,
                  colors = cols,
                  ncol = 5,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[1]]) # & ThemeLegendRight()
p2 <- MapFeatures(SCC, 
                  features = c("PI16","POSTN","ELN","IGHA1","IGLL5"), 
                  section_number = 2, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  colors = cols,
                  ncol = 5,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[2]])
dev.off()
p1
p2

p1 <- MapFeatures(SCC, 
                  features = c("PI16","POSTN","ELN"), 
                  section_number = 3, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[3]]) # & ThemeLegendRight()
p2 <- MapFeatures(SCC, 
                  features = c("PI16","POSTN","ELN"), 
                  section_number = 2, 
                  arrange_features = "row",
                  blend = TRUE,
                  # scale_alpha=TRUE,
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(name_list[[2]]) # & ThemeLegendRight()
p1
p2
