rm(list = ls())
gc()
library(patchwork)
library(SPOTlight)
library(semla)
library(Seurat)
library(qs)
library(harmony)

se <- qread('./semla_output/SCC.qs')
table(se$orig.ident)
name_list = list('CTB1NT','CTB1T','CTB3NT','CTB3T','CTB7NT','CTB7T')

se <- se |>
  NormalizeData() |>
  ScaleData() |>
  FindVariableFeatures() |>
  RunPCA() |>
  RunHarmony(group.by.vars = c("sample_id")) |>
  FindNeighbors(reduction = "harmony", dims = 1:10) |>
  FindClusters(resolution = 0.2) |>
  RunUMAP(reduction = "harmony", dims = 1:10)

MapLabels(se, 
          section_number = 1,
          #override_plot_dims = TRUE,
          column_name = "seurat_clusters", 
          ncol = 1) &
  theme(legend.position = "right")

# qsave(se, file = "./semla_output/se.qs")

se <- qread("./semla_output/se.qs")

cluster_colors <- setNames(RColorBrewer::brewer.pal(length(levels(se$seurat_clusters)), "Set3"), 
                           nm = levels(se$seurat_clusters))

p1 <- MapLabelsSummary(se, 
                       section_number = 5,
                       #override_plot_dims = TRUE,
                       label_by = "sample_id",
                       column_name = "seurat_clusters", 
                       ncol = 1, 
                       colors = cluster_colors)+
  ggtitle(name_list[[1]])
p2 <- MapLabelsSummary(se, 
                       section_number = 6,
                       #override_plot_dims = TRUE,
                       label_by = "sample_id",
                       column_name = "seurat_clusters", 
                       ncol = 1, 
                       colors = cluster_colors)+
  ggtitle(name_list[[2]])
p1
p2

table(se$celltype, se$seurat_clusters)

Idents(se) <- "orig.ident"
se@meta.data$orig.ident <- factor(se@meta.data$orig.ident, 
                                 levels = c('CTB1NT', 'CTB3NT', 'CTB7NT', 'CTB1T', 'CTB3T', 'CTB7T'))
se$orig.ident <- se@meta.data$sample_id 
levels(se$orig.ident)
levels(se$sample_id)
Idents(se) <- "seurat_clusters"
p1 <- DimPlot(se, 
              cols = cluster_colors,
              reduction = 'umap')
Idents(se) <- "celltype"
cluster_colors1 <- setNames(object = c("#FF9F1C", "#E71D36","#039BE5", "#FF6B6B", 
                              "#00A878", "#87CEEB","#6A4C93", "#FFD166", 
                              "#2EC4B6", "#F4511E"), 
                           nm = c('Keratinocyte','Fibroblast','Myeloid',
                                  "Endothelial",'T cell','Lymphatic Endothelial',
                                  'B cell','Smooth muscle','Melanocyte',
                                  'Mast cell'))
p2 <- DimPlot(se, 
              cols = cluster_colors1,
              reduction = 'umap')
p1+p2

new.cluster.ids <- c('unaffected dermis','granuloma region','epidermis',
                     'granuloma region','unaffected dermis','unaffected dermis',
                     'epidermis','unaffected dermis','granuloma region','granuloma region')
Idents(se) <- "seurat_clusters"
levels(se)
names(new.cluster.ids) <- levels(se)  
se <- RenameIdents(se, new.cluster.ids)
DimPlot(se, reduction = "umap", label = TRUE)
Idents(se)
se$region <- Idents(se)

cluster_colors2 <- setNames(
  object = c("#FFD700", "#FF4500", "#228B22"), 
  nm = c('epidermis', 'granuloma region', 'unaffected dermis')
)
p3 <- DimPlot(se, 
              cols = cluster_colors2,
              reduction = 'umap')
p1+p2+p3

cluster_colors <- setNames(RColorBrewer::brewer.pal(length(levels(se$seurat_clusters)), "Set3"), 
                           nm = levels(se$seurat_clusters))
se$region <- factor(se$region, levels = c('epidermis', 'granuloma region', 'unaffected dermis'))
MapLabelsSummary(se, 
                 section_number = 2,
                 #override_plot_dims = TRUE,
                 label_by = "sample_id",
                 column_name = "region", 
                 ncol = 1, 
                 colors = cluster_colors2)
MapLabels(se, 
          column_name = "region", 
          split_labels = TRUE, 
          section_number = 2, 
          ncol = 3,
          colors = cluster_colors2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top") & ggtitle(name_list[[2]])
k = 6
MapLabels(se, 
          column_name = "region", 
          # split_labels = TRUE, 
          # override_plot_dims = TRUE,
          section_number = k, 
          ncol = 1,
          colors = cluster_colors2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") & ggtitle(name_list[[k]])
ggsave(paste0('plots2/region_group_',name_list[k],'.jpeg'),height=8, width=7, dpi = 1000)

# Cell proportion
Cellratio <- prop.table(table(se@meta.data$region,se@meta.data$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('region','sample','Freq')
Cellratio$region <- factor(Cellratio$region,levels = c('epidermis', 'granuloma region', 'unaffected dermis'))
sample_colors <- c("#FFD700", "#FF4500", "#228B22")
ggplot(data = Cellratio, aes(x = sample, y = Freq, fill = region)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = sample_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 12),   
    plot.title = element_text(size = 12) 
  ) +
  labs(x = "", y = "Ratio", title = "Region Frequency by Sample in Spatial Data") +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )

# Cell proportion in region
se_gra <- subset(se, region=="granuloma region")
DimPlot(se_gra, 
        cols = cluster_colors2,
        reduction = 'umap')
table(se_gra@meta.data$celltype,se_gra@meta.data$orig.ident)
Cellratio <- prop.table(table(se_gra@meta.data$celltype,se_gra@meta.data$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('celltype','sample','Freq')
Cellratio$celltype <- factor(Cellratio$celltype,levels = c('Keratinocyte','Fibroblast','Myeloid',
                                                           "Endothelial",'T cell','Lymphatic Endothelial',
                                                           'B cell','Smooth muscle','Melanocyte',
                                                           'Mast cell'))
sample_colors <- c( "#FF9F1C",  "#E71D36",  "#039BE5",  "#FF6B6B",  
                    "#00A878",  "#87CEEB",  "#6A4C93",  "#FFD166",  
                    "#2EC4B6",  "#F4511E")
ggplot(data = Cellratio, aes(x = sample, y = Freq, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = sample_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 12),   
    plot.title = element_text(size = 12) 
  ) +
  labs(x = "", y = "Ratio", title = "Cell Type Frequency by Sample in Granuloma Region") +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
table(se_gra@meta.data$celltype,se_gra@meta.data$orig.ident)
table(se@meta.data$celltype,se@meta.data$orig.ident)
