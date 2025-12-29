rm(list = ls())
library(schard)
library(semla)
library(qs)
library(Seurat)
library(COSG)
library(tidyverse)
library(spacexr)
library(stringr)
library(scRNAtoolVis)
gc()
options(future.globals.maxSize = 1000 * 1024^2) 

result.RCTD.weights <- qread("qs_data/RCTD_results/result.RCTD.weights_full_subtype.qs")
result.RCTD <- qread("qs_data/RCTD_results/result.RCTD_full_subtype.qs")
head(result.RCTD.weights)[1:6]
head(result.RCTD)[1:6]

st_semla <- qread("qs_data/RCTD_results/RCTD_results_full_reference_subtype.qs")
table(st_semla$first_type, st_semla$sample_id)
sample_name <- c('CTB1_BL','CTB1_T','CTB3_BL','CTB3_T','CTB7_BL','CTB7_T')

?pheatmap::pheatmap
pheatmap::pheatmap(cor(result.RCTD.weights[,-ncol(result.RCTD.weights)]),
                   display_numbers = F,
                   RColorBrewer::brewer.pal(9,"Reds")) # viridis::mako(10)

cor_matrix <- cor(result.RCTD.weights[, -ncol(result.RCTD.weights)])
diag(cor_matrix) <- NA
library(RColorBrewer)
heatmap_colors <- brewer.pal(9, "Reds")  
na_color <- "grey90"                  
pheatmap::pheatmap(
  cor_matrix,
  color = c("white",heatmap_colors),
  na_col = na_color,          
  border_color = "white",     
  display_numbers = FALSE,     
  clustering_method = "complete", 
  fontsize_row = 8,
  fontsize_col = 8
)

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)  
cor_matrix <- cor(result.RCTD.weights[, -ncol(result.RCTD.weights)])
col_fun <- colorRamp2(
  breaks = seq(-0.1, 0.6, length.out = 8),# seq(min(cor_matrix), max(cor_matrix), length = 9),
  colors = brewer.pal(8, "Reds")
)
Heatmap(
  matrix = cor_matrix,
  col = col_fun,
  name = "Correlation",  
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "top",
  row_dend_side = "left",
  column_dend_side = "top",
  heatmap_legend_param = list(
    title_position = "leftcenter-rot"
  ),
  border = TRUE,  
  border_gp = gpar(col = "black", lwd = 1)
)
dev.off()

gc()
cols <- RColorBrewer::brewer.pal(11,"Spectral")|>rev()
MapFeatures(st_semla, 
            features = colnames(st_semla@meta.data)[29:30], 
            section_number = 1, 
            arrange_features = "row",
            colors = cols,
            ncol = 2,
            override_plot_dims = TRUE,
            pt_size =1)&ggtitle(sample_name[1])
MapFeatures(st_semla, 
            features = colnames(st_semla@meta.data)[35:36], 
            section_number = 1, 
            arrange_features = "row",
            colors = cols,
            ncol = 2,
            override_plot_dims = TRUE,
            pt_size =1)&ggtitle(sample_name[1])
MapFeatures(st_semla, 
            features = colnames(st_semla@meta.data)[55:56], 
            section_number = i, 
            arrange_features = "row",
            colors = cols,
            ncol = 2,
            override_plot_dims = TRUE,
            pt_size =1)&ggtitle(sample_name[i])
dev.off()
i = 3

library(ggrastr)
p <- SpatialFeaturePlot(st_semla, images = "CTB7T", 
                        features = colnames(st_semla@meta.data)[55])& 
  scale_y_reverse()&
  coord_fixed()&
  NoLegend()&
  theme_void()&
  DarkTheme()&
  # scale_fill_gradientn(colours=viridis_plasma_light_high, 
  #                      na.value = "black", 
  #                      limits=c(0,3))&
  ggtitle(sample_name[6])
p
rasterise(p,dpi = 300)
dev.off()
gc()


MapFeatures(st_semla, 
            features = colnames(st_semla@meta.data)[65:66], 
            section_number = i, 
            arrange_features = "row",
            colors = cols,
            ncol = 2,
            override_plot_dims = TRUE,
            pt_size =1)&ggtitle(sample_name[i])
factor_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                   '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                   '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3')
subtype_cols <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#FFE119", "#F58231", "#911EB4", "#46F0F0", "#F032E6",
  "#BCBD22", "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#AAFFC3", "#808000", "#FFD8B1",
  "#000075", "#A9A9A9", "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
  "#E377C2", "#7F7F7F", "#17BECF", "#FF9898", "#C5B0D5", "#C49C94", "#F7B6D2", "#C7C7C7",
  "#DBDB8D", "#9EDAE5", "#FFBB78", "#98DF8A", "#AEC7E8", "#FF9E9E", "#B5D0D9", "#D4B9D6",
  "#FDB462", "#80B1D3"
)
st_semla$first_type <- factor(st_semla$first_type, 
                      levels = c("Basal KC","Inflammatory KC",
                                 "Proliferating KC","Spinous KC",
                                 "Granular KC","Hair follicle",
                                 "Sweet Gland","Schwann Cell","Melanocyte",
                                 "CCL19+ FB","APCDD1+ FB","PI16+ FB",
                                 "POSTN+ FB","TNN+ FB","RAMP1+ FB",
                                 "Arterial EC","Capillary EC",
                                 "Venular EC1","Venular EC2",
                                 "Lymphatic Endothelial",
                                 "Smooth muscle","Mast cell",
                                 "Langerhans","cDC1", "cDC2A","cDC2B", "pDC",
                                 "M1-like Mac","M2-like Mac","SPP1 Mac",
                                 "TREM2 Mac", "CD8T", "Treg", "CD4 Naïve", 
                                 "Th17", "CD8T_NK", "Th17_Treg", "NK",
                                 "naïve B", "activated B", "memory B", "Plasma"))

desired_order <- c("Basal KC", "Inflammatory KC", 
                   "Proliferating KC", "Spinous KC",
                   "Granular KC", "Hair follicle",
                   "Sweet Gland", "Schwann Cell", "Melanocyte",
                   "CCL19+ FB", "APCDD1+ FB", "PI16+ FB",
                   "POSTN+ FB", "TNN+ FB", "RAMP1+ FB",
                   "Arterial EC", "Capillary EC",
                   "Venular EC1", "Venular EC2",
                   "Lymphatic Endothelial",
                   "Smooth muscle", "Mast cell",
                   "Langerhans", "cDC1", "cDC2A", "cDC2B", "pDC",
                   "M1-like Mac", "M2-like Mac", "SPP1 Mac",
                   "TREM2 Mac", "CD8T", "Treg", "CD4 Naïve", 
                   "Th17", "CD8T_NK", "Th17_Treg", "NK",
                   "naïve B", "activated B", "memory B", "Plasma")
current_cols <- colnames(st_semla@meta.data)
other_cols <- setdiff(current_cols, desired_order)
new_col_order <- c(other_cols, desired_order)
st_semla@meta.data <- st_semla@meta.data[, new_col_order]
# qsave(st_semla, file = "qs_data/RCTD_results/RCTD_results_full_reference_subtype.qs")

table(st_semla$first_type)
colnames(st_semla@meta.data)
MapLabels(st_semla, 
          section_number = 3,
          column_name = "first_type", 
          # image_use = "raw", 
          # pt_alpha = 0.6, 
          # override_plot_dims = TRUE,
          colors = subtype_cols,
          pt_size = 1) & 
  theme(legend.position = "right", legend.text = element_text(angle = 0)) &
  ggtitle(sample_name[5]) &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black")) &
  DarkTheme()

MapLabels(st_semla, 
          section_number = 5,
          column_name = "celltype", 
          colors = subtype_cols,
          pt_size = 1) & 
  theme(legend.position = "right", legend.text = element_text(angle = 0)) &
  ggtitle(sample_name[5]) &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))

cell_features <- colnames(st_semla@meta.data)[26:67]
MapMultipleFeatures(st_semla,  
                    section_number = 6,
                    features = cell_features[1:9], 
                    colors = subtype_cols[1:9], 
                    ncol = 1,
                    add_scalebar = TRUE, 
                    scalebar_gg = scalebar(x = 2000, 
                                           breaks = 11, 
                                           highlight_breaks = c(1, 6, 11)),
                    scalebar_position = c(0.55, 0.7), 
                    scalebar_height = 0.07,
                    pt_size = 1) & ggtitle(sample_name[1]) 
gc()
cell_features
MapMultipleFeatures(st_semla,  
                    section_number = 5,
                    features = cell_features[10:15], 
                    colors = subtype_cols[10:15], 
                    ncol = 1,
                    add_scalebar = TRUE, 
                    scalebar_gg = scalebar(x = 2000, breaks = 11, highlight_breaks = c(1, 6, 11)),
                    scalebar_position = c(0.55, 0.7), 
                    scalebar_height = 0.07,
                    pt_size = 1) & ggtitle(sample_name[1])
MapMultipleFeatures(st_semla,  
                    section_number = 5,
                    features = cell_features[28:31], 
                    colors = subtype_cols[28:31], 
                    ncol = 1,
                    add_scalebar = TRUE, 
                    scalebar_gg = scalebar(x = 2000, breaks = 11, highlight_breaks = c(1, 6, 11)),
                    scalebar_position = c(0.55, 0.7), 
                    scalebar_height = 0.07,
                    pt_size = 1) & ggtitle(sample_name[5])
MapFeatures(st_semla, 
            features = c("SPP1","TREM2"), 
            section_number = 1, 
            arrange_features = "row",
            colors = cols,
            ncol = 2,
            override_plot_dims = TRUE,
            pt_size =1)&ggtitle(sample_name[1])
MapFeatures(st_semla, 
            features = "SPP1", 
            section_number = 1, 
            arrange_features = "row",
            colors = cols,
            ncol = 1,
            min_cutoff = 0.3,
            # max_cutoff = 0.8,
            override_plot_dims = TRUE,
            pt_size =1
            )&ggtitle(sample_name[1])

# extract special cell identity
unique(st_semla$first_type)

filter_staffli_metadata <- function(seurat_obj, idents=c("SPP1 Mac", "TREM2 Mac")) {
  if (!"Staffli" %in% names(seurat_obj@tools)) {
    stop("对象中未找到Staffli工具")
  }
  
  # 过滤seurat
  cells_to_keep <- WhichCells(seurat_obj, idents = idents)
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  
  # 过滤Staffli元数据
  staffli_obj <- seurat_obj@tools$Staffli
  staffli_meta <- staffli_obj@meta_data
  staffli_meta_filtered <- staffli_meta[staffli_meta$barcode %in% cells_to_keep, ]
  
  # 更新并返回
  staffli_obj@meta_data <- staffli_meta_filtered
  seurat_obj@tools$Staffli <- staffli_obj
  return(seurat_obj)
}

se_small <- filter_staffli_metadata(st_semla, 
                                    idents=c("SPP1 Mac", "TREM2 Mac"))
Idents(st_semla) <- "first_type"
se_small <- filter_staffli_metadata(st_semla, 
                                    idents=cell_features[28:31])
se_small <- filter_staffli_metadata(st_semla, 
                                    idents=cell_features[c(2,10,30,31,35,37)])
dim(se_small)
table(se_small$first_type)
MapLabels(se_small, 
          section_number = 5,
          column_name = "first_type", 
          pt_alpha = 1, 
          colors = factor_colors[1:6],
          pt_size = 1) & 
  theme(legend.position = "right", legend.text = element_text(angle = 0)) &
  ggtitle(sample_name[3]) &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  DarkTheme()


## Create Spatial Reduction
# st_list <- split(st_semla, f = st_semla$sample_id)
st_list <- SplitObject(st_semla, split.by = "sample_id")

# rds <- st_list[[1]]
# spatial_corr <- rds@images$CTB1NT@coordinates[1:2]
# library(viridis)
# colnames(spatial_corr) <- c('s_1', 's_2')
# spatial_corr <- as.matrix(spatial_corr)
# rds[["spatial"]] <- CreateDimReducObject(embeddings = spatial_corr, key = "s_", assay = "Spatial")

st_list <- lapply(seq_along(st_list), FUN = function(x) {
  rds <- st_list[[x]]
  img_name <- names(rds@images)[1]
  spatial_corr <- rds@images[[img_name]]@coordinates[, 1:2]
  library(viridis)
  colnames(spatial_corr) <- c('s_1', 's_2')
  spatial_corr <- as.matrix(spatial_corr)
  rds[["spatial"]] <- CreateDimReducObject(
    embeddings = spatial_corr, 
    key = "s_", 
    assay = "Spatial"
  )
  return(rds)
})
gc()

# qsave(st_list, "qs_data/RCTD_results/RCTD_results_full_reference_subtype_split.qs")
st_list <- qread("qs_data/RCTD_results/RCTD_results_full_reference_subtype_split.qs")

rds <- st_list[[1]]
p1 <- FeaturePlot(rds, features = "nCount_Spatial", reduction = 'spatial') +
  scale_y_reverse() +
  coord_fixed() +
  scale_colour_viridis(option = "turbo",direction = 1) +
  # theme_void() +
  theme(legend.position = "right", plot.title = element_text(size = 12, face = "bold")) +
  ggtitle("Spatial UMI Count Distribution") &
  labs(subtitle = sample_name[6])
p1
dev.off()

rds <- st_list[[5]]
DimPlot(rds, group.by = "first_type", reduction = "spatial",pt.size = 1.5)&
  scale_y_reverse()&
  coord_fixed()&DarkTheme()& 
  scale_color_manual(values=subtype_cols)&
  ggtitle(sample_name[5])


gc()
library(ggrastr)

options(repr.plot.width=7, repr.plot.height=4)
p1 <- SpatialFeaturePlot(st_semla, images = "CTB1NT",features = "nCount_Spatial")
p2 <- p1 + theme_void() + xlab("") + ylab("") + scale_y_reverse() +
  scale_x_reverse() + coord_fixed() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        plot.background = element_rect(color = "black", size = 1), 
        legend.position = "top") 
p1 + p2

rds <- st_list[[1]]

SpatialDimPlot(rds, group.by = "first_type")& scale_y_reverse()&
  coord_fixed()&
  # theme_void()& 
  # theme(# axis.text = element_blank(), axis.ticks = element_blank(),
  #       plot.background = element_rect(color = "black", size = 1), 
  #       legend.position = "right")&
  DarkTheme()&
  scale_fill_manual(values=subtype_cols)&
  ggtitle(sample_name[5])&
  NoLegend()

table(rds$first_type)
dev.off()

## Plot multiple features
## Macrophages
library(patchwork)
Seurat::DefaultAssay(st_semla) <- "Spatial"
p1 <- MapMultipleFeatures(st_semla,
                          pt_size =1, 
                          section_number = 1, 
                          #max_cutoff =0.99,
                          override_plot_dims =TRUE,
                          features =c('CD68','TREM2','SPP1','CXCL9'))&
  ggtitle(sample_name[[1]])
p2 <- MapMultipleFeatures(st_semla,
                          pt_size =1, 
                          section_number = 2, 
                          #max_cutoff =0.99,
                          override_plot_dims =TRUE,
                          features =c('CD68','TREM2','SPP1','CXCL9'))&
  ggtitle(sample_name[[2]])
p <- wrap_plots(p1, p2, ncol = 2)
rasterise(p,dpi = 300)
dev.off()

p1 <- MapFeatures(st_semla, 
                  features = c("CD68","TREM2",'SPP1'), 
                  section_number = 1, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  # max_cutoff = 0.9, 
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(sample_name[[1]]) &
  ThemeLegendRight() &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))
p2 <- MapFeatures(st_semla, 
                  features = c("CD68","TREM2",'SPP1'), 
                  section_number = 2, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  # max_cutoff = 0.95, 
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(sample_name[[2]]) &
  ThemeLegendRight() &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))
p1
p2
p <- wrap_plots(p1, p2, ncol = 2)
rasterise(p2,dpi = 300)
dev.off()

## Tcells
p1 <- MapFeatures(st_semla, 
                  features = c("CD3D","G0S2",'RBPJ'), 
                  section_number = 1, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  # max_cutoff = 0.9, 
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(sample_name[[1]]) &
  ThemeLegendRight() &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))
p2 <- MapFeatures(st_semla, 
                  features = c("CD3D","G0S2",'RBPJ'), 
                  section_number = 2, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  # max_cutoff = 0.95, 
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(sample_name[[2]]) &
  ThemeLegendRight() &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))
p1
p2
rasterise(p2,dpi = 300)
dev.off()

## Fibroblasts
p1 <- MapFeatures(st_semla, 
                  features = c("COL1A2","FDCSP",'HLA-DRB5'), 
                  section_number = 1, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  # max_cutoff = 0.9, 
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(sample_name[[1]]) &
  ThemeLegendRight() &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))
p2 <- MapFeatures(st_semla, 
                  features = c("COL1A2","FDCSP",'HLA-DRB5'), 
                  section_number = 2, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  # max_cutoff = 0.9, 
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(sample_name[[2]]) &
  ThemeLegendRight() &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))
p1
p2
p <- wrap_plots(p1, p2, ncol = 2)
rasterise(p2,dpi = 300)
dev.off()
gc()


## new features
??semla::MapFeatures
p1 <- MapFeatures(st_semla, 
                  features = colnames(st_semla@meta.data)[55:56], 
                  section_number = 5, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  # max_cutoff = 0.9, 
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(sample_name[[5]]) &
  ThemeLegendRight() &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))



??semla:MapFeatures
p2 <- MapFeatures(st_semla, 
                  features = colnames(st_semla@meta.data)[55:56], 
                  section_number = 6, 
                  arrange_features = "row",
                  # scale_alpha=TRUE,
                  blend = TRUE,
                  # max_cutoff = 0.9, 
                  colors = cols,
                  ncol = 1,
                  override_plot_dims = TRUE,
                  pt_size =1)& ggtitle(sample_name[[6]]) &
  ThemeLegendRight() &
  guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))
p1
p2
p <- wrap_plots(p1, p2, ncol = 2)
rasterise(p2,dpi = 300)
dev.off()
gc()

??semla::MapMultipleFeatures
MapMultipleFeatures(st_semla,  
                    section_number = 1,
                    override_plot_dims = TRUE,
                    features = colnames(st_semla@meta.data)[55:56],  
                    colors = factor_colors[1:2], 
                    ncol = 1,
                    add_scalebar = TRUE, 
                    scalebar_gg = scalebar(x = 2000, 
                                           breaks = 11, 
                                           highlight_breaks = c(1, 6, 11)),
                    scalebar_position = c(0.55, 0.7), 
                    scalebar_height = 0.07,
                    pt_size = 1) & ggtitle(sample_name[[5]]) &
  ThemeLegendRight() &
  theme(plot.background = element_rect(fill = "grey98"),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "black"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))


??Seurat::SpatialFeaturePlot
p <- SpatialFeaturePlot(st_semla, 
                        images = "CTB1NT", 
                        blend = TRUE,
                        features = colnames(st_semla@meta.data)[55:56])& 
  scale_y_reverse()&
  coord_fixed()&
  theme_void()&
  DarkTheme()&
  ggtitle(sample_name[2])
rasterise(p,dpi = 300)
dev.off()

CellsByIdentities(object = st_semla)

st_semla$spec <- ifelse(grepl("TREM2",st_semla$first_type),"TREM2 Mac",
                        ifelse(grepl("SPP1",st_semla$first_type),
                               "SPP1 Mac","Other cells"))
table(st_semla$spec)
levels(st_semla$spec)
??Seurat::SpatialPlot

rds <- st_list[[1]]
DimPlot(rds, group.by = "spec", reduction = "spatial",pt.size = 0.0001)&
  scale_y_reverse()&
  coord_fixed()& 
  scale_color_manual(values=c("#E0EEEE","#DE2D26", '#911eb4'))&
  ggtitle(sample_name[1])
dev.off()

FeaturePlot(rds, features = c("TREM2 Mac", "SPP1 Mac"),
            blend = TRUE,
        reduction = "spatial",pt.size = 0.0001)&
  scale_y_reverse()&
  coord_fixed()& 
  # scale_fill_manual(values=c("#E0EEEE","#DE2D26"))&
  ggtitle(sample_name[1])
dev.off()


library(patchwork)
celltype_split <- function(sample){
  plist <- SpatialPlot(sample, 
                       cells.highlight = CellsByIdentities(object = sample),
                       cols.highlight = c("#DE2D26", "#E0EEEE"), 
                       facet.highlight = TRUE,
                       stroke = NA, 
                       combine = F, 
                       image.alpha = 0,
                       label.size=8, 
                       pt.size.factor=0.4)
  for (j in seq_len(length(plist))) {
    plist[[j]] <- plist[[j]] + 
      theme(plot.title = element_text(size=10,face="bold"))
  }
  patchwork::wrap_plots(plist, nrow = 1) + plot_annotation(title = sample_name[[i]])
}
Idents(rds) <- "spec"
celltype_split(rds)
dev.off()


#### Cell proportion
library(ggpubr)
fq <- prop.table(table(st_semla@meta.data$first_type, st_semla@meta.data[,"sample"]), margin=2) *100
colnames(fq) <- sample_name
df <- reshape2::melt(fq, value.name = "freq", varnames = c("subtype", "sample"))  #宽数据转长数据
df$Condition <- ifelse(grepl("BL",df$sample),"CTB_BL","CTB_T")
library(tidyr)
df_long <- pivot_wider(df,names_from=cell.type,values_from=freq)
library(RColorBrewer)
df$Condition <- factor(df$Condition,levels = c("CTB_BL","CTB_T"))
df$group <- df$Condition
str(df)
my_comparisons <- list( c("CTB_BL", "CTB_T") )
p <- ggboxplot(df, x = "subtype", y = "freq",
               color = "Condition", palette = c("#08519C","#F16913"),
               add = "jitter")
p + stat_compare_means(aes(group=Condition),
                       #label.y = c(50, 55, 80),
                       method="wilcox.test",
                       label="p.signif")+#
  theme_bw()+
  theme(axis.text.x = element_text(size=12,angle=45,hjust=1, color="black"),#x轴标签样式
        axis.text.y = element_text(size=12,colour = "black"),  
        axis.title.x = element_text(size=12,colour = 'black',vjust = -0.8,hjust = 0.5),#坐标轴标题
        axis.title.y = element_text(size=12,colour = 'black',vjust = -0.8,hjust = 0.5))+
  ylab("proportion")

# Cell proportion
Cellratio <- prop.table(table(st_semla@meta.data$first_type, st_semla@meta.data[,"sample"]), margin=2) 
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('celltype','sample','Freq')
dput(levels(Idents(sce)))
Cellratio$celltype <- factor(Cellratio$celltype,levels = c("CD8T", "Treg", "CD4 Naïve", "Th17", "CD8T/NK", "Th17/Treg", 
                                                           "NK"))
Cellratio$sample <- factor(Cellratio$sample, levels = c("N1","N2","N3","N4","N5","N6",
                                                        "CTB1_BL","CTB2_BL","CTB3_BL","CTB4_BL",
                                                        "CTB5_BL","CTB6_BL","CTB7_BL","CTB8_BL","CTB9_BL",
                                                        "CTB10_BL","CTB11_BL","CTB12_BL",
                                                        "CTB1_T","CTB2_T","CTB3_T","CTB5_T","CTB7_T",
                                                        "CTB8_T"))
sample_colors <- subtype_cols
p <- ggplot(data = Cellratio, aes(x = sample, y = Freq, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = sample_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    plot.title = element_text(size = 14) 
  ) +
  labs(x = "", y = "Ratio", title = "Cell Type Frequency by Sample in ST Data") +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
p
ggsave("./figures/cell_type_frequency.pdf", plot = p, width = 12, height = 6)