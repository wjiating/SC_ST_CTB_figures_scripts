# B10iocManager::install('saezlab/decoupleR')
# BiocManager::install('OmnipathR')

rm(list = ls())
library(decoupleR)
library(OmnipathR)
library(Seurat)
library(dplyr)
library(pheatmap)
library(semla)
library(qs)
gc()
options(future.globals.maxSize = 1000 * 1024^2)
st_semla <- qread("qs_data/RCTD_results/RCTD_results_full_reference_subtype.qs")
sample_name <- c('CTB1NT','CTB1T','CTB3NT','CTB3T','CTB7NT','CTB7T')


#### Progeny ####
net <- decoupleR::get_progeny(organism = 'human', top =500)
head(net)
table(net$source)

SpatialPlot(st_semla, images = "CTB1NT") & 
  guides(fill = guide_legend(override.aes = list(size = 4)))

mat <- as.matrix(st_semla@assays[["Spatial"]]@layers[["data"]])
dim(mat)
rownames(mat) <- rownames(st_semla)
colnames(mat) <- colnames(st_semla)
head(rownames(mat)); head(colnames(mat))

plan("multisession",workers = 10)
acts <- decoupleR::run_mlm(
  mat = mat,
  net = net,
  .source = "source",  
  .target = "target",  
  .mor = "weight",    
  minsize = 5         
)
head(acts)
table(acts$source)

shared_genes <- intersect(rownames(mat), net$target)

acts_pivot <- acts %>%
  tidyr::pivot_wider(
    id_cols = "source",
    names_from = "condition",
    values_from = "score"
) 
head(acts_pivot)

st_semla[["pathwaysmlm"]] <- acts %>%
  tidyr::pivot_wider(
    id_cols = "source",
    names_from = "condition",
    values_from = "score"
  ) %>%
  tibble::column_to_rownames(var = "source") %>%
  Seurat::CreateAssayObject()
Seurat::DefaultAssay(st_semla) <- "pathwaysmlm"
st_semla <- Seurat::ScaleData(st_semla)
st_semla@assays$pathwaysmlm@data <- st_semla@assays$pathwaysmlm@scale.data
rownames(st_semla)
qsave(st_semla, "qs_data/RCTD_results/RCTD_results_full_reference_subtype_progeny.qs")


#### pathways visualization ####
rm(list = ls())
gc()
library(RColorBrewer)
library(ggrastr)
st_semla <- qread("qs_data/RCTD_results/RCTD_results_full_reference_subtype_progeny.qs")
pathways <- dput(rownames(st_semla))
sample_name <- c('CTB1_BL','CTB1_T','CTB3_BL','CTB3_T','CTB7_BL','CTB7_T')
cols <- RColorBrewer::brewer.pal(9,"RdBu") %>% rev()
display.brewer.all() 
brewer.pal.info
Seurat::DefaultAssay(st_semla) <- "pathwaysmlm"
i = 1
MapFeatures(st_semla, 
            features = pathways[c(4)], 
            section_number = i, 
            arrange_features = "row",
            colors = cols,
            ncol = 1,
            override_plot_dims = TRUE,
            pt_size =1)&ggtitle(sample_name[i])
dev.off()

p <- SpatialFeaturePlot(st_semla, images = "CTB1NT", 
                        # max.cutoff = 'q99',
                        features = pathways[c(10)])& 
  scale_y_reverse()&
  coord_fixed()&
  theme_void()&
  DarkTheme()&
  scale_fill_viridis_c(option = "B")&
  ggtitle(sample_name[2])
rasterise(p,dpi = 300)
dev.off()

p <- SpatialDimPlot(st_semla, images = "CTB3NT", 
                    group.by = "first_type")& 
  scale_y_reverse()&
  coord_fixed()&
  theme_void()&
  DarkTheme()&
  NoLegend()&
  scale_fill_manual(values = cluster_colors)&
  guides(fill = guide_legend(override.aes = list(size = 6)))&
  ggtitle(sample_name[3])
rasterise(p,dpi = 300)
dev.off()

p <- SpatialFeaturePlot(st_semla, images = "CTB1NT", 
                        # max.cutoff = 'q99',
                        features = colnames(st_semla@meta.data)[55:56])& 
  scale_y_reverse()&
  coord_fixed()&
  NoLegend()&
  theme_void()&
  DarkTheme()&
  # scale_fill_gradientn(colours=viridis_plasma_light_high, 
  #                      na.value = "black", 
  #                      limits=c(0,3))&
  ggtitle(sample_name[1])
rasterise(p,dpi = 300)
dev.off()


Seurat::DefaultAssay(st_semla) <- "Spatial"
library(ggrastr)
cols <- RColorBrewer::brewer.pal(9,"Spectral") %>% rev()
p <- SpatialFeaturePlot(st_semla, images = "CTB1T", 
                        # max.cutoff = 'q99',
                        features = "POSTN+ FB")& 
  scale_y_reverse()&
  coord_fixed()&
  theme_void()&
  DarkTheme()&
  scale_fill_gradientn(colours=cols, 
                       na.value = "black")&
  ggtitle(sample_name[2])
rasterise(p,dpi = 300)
dev.off()

MapFeatures(st_semla, 
            features = "SPP1", 
            section_number = i, 
            arrange_features = "row",
            colors = viridis::rocket(11, direction = -1),
            ncol = 1,
            override_plot_dims = TRUE,
            pt_size =1)&ggtitle(sample_name[i])
dev.off()


MapFeatures(st_semla, 
            features = "POSTN+ FB", 
            section_number = i, 
            arrange_features = "row",
            colors = viridis::mako(11,direction = -1),
            ncol = 1,
            override_plot_dims = TRUE,
            pt_size =1)&ggtitle(sample_name[i])

pathway_df <- t(as.matrix(st_semla@assays$pathwaysmlm@data)) %>% 
  as.data.frame() %>% 
  dplyr::mutate(cluster = Seurat::Idents(st_semla)) %>% 
  tidyr::pivot_longer(
    cols = -cluster,  
    names_to = "source",  
    values_to = "score"
  ) %>% 
  dplyr::group_by(cluster, source) %>% 
  dplyr::summarise(mean = mean(score), .groups = "drop")

head(pathway_df)
588/42

pathway_matrix <- pathway_df %>% 
  tidyr::pivot_wider(
    id_cols = 'cluster',  
    names_from = 'source',  
    values_from = 'mean'
  ) %>% 
  tibble::column_to_rownames(var = 'cluster') %>% 
  as.matrix()

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
color_palette <- grDevices::colorRampPalette(colors = colors)(100)

breaks <- c(
  seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
  seq(0.05, 1.25, length.out = floor(100 / 2))
)

pheatmap::pheatmap(
  mat = pathway_matrix,
  color = color_palette,
  border_color = "white",
  breaks = breaks,
  cellwidth = 12,
  cellheight = 12,
  treeheight_row = 20,
  treeheight_col = 20
)

table(st_semla@meta.data[["sample_id"]])
st_semla$condition <- ifelse(grepl("NT",st_semla@meta.data[["sample_id"]]),"CTB_BL","CTB_T")
table(st_semla@meta.data[["condition"]],st_semla@meta.data[["sample_id"]])
pathway_df <- t(as.matrix(st_semla@assays$pathwaysmlm@data)) %>% 
  as.data.frame() %>% 
  dplyr::mutate(
    cluster = Seurat::Idents(st_semla), 
    condition = st_semla$condition       
  ) %>% 
  tidyr::pivot_longer(
    cols = -c(cluster, condition),      
    names_to = "source",  
    values_to = "score"
  ) %>% 
  dplyr::group_by(condition, cluster, source) %>%  
  dplyr::summarise(mean = mean(score), .groups = "drop")  
head(pathway_df)
gc()

pathway_matrix_BL <- pathway_df %>% 
  filter(condition == "CTB_BL") %>% 
  tidyr::pivot_wider(
    id_cols = "cluster",
    names_from = "source",
    values_from = "mean"
  ) %>% 
  tibble::column_to_rownames("cluster") %>% 
  as.matrix()

pathway_matrix_T <- pathway_df %>% 
  filter(condition == "CTB_T") %>% 
  tidyr::pivot_wider(
    id_cols = "cluster",
    names_from = "source",
    values_from = "mean"
  ) %>% 
  tibble::column_to_rownames("cluster") %>% 
  as.matrix()

pathway_matrix_combined <- cbind(pathway_matrix_BL, pathway_matrix_T)

colnames(pathway_matrix_combined) <- c(
  paste0(colnames(pathway_matrix_BL), "_CTB_BL"),
  paste0(colnames(pathway_matrix_T), "_CTB_T")
)

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))
color_palette <- colorRampPalette(colors)(100)
breaks <- c(seq(-1.25, 0, length.out = ceiling(100/2) + 1),
            seq(0.05, 1.25, length.out = floor(100/2)))
annotation_col <- data.frame(
  Condition = rep(c("CTB_BL", "CTB_T"), each = ncol(pathway_matrix_BL))
)
rownames(annotation_col) <- colnames(pathway_matrix_combined)

pheatmap::pheatmap(
  mat = pathway_matrix_combined,
  cluster_cols = FALSE,
  color = color_palette,
  breaks = breaks,
  border_color = "white",
  annotation_col = annotation_col,  
  annotation_colors = list(Condition = c(CTB_BL = "#F47252", CTB_T = "#00A0E9")),
  cellwidth = 10,
  cellheight = 10,
  treeheight_row = 10,
  treeheight_col = 10,
  main = "Pathway Activity by Condition and Subtype"
)



subtype_cols <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#FFE119", "#F58231", "#911EB4", "#46F0F0", "#F032E6",
  "#BCBD22", "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#AAFFC3", "#808000", "#FFD8B1",
  "#000075", "#A9A9A9", "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
  "#E377C2", "#7F7F7F", "#17BECF", "#FF9898", "#C5B0D5", "#C49C94", "#F7B6D2", "#C7C7C7",
  "#DBDB8D", "#9EDAE5", "#FFBB78", "#98DF8A", "#AEC7E8", "#FF9E9E", "#B5D0D9", "#D4B9D6",
  "#FDB462", "#80B1D3"
)
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
cluster_colors <- setNames(subtype_cols, desired_order)

existing_clusters <- rownames(pathway_matrix_combined)
annotation_row <- data.frame(
  Subtype = existing_clusters
)
rownames(annotation_row) <- existing_clusters

pheatmap::pheatmap(
  mat = pathway_matrix_combined,
  cluster_cols = FALSE,      
  cluster_rows = TRUE,      
  color = color_palette,
  breaks = breaks,
  border_color = "white",
  annotation_col = annotation_col,    
  annotation_row = annotation_row,     
  annotation_colors = list(
    Condition = c(CTB_BL = "#F47252", CTB_T = "#00A0E9"),
    Subtype = cluster_colors[existing_clusters]  
  ),
  cellwidth = 12,
  cellheight = 12,
  treeheight_row = 20,
  treeheight_col = 10,
  main = "Pathway Activity by Condition and Subtype"
)


library(ComplexHeatmap)

col_anno <- HeatmapAnnotation(
  Condition = annotation_col$Condition,
  col = list(Condition = c(CTB_BL = "#F47252", CTB_T = "#00A0E9"))
)
row_anno <- rowAnnotation(
  Subtype = annotation_row$Subtype,
  col = list(Subtype = cluster_colors),
  show_legend = FALSE
)

mat_width <- ncol(pathway_matrix_combined) * unit(4, "mm")
mat_height <- nrow(pathway_matrix_combined) * unit(4, "mm")

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
color_palette <- colorRampPalette(colors)(101)
breaks <- c(seq(-1.25, 0, length.out = ceiling(100/2) + 1),
            seq(0.05, 1.25, length.out = floor(100/2)))
Heatmap(
  pathway_matrix_combined,
  name = "Score", 
  col = circlize::colorRamp2(breaks, color_palette),
  rect_gp = gpar(col = "white", lwd = 0.5), 
  border_gp = gpar(col = "black", lwd = 1),   
  width = mat_width,
  height = mat_height,
  cluster_rows = TRUE,    
  cluster_columns = FALSE, 
  show_column_dend = FALSE, 
  top_annotation = col_anno,
  left_annotation = row_anno,
  row_names_side = "right",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 10),  
  column_names_gp = gpar(fontsize = 10), 
  row_dend_width = unit(15, "mm")  
)
dev.off()
