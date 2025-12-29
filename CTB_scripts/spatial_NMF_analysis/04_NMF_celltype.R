rm(list = ls())
gc()
library(patchwork)
library(SPOTlight)
library(semla)
library(Seurat)
library(qs)
library(singlet)
library(RcppML)


se_var <- qread(file = "NMF_output/se_var.qs")
name_list = list('CTB1NT','CTB1T','CTB3NT','CTB3T','CTB7NT','CTB7T')


factor_colors<- c('#ffe119','#3cb44b','#e6194b')
selected_factors <- c(1:3)

i=1
p <- MapMultipleFeatures(se_var,
                         section_number = i,
                         features = paste0("NMF_", selected_factors),
                         colors = factor_colors,
                         override_plot_dims = TRUE,
                         pt_size =1)  & ggtitle(name_list[[i]])
p

df <- FetchData(se_var, vars = paste0("NMF_", selected_factors))
nmf_embeddings <- Embeddings(se_var[["nmf"]])
nmf_embeddings <- as.data.frame(nmf_embeddings)
dim(nmf_embeddings)
head(df)[1:3]
celltype_df <- FetchData(se_var, vars = c("celltype","region"))

head(df)
head(celltype_df)
table(celltype_df$celltype)
dim(df); dim(celltype_df)


#==== NMF acitivity ====
merged_df <- merge(df, celltype_df, by = "row.names")
rownames(merged_df) <- merged_df$Row.names
merged_df$Row.names <- NULL


# across
library(dplyr)
mean_data <- merged_df %>%
  group_by(celltype) %>%
  summarise(across(starts_with("NMF"), mean, na.rm = TRUE)) %>%
  as.data.frame()

rownames(mean_data) <- mean_data$celltype
mean_matrix <- as.matrix(t(mean_data[, -1]))
head(mean_matrix); range(mean_matrix)
scaled_matrix <- t(scale(t(mean_matrix)))  
library(corrplot)
rev_col <- rev(corrplot::COL2('RdBu', 100))
corrplot(scaled_matrix,
         is.corr = FALSE,
         method = "color",
         col = rev_col,          
         tl.col = "black",
         addgrid.col = "gray",
         cl.pos = "b",
         addCoef.col = "black",  
         number.cex = 0.7,
         mar = c(0.5, 0.5, 1, 1),
         title = "Scaled NMF Activity by Cell Type",
         tl.cex = 0.8)
pheatmap::pheatmap(scaled_matrix,
                   cluster_rows = F,
                   cluster_cols = T,
                   color = RColorBrewer::brewer.pal(9,"Reds"))
pheatmap::pheatmap(
  scaled_matrix,
  cluster_rows = FALSE,  
  cluster_cols = TRUE,   
  color = RColorBrewer::brewer.pal(9, "Reds"),  
  main = "NMF Activity by Cell Type (Scaled)", 
  fontsize_row = 10,      
  fontsize_col = 10,     
  fontsize = 10,         
  angle_col = 90,        
  border_color = "gray"  
)

# mutate
mean_data <- merged_df %>%
  group_by(celltype) %>%
  mutate(NMF1 = mean(NMF_1)) %>% 
  mutate(NMF2 = mean(NMF_2)) %>% 
  mutate(NMF3 = mean(NMF_3)) %>% 
  ungroup() %>% 
  dplyr::select(4,6:8) %>% 
  distinct() %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames(var = "celltype")

scaled_matrix <- scale(mean_data)  
library(corrplot)
rev_col <- rev(corrplot::COL2('RdBu', 100))
corrplot(scaled_matrix,
         is.corr = FALSE,
         method = "color",
         col = rev_col,          
         tl.col = "black",
         addgrid.col = "gray",
         cl.pos = "b",
         addCoef.col = "black",  
         number.cex = 0.7,
         mar = c(0, 0, 1, 0),
         title = "Scaled NMF Activity by Cell Type",
         tl.cex = 0.8)
pheatmap::pheatmap(scaled_matrix,
                   cluster_rows = T,
                   cluster_cols = F,
                   main = "NMF Activity by Cell Type (Scaled)", 
                   color = RColorBrewer::brewer.pal(9,"Reds"))

#==== categroy ====
head(df)

dominant_topic <- apply(df, 1, function(x) {
  which.max(x)  
})
head(dominant_topic)

dominant_topic_label <- paste0("NMF_", dominant_topic)
head(dominant_topic_label)

se_var$NMF <- as.factor(dominant_topic_label)
table(se_var$NMF)

nmf_results <- as.data.frame(df) %>%
  tibble::rownames_to_column("CellID") %>%
  mutate(
    Dominant_Topic = names(df)[apply(df, 1, which.max)],
    Max_Score = apply(df, 1, max)
  ) %>%
  select(CellID, Dominant_Topic, Max_Score, everything())

head(nmf_results)
table(nmf_results$Dominant_Topic)

MapLabels(se_var, 
          section_number = 1,
          column_name = "NMF", 
          pt_alpha = 0.6, 
          colors = factor_colors,
          pt_size = 1) &
  guides(fill = guide_legend(override.aes = list(size = 4)))


#==== proportion ====
Cellratio <- prop.table(table(se_var@meta.data$NMF,se_var@meta.data$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('NMF','sample','Freq')
# Cellratio$region <- factor(Cellratio$region,levels = c('epidermis', 'granuloma region', 'unaffected dermis'))
sample_colors <- c('#ffe119','#3cb44b','#e6194b')
ggplot(data = Cellratio, aes(x = sample, y = Freq, fill = NMF)) +
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
  labs(x = "", y = "Ratio", title = "NMF Frequency by Sample in Spatial Data") +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
)

# Cell proportion in region
Idents(se_var) <- "NMF"
table(Idents(se_var))
se_gra <- subset(se_var, NMF =="NMF_3")
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
  labs(x = "", y = "Ratio", title = "Cell Type Frequency by Sample in NMF3") +
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
table(se_gra@meta.data$celltype,se_gra@meta.data$orig.ident)
table(se@meta.data$celltype,se@meta.data$orig.ident)



qsave(se_var, file = "NMF_output/se_var.qs")
