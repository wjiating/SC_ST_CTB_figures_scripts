rm(list = ls())
gc()
library(patchwork)
library(SPOTlight)
library(semla)
library(Seurat)
library(qs)
library(singlet)
library(RcppML)
packageVersion("Seurat")
packageVersion('spatstat.utils')
# install.packages("spatstat.utils")

se <- qread("./semla_output/se.qs")
name_list = list('CTB1_BL','CTB1_T','CTB3_BL','CTB3_T','CTB7_BL','CTB7_T')
se_var <- se[VariableFeatures(se), ]

# Set seed for reproducibility
set.seed(45)
table(se_var$sample_id)
se_var <- RunNMF(se_var, split.by = "sample_id")
k <- ncol(se_var@reductions$nmf@feature.loadings)
k
Idents(se_var) <- "sample_id"
DimPlot(se_var, dims = c(1, 2), reduction = "nmf", group.by = "sample_id",label = TRUE, pt.size = 0.5)&scale_color_startrek()
DimPlot(se_var, dims = c(2, 3), reduction = "nmf", group.by = "sample_id",label = TRUE, pt.size = 0.5)&scale_color_startrek()
DimPlot(se_var, dims = c(1, 3), reduction = "nmf", group.by = "sample_id",label = TRUE, pt.size = 0.5)&scale_color_startrek()

library(ggplot2)
library(ggsci)
library(readxl)
library(plotly)
nmf <- Embeddings(se_var, reduction = "nmf")
dim(nmf)
colnames(nmf)
nmf <- as.data.frame(nmf)
nmf$sample <- se_var$sample_id

ggplot(data = nmf, 
       aes(x = NMF_1, y = NMF_2, fill = sample)) +
  geom_point(
    color = "black", alpha = 0.7, 
    size = 1, shape = 21) +
  theme_light() +      
  scale_fill_startrek() +
  labs(x = "NMF 1",y = "NMF 2",fill = "Sample") +
  theme(aspect.ratio = 1,legend.position = "right",
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, color = "gray30", linewidth = 0.6))

library(plotly)
library(ggthemes)  # 用于pal_startrek()配色

# 3D NMF可视化
plotly::plot_ly(data = nmf,
                x = ~NMF_1,
                y = ~NMF_2, 
                z = ~NMF_3,
                marker = list(size = 1),
                color = ~sample,
                colors = pal_startrek()(7)) %>% 
  layout(title = '3D NMF',
         scene = list(xaxis = list(title = 'NMF_1'),
                      yaxis = list(title = 'NMF_2'),
                      zaxis = list(title = 'NMF_3')))

RankPlot(se_var)

qsave(se_var, file = "NMF_output/se_var.qs")
se_var <- qread(file = "NMF_output/se_var.qs")

ImageDimPlot(se_var,fov = "slice",group.by = 'celltype',dark.background = T, crop = T)+ggtitle(name_list[[i]])

MapFeatures(se_var, 
            section_number = 5,
            features = paste0("NMF_", 1:3), 
            override_plot_dims = TRUE, 
            ncol = 3,
            colors = viridis::magma(n = 11, direction = -1)) &
  theme(plot.title = element_blank())

factor_colors<- c('#ffe119','#3cb44b','#e6194b')
selected_factors <- c(1:3)

i=6
p <- MapMultipleFeatures(se_var,
                         section_number = i,
                         features = paste0("NMF_", selected_factors),
                         colors = factor_colors,
                         override_plot_dims = TRUE,
                         pt_size =1)& ggtitle(name_list[[i]]) &
  theme(legend.position = "right", legend.text = element_text(angle = 0)) &
  # guides(fill = guide_legend(override.aes = list(size = 6))) &
  theme(plot.background = element_rect(fill = "grey98"),
        # panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(colour = "grey98"),
        plot.subtitle = element_text(colour = "black"),
        legend.text = element_text(colour = "black"),
        legend.title = element_text(colour = "black"))
p
ggsave(paste0('figure_NMF/spatial_plot_crop_',name_list[i],'.jpeg'),height=8, width=7, dpi = 1000)


co_plist <- list()
for(i in seq_along(name_list)){
  p <- MapMultipleFeatures(se_var,
                           section_number = i,
                           features = paste0("NMF_", selected_factors),
                           colors = factor_colors,
                           override_plot_dims = TRUE,
                           pt_size =1)  & 
    ggtitle("") & theme(legend.position = "none")
  co_plist[[name_list[[i]]]] <- p
}
i = 6
co_plist[[i]]
ggsave(paste0('figure_NMF/spatial_plot_crop_',name_list[i],'.jpeg'),height=8, width=7, dpi = 1000)
wrap_plots(co_plist,ncol = 3,byrow = FALSE)


PlotFeatureLoadings(se_var,
                    dims = 1:3,
                    reduction ="nmf",
                    nfeatures = 30,
                    mode ="dotplot",
                    fill ="darkmagenta",
                    ncols = 3,
                    pt_size = 3)
PlotFeatureLoadings(se_var, 
                    dims = selected_factors, 
                    reduction = "nmf", 
                    nfeatures = 30, 
                    mode = "heatmap", 
                    gradient_colors = viridis::magma(n = 11, direction = -1))

# 功能富集分析
# 获取feature.loads
nmf_loadings<- se_var[["nmf"]]@feature.loadings

# 转换为长格式并按因子对数据进行分组
gene_loadings_sorted <- nmf_loadings |>
  as.data.frame() |>
  tibble::rownames_to_column(var ="gene") |>
  as_tibble() |>
  tidyr::pivot_longer(all_of(colnames(nmf_loadings)), names_to ="fctr", values_to ="loading") |>
  mutate(fctr = factor(fctr, colnames(nmf_loadings))) |>
  group_by(fctr) |>
  arrange(fctr, -loading)

# 每个nmf提取前100个基因
gene_loadings_sorted |>
  slice_head(n = 100)

# 获取 NMF_1-NMF_3的TOP100基因
gene_set_nmf <- gene_loadings_sorted |>
  group_by(fctr) |>
  slice_head(n =50)

