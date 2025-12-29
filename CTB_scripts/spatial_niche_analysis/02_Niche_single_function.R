


#### Calculate niches ----
Niche_Analysis <- function(object, 
                           image = "slice", 
                           celltype = "celltype",
                           neighbors.k = 20, 
                           niches.k = 6) {
  # 获取指定图像坐标
  coord_df <- object@images[[image]]@coordinates
  
  # 互换坐标 
  coord_df_swapped <- coord_df
  colnames(coord_df_swapped) <- c("image_y", "image_x") 
  coord_df_swapped$x <- coord_df$y  
  coord_df_swapped$y <- coord_df$x
  
  # 更新坐标到对象 
  object@images[[image]]@coordinates <- coord_df_swapped[, c("x", "y")]
  
  # 创建子集对象
  sce_ctb <- subset(object, cells = rownames(coord_df_swapped))
  
  # 获取互换后的坐标用于FOV
  fov_coord <- sce_ctb@images[[image]]@coordinates[,1:2]
  
  counts_new <- LayerData(sce_ctb, assay = "Spatial", layer = "counts")
  colData <- sce_ctb@meta.data
  
  # 创建Seurat对象
  SO <- CreateSeuratObject(
    counts = counts_new,
    assay = "Spatial",
    meta.data = colData 
  )
  
  # 创建FOV使用互换后的坐标
  fov <- SeuratObject::CreateFOV(
    coords = fov_coord,
    type = "centroids",
    assay = "Spatial"
  )
  
  SO@images$whole_image <- fov
  Idents(SO) <- celltype
  
  # 进行邻域分析
  SO <- BuildNicheAssay(
    object = SO,
    fov = "whole_image",
    group.by = celltype,
    assay = "niche",
    cluster.name = "niches",
    neighbors.k = neighbors.k,
    niches.k = niches.k
  )
  
  # 添加分析结果列和assay
  sce_ctb$niches <- SO$niches
  niche_assay <- SO[["niche"]]
  sce_ctb[["niche"]] <- niche_assay
  
  return(sce_ctb)
}

#### plot niches ----
Niche_split <- function(sample){
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
  patchwork::wrap_plots(plist, nrow = 3) + plot_annotation(title = names(sample@images))
}

#### celltype proportion in each niche ----
composition_niche <- function(object,
                              Niche = "niches",
                              Celltype = "celltype"){
  Cellratio <- prop.table(table(object@meta.data[[Niche]], object@meta.data[[Celltype]]), margin = 1)
  Cellratio <- as.data.frame(Cellratio)
  colnames(Cellratio) <- c('Niche','Celltype','Freq')
  Cellratio$Celltype <- factor(Cellratio$Celltype,levels = c('Keratinocyte','Fibroblast','Myeloid',
                                                       "Endothelial",'T cell','Lymphatic Endothelial',
                                                       'B cell','Smooth muscle','Melanocyte',
                                                       'Mast cell'))
  Cellratio$Niche <- factor(Cellratio$Niche, levels = as.character(sort(as.numeric(unique(Cellratio$Niche)))))
  sample_colors <- c( "#FF9F1C",  "#E71D36",  "#039BE5",  "#FF6B6B",  
                      "#00A878",  "#87CEEB",  "#6A4C93",  "#FFD166",  
                      "#2EC4B6",  "#F4511E")
  ggplot(data = Cellratio, aes(x = Niche, y = Freq, fill = Celltype)) +
    geom_bar(stat = "identity", width = 0.8, position = "fill") +
    scale_fill_manual(values = sample_colors) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 12, colour = "black"),
      legend.title = element_text(size = 12),  
      legend.text = element_text(size = 12),   
      plot.title = element_text(size = 12) 
    ) +
    labs(x = "Niche", y = "Ratio", title = paste0("Celltype Frequency by Spatial Niche in ", names(object@images))) 
}


