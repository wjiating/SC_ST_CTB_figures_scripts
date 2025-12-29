rm(list = ls())
gc()
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(dplyr)
library(ggraph)
library(igraph)
library(ggforce)
library(tidyverse)
library(ggplot2)
setwd("/home/data/gz0436/OneDrive_4_2025-3-8/Pyscenic")
source("R/compute_module_score.R")

# seu <- qs::qread("~/OneDrive_4_2025-3-8/data/ctb_scobj_v4.qs")
# Idents(seu) <- "disease"
# table(Idents(seu))
# seu <- subset(seu, idents = c('CTB_BL',"CTB_T"))
# table(seu$celltype, seu$disease)
# seu <- subset(seu, celltype %in% c("Fibroblast", "Myeloid", "T cell"))
# seu$celltype <- seu$subtype
# seu$group <- seu$orig.ident2
# qsave(seu, file = "~/OneDrive_4_2025-3-8/data/ctb_my_t_fb.qs")


seu <- qs::qread("~/OneDrive_4_2025-3-8/data/ctb_my_t_fb.qs")
DimPlot(seu, group.by = "subtype")&ggsci::scale_color_igv()
seu$celltype <- seu$subtype
seu$group <- seu$orig.ident2
## regulon (gene list)
regulons <- clusterProfiler::read.gmt("output/02-ctb_my_t_fb.regulons.gmt")
## data.frame -> list, list -- gene set
rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
# regulon.list <- split(regulons$gene, regulons$term)
# names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
names(regulon.list) <- rg.names

summary(sapply(regulon.list, length))
print(regulon.list[1])
saveRDS(regulon.list, "output/03-1.ctb_my_t_fb.regulons.rds")

regulon.list <- readRDS("output/03-1.ctb_my_t_fb.regulons.rds")
regulon.list[["MITF(58g)"]]
regulon.list[["MTF1(34g)"]]


## AUCell RAS matrix
## RAS = regulon activity score
seu <- ComputeModuleScore(seu, gene.sets = regulon.list, min.size = 10, cores = 10)
seu
rg.names
DefaultAssay(seu) <- "AUCell"

p1 <- FeaturePlot(seu, features = "MITF(58g)", split.by = "disease")
p2 <- FeaturePlot(seu, features = "MITF", split.by = "disease")
(p1 / p2) & scale_color_viridis_c()

p1 <- FeaturePlot(seu, features = "MTF1(34g)", split.by = "disease")
p2 <- FeaturePlot(seu, features = "MTF1", split.by = "disease")
(p1 / p2) & scale_color_viridis_c()

p1 <- FeaturePlot(seu, features = "ZNF226(15g)", split.by = "disease")
p2 <- FeaturePlot(seu, features = "ZNF226", split.by = "disease")
(p1 / p2) & scale_color_viridis_b()



seu <- qs::qread("~/OneDrive_4_2025-3-8/data/ctb_my_t_fb.qs")


## regulon module
rasMat <- seu[["AUCell"]]@data
dim(rasMat)
rasMat <- t(rasMat)
pccMat <- cor(rasMat) 

# Connection Specificity Index (CSI)
CSI <- function(r1, r2) {
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat), c(r1,r2))
  N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}

csiMat <- pbapply::pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)

## h for cut tree
library(dendextend)
library(ggsci)
h = 8
row_dend = as.dendrogram(hclust(dist(pccMat), method = "complete"))
clusters <- dendextend::cutree(row_dend, h = h) # dendextend::cutree()
table(clusters)
sum(table(clusters))
length(clusters)
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))
plot(row_dend)

## visualisation
library(ComplexHeatmap)
library(circlize)

col_range = c(0, 1)
col_fun <- colorRamp2(col_range, c("#FCF8DE", "#253177"))

for(i in 1:nrow(pccMat)) {
  pccMat[i, i] <- 0
}

ht <- Heatmap(
  matrix = pccMat,
  col = col_fun,
  name = "ht1",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  use_raster = TRUE
)

lgd <- Legend(
  col_fun = col_fun,
  title = "",
  at = col_range,
  labels = c("low", "high"),
  direction = "horizontal",
  legend_width = unit(1, "in"),
  border = FALSE
)

{
  draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))
  decorate_heatmap_body("ht1", {
    tree = column_dend(ht)
    ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    clusters <- names(table(ind))
    x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
    x2 = sapply(clusters, function(x) last_index(ind == x))
    x1 = x1/length(ind)
    x2 = x2/length(ind)
    grid.rect(x = x1, width = (x2 - x1), y = 1-x1, height = (x1 - x2),
              hjust = 0, vjust = 0, default.units = "npc",
              gp = gpar(fill=NA, col="#FCB800", lwd=2))
    grid.text(label = paste0("M",clusters),
              x = x2-length(clusters)/length(ind), y = 1-x1-(x2-x1)/2,
              default.units = "npc",
              hjust = 1, vjust = 0.5,
              gp = gpar(fontsize=12, fontface="bold"))
  })
  decorate_column_dend("ht1", {
    tree = column_dend(ht)
    ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    clusters <- names(table(ind))
    x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
    x2 = sapply(clusters, function(x) last_index(ind == x))
    grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
              default.units = "npc", gp = gpar(fill = pal_d3("category20")(20), alpha=.5, col = NA))
  })
}
dev.off()

tree <- column_dend(ht)  # return a dendrogram object
row_dend = color_branches(tree, h = 8, col = pal_d3("category20")(20))
plot(row_dend)
dev.off()

ind <- cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]

clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon=names(ind), cluster=paste0("M",ind))
table(regulon.clusters$cluster)
regulon.clusters
saveRDS(regulon.clusters, "output/03-2.ctb_my_t_fb.regulon_modules.rds")

regulon.clusters <- readRDS("output/03-2.ctb_my_t_fb.regulon_modules.rds")
table(regulon.clusters$cluster)
regulon.clusters %>% 
  dplyr::filter(cluster=="M1")

plot_module_network <- function(module_name, 
                                regulon_list, 
                                module_color = "#2E86AB", 
                                regulon_color = "#A8DADC",
                                module_size = 10, 
                                regulon_size = 6,
                                title = NULL,
                                bg_color = "#F9F2F2") {
  
  if(is.null(title)) title <- paste0(module_name, " Module - Regulons Network")
  edges <- data.frame(from = rep(module_name, length(regulon_list)), to = regulon_list)
  nodes <- data.frame(
    id = c(module_name, regulon_list),
    type = c("Module", rep("Regulon", length(regulon_list))),
    label = c(module_name, regulon_list)
  )
  graph <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  layout <- igraph::layout_as_star(graph, center = module_name)
  radius <- max(sqrt(layout[,1]^2 + layout[,2]^2)) * 1.2
  ggraph::ggraph(graph, layout = layout) +
    ggforce::geom_circle(
      aes(x0 = 0, y0 = 0, r = radius),
      fill = bg_color,
      color = NA,
      alpha = 0.3
    ) +
    ggraph::geom_edge_link(
      edge_colour = "black", 
      edge_width = 0.5, 
      alpha = 0.6
    ) +
    ggraph::geom_node_point(aes(color = type, size = type), alpha = 0.8) +
    ggraph::geom_node_text(aes(label = label), size = 2.5, repel = TRUE, max.overlaps = 20) +
    ggplot2::scale_color_manual(values = c("Module" = module_color, "Regulon" = regulon_color)) +
    ggplot2::scale_size_manual(values = c("Module" = module_size, "Regulon" = regulon_size)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::labs(title = title, color = "Node Type", size = "Node Type")
}
m1_regulons <- regulon.clusters %>% 
  dplyr::filter(cluster == "M1") %>% 
  dplyr::pull(regulon)
plot_module_network("M1", m1_regulons)


# regulon-module mean activity
k = 8 # length(clusters)
cell.info <- seu@meta.data
moduleRasMat <- lapply(paste0("M",1:k), function(x){
  regulon.use <- subset(regulon.clusters, cluster == x)$regulon
  rowMeans(rasMat[, regulon.use, drop=FALSE])
})
names(moduleRasMat) <- paste0("M",1:k)
moduleRasMat <- do.call(cbind, moduleRasMat)
dim(moduleRasMat)
cell.info <- cbind(cell.info, moduleRasMat[rownames(cell.info), ])

cell.info <- cbind(cell.info, FetchData(seu, vars = paste0("umap_", 1:2)))

p.list <- lapply(paste0("M",1:k), function(module){
  data.use <- cell.info
  # expression.color <- c("darkblue", "lightblue", "green", "yellow", "red")
  expression.color <- c('#5749a0', '#0f7ab0', '#00bbb1',
                        '#bef0b0', '#fdf4af', '#f9b64b',
                        '#ec840e', '#ca443d', '#a51a49')
  max.val <- quantile(data.use[, module], 0.99)
  low.val <- quantile(data.use[, module], 0.1)
  data.use[, module] <- ifelse(data.use[, module] > max.val, max.val, data.use[, module])
  ggplot(data.use, aes(umap_1, umap_2, color=get(module))) +
    geom_point(size=0.05) +
    theme_bw(base_size = 15) +
    ggtitle(module) +
    facet_wrap(~disease) +
    scale_color_gradientn(name = NULL, 
                          colors = expression.color#,
                          # guide = guide_colorbar(frame.colour = "black", 
                          #                       ticks.colour = NA),
                          # labels = c('low',"high"),
                          # breaks = c(low.val,max.val)
                          ) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold", size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          panel.border = element_rect(fill=NA,
                                      color="black", 
                                      size=1, 
                                      linetype="solid"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) 
})

library(ggrastr)
p <- cowplot::plot_grid(plotlist = p.list[1:3], ncol = 3)
rasterise(p,dpi=300)
dev.off()
p <- cowplot::plot_grid(plotlist = p.list[4:6], ncol = 3)
rasterise(p,dpi=300)
dev.off()
p <- cowplot::plot_grid(plotlist = p.list[7:8], ncol = 3)
rasterise(p,dpi=300)
dev.off()
cowplot::plot_grid(plotlist = p.list[10:12], ncol = 3)
dev.off()
cowplot::plot_grid(plotlist = p.list[13:15], ncol = 3)
dev.off()


df <- cell.info[,c("subtype","disease",paste0("M",seq(1,8,1)))]
table(df$subtype, df$disease)
head(df)

df <- df %>% 
  tidyr::pivot_longer(
    cols = -c("subtype","disease"),  
    names_to = "source",  
    values_to = "score"
) 
table(df$source)

library(ggpubr)
library(rstatix)
plot_module_comparison_simple <- function(df, 
                                          module, 
                                          color_palette = NULL, 
                                          figsize = c(14, 6)) {
  
  module_data <- df %>%
    dplyr::filter(source == module) %>%
    dplyr::filter(!is.na(score))
  
  stat_test <- module_data %>%
    group_by(subtype) %>%
    t_test(score ~ disease) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "subtype", dodge = 0.8)
  
  if (is.null(color_palette)) {
    color_palette <- c("CTB_BL" = "#ff7f0e", "CTB_T" = "#1f77b4")
  }
  
  p <- ggplot(module_data, aes(x = subtype, y = score, fill = disease)) +
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE, 
                color = NA, linewidth = 0, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA, 
                 position = position_dodge(0.9),
                 linewidth = 0.3) +
    scale_fill_manual(values = color_palette) +
    labs(
      title = paste0(module, " - Module Activity"),
      x = "Cell Types",
      y = "Module Activity Score"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color="black"),
      axis.text.y = element_text(size = 10, color="black"),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  if (nrow(stat_test) > 0) {
    p <- p + 
      stat_pvalue_manual(
        stat_test, 
        label = "p.adj.signif",
        tip.length = 0.01,
        bracket.shorten = 0.05,
        hide.ns = FALSE,
        inherit.aes = FALSE
      )
  }
  
  return(list(plot = p, statistics = stat_test))
}


generate_all_module_plots_simple <- function(df, modules = paste0("M", 1:8), 
                                             color_palette = NULL, 
                                             figsize = c(14, 6)) {
  
  plot_list <- list()
  stat_list <- list()
  
  for (module in modules) {
    if (module %in% unique(df$source)) {
      result <- plot_module_comparison_simple(df, module, color_palette, figsize)
      plot_list[[module]] <- result$plot
      stat_list[[module]] <- result$statistics
      
      cat("\n=== Statistics for", module, "===\n")
      print(result$statistics)
    } else {
      warning(paste("Module", module, "not found in dataframe"))
    }
  }
  
  return(list(plots = plot_list, statistics = stat_list))
}

results <- generate_all_module_plots_simple(df)

print(results$plots$M1)
print(results$plots$M2)

for (module_name in names(results$plots)) {
  filename <- paste0("figs/","module_comparison_", module_name, ".pdf")
  ggsave(filename, plot = results$plots[[module_name]], 
         width = 16, height = 6, dpi = 300)
}


## TF module: a list of regulons
## regulon:   a list of genes
## genes?? enrichment !!

regulon.clusters <- readRDS("output/03-2.ctb_my_t_fb.regulon_modules.rds")

rms <- unique(regulon.clusters$cluster)
module.list <- lapply(rms, function(xx) subset(regulon.clusters, cluster == xx)$regulon)
names(module.list) <- rms

regulon.list <- readRDS("output/03-1.ctb_my_t_fb.regulons.rds")

regulon.module.genes <- lapply(seq_along(module.list), function(xx) {
  unique(unlist(regulon.list[module.list[[xx]] ]))
})
names(regulon.module.genes) <- names(module.list)

sapply(regulon.module.genes, length)

## enrichment
library(clusterProfiler)
library(msigdbr)

hs_go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
head(hs_go_bp)
go_bp_terms <- split(hs_go_bp$gene_symbol, hs_go_bp$gs_name)
hm <- clusterProfiler::read.gmt("resource/h.all.v2022.1.Hs.symbols.gmt")
gb <- clusterProfiler::read.gmt("resource/c5.go.bp.v2023.2.Hs.symbols.gmt")

# t2g <- readRDS("resource/msigdb.all.human_symbol.rds")

eres.list <- lapply(sort(names(regulon.module.genes)), function(clu) {
  message(glue::glue("processing {clu} ..."))
  genes <- regulon.module.genes[[clu]]
  ## downsample to speed up the calculation
  # if (length(genes) >= 200) {
  #   genes <- sample(genes, size = 200)
  # }
  clusterProfiler::enricher(gene = genes, 
                            TERM2GENE = hm,  # gb
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.05)
})
names(eres.list) <- sort(names(regulon.module.genes))

## plot data
data.use <- lapply(names(eres.list), function(xx) {
  res <- subset(eres.list[[xx]]@result, grepl("^HALLMARK_", ID)) # only show HALLMARK results
  # res <- subset(eres.list[[xx]]@result, grepl("^GOBP", ID)) # only show HALLMARK results
  res <- head(res, 5)
  res$group <- xx
  res
}) %>% Reduce(rbind, .)

term.levels <- unique(data.use$Description)

data.use <- lapply(names(eres.list), function(xx) {
  res <- subset(eres.list[[xx]]@result, Description %in% term.levels)
  res$group <- xx
  res
}) %>% Reduce(rbind, .)

data.use$Description <- factor(data.use$Description, levels = rev(term.levels))
data.use <- subset(data.use, p.adjust < 0.05 & Count > 2)

## plot
ggplot(data.use, aes(group, Description)) +
  geom_point(aes(size = Count, fill = -log10(p.adjust)), shape=21, color = "black") +
  scale_fill_gradientn(colors = c("white","red")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "gene cluster") +
  theme_bw(base_size = 13) +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color = "black"))

qs::qsave(eres.list, "output/03-3.ctb_my_t_fb.regulon_modules.enriched_terms.qs")


## regulon-TF-Targets
plot_tf_network <- function(tf_name, target_genes, 
                            tf_color = "#6D597A",        
                            target_color = "#B56576",    
                            tf_size = 8, 
                            target_size = 4,
                            title = NULL,
                            bg_color = "#F8EDEB") {
  
  if(is.null(title)) title <- paste0(tf_name, "-mediated gene network")
  edges <- data.frame(from = rep(tf_name, length(target_genes)), to = target_genes)
  nodes <- data.frame(
    id = c(tf_name, target_genes),
    type = c("TF", rep("Target", length(target_genes))),
    label = c(tf_name, target_genes)
  )
  graph <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
  layout <- igraph::layout_as_star(graph, center = tf_name)
  
  radius <- max(sqrt(layout[,1]^2 + layout[,2]^2)) * 1.2
  
  ggraph::ggraph(graph, layout = layout) +
    ggforce::geom_circle(
      aes(x0 = 0, y0 = 0, r = radius),
      fill = bg_color,
      color = NA,
      alpha = 0.3
    ) +
    ggraph::geom_edge_link(
      edge_colour = "black", 
      edge_width = 0.5, 
      alpha = 0.6,
      arrow = ggplot2::arrow(length = unit(2, "mm"), type = "closed"),
      end_cap = ggraph::circle(2, "mm")  
    ) +
    ggraph::geom_node_point(aes(color = type, size = type), alpha = 0.8) +
    ggraph::geom_node_text(aes(label = label), size = 2.5, repel = TRUE, max.overlaps = 20) +
    ggplot2::scale_color_manual(values = c("TF" = tf_color, "Target" = target_color)) +
    ggplot2::scale_size_manual(values = c("TF" = tf_size, "Target" = target_size)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::labs(title = title, color = "Node Type", size = "Node Type")
}

plot_tf_network("MITF", regulon.list[["MITF(58g)"]])
plot_tf_network("MTF1", regulon.list[["MTF1(34g)"]])

## enrichment
library(clusterProfiler)
library(msigdbr)

hm <- clusterProfiler::read.gmt("resource/h.all.v2022.1.Hs.symbols.gmt")
gb <- clusterProfiler::read.gmt("resource/c5.go.bp.v2023.2.Hs.symbols.gmt")

# t2g <- readRDS("resource/msigdb.all.human_symbol.rds")


regulon.list[["MTF1(34g)"]]
eres.list <- lapply("MTF1(34g)", function(clu) {
  message(glue::glue("processing {clu} ..."))
  genes <- regulon.list[[clu]]
  ## downsample to speed up the calculation
  # if (length(genes) >= 200) {
  #   genes <- sample(genes, size = 200)
  # }
  clusterProfiler::enricher(gene = genes, 
                            TERM2GENE = gb,  # gb
                            minGSSize = 5,
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.1)
})
names(eres.list) <- "MTF1"

## plot data
data.use <- lapply(names(eres.list), function(xx) {
  res <- subset(eres.list[[xx]]@result, grepl("^HALLMARK_", ID)) # only show HALLMARK results
  # res <- subset(eres.list[[xx]]@result, grepl("^GOBP", ID)) # only show HALLMARK results
  res <- head(res, 5)
  res$group <- xx
  res
}) %>% Reduce(rbind, .)

term.levels <- unique(data.use$Description)

data.use <- lapply(names(eres.list), function(xx) {
  res <- subset(eres.list[[xx]]@result, Description %in% term.levels)
  res$group <- xx
  res
}) %>% Reduce(rbind, .)

data.use$Description <- factor(data.use$Description, levels = rev(term.levels))
data.use <- subset(data.use, pvalue < 0.05 & Count > 2)

## plot
data.use <- data.use[order(data.use$pvalue, decreasing = FALSE), ]

ggplot(data.use, aes(x = Count, y = reorder(Description, Count))) +
  geom_col(fill = "#4A90A4", width = 0.7, alpha = 0.8) +  
  labs(
    x = "Gene Count",
    y = "Pathway",
    title = "MTF1 Pathway Enrichment Analysis"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )



GetTopRegulons <- function(rssMat, topn = 10) {
  top_regulons_list <- list()
  for (cell.type in colnames(rssMat)) {
    rss_values <- rssMat[, cell.type]
    sorted_regulons <- names(sort(rss_values, decreasing = TRUE))
    top_regulons <- head(sorted_regulons, n = topn)
    top_regulons_list[[cell.type]] <- top_regulons
  }
  
  return(top_regulons_list)
}

top_regulons_list <- GetTopRegulons(rssMat)

GetTopRegulons <- function(rssMat, topn = 5) {
  top_regulons_list <- list()
  for (cell.type in colnames(rssMat)) {
    rss_values <- rssMat[, cell.type]
    sorted_regulons <- names(sort(rss_values, decreasing = TRUE))
    top_regulons <- head(sorted_regulons, n = topn)
    top_regulons_list[[cell.type]] <- top_regulons
  }
  return(top_regulons_list)
}

top_regulons_list <- GetTopRegulons(rssMat)
rasMat <- seu[["AUCell"]]@data
dim(rasMat)
rasMat <- t(rasMat)

ct_mat <- rasMat[, top_regulons_list[["SPP1 Mac"]]]
data <- FetchData(seu, vars = c("celltype","disease"))
identical(rownames(rasMat), rownames(data))

df <- cbind(ct_mat, data)
head(df)


df$is_SPP1_Mac <- ifelse(df$celltype == "SPP1 Mac", 1, 0)
regulon_names <- c("ZNF226(15g)", "NR1H3(34g)", "MTF1(34g)", "IRF5(375g)", 
                   "MITF(58g)")

colnames(df)

cor_results <- data.frame(
  Regulon = character(),
  Correlation = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for (regulon in regulon_names) {
  cor_test <- cor.test(df[[regulon]], df$is_SPP1_Mac, method = "pearson")
  cor_results <- rbind(cor_results, data.frame(
    Regulon = regulon,
    Correlation = cor_test$estimate,
    P_value = cor_test$p.value
  ))
}
cor_results <- cor_results[order(-cor_results$Correlation), ]

print(cor_results)


plot_data <- df %>%
  select(all_of(regulon_names), celltype) %>%
  mutate(Group = ifelse(celltype == "SPP1 Mac", "SPP1 Mac", "Other Cells")) %>%
  select(-celltype) %>%
  tidyr::pivot_longer(cols = all_of(regulon_names), 
                      names_to = "Regulon", 
                      values_to = "Activity")
regulon_order <- cor_results$Regulon[order(cor_results$Correlation)]
plot_data$Regulon <- factor(plot_data$Regulon, levels = regulon_order)

p2 <- ggplot(plot_data, aes(x = Regulon, y = Activity, fill = Group)) +
  geom_violin(alpha = 0.7, scale = "width", width = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.size = 0.5, 
               position = position_dodge(0.8)) +
  scale_fill_manual(values = c("SPP1 Mac" = "#E64B35", "Other Cells" = "#4DBBD5")) +
  labs(
    title = "Regulon Activity Distribution",
    subtitle = "Comparison between SPP1 Mac and Other Cell Types",
    x = "Regulon",
    y = "Regulon Activity Score",
    fill = "Cell Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  ) +
  coord_flip()
print(p2)
library(ggrastr)
library(rstatix)
library(ggpubr)
rasterise(p2, dpi=300)

summary_data <- plot_data %>%
  group_by(Regulon, Group) %>%
  summarise(Mean_Activity = mean(Activity, na.rm = TRUE),
            .groups = 'drop')

p3 <- ggplot(cor_results, aes(x = reorder(Regulon, Correlation), y = Correlation, 
                              color = Correlation, size = abs(Correlation))) +
  geom_point(alpha = 0.8) +
  geom_segment(aes(xend = reorder(Regulon, Correlation), yend = 0), 
               linewidth = 1, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                        midpoint = 0, name = "Correlation") +
  scale_size_continuous(range = c(3, 8), name = "|Correlation|") +
  geom_text(aes(label = sprintf("r = %.3f\n%s", Correlation, Significance)), 
            nudge_y = ifelse(cor_results$Correlation > 0, 0.03, -0.03),
            size = 3.5) +
  labs(
    title = "Regulon-SPP1 Mac Correlation Dot Plot",
    x = "Regulon",
    y = "Correlation Coefficient"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold")
  ) +
  coord_flip()
print(p3)


