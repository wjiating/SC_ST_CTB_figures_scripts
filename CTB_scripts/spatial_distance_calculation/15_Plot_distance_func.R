
#### calculate_sample_distances #### 
calculate_sample_distances <- function(data, sample_name) {
  sample_data <- data %>% 
    dplyr::filter(Sample == sample_name)
  distances <- phenoptr::find_nearest_distance(sample_data)
  result <- dplyr::bind_cols(sample_data, distances)
  return(result)
}

#### named colors ####
subtype_cols <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#FFE119", "#F58231", "#911EB4", "#46F0F0", "#F032E6",
  "#BCBD22", "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#AAFFC3", "#808000", "#FFD8B1",
  "#000075", "#A9A9A9", "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
  "#E377C2", "#7F7F7F", "#17BECF", "#FF9898", "#C5B0D5", "#C49C94", "#F7B6D2", "#C7C7C7",
  "#DBDB8D", "#9EDAE5", "#FFBB78", "#98DF8A", "#AEC7E8", "#FF9E9E", "#B5D0D9", "#D4B9D6",
  "#FDB462", "#80B1D3"
)

phenotypes <- c(
  "Basal KC", "Inflammatory KC", "Proliferating KC", "Spinous KC",
  "Granular KC", "Hair follicle", "Sweet Gland", "Schwann Cell", "Melanocyte",
  "CCL19+ FB", "APCDD1+ FB", "PI16+ FB", "POSTN+ FB", "TNN+ FB", "RAMP1+ FB",
  "Arterial EC", "Capillary EC", "Venular EC1", "Venular EC2", "Lymphatic Endothelial",
  "Smooth muscle", "Mast cell", "Langerhans", "cDC1", "cDC2A", "cDC2B", "pDC",
  "M1-like Mac", "M2-like Mac", "SPP1 Mac", "TREM2 Mac", "CD8T", "Treg", "CD4 Naïve",
  "Th17", "CD8T_NK", "Th17_Treg", "NK", "naïve B", "activated B", "memory B", "Plasma"
)
named_colors <- setNames(subtype_cols, phenotypes)


#### plot cell distance (boxplot and jitter) ####
plot_cell_distance <- function(csd_data, target_cell_type, distance_col = NULL, metric = c("median", "mean")) {
  if (!target_cell_type %in% csd_data$Phenotype) {
    stop("Target cell type not found in Phenotype column")
  }
  metric <- match.arg(metric) 
  
  if (is.null(distance_col)) {
    distance_col <- paste("Distance to", target_cell_type)
    if (!distance_col %in% names(csd_data)) {
      stop(paste("Distance column", distance_col, "not found in data"))
    }
  }
  
  distance_metric <- csd_data %>%
    group_by(Phenotype) %>%
    summarise(
      Median_Distance = median(.data[[distance_col]], na.rm = TRUE),
      Mean_Distance = mean(.data[[distance_col]], na.rm = TRUE)
    ) %>%
    arrange(if (metric == "median") Median_Distance else Mean_Distance)
  
  csd_data$Phenotype <- factor(
    csd_data$Phenotype,
    levels = distance_metric$Phenotype
  )
  
  p <- ggplot(csd_data, aes(x = Phenotype, y = .data[[distance_col]], color = Phenotype)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
    geom_point(
      data = distance_metric,
      aes(y = if (metric == "median") Median_Distance else Mean_Distance, 
          x = Phenotype),
      color = "black", size = 2, shape = 19, inherit.aes = FALSE
    ) +
    scale_color_manual(values = named_colors) +
    labs(
      title = paste("Distance to", target_cell_type, "Cells Across Phenotypes"),
      subtitle = paste("Sorted by", metric, "distance (black dots show", metric, "distance)"),
      x = NULL, 
      y = paste("Distance to", target_cell_type, "(μm)")
    ) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "right")
  
  return(p)
}


#### plot distance density #### 
plot_distance_density <- function(
    data, 
    target_cell_type, 
    highlight_types = NULL,
    distance_col = NULL,
    color_palette = named_colors,
    other_color = "black",
    line_size = 1,
    alpha = 0.7,
    title = NULL,
    x_lab = NULL,
    y_lab = "Density",
    legend_title = "Cell Phenotype",
    legend_position = "right",
    theme = theme_classic(),
    show_plot = TRUE,
    ...
) {
  if (is.null(distance_col)) {
    distance_col <- paste("Distance to", target_cell_type)
    if (!distance_col %in% names(data)) {
      stop(paste("Distance column", distance_col, "not found in data"))
    }
  }
  
  if (!is.null(highlight_types)) {
    if (!all(highlight_types %in% data$Phenotype)) {
      missing_types <- setdiff(highlight_types, unique(data$Phenotype))
      warning("Some highlight_types not found in data: ", paste(missing_types, collapse = ", "))
      highlight_types <- intersect(highlight_types, unique(data$Phenotype))
    }
  } else {
    highlight_types <- unique(data$Phenotype)
  }
  
  if (is.null(title)) {
    title <- paste("Distance Distribution to", target_cell_type)
  }
  if (is.null(x_lab)) {
    x_lab <- paste("Distance to", target_cell_type, "(μm)")
  }
  
  plot_data <- data %>%
    mutate(
      plot_group = ifelse(
        .data$Phenotype %in% highlight_types,
        as.character(.data$Phenotype),
        "Other"
      )
    )
  
  plot_colors <- c(color_palette[highlight_types], "Other" = other_color)
  
  p <- ggplot(plot_data, aes(.data[[distance_col]])) +
    geom_density(
      linewidth = line_size,
      aes(color = plot_group),
      alpha = alpha,
      ...
    ) +
    scale_color_manual(
      name = legend_title,
      values = plot_colors,
      guide = guide_legend(
        override.aes = list(alpha = 1, linewidth = 1.5)
      )
    ) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    theme +
    theme(
      legend.position = legend_position,
      legend.title = element_text(face = "bold")
    )
  
  if (show_plot) {
    print(p)
  }
  
  invisible(p)
}
#### plot relative abundance ####
plot_relative_abundance <- function(
    data,
    target_cell_type,
    distance_col = NULL,
    fill_palette = named_colors,
    alpha = 0.8,
    position = "fill",
    title = NULL,
    x_lab = NULL,
    y_lab = "Proportion of Phenotype",
    legend_title = "Cell Phenotype",
    legend_position = "right",
    percent_labels = TRUE,
    theme = theme_classic(),
    show_plot = TRUE,
    ...
) {
  if (is.null(distance_col)) {
    distance_col <- paste("Distance to", target_cell_type)
    if (!distance_col %in% names(data)) {
      stop(paste("Distance column", distance_col, "not found in data"))
    }
  }
  
  if (is.null(title)) {
    title <- paste("Relative Abundance by Distance to", target_cell_type)
  }
  if (is.null(x_lab)) {
    x_lab <- paste("Distance to", target_cell_type, "(μm)")
  }
  
  missing_fills <- setdiff(unique(data$Phenotype), names(fill_palette))
  if (length(missing_fills) > 0) {
    warning("Missing colors for phenotypes: ", paste(missing_fills, collapse = ", "),
            "\nAssigning random colors.")
    fill_palette <- c(fill_palette, randomColors(length(missing_fills)))
    names(fill_palette)[(length(fill_palette)-length(missing_fills)+1):length(fill_palette)] <- missing_fills
  }
  
  p <- ggplot(data, aes(.data[[distance_col]], after_stat(count), fill = Phenotype)) +
    geom_density(
      position = position,
      alpha = alpha,
      ...
    ) +
    scale_fill_manual(
      name = legend_title,
      values = fill_palette,
      guide = guide_legend(
      override.aes = list(alpha = 1, size = 0.5)
      )
    ) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    theme +
    theme(
      legend.position = legend_position,
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(12, "pt")
    )
  
  if (percent_labels) {
    p <- p + scale_y_continuous(labels = scales::percent)
  }
  
  # if (length(unique(data$Phenotype)) > 15) {
  #   p <- p + theme(
  #     legend.key.size = unit(8, "pt"),
  #     legend.text = element_text(size = 6)
  #   )
  # }
  
  if (show_plot) {
    print(p)
  }
  
  invisible(p)
}

randomColors <- function(n) {
  hues <- seq(15, 375, length.out = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#### plot spatial distribution ####
plot_spatial_distribution <- function(
    data,
    phenotype_col = "Phenotype",
    x_col = "Cell X Position",
    y_col = "Cell Y Position",
    color_palette = named_colors,
    point_size = 0.5,
    alpha = 0.7,
    title = "Spatial Distribution of Cell Phenotypes",
    legend_title = "Phenotype",
    reverse_y = TRUE,  
    show_tissue_edge = FALSE,  
    tissue_edge_color = "black", 
    theme = theme_classic(),
    show_plot = TRUE,
    ...
) {
  if(!is.factor(data[[phenotype_col]])) {
    data[[phenotype_col]] <- factor(
      data[[phenotype_col]],
      levels = phenotypes 
    )
  }
  
  required_cols <- c(phenotype_col, x_col, y_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  p <- ggplot(data, aes(x = .data[[x_col]], 
                        y = if(reverse_y) -1*.data[[y_col]] else .data[[y_col]],  
                        color = .data[[phenotype_col]])) +
    geom_point(size = point_size, alpha = alpha, ...) +
    scale_color_manual(
      name = legend_title,
      values = color_palette,
      na.value = "grey70"
    ) +
    labs(
      title = title,
      x = "X Position (μm)",
      y = "Y Position (μm)"
    ) +
    coord_fixed(ratio = 1) + 
    theme +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.background = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  if(reverse_y){
    p <- p + scale_y_continuous(labels = function(x) abs(x))
  }
  
  if(show_tissue_edge){
    if(exists("tissue_boundary_data")){
      p <- p + 
        geom_path(
          data = tissue_boundary_data,
          aes(x = x, y = if(reverse_y) -1*y else y),
          color = tissue_edge_color,
          size = 0.5,
          inherit.aes = FALSE
        )
    } else {
      message("Note: tissue_boundary_data not found, skipping edge plotting")
    }
  }
  
  n_phenotypes <- length(unique(data[[phenotype_col]]))
  if(n_phenotypes > 12){
    p <- p + guides(color = guide_legend(
      override.aes = list(size = 3, alpha = 1),
      ncol = ceiling(n_phenotypes/20)
    ))
  }
  
  if(show_plot) print(p)
  invisible(p)
}

#### plot nearest neighbors ####
# plot_nearest_neighbors <- function(
#     data, 
#     target_type, 
#     neighbor_type,
#     target_color = "cyan2",
#     neighbor_color = "red",
#     line_color = "purple",
#     point_size = 1,
#     line_width = 0.5
# ) {
#   target_cells <- data %>% filter(select_rows(data, target_type))
#   neighbor_cells <- data %>% filter(select_rows(data, neighbor_type))
#   
#   neighbor_col <- paste("Cell ID", neighbor_type)
#   joined_data <- target_cells %>% 
#     left_join(neighbor_cells, 
#               by = setNames("Cell ID", neighbor_col),
#               suffix = c("", paste0(".", neighbor_type)))
#   
#   ggplot() +
#     geom_segment(
#       data = joined_data,
#       aes(x = `Cell X Position`, y = `Cell Y Position`,
#           xend = get(paste0("Cell X Position.", neighbor_type)),
#           yend = get(paste0("Cell Y Position.", neighbor_type))),
#       color = line_color, linewidth = line_width
#     ) +
#     geom_point(
#       data = target_cells,
#       aes(`Cell X Position`, `Cell Y Position`, color = "Target"),
#       size = point_size
#     ) +
#     geom_point(
#       data = neighbor_cells,
#       aes(`Cell X Position`, `Cell Y Position`, color = "Neighbor"),
#       size = point_size
#     ) +
#     scale_color_manual(
#       name = "Cell Type",
#       values = c("Target" = target_color, "Neighbor" = neighbor_color),
#       breaks = c("Target", "Neighbor"),
#       labels = c(target_type, neighbor_type)
#     ) +
#     coord_fixed() +
#     scale_y_reverse() +
#     labs(
#       title = paste("Spatial Relationship Between", target_type, "and", neighbor_type),
#       x = "X Position (μm)", 
#       y = "Y Position (μm)"
#     ) +
#     theme_bw() +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       panel.grid = element_blank(),
#       legend.position = "right"
#     )
# }
plot_nearest_neighbors <- function(
    data, 
    target_type, 
    neighbor_type,
    target_color = "cyan2",
    neighbor_color = "red",
    line_color = "purple",
    point_size = 1,
    line_width = 0.3,
    show_arrows = TRUE,     
    arrow_length = 0.02,     
    label_lines = FALSE      
) {
  target_cells <- data %>% filter(select_rows(data, target_type))
  neighbor_cells <- data %>% filter(select_rows(data, neighbor_type))
  
  neighbor_col <- paste("Cell ID", neighbor_type)
  joined_data <- target_cells %>% 
    left_join(neighbor_cells, 
              by = setNames("Cell ID", neighbor_col),
              suffix = c("", paste0(".", neighbor_type))) %>%
    mutate(
      Distance = sqrt(
        (`Cell X Position` - get(paste0("Cell X Position.", neighbor_type)))^2 +
          (`Cell Y Position` - get(paste0("Cell Y Position.", neighbor_type)))^2
      )
    )
  
  p <- ggplot() +
    geom_segment(
      data = joined_data,
      aes(x = `Cell X Position`, y = `Cell Y Position`,
          xend = get(paste0("Cell X Position.", neighbor_type)),
          yend = get(paste0("Cell Y Position.", neighbor_type))),
      color = line_color, 
      linewidth = line_width,
      arrow = if (show_arrows) arrow(length = unit(arrow_length, "npc")) else NULL
    ) +
    geom_point(
      data = target_cells,
      aes(`Cell X Position`, `Cell Y Position`, color = "Target"),
      size = point_size
    ) +
    geom_point(
      data = neighbor_cells,
      aes(`Cell X Position`, `Cell Y Position`, color = "Neighbor"),
      size = point_size
    ) +
    scale_color_manual(
      name = "Direction",
      values = c("Target" = target_color, "Neighbor" = neighbor_color),
      breaks = c("Target", "Neighbor"),
      labels = c(
        paste(target_type, "(Source)"), 
        paste(neighbor_type, "(Target)")
      )
    ) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = paste(target_type, "to Nearest", neighbor_type),
      x = "X Position (μm)", 
      y = "Y Position (μm)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  if (label_lines) {
    p <- p + 
      geom_text(
        data = joined_data,
        aes(
          x = (`Cell X Position` + get(paste0("Cell X Position.", neighbor_type))) / 2,
          y = (`Cell Y Position` + get(paste0("Cell Y Position.", neighbor_type))) / 2,
          label = sprintf("%.1f", Distance)
        ),
        color = "black", size = 3, vjust = -0.5
      )
  }
  
  return(p)
}
#### compare_cell_distance ----
compare_cell_distance <- function(
    data, 
    target_cell_type, 
    source_cell_type,
    group_col = "Group",
    sample_col = "Sample",
    phenotype_col = "Phenotype",
    metric = c("median", "mean"),
    test_method = c("wilcox.test", "t.test")
) {
  # 加载所需包
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(dplyr)) install.packages("dplyr")
  if (!require(tidyr)) install.packages("tidyr")
  if (!require(ggpubr)) install.packages("ggpubr")
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  
  # 参数验证
  if (!target_cell_type %in% data[[phenotype_col]]) {
    stop("Target cell type not found in Phenotype column")
  }
  
  if (!source_cell_type %in% data[[phenotype_col]]) {
    stop("Source cell type not found in Phenotype column")
  }
  
  if (!group_col %in% names(data)) {
    stop("Group column not found in data")
  }
  
  if (!sample_col %in% names(data)) {
    stop("Sample column not found in data")
  }
  
  metric <- match.arg(metric)
  test_method <- match.arg(test_method)
  
  # 自动确定距离列名
  distance_col <- paste("Distance to", target_cell_type)
  if (!distance_col %in% names(data)) {
    stop(paste("Distance column", distance_col, "not found in data"))
  }
  
  # 筛选特定源细胞类型的数据
  source_data <- data[data[[phenotype_col]] == source_cell_type, ]
  
  # 计算每个样本的汇总统计量
  sample_stats <- source_data %>%
    group_by(.data[[group_col]], .data[[sample_col]]) %>%
    summarise(
      Median_Distance = median(.data[[distance_col]], na.rm = TRUE),
      Mean_Distance = mean(.data[[distance_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # 准备统计检验
  groups <- unique(sample_stats[[group_col]])
  if (length(groups) != 2) {
    stop("Exactly two groups are required for comparison")
  }
  
  # 提取配对样本标识
  sample_stats$Pair_ID <- gsub("_.*$", "", sample_stats[[sample_col]])
  
  # 提取各组数据
  group1_data <- sample_stats[sample_stats[[group_col]] == groups[1], ]
  group2_data <- sample_stats[sample_stats[[group_col]] == groups[2], ]
  
  if (metric == "median") {
    group1_values <- group1_data$Median_Distance
    group2_values <- group2_data$Median_Distance
  } else {
    group1_values <- group1_data$Mean_Distance
    group2_values <- group2_data$Mean_Distance
  }
  
  # 执行统计检验
  if (test_method == "t.test") {
    test_result <- t.test(group1_values, group2_values, paired = TRUE)
    test_name <- "Paired t-test"
  } else {
    test_result <- wilcox.test(group1_values, group2_values, paired = TRUE)
    test_name <- "Paired Wilcoxon test"
  }
  
  # 创建结果数据框
  result_df <- data.frame(
    Source_Cell_Type = source_cell_type,
    Target_Cell_Type = target_cell_type,
    Group1 = groups[1],
    Group2 = groups[2],
    Group1_Value = if (metric == "median") median(group1_values) else mean(group1_values),
    Group2_Value = if (metric == "median") median(group2_values) else mean(group2_values),
    Group1_SD = sd(group1_values),
    Group2_SD = sd(group2_values),
    Statistic = test_result$statistic,
    P_value = test_result$p.value,
    Method = test_name,
    Metric = metric,
    stringsAsFactors = FALSE
  )
  
  # 创建可视化
  # 将分组转换为数值坐标
  sample_stats$group_numeric <- ifelse(sample_stats[[group_col]] == groups[1], 1, 2)
  
  # 创建 y_value 列
  sample_stats$y_value <- if (metric == "median") {
    sample_stats$Median_Distance
  } else {
    sample_stats$Mean_Distance
  }
  
  # 创建配对数据框用于连线 - 修复列名问题
  paired_data <- sample_stats %>%
    select(Pair_ID, group_numeric, y_value) %>%
    pivot_wider(
      names_from = group_numeric, 
      values_from = c(y_value),
      names_prefix = "y_value_"
    )
  
  # 确保配对数据框包含所有配对样本
  paired_data <- paired_data[complete.cases(paired_data), ]
  
  # 计算显著性标注的位置
  y_max <- max(sample_stats$y_value) * 1.1
  signif_label <- ifelse(test_result$p.value < 0.001, "***",
                         ifelse(test_result$p.value < 0.01, "**",
                                ifelse(test_result$p.value < 0.05, "*", "ns")))
  
  # 创建基础绘图对象
  p <- ggplot(sample_stats, aes(x = group_numeric)) +
    # 小提琴图展示分布
    geom_violin(
      aes(y = y_value, fill = .data[[group_col]]),
      width = 0.6, alpha = 0.3, trim = TRUE, scale = "width"
    ) +
    # 箱线图展示统计量
    geom_boxplot(
      aes(y = y_value, fill = .data[[group_col]]),
      width = 0.3, alpha = 0.7, outlier.shape = NA
    ) +
    # 添加配对连线 - 所有点位于组中心
    geom_segment(
      data = paired_data,
      aes(x = 1, xend = 2,
          y = y_value_1, yend = y_value_2,
          color = Pair_ID),
      size = 1, alpha = 0.8
    ) +
    # 添加样本点 - 所有点位于组中心
    geom_point(
      aes(y = y_value, x = group_numeric, color = Pair_ID, shape = Pair_ID),
      size = 4#, position = position_jitter(width = 0.05, height = 0)
    ) +
    # 美化设置
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF")) + 
    scale_color_manual(values = c("#00A087FF", "#3C5488FF", "#F39B7FFF")) + 
    scale_shape_manual(values = c(16, 17, 15)) +
    scale_x_continuous(
      breaks = c(1, 2),
      labels = groups,
      limits = c(0.5, 2.5)
    ) +
    labs(
      title = paste("Distance from", source_cell_type, "to", target_cell_type),
      subtitle = paste("Comparison between", groups[1], "and", groups[2], "\n",
                       test_name, "p-value =", format.pval(test_result$p.value, digits = 3)),
      x = "Group",
      y = paste("Distance (μm) -", metric),
      color = "Sample Pair",
      shape = "Sample Pair",
      fill = "Group"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.grid.major.y = element_line(color = "grey90", size = 0.2),
      panel.grid.minor.y = element_blank(),
      legend.box = "horizontal",
      legend.spacing.x = unit(0.2, "cm"),
      legend.key.size = unit(0.8, "cm"),
      aspect.ratio = 1 
    ) +
    # 添加统计显著性标记
    annotate("segment", 
             x = 1, xend = 2, 
             y = y_max, yend = y_max,
             size = 0.5) +
    annotate("text", 
             x = 1.5, y = y_max * 1.05, 
             label = signif_label, 
             size = 5)
  
  # 返回结果和图表
  return(list(
    statistics = result_df,
    plot = p,
    sample_data = sample_stats,
    paired_data = paired_data
  ))
}
#### compare_cell_distance_all_cells ----
compare_cell_distance_all_cells <- function(
    data, 
    target_cell_type, 
    source_cell_type,
    group_col = "Group",
    sample_col = "Sample",
    phenotype_col = "Phenotype",
    test_method = c("wilcox.test", "t.test")
) {
  # 加载所需包
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(dplyr)) install.packages("dplyr")
  if (!require(ggpubr)) install.packages("ggpubr")
  if (!require(ggsignif)) install.packages("ggsignif")
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(ggsignif)
  
  # 参数验证
  if (!target_cell_type %in% data[[phenotype_col]]) {
    stop("Target cell type not found in Phenotype column")
  }
  
  if (!source_cell_type %in% data[[phenotype_col]]) {
    stop("Source cell type not found in Phenotype column")
  }
  
  if (!group_col %in% names(data)) {
    stop("Group column not found in data")
  }
  
  if (!sample_col %in% names(data)) {
    stop("Sample column not found in data")
  }
  
  test_method <- match.arg(test_method)
  
  # 自动确定距离列名
  distance_col <- paste("Distance to", target_cell_type)
  if (!distance_col %in% names(data)) {
    stop(paste("Distance column", distance_col, "not found in data"))
  }
  
  # 筛选特定源细胞类型的数据
  source_data <- data[data[[phenotype_col]] == source_cell_type, ]
  
  # 准备统计检验
  groups <- unique(source_data[[group_col]])
  if (length(groups) != 2) {
    stop("Exactly two groups are required for comparison")
  }
  
  # 提取配对样本标识
  source_data$Pair_ID <- gsub("_.*$", "", source_data[[sample_col]])
  
  # 提取各组数据
  group1_data <- source_data[source_data[[group_col]] == groups[1], ]
  group2_data <- source_data[source_data[[group_col]] == groups[2], ]
  
  # 执行统计检验 - 直接使用所有细胞
  if (test_method == "t.test") {
    test_result <- t.test(group1_data[[distance_col]], group2_data[[distance_col]])
    test_name <- "Student's t-test"
  } else {
    test_result <- wilcox.test(group1_data[[distance_col]], group2_data[[distance_col]])
    test_name <- "Wilcoxon rank-sum test"
  }
  
  # 创建结果数据框
  result_df <- data.frame(
    Source_Cell_Type = source_cell_type,
    Target_Cell_Type = target_cell_type,
    Group1 = groups[1],
    Group2 = groups[2],
    Group1_Mean = mean(group1_data[[distance_col]], na.rm = TRUE),
    Group2_Mean = mean(group2_data[[distance_col]], na.rm = TRUE),
    Group1_Median = median(group1_data[[distance_col]], na.rm = TRUE),
    Group2_Median = median(group2_data[[distance_col]], na.rm = TRUE),
    Statistic = test_result$statistic,
    P_value = test_result$p.value,
    Method = test_name,
    stringsAsFactors = FALSE
  )
  
  # 计算显著性标注的位置
  y_max <- max(source_data[[distance_col]], na.rm = TRUE) * 1.1
  signif_label <- ifelse(test_result$p.value < 0.001, "***",
                         ifelse(test_result$p.value < 0.01, "**",
                                ifelse(test_result$p.value < 0.05, "*", "ns")))
  
  # 创建箱线图
  p_box <- ggplot(source_data, aes(x = .data[[group_col]], y = .data[[distance_col]])) +
    geom_jitter(aes(color = Pair_ID), width = 0.2, alpha = 0.3, size = 2) +
    geom_boxplot(aes(fill = .data[[group_col]]), width = 0.6, alpha = 0.7, outlier.shape = NA) +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
    scale_color_manual(values = c("#00A087FF", "#3C5488FF", "#F39B7FFF")) +
    labs(
      title = paste("Distance from", source_cell_type, "to", target_cell_type),
      subtitle = paste(test_name, "p-value =", format.pval(test_result$p.value, digits = 3)),
      x = "Group",
      y = paste("Distance to", target_cell_type, "(μm)"),
      fill = "Group",
      color = "Sample Pair"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    geom_signif(
      comparisons = list(groups),
      map_signif_level = TRUE,
      textsize = 5,
      vjust = 0.2
    )
  
  # 创建密度图
  p_density <- ggplot(source_data, aes(x = .data[[distance_col]], fill = .data[[group_col]])) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Pair_ID, ncol = 1) +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
    labs(
      title = paste("Distance Distribution from", source_cell_type, "to", target_cell_type),
      x = paste("Distance to", target_cell_type, "(μm)"),
      y = "Density",
      fill = "Group"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
  
  # 返回结果和图表
  return(list(
    statistics = result_df,
    boxplot = p_box,
    density_plot = p_density,
    cell_data = source_data
  ))
}