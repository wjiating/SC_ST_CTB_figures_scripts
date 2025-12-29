rm(list = ls())
gc()
library(patchwork)
library(SPOTlight)
library(semla)
library(Seurat)
library(qs)
library(harmony)


SCC <- qread('./semla_output/SCC.qs')
SCC

se <- qread('./semla_output/st_semla_clean.qs')
mt1 <- se@meta.data %>% tibble::rownames_to_column("cell_id")
mt2 <- SCC@meta.data %>% tibble::rownames_to_column("cell_id")
table(identical(rownames(mt1), rownames(mt2)))
mt_new <- left_join(mt2, mt1, by = "cell_id")%>% tibble::column_to_rownames("cell_id")
SCC@meta.data <-mt_new

if (!"spatial_x" %in% colnames(SCC@meta.data)) {
  SCC@meta.data$spatial_x <- NA_real_
}
if (!"spatial_y" %in% colnames(SCC@meta.data)) {
  SCC@meta.data$spatial_y <- NA_real_
}
colnames(SCC@meta.data)

for (img_name in names(SCC@images)) {
  img_coords <- SCC@images[[img_name]]@coordinates
  img_cells <- rownames(img_coords)
  missing_cells <- setdiff(img_cells, colnames(SCC))
  if (length(missing_cells) > 0) {
    warning(paste(length(missing_cells), "cells in image", img_name, "not found in object"))
    valid_idx <- img_cells %in% colnames(SCC)
    img_cells <- img_cells[valid_idx]
    img_coords <- img_coords[valid_idx, ]
  }
  SCC@meta.data[img_cells, "spatial_x"] <- img_coords$x
  SCC@meta.data[img_cells, "spatial_y"] <- img_coords$y
  cat("Added coordinates for", length(img_cells), "cells in", img_name, "\n")
}

mt <- SCC@meta.data
colnames(mt)

ctb_csd <- mt[,c("orig.ident.x","celltype.x","first_type","spatial_x","spatial_y" )]
ctb_csd <- ctb_csd %>% 
  tibble::rownames_to_column("CellID")

colnames(ctb_csd) <- c("Cell ID", "Sample","Celltype", "Phenotype", 
                       "Cell X Position", "Cell Y Position")
str(ctb_csd)
ctb_csd$Phenotype <- as.character(ctb_csd$Phenotype)
write.csv(ctb_csd, file = "./semla_output/ctb_csd_data.csv", row.names = FALSE)


table(ctb_csd$Sample)
csd_1bl <- ctb_csd %>% 
  dplyr::filter(Sample == "CTB1_NT")
head(csd_1bl)
library(phenoptr)
names(csd_1bl) <- c("Cell ID", "Sample", "Celltype", "Phenotype", 
                    "Cell X Position", "Cell Y Position")
csd_1bl <- csd_1bl %>% 
  dplyr::select(c("Phenotype","Cell ID","Cell X Position", "Cell Y Position",
                  "Sample", "Celltype"))
csd_1bl$Phenotype <- as.character(csd_1bl$Phenotype)
str(csd_1bl)

# Get unique phenotypes first
phenotypes <- unique(csd_1bl$Phenotype)
distances <- find_nearest_distance(csd_1bl, phenotypes = phenotypes)

distances <- find_nearest_distance(csd_1bl)
glimpse(distances)

csd_with_distance <- bind_cols(csd_1bl, distances)


ggplot(csd_with_distance, aes(`Distance to SPP1 Mac`, color=Phenotype)) +
  geom_density(linewidth=1)+
  scale_color_manual(values = named_colors)+
  theme_classic()

highlight_types <- c("SPP1 Mac", "Th17_Treg", "CCL19+ FB", "Langerhans")
ggplot(csd_with_distance, aes(`Distance to SPP1 Mac`)) +
  geom_density(
    linewidth = 1,
    aes(color = ifelse(.data$Phenotype %in% highlight_types, 
                       as.character(.data$Phenotype), 
                       "Other")),
    alpha = 0.7
  ) +
  scale_color_manual(
    name = "Cell Phenotype",
    values = c(named_colors[highlight_types], "Other" = "grey80")
  ) +
  theme_classic()

library(ggridges)
ggplot(csd_with_distance, aes(x = `Distance to SPP1 Mac`, y = Phenotype, fill = Phenotype)) +
  geom_density_ridges(scale = 0.9, alpha = 0.7) +
  scale_fill_manual(values = named_colors) +
  theme_ridges() +
  theme(legend.position = "none")

ggplot(csd_with_distance, aes(`Distance to SPP1 Mac`, after_stat(count), fill=Phenotype)) +
  geom_density(position="fill", alpha=0.8) +
  scale_fill_manual(values = named_colors) +
  labs(title="Relative Abundance by Distance to SPP1 Mac",
       x="Distance to SPP1 Mac (μm)", y="Proportion of Phenotype") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic()

ggplot(csd_with_distance, aes(`Distance to SPP1 Mac`, after_stat(count), fill=Phenotype)) +
  geom_density(position="fill", alpha=0.8) +
  scale_fill_manual(values = named_colors) +
  labs(title="Relative Abundance by Distance to SPP1 Mac",
       x="Distance to SPP1 Mac (μm)", y="Proportion of Phenotype") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic()

df <- csd_with_distance %>% group_by(Phenotype) %>% 
  select(Phenotype, starts_with('Distance to')) %>% 
  summarize_all(~round(mean(.), 1)) %>% 
  tibble::column_to_rownames("Phenotype") %>% 
  t()


ggplot(csd_with_distance, aes(Phenotype, `Distance to SPP1 Mac`, color=Phenotype)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.3, size=1.5) +
  scale_color_manual(values=c(named_colors)) +
  labs(title="Distance to SPP1 Mac Cells Across Phenotypes",
       x=NULL, y="Distance (μm)") +
  coord_flip() +  
  theme_classic()


source("scripts/15_Plot_distance_func.R")

plot_cell_distance(csd_with_distance, "SPP1 Mac", metric = "median")
plot_cell_distance(csd_with_distance, "Th17_Treg", metric = "median")
plot_cell_distance(csd_with_distance, "TREM2 Mac", metric = "median")
plot_cell_distance(csd_with_distance, "CCL19+ FB", metric = "median")

plot_distance_density(csd_with_distance,"SPP1 Mac")
plot_distance_density(csd_with_distance,"SPP1 Mac",
  highlight_types = c("SPP1 Mac", "Th17_Treg", "CCL19+ FB", "Langerhans"))

plot_relative_abundance(csd_with_distance, "SPP1 Mac")
plot_relative_abundance(csd_with_distance,"Th17_Treg")

head(csd_1bl)
plot_spatial_distribution(csd_1bl)

plot_nearest_neighbors(csd_with_distance, "SPP1 Mac", "CCL19+ FB")
plot_nearest_neighbors(csd_with_distance, "SPP1 Mac", "Th17_Treg")
plot_nearest_neighbors(csd_with_distance, "SPP1 Mac", "Th17_Treg", show_arrows = F)
plot_nearest_neighbors(csd_with_distance, "SPP1 Mac", "cDC2A")
plot_nearest_neighbors(csd_with_distance, "Th17_Treg", "SPP1 Mac")

count_within(csd, from='CD68+', to='CK+', radius=25)

table(csd_with_distance$Phenotype)

count_within(csd_1bl, from='SPP1 Mac', to='Th17_Treg', radius=c(25, 50))
count_within(csd_1bl, from='Th17_Treg', to='SPP1 Mac', radius=c(25, 50))

