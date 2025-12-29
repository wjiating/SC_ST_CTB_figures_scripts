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
# library(scCustomize)
gc()
options(future.globals.maxSize = 1000 * 1024^2) 


scRNA <- qread(file = "qs_data/CTB_all_subtype.qs")
Idents(scRNA) <- "subtype"
subtype_cols <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#FFE119", "#F58231", "#911EB4", "#46F0F0", "#F032E6",
  "#BCBD22", "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#AAFFC3", "#808000", "#FFD8B1",
  "#000075", "#A9A9A9", "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
  "#E377C2", "#7F7F7F", "#17BECF", "#FF9898", "#C5B0D5", "#C49C94", "#F7B6D2", "#C7C7C7",
  "#DBDB8D", "#9EDAE5", "#FFBB78", "#98DF8A", "#AEC7E8", "#FF9E9E", "#B5D0D9", "#D4B9D6",
  "#FDB462", "#80B1D3"
)
plot <- DimPlot(scRNA, cols = subtype_cols, reduction = "umap")
LabelClusters(plot = plot, id = "ident")
dev.off()

table(scRNA$subtype) %>% sort()

ref.ds <- subset(scRNA, downsample = 300)
table(ref.ds$subtype) %>% sort()
dim(ref.ds)
ref.ds$subtype <- gsub("/", "_", ref.ds$subtype)
table(ref.ds$subtype) %>% sort()

# extract information to pass to the RCTD Reference function
counts <- as.matrix(ref.ds@assays$RNA@layers$counts)
dim(counts)
head(counts)[1:5, 1:5]
rownames(counts) <- rownames(ref.ds)
colnames(counts) <- colnames(ref.ds)
head(counts)[1:5, 1:5]

nUMI <- colSums(counts)

cellType <- data.frame(barcode = colnames(ref.ds), celltype = ref.ds$subtype)
names(cellType) <- c('barcode', 'cell_type')
cell_types <- cellType$cell_type
names(cell_types) <- cellType$barcode
cell_types <- as.factor(cell_types)

reference <- Reference(counts, cell_types, nUMI)
class(reference)
rm(scRNA);gc()


stRNA <- qread(file="./qs_data/SCC_v5.qs")
SpatialDimPlot(stRNA, group.by="celltype", images = "CTB1NT")
MapLabels(stRNA, section_number = 2, column_name = "celltype") & 
  theme(legend.position = "right", legend.text = element_text(angle = 0)) &
  ggtitle("CTB1NT") &
  guides(fill = guide_legend(override.aes = list(size = 6)))
# sample_name <- c('CTB1NT','CTB1T','CTB3NT','CTB3T','CTB7NT','CTB7T')

table(stRNA$sample_id)
sto.list <- Seurat::SplitObject(stRNA, split.by = "sample_id")
slice <- names(sto.list)
str(stRNA@assays$Spatial)
# Formal class 'Assay5' [package "SeuratObject"] with 8 slot

result.RCTD <- data.frame()
result.RCTD.weights<- data.frame()

# counts = GetAssayData(object = sto.list[["CTB1_NT"]],
#                       assay = "Spatial",
#                       layer = "counts")
# head(counts)[1:5,1:5]
# counts@x
# coords <- GetTissueCoordinates(sto.list[["CTB1_NT"]])[,1:2]
# head(coords)
# ggplot(coords,aes(x,-y))+geom_point(size=1)

for(i in slice){
  cat("Processing sample:", i, "\n")

  current_sample <- sto.list[[i]]

  coords <- GetTissueCoordinates(current_sample)[,1:2]

  current_counts <- GetAssayData(object = current_sample,
                                 assay = "Spatial",
                                 layer = "counts")

  nUMI_vector <- colSums(current_counts) 
  names(nUMI_vector) <- colnames(current_sample)

  query <- SpatialRNA(
    coords = coords,
    counts = current_counts,
    nUMI = nUMI_vector
  )
  
  RCTD <- create.RCTD(spatialRNA = query, reference = reference, max_cores = 20)
  RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet')

  result_df <- as.data.frame(RCTD@results$results_df)
  result_df$sample <- i  
  result.RCTD <- rbind(result.RCTD, result_df)
  
  result_weights <- as.data.frame(as.matrix(RCTD@results$weights))
  result_weights$sample <- i  
  result.RCTD.weights <- rbind(result.RCTD.weights, result_weights)
}

head(result.RCTD)[1:5,1:5]
table(result.RCTD$spot_class,result.RCTD$first_class)

gc()

dim(stRNA)
identical(rownames(result.RCTD), colnames(stRNA))

stRNA1 <- stRNA[, rownames(result.RCTD)]

stRNA1 <- AddMetaData(stRNA1, metadata = result.RCTD)
stRNA1 <- AddMetaData(stRNA1, metadata = result.RCTD.weights)

Idents(stRNA1) <- stRNA1$first_type

table(stRNA1$first_type, stRNA1$celltype)
qsave(result.RCTD, file = "qs_data/RCTD_results/result.RCTD_full_subtype.qs")
qsave(result.RCTD.weights, file = "qs_data/RCTD_results/result.RCTD.weights_full_subtype.qs")
qsave(stRNA1, file = "qs_data/RCTD_results/RCTD_results_full_reference_subtype.qs")

