
rm(list = ls())
library(Seurat)
library(tidyverse)
setwd(here::here())
setwd("~/OneDrive_4_2025-3-8/Pyscenic/")
source("R/makeMetaCells.R")
if (!dir.exists("output")) {
  dir.create("output")
}


## clusters were annotated after batch effect removal via harmony.
seu <- qs::qread("~/OneDrive_4_2025-3-8/data/ctb_scobj_v4.qs")
Idents(seu) <- "disease"
table(Idents(seu))
seu <- subset(seu, idents = c('CTB_BL',"CTB_T"))
table(seu$celltype, seu$disease)
seu <- subset(seu, celltype %in% c("Fibroblast", "Myeloid", "T cell"))
table(seu$subtype, seu$disease)
seu$celltype <- seu$subtype
seu$group <- seu$orig.ident2
DimPlot(seu, group.by = "celltype", split.by = "disease") & ggsci::scale_color_d3("category20")

seu.list <- SplitObject(seu, split.by = "group")
seu.list <- lapply(seu.list, function(object) {
  object@project.name <- unique(object$group)
  return(object)
})


metacells.list <- lapply(seq_along(seu.list), function(ii) {
  makeMetaCells(
    seu       = seu.list[[ii]],
    min.cells = 10,
    reduction = "umap",
    dims      = 1:2,
    k.param   = 10,
    resolution = 30,
    cores     = 10)
})

mc.mat <- lapply(metacells.list, function(mc) mc$mat) %>% Reduce(cbind, .)
mc.cellmeta <- lapply(metacells.list, function(mc) mc$metadata) %>% Reduce(rbind, .)
dim(mc.mat)
summary(mc.cellmeta$CELL_COUNT)

saveRDS(mc.mat, "output/00-1.mc.mat.rds")


## 准备pySCENIC的输入文件

## (1) TF list文件(optional)，可以使用预定义的TF list，例如pySCENIC官方提供的，或者animalTFDB提供的。
motif2tfs <- data.table::fread("cisTarget_db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
TFs <- sort(unique(motif2tfs$gene_name))
writeLines(TFs, "cisTarget_db/hsa_hgnc_tfs.motifs-v10.txt")

## (2) meta cell matrix (for step1): *.csv or *.loom
mc.mat <- readRDS("output/00-1.mc.mat.rds")
## (2.1) 过滤低表达基因
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5, ]
dim(mc.mat)
## (2.2) 过滤不在cisTargetDB中的基因
cisdb <- arrow::read_feather("cisTarget_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
length(genes.use)
dim(mc.mat)
mc.mat <- mc.mat[genes.use, ]
dim(mc.mat)

loom <- SCopeLoomR::build_loom(
  file.name         = "output/00-2.mc_mat_for_step1.loom",
  dgem              = mc.mat,
  default.embedding = NULL
)
loom$close()

rm(loom)
gc()

# install.packages("devtools")
# devtools::install_github("aertslab/SCopeLoomR")
