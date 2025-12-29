rm(list = ls())
library(schard)
library(semla)
library(qs)
library(Seurat)
library(ClusterR)
library(COSG)
library(tidyverse)
library(reticulate)
library(patchwork)
# library(SCP)
gc()
options(future.globals.maxSize = 1000 * 1024^2)

sample_name <- c('CTB1NT','CTB1T','CTB3NT','CTB3T','CTB7NT','CTB7T')

source("scripts/02_Niche_single_function.R")

st_semla <- qread("qs_data/RCTD_results/RCTD_results_full_reference_subtype.qs")
mt <- st_semla@meta.data
colnames(mt)
mt <- mt[, !grepl("niche", colnames(mt), ignore.case = TRUE)]
colnames(mt)
st_semla@meta.data <- mt 
table(st_semla$first_type)
table(st_semla$celltype)

# neighbors.k = 30, niches.k = 10
CTB1NT <- Niche_Analysis(st_semla, image = "CTB1NT",celltype="first_type", neighbors.k = 30, niches.k = 10)
CTB1T <- Niche_Analysis(st_semla, image = "CTB1T",celltype="first_type", neighbors.k = 30, niches.k = 10)
CTB3NT <- Niche_Analysis(st_semla, image = "CTB3NT", celltype="first_type",neighbors.k = 30, niches.k = 10)
CTB3T <- Niche_Analysis(st_semla, image = "CTB3T", celltype="first_type",neighbors.k = 30, niches.k = 10)
CTB7NT <- Niche_Analysis(st_semla, image = "CTB7NT", celltype="first_type", neighbors.k = 30, niches.k = 10)
CTB7T <- Niche_Analysis(st_semla, image = "CTB7T", celltype="first_type", neighbors.k = 30, niches.k = 10)

CTB_list <- lapply(sample_name,FUN=function(x){
  niche <- Niche_Analysis(st_semla, 
                          image = x,
                          celltype="first_type", 
                          neighbors.k = 30, 
                          niches.k = 12)
  return(niche)
})

niches_data <- lapply(seq_along(CTB_list), function(x) {
  data.frame(
    Sample = sample_name[x],
    Cell_Barcode = rownames(CTB_list[[x]]@meta.data), 
    niches = CTB_list[[x]]@meta.data$niches,
    stringsAsFactors = FALSE
  )
})
combined_niches <- do.call(rbind, niches_data)
head(combined_niches)
tail(combined_niches)


## function
run_multiple_niche_analysis <- function(st_semla, 
                                        sample_name, 
                                        k_values = c(10, 12, 15, 18, 20)) {
  all_results <- list()
  for (k in k_values) {
    message("Processing niches.k = ", k, "...")
    CTB_list <- lapply(sample_name, function(x) {
      niche_result <- Niche_Analysis(
        object = st_semla, 
        image = x,
        celltype = "first_type", 
        neighbors.k = 30, 
        niches.k = k
      )
      return(niche_result)
    })
    
    niches_data <- lapply(seq_along(CTB_list), function(i) {
      df <- data.frame(
        Sample = sample_name[i],
        Cell_Barcode = rownames(CTB_list[[i]]@meta.data),
        stringsAsFactors = FALSE
      )
      df[[paste0("niches_", k)]] <- CTB_list[[i]]@meta.data$niches
      return(df)
    })
    combined_data <- do.call(rbind, niches_data)
    all_results[[as.character(k)]] <- combined_data
  }
  
  final_result <- Reduce(
    function(x, y) {
      merge(x, y, by = c("Sample", "Cell_Barcode"), all = TRUE)
    }, 
    all_results
  )
  
  return(final_result)
}

final_dataframe <- run_multiple_niche_analysis(
  st_semla = st_semla,
  sample_name = sample_name,
  k_values = c(10, 12, 15, 18, 20) 
)
head(final_dataframe)

write.csv(final_dataframe,file = "./qs_data/RCTD_results/all_niche_single.csv")





Idents(CTB3NT) <- "niches"
Niche_split(CTB3NT)
Idents(CTB1T) <- "niches"
Niche_split(CTB1T)
dev.off()
colnames(CTB1NT@meta.data)
table(CTB1NT$niches)

CTB1NT$niches <- factor(CTB1NT$niches)
CTB1NT$niches <- factor(CTB1NT$niches)
library(ggrastr)
p <- SpatialDimPlot(CTB1NT, images = "CTB1NT", 
                        group.by = "niches")& 
  # scale_y_reverse()&
  coord_fixed()&
  ggsci::scale_fill_d3("category20")&
  NoLegend()&
  theme_void()&
  guides(fill = guide_legend(override.aes = list(size = 6)))&
  # DarkTheme()&
  # scale_fill_gradientn(colours=viridis_plasma_light_high, 
  #                      na.value = "black", 
  #                      limits=c(0,3))&
  ggtitle(sample_name[1])
p
rasterise(p,dpi = 300)
dev.off()
gc()


