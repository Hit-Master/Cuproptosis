rm(list = ls())

library(Seurat)
library(dplyr)
library(readxl)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)

nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 30
pHB_lower <- 0
pHB_upper <- 5

theme_set(theme_cowplot())

use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")

sample_list <- read_excel("data/metadata/patients_metadata.xlsx", range = cell_cols("A:A")) %>% .$sample_id

sc_list1 <- CreateSeuratObject(counts = Read10X(data.dir = paste0("data/cellranger/", sample_list[1], "/filtered_feature_bc_matrix")), project = sample_list[1], min.cells = 3, min.features = 200)
sc_list2 <- lapply(sample_list[-1], function(sample_id){ 
  CreateSeuratObject(counts = Read10X(data.dir = paste0("data/cellranger/", sample_id, "/filtered_feature_bc_matrix")), 
                     project = sample_id, min.cells = 3, min.features = 200)
})

seu_obj <- merge(sc_list1, y = sc_list2, add.cell.ids = sample_list, project = "lung")

rm(sc_list1)
rm(sc_list2)

seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL", col.name = "pRP")

qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")

for (i in seq_along(qcparams)){
  print(VlnPlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident", pt.size = 0))
}
for (i in seq_along(qcparams)){
  print(RidgePlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident"))
}

VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "pMT"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
ggsave2("VlnPlot.pdf", path = "figure/Preprocessing/QC/", width = 20, height = 20, units = "cm")

RidgePlot(object = seu_obj, features = c("nFeature_RNA", "nCount_RNA", "pMT"), group.by = "orig.ident", ncol = 1, log = T)
ggsave2("RidgePlot.pdf", path = "figure/Preprocessing/QC/", width = 20, height = 20, units = "cm")

qc_std_plot_helper <- function(x) x + 
  scale_color_viridis() +
  geom_point(size = 0.01, alpha = 0.3)

qc_std_plot <- function(seu_obj) {
  qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))
  plot_grid(
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pHB, color = nFeature_RNA))) + 
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pRP, color = nFeature_RNA))) + 
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pHB, color = nCount_RNA))) + 
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pRP, color = nCount_RNA))) + 
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    
    
    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 5),
    
    ncol = 3, align = "hv"
  )
}


# Before filtering
seu_obj_unfiltered <- seu_obj
qc_std_plot(seu_obj_unfiltered)
ggsave2("before_filtering.pdf", path = "figure/Preprocessing/QC/", width = 30, height = 30, units = "cm")

# After filtering
seu_obj <- subset(seu_obj_unfiltered, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & pMT < pMT_upper & pHB < pHB_upper)
qc_std_plot(seu_obj)
ggsave2("after_filtering.pdf", path = "figure/Preprocessing/QC/", width = 30, height = 30, units = "cm")

seu_obj_unfiltered
seu_obj

# Data normalization
seu_obj <- SCTransform(seu_obj, verbose = T, vars.to.regress = c("nCount_RNA", "pMT"), conserve.memory = T)
# saveRDS(seu_obj, file = "seurat_object/Preprocessing/SCTransform.RDS")
