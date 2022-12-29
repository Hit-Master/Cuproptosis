library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(readr)
library(stringr)
library(progeny)
library(scales)

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

epi <- readRDS("seurat_object/Preprocessing/epi.RDS")
imm <- readRDS("seurat_object/Preprocessing/imm.RDS")
str <- readRDS("seurat_object/Preprocessing/str.RDS")

epi <- RunPCA(epi)
ElbowPlot(epi,  ndims = 50)
ggsave2("epi_Elbow.pdf", path = "figure/Recluster/epithelial/resolution/", width = 10, height = 5)

# for (i in c(10, 15, 20, 25)){
#   umaptest <- RunUMAP(epi, dims = 1:i, verbose = F)
#   print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", split.by = "tissue_type") + labs(title = paste0(i, "dimensions")))
#   remove(umaptest)
# }

epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  epi <- FindClusters(epi, resolution = i)
  DimPlot(epi, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("resolution", i,".pdf"), path = "figure/Recluster/epithelial/resolution/", width = 10, height = 10, units = "cm")
}

Idents(epi) <- epi@meta.data$SCT_snn_res.1


imm <- RunPCA(imm)
ElbowPlot(imm,  ndims = 50)
ggsave2("imm_Elbow.pdf", path = "figure/Recluster/immune/resolution/", width = 10, height = 5)

# for (i in c(10, 15, 20, 25)){
#  umaptest <- RunUMAP(imm, dims = 1:i, verbose = F)
#  print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", split.by = "tissue_type") + labs(title = paste0(i, " dimensions")))
#  remove(umaptest)
# }

imm <- RunUMAP(imm, dims = 1:20)
imm <- FindNeighbors(imm, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  imm <- FindClusters(imm, resolution = i)
  DimPlot(imm, reduction = "umap") + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("resolution", i,".pdf"), path = "figure/Recluster/immune/resolution/", width = 10, height = 10, units = "cm")
}

Idents(imm) <- imm@meta.data$SCT_snn_res.0.5


str <- RunPCA(str)
ElbowPlot(str, ndims = 50)
ggsave2("str_Elbow.pdf", path = "figure/Recluster/stromal/resolution/", width = 10, height = 5)

str <- RunUMAP(str, dims = 1:20)
str <- FindNeighbors(str, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  str <- FindClusters(str, resolution = i)
  DimPlot(str, reduction = "umap") + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("resolution", i,".pdf"), path = "figure/Recluster/stromal/resolution/", width = 10, height = 10, units = "cm")
}

Idents(str) <- str@meta.data$SCT_snn_res.1


DotPlot(imm, features = myeloid_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_myeloid_markers.png", path = "figure/Recluster/immune/annotation/", width = 20, height = 30, units = "cm")

DotPlot(imm, features = tcell_nk_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_T_NK_markers.png", path = "figure/Recluster/immune/annotation/", width = 20, height = 20, units = "cm")

DotPlot(imm, features = bcell_plasma_mast_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_B_Plasma_markers.png", path = "figure/Recluster/immune/annotation/", width = 20, height = 20, units = "cm")

DimPlot(imm, group.by = "SCT_snn_res.0.5", label = T, split.by = "tissue_type")
ggsave2("DimPlot_immune_clusters.png", path = "figure/Recluster/immune/annotation/", width =10, height =5)

DimPlot(imm, group.by = "patient_id", cols = use_colors)
ggsave2("DimPlot_immune_patients.png", path = "figure/Recluster/immune/annotation/", width = 10, height = 5)

