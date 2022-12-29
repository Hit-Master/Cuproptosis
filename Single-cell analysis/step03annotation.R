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

epi <- readRDS("seurat_objects/epi.RDS")
imm <- readRDS("seurat_objects/imm.RDS")
str <- readRDS("seurat_objects/str.RDS")

### epithelial subclustering
epi <- RunPCA(epi)
ElbowPlot(epi,  ndims = 50)

epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  epi <- FindClusters(epi, resolution = i)
  print(DimPlot(epi, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

Idents(epi) <- epi@meta.data$SCT_snn_res.1


### immune subclustering
imm <- RunPCA(imm)
ElbowPlot(imm,  ndims = 50)

imm <- RunUMAP(imm, dims = 1:20)
imm <- FindNeighbors(imm, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  imm <- FindClusters(imm, resolution = i)
  print(DimPlot(imm, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}

Idents(imm) <- imm@meta.data$SCT_snn_res.0.5


### stromal sublustering
str <- RunPCA(str)
ElbowPlot(str, ndims = 50)

str <- RunUMAP(str, dims = 1:20)
str <- FindNeighbors(str, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  str <- FindClusters(str, resolution = i)
  print(DimPlot(str, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}

Idents(str) <- str@meta.data$SCT_snn_res.1
