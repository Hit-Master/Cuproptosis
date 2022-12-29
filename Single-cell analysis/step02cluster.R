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
library(progeny)
library(readr)
library(stringr)
seu_obj=readRDS('scRNA_SCTransform.RDS')

seu_obj <- RunPCA(seu_obj)
ElbowPlot(seu_obj, ndims = 50)

seu_obj <- RunUMAP(seu_obj, dims = 1:15, verbose = T)
seu_obj <- RunTSNE(seu_obj, dims = 1:15, verbose = T)
seu_obj <- FindNeighbors(seu_obj, dims = 1:15)
library(ggplot2)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  seu_obj <- FindClusters(seu_obj, resolution = i)
  print(DimPlot(seu_obj, reduction = "tsne") + labs(title = paste0("resolution: ", i)))
}

DimPlot(seu_obj, group.by = "orig.ident",reduction = 'tsne')
DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T,reduction = 'umap')

mainmarkers <- c("PECAM1", "VWF", "ACTA2", "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", "EPCAM", "CDH1", "KRT7", "KRT19")
dir.create('annotation')
dev.off()

for (i in seq_along(mainmarkers)) {
  FeaturePlot(seu_obj, features = mainmarkers[i], coord.fixed = T, order = T, cols = cividis(10))
  ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"), path = "./annotation", width = 10, height = 10, units = "cm")
}

DotPlot(seu_obj, features = mainmarkers, group.by = "SCT_snn_res.0.2") + 
  coord_flip() + 
  scale_color_binned()
ggsave2("DotPlot_mainmarkers.png", path = "./annotation", width = 10, height = 5)

DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)
ggsave2("DimPlot_all_clusters.png", path = "./annotation", width = 5, height = 5)

Idents(seu_obj) <- seu_obj$SCT_snn_res.0.2
annotation_curated_main <- read_excel("./annotation/annotation_main.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new_ids_main)
seu_obj@meta.data$Main_cell_type <- Idents(seu_obj)

DimPlot(seu_obj, label = T, label.size = 5,group.by="Main_cell_type" ,
       cols= c("seagreen","darkgoldenrod2","steelblue"))
ggsave2("DimPlot_main_cell_type.pdf", path = "./annotation", width = 7, height = 6)
saveRDS(seu_obj, file = "scRNA_main_annotated.RDS")
