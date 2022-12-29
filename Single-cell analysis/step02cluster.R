rm(list = ls())

library(Seurat)
library(cowplot)
library(ggplot2)
library(viridis)
library(readxl)
library(clustree)
library(ROGUE)

seu_obj <- readRDS('seurat_object/Preprocessing/SCTransform.RDS')

seu_obj <- RunPCA(seu_obj)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
ggsave2("Elbow.pdf", path = "figure/Preprocessing/resolution/", width = 10, height = 5)
      
seu_obj <- RunUMAP(seu_obj, dims = 1:15, verbose = T)
# seu_obj <- RunTSNE(seu_obj, dims = 1:15, verbose = T)
seu_obj <- FindNeighbors(seu_obj, dims = 1:15)

for (i in c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2)) {
  seu_obj <- FindClusters(seu_obj, resolution = i)
  DimPlot(seu_obj, reduction = "umap") + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("umap_resolution", i,".pdf"), path = "figure/Preprocessing/resolution/", width = 10, height = 10, units = "cm")
}
  
mainmarkers <- c("PECAM1", "VWF", "ACTA2", "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", "EPCAM", "CDH1", "KRT7", "KRT19")

for (i in seq_along(mainmarkers)) {
  FeaturePlot(seu_obj, features = mainmarkers[i], coord.fixed = T, order = T, cols = cividis(10))
  ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"), path = "figure/Preprocessing/annotation/", width = 10, height = 10, units = "cm")
}

DotPlot(seu_obj, features = mainmarkers, group.by = "SCT_snn_res.0.2") + 
  coord_flip() + 
  scale_color_binned()
ggsave2("DotPlot_mainmarkers.png", path = "figure/Preprocessing/annotation/", width = 10, height = 5)

DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)
ggsave2("DimPlot_all_clusters.png", path = "figure/Preprocessing/annotation/", width = 5, height = 5)

Idents(seu_obj) <- seu_obj$SCT_snn_res.0.2
annotation_curated_main <- read_excel("data/annotation/annotation_main.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seu_obj)

seu_obj <- RenameIdents(seu_obj, new_ids_main)
seu_obj@meta.data$main_cell_type <- Idents(seu_obj)

DimPlot(seu_obj, label = T, label.size = 5, group.by="main_cell_type",
        cols= c("#459943", "#db6968", "#88c4e8"))
ggsave2("DimPlot_main_cell_type.pdf", path = "figure/Preprocessing/annotation/", width = 7, height = 6)

# saveRDS(seu_obj, file = "seurat_object/Preprocessing/scRNA_main_annotation.RDS")
