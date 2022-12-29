rm(list = ls())

library(Seurat)
library(cowplot)
library(ggplot2)
library(viridis)
library(readxl)
library(clustree)
library(ROGUE)

seu_obj <- readRDS('seurat_object/Preprocessing/SCTransform.RDS')

# PCA降维
seu_obj <- RunPCA(seu_obj)

# Determine percent of variation associated with each PC
pct <- seu_obj[["pca"]]@stdev / sum( seu_obj [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
ggsave2("Elbow.pdf", path = "figure/Preprocessing/resolution/", width = 10, height = 5)

# PCA热图
DimHeatmap(seu_obj, dims = 13:18, cells = 500, balanced = TRUE)


# for (i in c(15, 20)) {
#  umaptest <- RunUMAP(seu_obj, dims = 1:i, verbose = F)
#  print(DimPlot(umaptest, reduction = "umap", group.by = "orig.ident") + labs(title = paste0(i, " dimensions")))
#  print(FeaturePlot(umaptest, features = c("EPCAM", "PTPRC"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("MARCO", "KIT"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("FOXJ1", "AGER"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("JCHAIN", "VWF"), sort.cell = T))
#  remove(umaptest)
# }

seu_obj <- RunUMAP(seu_obj, dims = 1:15, verbose = T)
# seu_obj <- RunTSNE(seu_obj, dims = 1:15, verbose = T)
seu_obj <- FindNeighbors(seu_obj, dims = 1:15)

for (i in c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2)) {
  seu_obj <- FindClusters(seu_obj, resolution = i)
  DimPlot(seu_obj, reduction = "umap") + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("umap_resolution", i,".pdf"), path = "figure/Preprocessing/resolution/", width = 10, height = 10, units = "cm")
}

# 判断resolution值
clustree(seu_obj, prefix = 'SCT_snn_res.') + coord_flip()
ggsave2("clustree.pdf", path = "figure/Preprocessing/resolution/", width = 15, height = 20, units = "cm")

# ROGUE评分
sce_sub <- seu_obj %>% subset(., downsample = 100)
expr <- GetAssayData(sce_sub, slot = 'counts') %>% as.matrix()
expr <- expr[1:9000,]
meta <- sce_sub@meta.data
rogue(expr = expr, 
      labels = meta$SCT_snn_res.0.2, 
      samples = meta$orig.ident, 
      platform = 'UMI',
      span = 0.6) %>% rogue.boxplot()
ggsave2("ROGUE.pdf", path = "figure/Preprocessing/resolution/", width = 15, height = 10, units = "cm")


# 主要的细胞类型注释——区分上皮、间质、免疫细胞

# 选定marker,上皮marker如EPCAM,间质如VWF,免疫如PTPRC,KIT是间叶特征marker
mainmarkers <- c("PECAM1", "VWF", "ACTA2", "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", "EPCAM", "CDH1", "KRT7", "KRT19")

# 对每个marker基因作图并观察,这里建议画umap更清楚一些
for (i in seq_along(mainmarkers)) {
  FeaturePlot(seu_obj, features = mainmarkers[i], coord.fixed = T, order = T, cols = cividis(10))
  ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"), path = "figure/Preprocessing/annotation/", width = 10, height = 10, units = "cm")
}

# 我们就选择0.2,共19个cluster,resolution值选择应该稍大一点

DotPlot(seu_obj, features = mainmarkers, group.by = "SCT_snn_res.0.2") + 
  coord_flip() + 
  scale_color_binned()
ggsave2("DotPlot_mainmarkers.png", path = "figure/Preprocessing/annotation/", width = 10, height = 5)

DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)
ggsave2("DimPlot_all_clusters.png", path = "figure/Preprocessing/annotation/", width = 5, height = 5)


# 同时考虑到Umap的接近程度和marker的表达,确定细胞类型

Idents(seu_obj) <- seu_obj$SCT_snn_res.0.2
annotation_curated_main <- read_excel("data/annotation/annotation_main.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seu_obj)
# 加入seurat对象中
seu_obj <- RenameIdents(seu_obj, new_ids_main)
seu_obj@meta.data$main_cell_type <- Idents(seu_obj)

DimPlot(seu_obj, label = T, label.size = 5, group.by="main_cell_type",
        cols= c("#459943", "#db6968", "#88c4e8"))
ggsave2("DimPlot_main_cell_type.pdf", path = "figure/Preprocessing/annotation/", width = 7, height = 6)

# saveRDS(seu_obj, file = "seurat_object/Preprocessing/scRNA_main_annotation.RDS")

