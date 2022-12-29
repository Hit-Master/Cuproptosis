seu_obj <- RunPCA(seu_obj)
ElbowPlot(seu_obj, ndims = 50)

seu_obj <- RunUMAP(seu_obj, dims = 1:15, verbose = T)
seu_obj <- FindNeighbors(seu_obj, dims = 1:15)

for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  seu_obj <- FindClusters(seu_obj, resolution = i)
  print(DimPlot(seu_obj, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}

for (i in c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")) {
  print(FeaturePlot(seu_obj, features = i, coord.fixed = T, sort.cell = T))
}

mainmarkers <- c("PECAM1", "VWF", "ACTA2", "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", "EPCAM", "CDH1", "KRT7", "KRT19")

for (i in seq_along(mainmarkers)) {
  FeaturePlot(seu_obj, features = mainmarkers[i], coord.fixed = T, order = T, cols = viridis(10))
  #ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"), path = "output/annotation", width = 10, height = 10, units = "cm")
}

DotPlot(seu_obj, features = mainmarkers, group.by = "SCT_snn_res.0.2") + 
  coord_flip() + 
  scale_color_viridis()

DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)

Idents(seu_obj) <- seu_obj$SCT_snn_res.0.2
annotation_curated_main <- read_excel("../data/curated_annotation/curated_annotation_main.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new_ids_main)
seu_obj@meta.data$main_cell_type <- Idents(seu_obj)

metatable <- read_excel("../data/metadata/patients_metadata.xlsx")

metadata <- FetchData(seu_obj, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$sample_id <- metadata$orig.ident
metadata <- left_join(x = metadata, y = metatable, by = "sample_id")
rownames(metadata) <- metadata$cell_id

seu_obj <- AddMetaData(seu_obj, metadata = metadata)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

score_cc <- function(seu_obj) {
  seu_obj <- CellCycleScoring(seu_obj, s.genes, g2m.genes)
  seu_obj@meta.data$CC.Diff <- seu_obj@meta.data$S.Score - seu_obj@meta.data$G2M.Score
  return(seu_obj)
}

seu_obj <- score_cc(seu_obj)

FeatureScatter(seu_obj, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
  coord_fixed(ratio = 1)

Idents(seu_obj) <- seu_obj@meta.data$main_cell_type

epi <- subset(seu_obj, idents = "Epithelial")
imm <- subset(seu_obj, idents = "Immune")
str <- subset(seu_obj, idents = "Stromal")

epi <- ScaleData(epi)
imm <- ScaleData(imm)
str <- ScaleData(str)

saveRDS(epi, file = "seurat_objects/epi.RDS")
saveRDS(imm, file = "seurat_objects/imm.RDS")
saveRDS(str, file = "seurat_objects/str.RDS")
