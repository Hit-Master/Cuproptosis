# Habermann et al.
# https://www.biorxiv.org/content/10.1101/753806v1

habermann_epi <- c("ABCA3", "SFTPB", "SFTPC", "AGER", "PDPN",  "KRT5", "TRP63", "NGFR", "SCGB1A1", "MUC5B", "KRT17", "FOXJ1", "TMEM190", "CAPS", "CHGA", "CALCA", "ASCL1", "PTPRC", "EPCAM")

habermann_imm <- c("CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", "LILRA4", "CLEC4C", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", "CD79A", "CPA3", "KIT", "MKI67", "CDK1", "EPCAM")

habermann_oth <- c("VWF", "PECAM1", "CCL21", "PROX1", "ACTA2", "MYH11", "PDGFRB", "WT1", "UPK3B", "LUM", "PDGFRA", "MYLK", "HAS1", "PLIN2", "FAP", "PTPRC", "EPCAM")


# Epithelial genes according to Habermann et al.
# manual annotation of normal cell types & detection of immune cell contaminated clusters
DotPlot(subset(epi, subset = cluster_type == "Normal"), features = habermann_epi, group.by = "SCT_snn_res.1") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave2(filename = 'DotPlot_Normal_epi_markers.pdf',path = "figure/Recluster/epithelial/annotation", width =10, height =15)
# manual detection of immune cell contaminated clusters
DotPlot(subset(epi, subset = cluster_type == "Tumor"), features = habermann_epi, group.by = "SCT_snn_res.1") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave2(filename = 'DotPlot_Tumor_epi_markers.pdf',path = "figure/Recluster/epithelial/annotation", width =10, height =15)

# Immune genes according to Habermann et al.
DotPlot(imm, features = habermann_imm, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2(filename = 'DotPlot_imm_markers.pdf',path = "figure/Recluster/immune/annotation", width =10, height =15)

# for (i in seq_along(habermann_imm)) {
#  plotlist <- list()
#  plotlist[1] <- FeaturePlot(imm, features = habermann_imm[i], sort.cell = T, combine = F)
#  plotlist[2] <- VlnPlot(imm, features = habermann_imm[i], pt.size = 0, combine = F)
#  print(CombinePlots(plots = plotlist))
# }

# Stromal genes according to Habermann et al.
DotPlot(str, features = habermann_oth, group.by = "SCT_snn_res.1") + coord_flip()
ggsave2(filename = 'DotPlot_str_markers.pdf',path = "figure/Recluster/stromal/annotation", width =10, height =15)

# for (i in seq_along(habermann_oth)) {
#  plotlist <- list()
#  plotlist[1] <- FeaturePlot(str, features = habermann_oth[i], sort.cell = T, combine = F)
#  plotlist[2] <- VlnPlot(str, features = habermann_oth[i], pt.size = 0, combine = F)
#  print(CombinePlots(plots = plotlist, ncol = 3))
# }

# cell type annotation and subsetting
# epithelial
annotation_curated_epi <- read_excel("data/annotation/annotation_epi.xlsx")
epi_anno <- epi
new_ids_epi <- annotation_curated_epi$cell_type_epi
names(new_ids_epi) <- levels(epi_anno)
epi_anno <- RenameIdents(epi_anno, new_ids_epi)
epi_anno@meta.data$cell_type_epi <- Idents(epi_anno)

epi_anno <- subset(epi_anno, subset = cell_type_epi != "Immune_contamination")
epi_anno <- ScaleData(epi_anno)

# immune
annotation_curated_imm <- read_excel("data/annotation/annotation_imm.xlsx")
imm_anno <- imm
new_ids_imm <- annotation_curated_imm$cell_type_imm
names(new_ids_imm) <- levels(imm_anno)
imm_anno <- RenameIdents(imm_anno, new_ids_imm)
imm_anno@meta.data$cell_type_imm <- Idents(imm_anno)

imm_anno <- subset(imm_anno, subset = cell_type_imm != "Epithelial_contamination")
imm_anno <- ScaleData(imm_anno)

# stromal
annotation_curated_str <- read_excel("data/annotation/annotation_str.xlsx")
str_anno <- str
new_ids_str <- annotation_curated_str$cell_type_str
names(new_ids_str) <- levels(str_anno)
str_anno <- RenameIdents(str_anno, new_ids_str)
str_anno@meta.data$cell_type_str <- Idents(str_anno)

str_anno <- subset(str_anno, subset = cell_type_str != "Immune/Epithelial contamination")
str_anno <- ScaleData(str_anno)

saveRDS(epi_anno, file = "seurat_object/epi_anno.RDS")
saveRDS(imm_anno, file = "seurat_object/imm_anno.RDS")
saveRDS(str_anno, file = "seurat_object/str_anno.RDS")
