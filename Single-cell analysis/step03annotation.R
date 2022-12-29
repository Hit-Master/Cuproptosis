seu_obj=readRDS('scRNA_main_annotated.RDS')
metatable <- read_excel("./metadata/patients_metadata.xlsx")

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


Idents(seu_obj) <- seu_obj@meta.data$Main_cell_type
epi <- subset(seu_obj, idents = "Epithelial")
imm <- subset(seu_obj, idents = "Immune")
str <- subset(seu_obj, idents = "Stromal")

epi <- ScaleData(epi)
imm <- ScaleData(imm)
str <- ScaleData(str)

saveRDS(epi, file = "epi.RDS")
saveRDS(imm, file = "imm.RDS")
saveRDS(str, file = "str.RDS")
