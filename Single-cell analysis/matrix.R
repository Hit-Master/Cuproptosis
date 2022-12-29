library(Seurat)
library(dplyr)
library(readxl)
library(data.table)


# epi_t <- subset(epi, subset = tissue_type == 'Tumor')

epi_matrix <- epi@assays$SCT@data
# epi_matrix_1 <- epi@assays$SCT
# epi_matrix_1 <- GetAssayData(object = epi_matrix_1, slot = "data")
# identical(epi_matrix, epi_matrix_1)

epi_matrix <- as.matrix(epi_matrix)
epi_matrix[1:10, 1:10]

rm(epi)

brown_list <- brown$symbol
brown_matrix <- epi_matrix[row.names(epi_matrix) %in% brown_list,]
brown_matrix <- as.data.frame(t(brown_matrix))

lightcyan_list <- lightcyan$symbol
lightcyan_matrix <- epi_matrix[row.names(epi_matrix) %in% lightcyan_list,]
lightcyan_matrix <- as.data.frame(t(lightcyan_matrix))

turquoise_list <- turquoise$symbol
turquoise_matrix <- epi_matrix[row.names(epi_matrix) %in% turquoise_list,]
turquoise_matrix <- as.data.frame(t(turquoise_matrix))

cnv <- fread('data/inferCNV_output/infercnv_clone_scores_nsclc.tsv')

lightcyan_matrix$group <- ifelse(row.names(lightcyan_matrix) %in% cnv$cell_id, 1, 2) # 1 T 2 N
brown_matrix$group <- ifelse(row.names(brown_matrix) %in% cnv$cell_id, 1, 2)
turquoise_matrix$group <- ifelse(row.names(turquoise_matrix) %in% cnv$cell_id, 1, 2)

prob <- data.frame(cell_id = row.names(brown_matrix), sample_id = 1:20582, group = brown_matrix$group)
row.names(lightcyan_matrix) <- 1:20582
row.names(brown_matrix) <- 1:20582
row.names(turquoise_matrix) <- 1:20582

brown_matrix <- brown_matrix[rowSums(brown_matrix[, -257])>0, ]
lightcyan_matrix <- lightcyan_matrix[rowSums(lightcyan_matrix[, -16])>0, ]
turquoise_matrix <- turquoise_matrix[rowSums(turquoise_matrix[, -143])>0, ]

save(prob, brown_matrix, lightcyan_matrix, turquoise_matrix, file = 'seurat_object/matrix/matrix.Rdata')

fwrite(brown_matrix, file = 'seurat_object/matrix/brown_matrix.txt', sep = '\t')
fwrite(lightcyan_matrix, file = 'seurat_object/matrix/lightcyan_matrix.txt', sep = '\t')
fwrite(turquoise_matrix, file = 'seurat_object/matrix/turquoise_matrix.txt', sep = '\t')
fwrite(prob, file = 'seurat_object/matrix/prob.txt', sep = '\t')
