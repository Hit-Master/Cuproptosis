# (1)计算模块特征值
# MEs = moduleEigengenes(exp_dat, net$colors)$eigengenes
MEs = net$MEs
MEs = orderMEs(MEs)

# (2)计算18个module与25个表型的相关性以及对应的P值
moduleTraitCor = cor(MEs, phenotype, use = "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs))


# (3)可视化相关性与P值
# sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(phenotype),
               yLabels = gsub('ME', '',names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

library(stringr)
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

# 不需要重新计算,改下列名就行
MEs_col = MEs
colnames(MEs_col) = labels2colors(
  str_replace_all(colnames(MEs),"ME",""))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(1,6,2,6),
                      marHeatmap = c(6,6,1,2), plotDendrograms = T,
                      xLabelsAngle = 90)

