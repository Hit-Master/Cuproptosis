# 计算指定表型的相关分析
pheno = 'stage'
pheno_type = phenotype[ , pheno, drop=F]
colnames(pheno_type) = pheno
# (1) Gene significance,GS: 即比较样本某个基因与对应表型的相关性
GS = as.data.frame(cor(wgcna_matrix, pheno_type, use = "p"))
colnames(GS) = "GS"

GS.p = as.data.frame(corPvalueStudent(as.matrix(GS), nrow(wgcna_matrix)))


# (2) Module Membership: 模块内基因表达与模块特征值的相关性
modNames = substring(names(MEs), 3)
# 计算3600个基因与18个模块的相关性
MM = as.data.frame(cor(wgcna_matrix, MEs, use = "p"))
colnames(MM) = paste("MM", modNames, sep="")

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(MM), nrow(wgcna_matrix)))
colnames(MMPvalue) = paste("p.MM", modNames, sep="");
MMPvalue[1:4,1:4]

# (3) 可视化blue模块基因特征
module = "turquoise"
moduleGenes = names(net$colors)[net$colors==module]
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(MM[moduleGenes, paste("MM", module, sep="")]),
                   abs(GS[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

