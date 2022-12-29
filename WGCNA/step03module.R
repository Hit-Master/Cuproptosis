power = sft$powerEstimate

net = blockwiseModules(wgcna_matrix, power = 2,
                       corType = "pearson",
                       networkType="unsigned",
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       mergeCutHeight = 0.25,
                       verbose = 3)

# 最终得到的网络模块(合并之后)
unique(net$colors)
table(net$colors)
# grey module是 unassigned gene

# 合并之前的网络模块
unique(net$unmergedColors)
table(net$unmergedColors)

# 样本对于模块的特征值
dim(net$MEs)
net$MEs[1:4,1:4]
net$MEsOK


# 灰色的为**未分类**到模块的基因
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意,还可以recutBlockwiseTrees
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

color <- as.data.frame(net$colors)
color <- cbind(color, row.names(color))
colnames(color) <- c('color', 'symbol')
gene_set <- color[gene_set, ]

brown <- subset(color, subset = color == 'brown')
turquoise <- subset(color, subset = color == 'turquoise')
lightcyan <- subset(color, subset = color == 'lightcyan')

save(color, brown, turquoise, lightcyan, gene_set, file = 'Gene_set.Rdata')

save(net, file = 'Net.Rdata')
