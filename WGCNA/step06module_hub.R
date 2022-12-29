## (1) 模块内基因连接度
adjacency = adjacency(wgcna_matrix, power = 2)
TOM = TOMsimilarity(adjacency)
TOM[1:4,1:4]
Alldegrees <- intramodularConnectivity(adjacency, net$colors)
Alldegrees <- Alldegrees[order(Alldegrees$kTotal, Alldegrees$kWithin, decreasing = T),]
head(Alldegrees)
#                 kTotal    kWithin       kOut      kDiff
# MMT00000044  0.4092743  0.2862358  0.1230385  0.1631973
# MMT00000046 37.8927830 24.9652317 12.9275513 12.0376805
# MMT00000051 28.3866248 17.2076759 11.1789488  6.0287271
# MMT00000076  1.3015473  1.1992420  0.1023053  1.0969366
# MMT00000080 25.9713107 16.3954194  9.5758914  6.8195280
# MMT00000102 10.5051504  2.4713718  8.0337786 -5.5624067

# kTotal:基因在整个网络中的连接度
# kWithin: 基因在所属模块中的连接度,即Intramodular connectivity
# kOut: kTotal-kWithin
# kDiff: kIn-kOut

save(TOM, adjacency, Alldegrees, file = 'module_hub.Rdata')
#也可以绘制一个模块基因的Intramodular connectivity与Gene significance的散点图


## (2) Module Membership: 即上一节计算的模块内基因表达与模块特征值的相关性

## (3) Gene significance,GS: 即比较样本某个基因与对应表型的相关性

## (4) 综合上述指标,自定义合适阈值,选取一定数量的模块Hub基因