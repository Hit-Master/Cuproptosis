rm(list = ls())

# 指定候选值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(wgcna_matrix, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# (1)是否符合幂律分布
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
# (2)节点的平均连接度
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

par(mfrow = c(1,1))

save(sft, file = 'softThres.Rdata')
