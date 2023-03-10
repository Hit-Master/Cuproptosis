rm(list = ls())  

library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

load('seurat_object/Drug_Predict.Rdata')
tpm <- survival[, -c(1, 2)]
tpm$score <- NULL
tpm$group <- NULL
testExpr <- t(tpm)

th <- theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir <- 'data/DataFiles/DataFiles/Training Data/'
GDSC2_Expr <- readRDS(file = file.path(dir, 'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res <- readRDS(file = file.path(dir, "GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              printOutput = TRUE,
              removeLowVaringGenesFrom = 'rawData')

library(data.table)
testPtype <- fread('./Drug/GDSC2/DrugPredictions-low.csv', data.table = F)
testPtype1 <- fread('./Drug/GDSC2/DrugPredictions-high.csv', data.table = F)
testPtype$group <- 'low'
testPtype1$group <- 'high'
testPtype <- rbind(testPtype, testPtype1)
row.names(testPtype) <- testPtype$V1
testPtype <- testPtype[, -1]

# fwrite(testPtype, file = 'Drug/TCGA-drug-predict.csv', row.names = T)

drug <- testPtype
drug <- t(drug)
drug <- drug[complete.cases(drug),]
drug <- as.data.frame(t(drug))

count_norm <- drug
count_norm$group <- NULL
count_norm <- apply(count_norm, 2 ,as.numeric)
row.names(count_norm) <- row.names(drug)
count_norm <- t(count_norm)
expr <- count_norm

col_vector <- colSums(expr)*length(colnames(expr))/sum(expr)
expr <- sweep(expr, MARGIN = 2, STATS = col_vector, FUN = "/")
count_norm <- expr

library(stringr)

group_list <- drug$group
group_list <- factor(t(group_list), levels = c('low', 'high'))

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])), group_list)
  p=wilcox.test(gene~group_list, data)$p.value
  return(p)
})
fdr <- p.adjust(pvalues, method = "fdr")

conditionsLevel <- levels(group_list)
dataCon1 <- count_norm[,c(which(group_list==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(group_list==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))

out <- data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
out <- na.omit(out)

degs <- out[out$FDR<0.05, ]
degs <- degs[degs$pValues<0.05, ]

fwrite(degs, file = 'Drug/GDSC2/Drug-Predict-DEG.csv', row.names = T)

# plot
degs <- degs[order(degs$log2foldChange, decreasing = T), ]
df <- count_norm[c(row.names(degs)[1:6]), ]
df <- as.data.frame(t(df))
df <- melt(df)
df$group <- rep(drug$group, 6)

ggplot(df, aes(group, log10(value+1)))+
  geom_boxplot(aes(color=group), fill="white", outlier.size = 0.1) +
  scale_color_manual(values = c("#982b2b", "#0074b3")) +
  scale_fill_manual(values = c("#982b2b","#0074b3")) +
  xlab("") +
  ylab("") +
  ylim(7, 9) +
  theme(legend.position = 'none') +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  facet_wrap(.~variable, nrow = 2) +
  stat_compare_means(label = 'p.signif', label.x = 1.4, label.y = 14, color = '#DE6757', hide.ns = TRUE)


# pearson
load('seurat_object/score.Rdata')
drug <- as.data.frame(t(count_norm))

score$num <- 1
score <- score[row.names(drug), ]
drug$score <- score$score

cor_score <- cor(drug, method = 'pearson')
cor_score <- as.data.frame(cor_score)
cor_degs <- subset(cor_score, subset = cor_score$score < 0)
cor_degs <- cor_degs[, c('score', 'Cytarabine_1006')]
cor_degs <- cor_degs[order(cor_degs$score),]
