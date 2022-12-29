library(data.table)

fpkm <- fread('TCGA-LUAD.htseq_fpkm.tsv', data.table = F)
annotation <- fread('gencode.v22.annotation.gene.probeMap', data.table = F)
row.names(fpkm) <- fpkm[,1]
fpkm <- fpkm[,-1]
fpkm <- 2^fpkm - 1
fpkm <- fpkm[rownames(fpkm) %in% annotation$id,] 

ids <- annotation[, 1:2]
ids <- ids[match(row.names(fpkm), ids$id),]
colnames(ids) <- c('probe_id', 'symbol')
ids$median <- apply(fpkm,1,median) 
ids <- ids[order(ids$symbol,ids$median,decreasing = T),] 
ids <- ids[!duplicated(ids$symbol),] 
fpkm <- fpkm[ids$probe_id,] 
rownames(fpkm) <- ids$symbol 

fpkm_to_tpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
tpm <- apply(fpkm,2,fpkm_to_tpm)

tpm <- as.data.frame(tpm)
 
count_norm <- tpm


library(stringr)

FDR <- 0.05 
FC <- 2
p_value <- 0.05

group_list <- ifelse(as.numeric(str_sub(colnames(count_norm), 14, 15))<10, 'tumor', 'normal')
group_list <- factor(t(group_list))

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
degs <- out[out$FDR<FDR, ]
degs <- degs[degs$pValues<p_value, ]
degs <- degs[abs(degs$log2foldChange)>FC, ]
wilcox_report <- list(out=out, degs=degs)


fwrite(count_norm, file = 'TPM.txt', sep = '\t', row.names = T)
fwrite(wilcox_report$out, file = 'Wilcoxon_test.txt', sep = '\t', row.names = T)
fwrite(wilcox_report$degs, file = 'DEG.txt', sep = '\t', row.names = T)
