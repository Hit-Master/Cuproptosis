count_max_tpm <- function(fpkm_file='', annotation_file=''){
  library(data.table)
  
  # 读取 UCSC xena 下载的FPKM矩阵,反log转换成FPKM值
  # 读取注释文件(Ensembl ID, Symbol, 染色体起始终止位置及正负链)
  # 删除FPKM矩阵中没有注释的探针和其余信息
  fpkm <- fread(fpkm_file, data.table = F)
  annotation <- fread(annotation_file, data.table = F)
  row.names(fpkm) <- fpkm[,1]
  fpkm <- fpkm[,-1]
  fpkm <- 2^fpkm - 1
  fpkm <- fpkm[rownames(fpkm) %in% annotation$id,] # 删除FPKM矩阵中没有注释的探针和其余信息
  
  # 读取注释文件(Ensembl ID, Symbol, 染色体起始终止位置及正负链)
  # 选取同一symbol最大表达的FPKM矩阵,并注释到基因名
  ids <- annotation[, 1:2]
  ids <- ids[match(row.names(fpkm), ids$id),]
  colnames(ids) <- c('probe_id', 'symbol')
  ids$median <- apply(fpkm,1,median) # ids新建median这一列，列名为median，同时对FPKM这个矩阵按行操作，取每一行的中位数,将结果给到median这一列的每一行
  ids <- ids[order(ids$symbol,ids$median,decreasing = T),] # 对ids$symbol按照ids$median中位数从大到小排列的顺序排序,将对应的行赋值为一个新的ids
  ids <- ids[!duplicated(ids$symbol),] # 将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ,保留每个基因最大表达量结果
  fpkm <- fpkm[ids$probe_id,] # 新的ids取出probe_id这一列,将FPKM按照取出的这一列中的每一行组成一个新的FPKM
  rownames(fpkm) <- ids$symbol # 把ids的symbol这一列中的每一行给FPKM作为FPKM的行名
  
  # 转换成TPM标准化的count矩阵
  fpkm_to_tpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
  tpm <- apply(fpkm,2,fpkm_to_tpm)
  
  tpm <- as.data.frame(tpm)
  
  return(tpm)
}


count_mean_tpm <- function(fpkm_file='', annotation_file=''){
  library(data.table)
  library(tibble)
  library(dplyr)
  library(tidyr)
  
  # 从UCSC Xena下载的fpkm数据格式是log2(FPKM+1)
  fpkm_matrix <- fread(fpkm_file, data.table = F)
  fpkm <- fread(annotation_file, data.table = F)%>%
    select(id, gene)%>%
    inner_join(fpkm_matrix, by=c('id'='Ensembl_ID'))%>%
    select(-id)%>%
    group_by(gene)%>% 
    summarise_all(mean)%>% # 重复基因取均值
    column_to_rownames('gene')
  
  # 反log转换成FPKM,再转换成TPM
  fpkm <- 2^fpkm-1
  fpkm_to_tpm <- function(fpkm){exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
  tpm <- apply(fpkm,2,fpkm_to_tpm) 
  
  tpm <- as.data.frame(tpm)
  
  return(tpm)
}


count_Wilcoxon_test <- function(count_norm, FDR=0.05, FC=2, p_value=0.05){
  library(stringr)
  
  # 确定分组信息
  group_list <- ifelse(as.numeric(str_sub(colnames(count_norm), 14, 15))<10, 'tumor', 'normal')
  group_list <- factor(t(group_list))
  
  # Wilcoxon rank-sum test 
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])), group_list)
    p=wilcox.test(gene~group_list, data)$p.value
    return(p)
  })
  fdr <- p.adjust(pvalues, method = "fdr")
  
  # 计算Fold Change
  conditionsLevel <- levels(group_list)
  dataCon1 <- count_norm[,c(which(group_list==conditionsLevel[1]))]
  dataCon2 <- count_norm[,c(which(group_list==conditionsLevel[2]))]
  foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  
  # 输出FDR阈值的结果
  out <- data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  out <- na.omit(out)
  degs <- out[out$FDR<FDR, ]
  degs <- degs[degs$pValues<p_value, ]
  degs <- degs[abs(degs$log2foldChange)>FC, ]
  wilcox_report <- list(out=out, degs=degs)
  
  return(wilcox_report)
}


count_norm <- count_max_tpm(fpkm_file='TCGA-LUAD.htseq_fpkm.tsv', annotation_file='gencode.v22.annotation.gene.probeMap')
wilcox_report <- count_Wilcoxon_test(count_norm)

fwrite(count_norm, file = 'TPM.txt', sep = '\t', row.names = T)
fwrite(wilcox_report$out, file = 'Wilcoxon_test.txt', sep = '\t', row.names = T)
fwrite(wilcox_report$degs, file = 'DEG.txt', sep = '\t', row.names = T)
