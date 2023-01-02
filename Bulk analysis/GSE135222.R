rm(list = ls())

library(GEOquery)
library(stringr)
library(data.table)
library(devtools)
library(org.Hs.eg.db)
library(clusterProfiler)

expr <- fread('GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz', data.table = F)

expr$gene_id <- gsub("\\..*", "",  expr$gene_id)
id <- bitr(expr$gene_id, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")

expr <- expr[expr$gene_id %in% id$ENSEMBL, ]
id <- id[id$ENSEMBL %in% expr$gene_id, ]
id <- id[order(id$ENSEMBL, decreasing = T), ] 
mutli <- id[duplicated(id$ENSEMBL), ]
for (i in unique(mutli$ENSEMBL)){
  id <- subset(id, subset = id$ENSEMBL != i)
}

expr <- expr[expr$gene_id %in% id$ENSEMBL, ]
id <- id[id$ENSEMBL %in% expr$gene_id, ]

row.names(expr) <- expr$gene_id
expr <- expr[, -1]
expr <- as.data.frame(expr)
id$median <- apply(expr, 1, median)
id <- id[order(id$SYMBOL, id$median, decreasing = T),] 
id <- id[!duplicated(id$SYMBOL),] 
row.names(id) <- id$ENSEMBL

expr <- expr[row.names(expr) %in% row.names(id), ]
id <- id[row.names(id) %in% row.names(expr), ]
id <- id[row.names(expr), ]

row.names(expr) <- id$SYMBOL

# clinic
gset <- getGEO('GSE135222', destdir=".", AnnotGPL = F, getGPL = F)
gset <- gset[[1]]
pheno <- pData(gset)

survival <- pheno[, 41:44]

colnames(survival) <- c('Age', 'Sex', 'time', 'status')
row.names(survival) <- gsub(' ', '', pheno$title)

expr <- as.data.frame(t(expr))
expr <- expr[row.names(expr) %in% row.names(survival), ]
survival <- survival[row.names(survival) %in% row.names(expr), ]
survival <- survival[row.names(expr), ]

features <- list(importance$gene)
tpm <- t(expr)

AddModuleScore_bulk <- function (object, features, nbin = 24, ctrl = 100, 
                                 k = FALSE, assay = NULL, name = "Cluster", seed = 1) 
{
  set.seed(seed = 1)
  
  features <- lapply(X = features, FUN = function(x) {
    missing.features <- setdiff(x = x, y = rownames(x = object))
    if (length(x = missing.features) > 0) {
      warning("The following features are not present in the object: ", 
              paste(missing.features, collapse = ", "), 
              ifelse(test = FALSE, yes = ", attempting to find updated synonyms", 
                     no = ", not searching for symbol synonyms"), 
              call. = FALSE, immediate. = TRUE)
    }
    return(intersect(x = x, y = rownames(x = object)))
  })
  
  cluster.length <- length(x = features)
  
  pool <- rownames(x = object)
  data.avg <- apply(object,1, mean)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- apply(object[features.use,], 2, mean)
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- object[features.use, , drop = FALSE]
    features.scores[i, ] <- apply(data.use, 2, mean)
  }
  features.scores.use <- features.scores - ctrl.scores
  
  colnames(features.scores.use) <- colnames(object)
  rownames(features.scores.use) <- name
  return(features.scores.use)
}

score <- AddModuleScore_bulk(tpm, features, name='score')
score <- as.data.frame(t(score))
identical(row.names(score), row.names(survival))

survival$score <- score$score
survival_matrix <- cbind(survival, expr)

# fwrite(survival_matrix, file = 'survival_matrix.txt', row.names = T)

# survival analysis
library(survival)
library(dplyr)
library(survminer)
library(ggplot2)

gene <- as.character(importance$gene)
survival <- survival_matrix[, c('Age', 'Sex', 'time', 'status', 'score', 'CD274', gene)]

survival$status <- as.numeric(survival$status)
survival$time <- as.numeric(survival$time)

num <- mean(survival$score)
survival$group <- ifelse(survival$score >= num, 'high', 'low')

# 创建生存模型
# status:1表示数据删失,2表示患者死亡
# sex:1表示男性,2表示女性
# 绘制生存曲线需要三个变量:time,status和分组变量
# time对应生存时间,status对应状态,如果用性别来分组,~sex是分组
# 把性别那一列的列名换掉group,data就是数据名
fit <- survfit(Surv(time, status) ~ group, data = survival)
fit

# pdf("1.pdf") # 保存为pdf
ggsurvplot(fit,
           pval = T, # 显示差异的p值
           conf.int = T, # 加置信区间
           conf.int.style = "ribbon", # 设置置信区间的样式
           conf.int.alpha = 0.1, # 置信区间透明度
           surv.median.line = "hv", # 增加中位生存时间
           palette = "jco", # 设置颜色模式,可选调色板有:"grey","npg","aaas","lancet","jco"等
           risk.table = T, # 添加风险表
)



