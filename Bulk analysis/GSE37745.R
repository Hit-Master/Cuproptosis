rm(list = ls())

library(GEOquery)
library(stringr)
library(data.table)
library(ggplot2)

gset <- getGEO('GSE37745', destdir=".", AnnotGPL = F, getGPL = F)
gset <- gset[[1]]
expr <- exprs(gset)

anno <- fread('GPL570.annot.gz',sep = '\t', header = T, data.table = F)
anno <- anno[, c('ID', 'Gene symbol')]
gene <- strsplit(anno$`Gene symbol`, split = '///', fixed = T)  
gene = sapply(gene, function(x){x[1]})
id  <- data.frame(anno$ID, gene)
row.names(id) <- id[, 1]

expr <- expr[row.names(expr) %in% row.names(id), ]
id <- id[row.names(id) %in% row.names(expr), ]
id <- id[row.names(expr), ]

id$median <- apply(expr, 1, median)
id <- id[order(id$gene, id$median, decreasing = T),] 
id <- id[!duplicated(id$gene),] 

expr <- expr[row.names(expr) %in% row.names(id), ]
id <- id[row.names(id) %in% row.names(expr), ]
id <- id[row.names(expr), ]

row.names(expr) <- id$gene

load("~/RstudioProjects/铜死亡/importance.Rdata")
features <- list(importance$gene)
tpm <- as.matrix(expr)

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
# Warning: The following features are not present in the object: AK4, not searching for symbol synonyms

# clinic
pheno <- pData(gset)
survival <- pheno[, c(10, 11, 13, 15, 16)]
colnames(survival) <- c('status', 'time', 'Stage', 'Age', 'Sex')

# fwrite(survival, file = 'survival.csv', sep = ',', row.names = T)

survival <- fread('survival.csv', header = T, data.table = F)
row.names(survival) <- survival$V1
survival <- survival[, -1]

survival$Sex <- gsub(survival$Sex, pattern = 'gender: ', replacement = '')
survival$Stage <- gsub(survival$Stage, pattern = 'tumor stage: ', replacement = '')
survival$status <- gsub(survival$status, pattern = 'dead: ', replacement = '')
survival$time <- gsub(survival$time, pattern = 'days to determined death status: ', replacement = '')
survival$Age <- gsub(survival$Age, pattern = 'age: ', replacement = '')

survival$Sex <- ifelse(survival$Sex == 'female', 'F', 'M')
survival$status <- ifelse(survival$status == 'no', 0, 1)

expr <- as.data.frame(t(expr))
survival <- survival[row.names(expr), ]

identical(row.names(expr), row.names(survival))
identical(row.names(score), row.names(survival))

survival$score <- score$score
survival <- cbind(survival, expr)

# fwrite(survival, file = 'survival_matrix.txt', row.names = T)

# survival analysis
library(survival)
library(dplyr)
library(survminer)

survival$time <- as.numeric(survival$time)

num <- median(survival$CDKN2A)
survival$CDKN2A <- ifelse(survival$CDKN2A >= num, 'high', 'low')

num <- median(survival$DLD)
survival$DLD <- ifelse(survival$DLD >= num, 'high', 'low')

num <- median(survival$FDX1)
survival$FDX1 <- ifelse(survival$FDX1 >= num, 'high', 'low')

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


# estimate
library(estimate)

tpm <- as.matrix(tpm)
estimate <- function(dat, project){
  library(estimate)
  
  input.f <- paste0(project,'_estimate_input.txt')
  output.f <- paste0(project,'_estimate_gene.gct')
  output.ds <- paste0(project,'_estimate_score.gct')
  
  write.table(dat, file = input.f, sep = '\t', quote = F)
  
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds=output.f,
                output.ds=output.ds,
                platform="affymetrix"   ## 注意这个platform参数
  )
  scores <- read.table(output.ds, skip = 2, header = T)
  rownames(scores) <- scores[,1]
  scores <- t(scores[,3:ncol(scores)])
  
  return(scores)
}
project <- 'lung'
scores <- estimate(tpm, project) 

scores <- as.data.frame(scores)
identical(row.names(score), row.names(scores))

scores$Cu_score <- score$score
scores$group <- survival$group

fwrite(scores, file = 'estimate_score.txt', row.names = T)

ggplot(scores) +
  geom_point(size=1.5, aes(x=Cu_score, y=StromalScore, color=group)) +
  scale_colour_manual(values=c('#5B9BD5', 'grey60')) + 
  theme_bw() +  
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  labs(x='Cu_score', y='StromalScore', title = '') + 
  theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +
  theme(legend.position = 'none') +
  geom_smooth(method = lm, aes(x=Cu_score, y=StromalScore), level=0.99 , colour="#DE6757") +
  stat_cor(method = 'pearson', aes(x=Cu_score, y=StromalScore), color="#5B9BD5", label.x = 0.2)


ggplot(scores) +
  geom_point(size=1.5, aes(x=Cu_score, y=ImmuneScore, color=group)) +
  scale_colour_manual(values=c('#5B9BD5', 'grey60')) + 
  theme_bw() +  
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  labs(x='Cu_score', y='ImmuneScore', title = '') + 
  theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +
  theme(legend.position = 'none') +
  geom_smooth(method = lm, aes(x=Cu_score, y=ImmuneScore), level=0.99 , colour="#DE6757") +
  stat_cor(method = 'pearson', aes(x=Cu_score, y=ImmuneScore), color="#5B9BD5", label.x = 0.2)

ggplot(scores) +
  geom_point(size=1.5, aes(x=Cu_score, y=TumorPurity, color=group)) +
  scale_colour_manual(values=c('#5B9BD5', 'grey60')) + 
  theme_bw() +  
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  labs(x='Cu_score', y='TumorPurity', title = '') + 
  theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +
  theme(legend.position = 'none') +
  geom_smooth(method = lm, aes(x=Cu_score, y=TumorPurity), level=0.99 , colour="#DE6757") +
  stat_cor(method = 'pearson', aes(x=Cu_score, y=TumorPurity), color="#5B9BD5", label.x = 0.2)

# PD-L1 CD274
gene <- c('score', 'group', 'CD274')
df <- survival[, gene]

ggplot(df) +
  geom_point(size=1.5, aes(x=score, y=CD274, color=group)) +
  scale_colour_manual(values=c('#5B9BD5', 'grey60')) + 
  theme_bw() +  
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  labs(x='Cu_score', y='PD-L1', title = '') + 
  theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +
  # theme(legend.position = 'none') +
  geom_smooth(method = lm, aes(x=score, y=CD274), level=0.99 , colour="#DE6757") +
  stat_cor(method = 'pearson', aes(x=score, y=CD274), color="#5B9BD5", label.x = 0)

