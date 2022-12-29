rm(list = ls())

library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(org.Hs.eg.db)
library(enrichplot)

# load scRNA data
scRNA <- readRDS("~/Cu/codeOcean/seurat_object/Preprocessing/SCTransform.RDS")
scRNA <- subset(scRNA, subset = orig.ident == 'p018t' | orig.ident == 'p019t' | orig.ident == 'p023t' | 
                                orig.ident == 'p024t' | orig.ident == 'p027t' | orig.ident == 'p030t' |
                                orig.ident == 'p031t' | orig.ident == 'p032t' | orig.ident == 'p033t' |
                                orig.ident == 'p034t')

scRNA$group <- recode(scRNA$orig.ident, 
                      'p018t' = 'low',
                      'p024t' = 'low',
                      'p030t' = 'low',
                      'p031t' = 'low',
                      'p032t' = 'low',
                      'p019t' = 'high',
                      'p023t' = 'high',
                      'p027t' = 'high',
                      'p033t' = 'high',
                      'p034t' = 'high',)

gene <- FindMarkers(scRNA, ident.1 = 'high', ident.2 = 'low',
                   group.by = 'group', logfc.threshold = 0, min.pct = 0)

deg <- subset(gene, p_val_adj<0.05&abs(avg_log2FC)>0.15)

# GO
ego <- enrichGO(gene          = row.names(deg),
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

ego_results <- as.data.frame(ego)

# KEGG
eg = bitr(row.names(deg), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

kegg <- enrichKEGG(
  gene = eg$ENTREZID,
  keyType = "kegg",
  organism  = 'human',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05
)

kegg_results <- as.data.frame(kegg)

save(gene, ego_results, kegg_results, file = 'seurat_object/GO_KEGG.Rdata')


# plot
go_list <- by(ego_results, ego_results$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
go_enrich_df <- rbind(go_list$BP, go_list$CC, go_list$MF)

shorten_names <- function(x, n_word=10, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char))
  {
    if (nchar(x) > n_char) x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
labels=(sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names))

# 排序
names(labels) = rev(1:nrow(go_enrich_df))
CPCOLS <- c("#66C2A5", "#FC8D62", "#8DA0CB")

library(ggplot2)

ggplot(data=go_enrich_df, aes(x=number, y=-log2(qvalue), fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_test() + 
  scale_x_discrete(labels=labels) +
  xlab("") + ylab("-log2(qvalue)") + 
  theme(axis.text=element_text(face = "plain"))+
  theme(axis.text.y=element_text(color=rep(rev(CPCOLS),rep(10,3)))) + 
  theme(legend.title=element_blank())

kegg_enrich_df <- kegg_results[c('hsa03010', 'hsa04612', 'hsa05012', 'hsa04210', 'hsa05010',
                                 'hsa04660', 'hsa04650', 'hsa05235', 'hsa04672', 'hsa04218'), ]

kegg_enrich_df <- kegg_enrich_df[order(kegg_enrich_df$qvalue), ]
kegg_enrich_df$Description <- factor(kegg_enrich_df$Description, levels = rev(kegg_enrich_df$Description))

shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char))
  {
    if (nchar(x) > n_char) x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

kegg_enrich_df$number <- factor(rev(1:nrow(kegg_enrich_df))) # 生成倒序的排序
labels=(sapply(
  levels(kegg_enrich_df$Description)[as.numeric(kegg_enrich_df$Description)],
  shorten_names))

# 排序
names(labels) = rev(1:nrow(kegg_enrich_df))

#绘制KEGG气泡图
ggplot(kegg_enrich_df, aes(x=GeneRatio, y=number, colour=-1*log10(qvalue), size=Count))+
  geom_point()+
  # scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "#8DA0CB",high = "#66C2A5")+
  scale_y_discrete(labels=labels) +
  theme_bw()+
  ylab("KEGG Pathway")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)))
