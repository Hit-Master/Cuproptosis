gene_set_xls <- function(DEG_file='', gene_file='', merge_method = 'intersect'){
  # 差异基因 Gene Symbol
  degs <- data.table::fread(DEG_file, data.table = F)
  degs <- degs[[1]]
  # degs <- gsub(degs, pattern = '-', replacement = '_')
  
  # 特定基因集 Gene Symbol
  if (gene_file != ''){
    gene_set <- readxl::read_excel(gene_file)
    gene_set <- gene_set[[1]]
    # gene_set <- gsub(gene_set, pattern = '-', replacement = '_')
  } else{
    gene_set <- NULL
  }
  
  # gene_list
  if (merge_method == 'intersect') {
    gene_list <- intersect(degs, gene_set)
  } else if (merge_method == 'union') {
    gene_list <- union(degs, gene_set)
  } else if (merge_method == 'setdiff') {
    gene_list <- setdiff(degs, gene_set)
  } else {
    print('Please select valid merge_method: intersect, union or setdiff')
  }
  
  return(gene_list)
}


wgcna_tpm <- function(tpm_file, gene_list){
  library(data.table)
  library(WGCNA)
  library(stringr)
  
  # 构建gene list的矩阵
  tpm <- fread(tpm_file, data.table = F)
  row.names(tpm) <- tpm[, 1]
  tpm <- tpm[, -1]
  t_tpm <- as.data.frame(t(tpm))
  t_tpm <- t_tpm[, gene_list]
  
  # 缺失值控制
  sample_qc = goodSamplesGenes(t_tpm)
  
  # 离群样本控制
  sampleTree = hclust(dist(t_tpm), method = "average")
  sizeGrWindow(12,9)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", 
       sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  
  wgcna_tpm_report <- list(tpm = t_tpm, sample_qc = sample_qc, sampleTree = sampleTree)
  
  return(wgcna_tpm_report)
}


wgcna_sample_filter <- function(tpm, sampleTree, cutHeight, clust_num=1){
  library(WGCNA)
  
  # 过滤显著离群样本
  clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
  print(table(clust))

  keepSamples = (clust==clust_num)
  tpm = tpm[keepSamples, ]
  
  return(tpm)
}


phenotype_filter <- function(phenotype_file, survival_file, wgcna_matrix, phenotype_format, phenotype_list=''){
  library(data.table)
  
  phenotype <- fread(phenotype_file, data.table = F)
  survival <- fread(survival_file, data.table = F)
  
  # 构建phenotype list的表型矩阵
  if (phenotype_format == 'all') {
    if (phenotype_list != '') {
      phenotype <- phenotype[, phenotype_list]
    } else {
      phenotype <- phenotype[, c('submitter_id.samples', 'age_at_initial_pathologic_diagnosis', 'pathologic_M', 
                                 'pathologic_N', 'pathologic_T','gender.demographic', 'tumor_stage.diagnoses')]
    }
  } else if (phenotype_format == 'filter') {
    print('Please ensure you have filtered phenotype information.')
  } else {
    print('Please assign valid phenotype file format: all or filter.')
    print('Meanwhile, all phenotype file needs to submit a list of phenotype information of interest.')
  }
  phenotype <- phenotype[match(row.names(wgcna_matrix), phenotype$submitter_id.samples),]
  row.names(phenotype) <- phenotype[, 1]
  phenotype <- phenotype[, -1]
  
  # 表型矩阵加生存信息
  survival <- survival[, c(1, 2, 4)]
  survival <- survival[match(row.names(wgcna_matrix), survival$sample),]
  phenotype <- cbind(phenotype, survival[, c(2, 3)])
  
  return(phenotype)
}


gene_list <- gene_set_xls('DEG.txt', 'Gene_set.xlsx', merge_method = 'union')

wgcna_tpm_report <- wgcna_tpm(tpm_file = 'tpm.txt', gene_list)
tpm <- wgcna_tpm_report$tpm
sample_qc <- wgcna_tpm_report$sample_qc
sampleTree <- wgcna_tpm_report$sampleTree
abline(h = 60000, col = "red") #根据实际情况而定

wgcna_matrix <- wgcna_sample_filter(tpm, sampleTree, cutHeight = 60000)
phenotype <- phenotype_filter(phenotype_file = 'TCGA-LUAD.GDC_phenotype.tsv', survival_file = 'TCGA-LUAD.survival.tsv', wgcna_matrix, phenotype_format = 'all')

phenotype$gender.demographic <- ifelse(phenotype$gender.demographic == 'female', 1, 2)

phenotype$pathologic_M <- gsub('a', '',phenotype$pathologic_M)
phenotype$pathologic_M <- gsub('b', '',phenotype$pathologic_M)
phenotype$pathologic_M <- ifelse(phenotype$pathologic_M == 'M0', 1, 
                          ifelse(phenotype$pathologic_M == 'M1', 2, NA))

phenotype$pathologic_N <- ifelse(phenotype$pathologic_N == 'N0', 1, 
                          ifelse(phenotype$pathologic_N == 'N1', 2, 
                          ifelse(phenotype$pathologic_N == 'N2', 3, 
                          ifelse(phenotype$pathologic_N == 'N3', 4, NA))))

phenotype$pathologic_T <- gsub('a', '',phenotype$pathologic_T)
phenotype$pathologic_T <- gsub('b', '',phenotype$pathologic_T)
phenotype$pathologic_T <- ifelse(phenotype$pathologic_T == 'T1', 1, 
                          ifelse(phenotype$pathologic_T == 'T2', 2, 
                          ifelse(phenotype$pathologic_T == 'T3', 3, 
                          ifelse(phenotype$pathologic_T == 'T4', 4, NA))))

phenotype$tumor_stage.diagnoses <- gsub('a', '',phenotype$tumor_stage.diagnoses)
phenotype$tumor_stage.diagnoses <- gsub('b', '',phenotype$tumor_stage.diagnoses)
phenotype$tumor_stage.diagnoses <- ifelse(phenotype$tumor_stage.diagnoses == 'stge i', 1, 
                                   ifelse(phenotype$tumor_stage.diagnoses == 'stge ii', 2, 
                                   ifelse(phenotype$tumor_stage.diagnoses == 'stge iii', 3, 
                                   ifelse(phenotype$tumor_stage.diagnoses == 'stge iv', 4, NA))))

colnames(phenotype) <- c('age', 'M', 'N', 'T', 'gender', 'stage', 'OS', 'OS.time')

save(wgcna_matrix, phenotype, file = 'step01-preprocessing.Rdata')
