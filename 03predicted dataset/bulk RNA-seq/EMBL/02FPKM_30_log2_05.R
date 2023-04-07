#===========================EMBL_expression_atlas===============================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/EMBL_expression_atlas/")
####基因长度(37)####
gene_length <- read.table("I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/gene_length.txt",sep='\t',header = F,as.is = T)
colnames(gene_length) <- c("gene","length")
gene_length$gene <- gsub("\\..*","",gene_length$gene)
head(gene_length)

dir.create("FPKM")
dir.create("30_log_05")
library(tibble)

cancer <- c("fibrosarcoma","nasopharyngeal carcinoma","neuroblastoma","non-small-cell lung cancer",
            "ovarian serous adenocarcinoma","pancreatic ductal carcinoma(PDAC)")

nrows <- c()
for(i in 1:length(cancer)){
  file <- list.files(paste0(getwd(),"/",cancer[i]))
  count_file <- file[grepl("raw-counts",file)]
  count <- read.table(paste0(getwd(),"/",cancer[i],"/",count_file),sep='\t',header = T,as.is = T,check.names = F)
  count2 <- count[,-2]
  colnames(count2)[1] <- "gene"
  count2[1:4,1:4]
  ####FPKM####
  #匹配gene_count与gene_length
  count_merge <- merge(count2,gene_length,by="gene")
  count_sam <- setdiff(colnames(count_merge),c("gene","length"))
  #计算每个样本的mapped reads数
  mapped_reads <- colSums(count_merge[,count_sam])
  #计算FPKM值
  FPKM <- count_merge[,count_sam]/(10^-9*matrix(as.numeric(count_merge[,"length"]),ncol=1) %*% matrix(mapped_reads,nrow=1))
  #相同ENSG合并
  frenq <- table(count_merge$gene)
  cat(paste0(max(frenq),"\n"))
  rownames(FPKM) <- count_merge$gene
  
  save(FPKM,file = paste0(getwd(),"/FPKM/",cancer[i],"_cancer_FPKM.RData"))
  ####移除FPKM在超过30%（包括30%）样本为0的基因，log2(FPKM+0.05)转化####
  len <- apply(FPKM,1,function(x) length(which(x==0)))
  FPKM2 <- FPKM[len < round(ncol(FPKM)*0.3),,drop=F]
  FPKM3 <- log2(FPKM2+0.05)
  
  save(FPKM3, file = paste0(getwd(),"/30_log_05/",cancer[i],"_cancer_FPKM_30_log_05.RData"))
  cat(paste0(cancer[i],"\n"))
  nrows <- rbind(nrows,cbind(cancer[i],nrow(FPKM),ncol(FPKM),nrow(FPKM3)))
  
}

colnames(nrows) <- c("cancer","gene_num","sample_num","gene_30_log_05_num")
write.table(nrows,"cancer_count_info.txt",sep='\t',quote = F,row.names = F) 

####cancer-normal####
cancer <- c("papillary thyroid carcinoma","Small Cell Lung Cancer","uterine leiomyosarcoma")

dir.create("normal/FPKM",recursive =T)
dir.create("normal/30_log_05",recursive =T)

nrows <- c()
for(i in 1:length(cancer)){
  file <- list.files(paste0(getwd(),"/",cancer[i]))
  cancer_count_file <- file[grepl("cancer_raw-counts",file)]
  count <- read.table(paste0(getwd(),"/",cancer[i],"/",cancer_count_file),sep='\t',header = T,as.is = T,check.names = F)
  count2 <- count[,-2]
  colnames(count2)[1] <- "gene"
  count2[1:4,1:4]
  ####FPKM####
  #匹配gene_count与gene_length
  count_merge <- merge(count2,gene_length,by="gene")
  count_sam <- setdiff(colnames(count_merge),c("gene","length"))
  #计算每个样本的mapped reads数
  mapped_reads <- colSums(count_merge[,count_sam])
  #计算FPKM值
  FPKM <- count_merge[,count_sam]/(10^-9*matrix(as.numeric(count_merge[,"length"]),ncol=1) %*% matrix(mapped_reads,nrow=1))
  #相同ENSG合并
  frenq <- table(count_merge$gene)
  cat(paste0(max(frenq),"\n"))
  rownames(FPKM) <- count_merge$gene

  save(FPKM,file = paste0(getwd(),"/FPKM/",cancer[i],"_cancer_FPKM.RData"))
  ####移除FPKM在超过30%（包括30%）样本为0的基因，log2(FPKM+0.05)转化####
  len <- apply(FPKM,1,function(x) length(which(x==0)))
  FPKM2 <- FPKM[len < round(ncol(FPKM)*0.3),,drop=F]
  FPKM3 <- log2(FPKM2+0.05)
  
  save(FPKM3, file = paste0(getwd(),"/30_log_05/",cancer[i],"_cancer_FPKM_30_log_05.RData"))
  cat(paste0(cancer[i],"\n"))
  ####normal####
  normal_count_file <- file[grepl("normal_raw-counts",file)]
  count_n <- read.table(paste0(getwd(),"/",cancer[i],"/",normal_count_file),sep='\t',header = T,as.is = T,check.names = F)
  count_n2 <- count_n[,-2]
  colnames(count_n2)[1] <- "gene"
  count_n2[1:4,1:4]
  ####FPKM####
  #匹配gene_count与gene_length
  count_merge_n <- merge(count_n2,gene_length,by="gene")
  count_sam_n <- setdiff(colnames(count_merge_n),c("gene","length"))
  #计算每个样本的mapped reads数
  mapped_reads_n <- colSums(count_merge_n[,count_sam_n])
  #计算FPKM值
  FPKM_n <- count_merge_n[,count_sam_n]/(10^-9*matrix(as.numeric(count_merge_n[,"length"]),ncol=1) %*% matrix(mapped_reads_n,nrow=1))
  #相同ENSG合并
  frenq <- table(count_merge_n$gene)
  cat(paste0(max(frenq),"\n"))
  rownames(FPKM_n) <- count_merge_n$gene
  
  save(FPKM_n,file = paste0(getwd(),"/normal/FPKM/",cancer[i],"_normal_FPKM.RData"))
  ####移除FPKM在超过30%（包括30%）样本为0的基因，log2(FPKM+0.05)转化####
  len <- apply(FPKM_n,1,function(x) length(which(x==0)))
  FPKM_n2 <- FPKM_n[len < round(ncol(FPKM_n)*0.3),,drop=F]
  FPKM_n3 <- log2(FPKM_n2+0.05)
  
  save(FPKM_n3, file = paste0(getwd(),"/normal/30_log_05/",cancer[i],"_normal_FPKM_30_log_05.RData"))
  cat(paste0(cancer[i],"\n"))
  nrows <- rbind(nrows,cbind(cancer[i],nrow(FPKM),ncol(FPKM),nrow(FPKM3),nrow(FPKM_n),ncol(FPKM_n),nrow(FPKM_n3)))
  
}

colnames(nrows) <- c("cancer","gene_num","sample_num","gene_30_log_05_num",
                     "normal_gene_num","normal_sample_num","normal_gene_30_log_05_num")
write.table(nrows,"cancer-normal_count_info.txt",sep='\t',quote = F,row.names = F) 
