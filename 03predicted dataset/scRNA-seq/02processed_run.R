#=====================================预测单细胞细胞互作（放进TICCom browse里）======================================
####GSE72056_melanoma####
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
path <- "GSE72056_melanoma"

out <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"

meta <- read.table(paste(d,path,"GSE72056_melanoma_meta_filter.txt",sep='/'),sep='\t',header = T,as.is = T)
celllabel <- meta[,c("Cell","celltype")]
colnames(celllabel) <- c("cell","cell_type")
dim(celllabel)
out_name <- paste0(path,"_cellLabel.txt")
dir.create(paste(out,path,sep='/'))
write.table(celllabel,paste(out,path,out_name,sep='/'),sep='\t',quote = F,row.names = F)

count <- readRDS(paste(d,path,"GSE72056_melanoma_count.rds",sep='/'))
count[1:4,1:4]

rm(list=ls())
####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE72056_melanoma"

dataFile <-paste(d2,path,"GSE72056_melanoma_count.rds",sep='/')
metaFile <- paste(d2,path,"GSE72056_melanoma_cellLabel.txt",sep='/')
referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")

start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,organism = "human",stats = "mean",outPath,
                                   databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                            outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","icellnet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE75688_breast cancer####
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
path <- "GSE75688_breast cancer"

out <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"

meta <- read.table(paste(d,path,"GSE75688_GEO_processed_Breast_Cancer_meta_filter.txt",sep='/'),sep='\t',header = T,as.is = T)
dim(meta)

count <- readRDS(paste(d,path,"GSE75688_GEO_processed_Breast_Cancer_count.rds",sep='/'))
count[1:4,1:4]

rm(list=ls())
####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE75688_breast cancer"

dataFile <-paste(d2,path,"GSE75688_GEO_processed_Breast_Cancer_count.rds",sep='/')
metaFile <- paste(d2,path,"GSE75688_GEO_processed_Breast_Cancer_meta_filter.txt",sep='/')
referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")

start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","icellnet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE84465_GBM####
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
path <- "GSE84465_GBM"

out <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"

meta <- readxl::read_xlsx(paste(d,path,"GSE84465_GBM_umap.xlsx",sep='/'))
dim(meta)
table(meta$cell_sample)

meta$compare_group <- NA
meta[grep("Periphery",meta$cell_sample),]$compare_group <- "Periphery"
meta[grep("Tumor",meta$cell_sample),]$compare_group <- "Tumor"

table(meta$cell_type)
num <- as.matrix(table(meta$cell_type,meta$compare_group))
cells <- num[!(num[,1]<2|num[,2]<2),]
celltypes <- rownames(cells)
#去掉不在两组样本中出现的细胞类型，且每组细胞个数小于2个的
meta_filter <- meta[meta$cell_type %in% celltypes,]

#each replicate has the same cell clusters
num2 <- as.matrix(table(meta_filter$cell_type,meta_filter$cell_sample))
cells2_index <- apply(num2,1,function(x){ length(which(x==0))})
celltypes2 <- names(cells2_index[cells2_index==0])
meta_filter2 <- meta_filter[meta_filter$cell_type %in% celltypes2,]
dim(meta_filter2)
#2345    6
library(dplyr)
celllabel <- meta_filter2 %>% select(cell = cell_id,cell_type)

group <-  meta_filter2 %>% select(cell = cell_id,compare_group = compare_group,
                          replicate = cell_sample)

dir.create(paste(out,path,sep = "/"))
write.table(celllabel,paste(out,path,"GSE84465_GBM_cellLabel.txt",sep='/'),sep = "\t",quote = F,row.names = F)
write.table(group,paste(out,path,"GSE84465_GBM_compare_replicate_group.txt",sep='/'),sep = "\t",quote = F,row.names = F)

count <- read.table(paste(d,path,"GSE84465_GBM_All_data.csv",sep='/'),sep = " ",header = T,as.is = T,check.names = F)
dim(count)
count[1:4,1:4]
saveRDS(count,paste(out,path,"GSE84465_GBM_expression.rds",sep='/'))
rm(list=ls())

####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE84465_GBM"

dataFile <-paste(d2,path,"GSE84465_GBM_expression.rds",sep='/')
metaFile <- paste(d2,path,"GSE84465_GBM_cellLabel.txt",sep='/')
groupFile <- paste(d2,path,"GSE84465_GBM_compare_replicate_group.txt",sep='/')

referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_DEG.R")
start <- Sys.time()
iTALK_method_DEG(dataFile,metaFile,groupFile,referencePath,gene_names = "symbol",organism = "human",
                  min_valid_cells = 0,min_gene_expressed = 0,
                             DEG_method = "Wilcox",q_cutoff = 0.05,databaseFile= "union_LR_database",outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/celltalker_method.R")
start <- Sys.time()
celltalker_method(dataFile,metaFile,group_replicateFile = groupFile,referencePath,gene_names = "symbol",
                  organism = "human",cells.reqd = 0,
                             freq.pos.reqd = 0,freq.group.in.cluster = 0,
                  databaseFile = "union_LR_database",outPath)

end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/nichenet_method.R")
start <- Sys.time()
NicheNet_method(dataFile,metaFile,defined_receiver_cell = "false",genesetFile="empty",
                groupFile,referencePath,
                            gene_names = "symbol",organism = "human",q_cutoff=0.05,
                databaseFile = "union_LR_database",outPath)

end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")
start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,
                       organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start
#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","italk_deg","celltalker","icellnet","nichenet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE116256_AML(除italk外其他没运行，太大)####
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
path <- "GSE116256_AML"

out <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"

count1 <- readRDS(paste(d,path,"GSE116256_AML_D0_count.rds",sep='/'))
count2 <- readRDS(paste(d,path,"GSE116256_AML_D14after_count.rds",sep='/'))
identical(rownames(count1),rownames(count2))
count <- cbind(count1,count2)
saveRDS(count,paste(d,path,"GSE116256_AML_D0_D14after_count.rds",sep='/'))

meta1 <- read.table(paste(d,path,"GSE116256_AML_D0_meta.txt",sep='/'),sep='\t',header = T,as.is = T)
meta2 <- read.table(paste(d,path,"GSE116256_AML_D14after_metafilter.txt",sep='/'),sep='\t',header = T,as.is = T)
meta <- rbind(meta1,meta2)
meta$group <- NA
meta[meta$Cell %in% meta1$Cell,]$group <- "D0"
meta[meta$Cell %in% meta2$Cell,]$group <- "D14after"

####太多，筛选一下####
meta2 <- meta[meta$PredictionRF2 == meta$PredictionRefined,]
dim(meta2)

library(dplyr)
meta2_filter <- ddply(meta2,.(cell_type,group),function(x){
  len <- round(nrow(x) * 0.5)
  x[1:len,]
})

count <- readRDS(paste(d,path,"GSE116256_AML_D0_D14after_count.rds",sep='/'))
count_filter <- count[,meta2_filter$Cell]
dim(count_filter)
saveRDS(count_filter,paste(out,path,"GSE116256_AML_D0_D14after_count_filter.rds",sep='/'))

celllabel <- meta2_filter %>% select(cell = Cell,cell_type)
group <- meta2_filter %>% select(cell = Cell,compare_group = group)
dir.create(paste(out,path,sep='/'))
write.table(celllabel,paste(out,path,"GSE116256_AML_D0_D14after_cellLabel.txt",sep='/'),sep='\t',quote = F,row.names = F)
write.table(group,paste(out,path,"GSE116256_AML_D0_D14after_compare_group.txt",sep='/'),sep='\t',quote = F,row.names = F)

####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE116256_AML"

dataFile <-paste(d2,path,"GSE116256_AML_D0_D14after_count_filter.rds",sep='/')
metaFile <- paste(d2,path,"GSE116256_AML_D0_D14after_cellLabel.txt",sep='/')
groupFile <- paste(d2,path,"GSE116256_AML_D0_D14after_compare_group.txt",sep='/')

referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_DEG.R")
start <- Sys.time()
iTALK_method_DEG(dataFile,metaFile,groupFile,referencePath,gene_names = "symbol",organism = "human",
                 min_valid_cells = 0,min_gene_expressed = 0,
                 DEG_method = "Wilcox",q_cutoff = 0.05,databaseFile= "union_LR_database",outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/nichenet_method.R")
start <- Sys.time()
NicheNet_method(dataFile,metaFile,defined_receiver_cell = "false",genesetFile="empty",
                groupFile,referencePath,
                gene_names = "symbol",organism = "human",q_cutoff=0.05,
                databaseFile = "union_LR_database",outPath)

end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")
start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,
                       organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start
#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","italk_deg","icellnet","nichenet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE117988_Merkel cell carcinoma####
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
path <- "GSE117988_Merkel cell carcinoma"

out <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"

meta <- read.table(paste(d,path,"GSE117988_Merkel cell carcinoma_meta_filter.txt",sep='/'),sep='\t',header = T,as.is = T)
meta$compare_group <- NA
meta[grepl("PBMC",meta$cell_sample),]$compare_group <- "PBMC"
meta[grepl("Tumor",meta$cell_sample),]$compare_group <- "Tumor"
colnames(meta)[3] <- "replicate"
dir.create(paste(out,path,sep='/'))

table(meta$cell_type,meta$compare_group)

library(dplyr)
celllabels <- meta %>% select(cell,cell_type)
group <- meta %>% select(cell,compare_group,replicate)

write.table(celllabels,paste(out,path,"GSE117988_Merkel_cell_carcinoma_cellLabels.txt",sep='/'),sep='\t',quote = F,row.names = F)
#write.table(group,paste(out,path,"GSE117988_Merkel_cell_carcinoma_compare_group.txt",sep='/'),sep='\t',quote = F,row.names = F)
####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE117988_Merkel cell carcinoma"

dataFile <-paste(d2,path,"GSE117988_Merkel cell carcinoma_count.rds",sep='/')
metaFile <- paste(d2,path,"GSE117988_Merkel_cell_carcinoma_cellLabels.txt",sep='/')

referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")
start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,
                       organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","icellnet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE125449_liver cancer####
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
path <- "GSE125449_liver cancer"
out <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"

meta <- read.table(paste(d,path,"GSE125449_liver_cancer_cellLabel.txt",sep='/'),sep='\t',header = T,as.is = T)
meta <- meta[meta$cell_type!="unclassified",]
dir.create(paste(out,path,sep="/"))
write.table(meta,paste(out,path,"GSE125449_liver_cancer_cellLabel_filter.txt",sep="/"),sep='\t',quote = F,row.names = F)

count <- readRDS(paste(d,path,"GSE125449_Hepatocellular_carcinoma_count.rds",sep='/'))

####as.matrix####
count <- as.matrix(count)
saveRDS(count,paste(out,path,"GSE125449_Hepatocellular_carcinoma_denseM_count.rds",sep='/'))
####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE125449_liver cancer"

dataFile <-paste(d2,path,"ExpressionData.rds",sep='/')
metaFile <- paste(d2,path,"cellLabel.txt",sep='/')

referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")
start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,
                       organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","icellnet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE145137_bladder_cancer####

####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE145137_bladder_cancer"

dataFile <-paste(d2,path,"ExpressionData.rds",sep='/')
metaFile <- paste(d2,path,"cellLabel.txt",sep='/')

referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")
start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,
                       organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","icellnet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE141982_glioma####

####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE141982_glioma"

dataFile <-paste(d2,path,"ExpressionData.rds",sep='/')
metaFile <- paste(d2,path,"cellLabel.txt",sep='/')

referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")
start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,
                       organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","icellnet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE134520_Early Gastric Cancer####
####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE134520_Early Gastric Cancer"

dataFile <-paste(d2,path,"ExpressionData.rds",sep='/')
metaFile <- paste(d2,path,"cellLabel.txt",sep='/')

referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")
start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,
                       organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","icellnet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE103322_HNSCC\multiple####
####run multiple####
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
path <- "GSE103322_HNSCC"
outPath <- paste(out_d,path,sep='/')
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","italk_deg","icellnet","nichenet")
path <- outPath
referencePath <-"I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
multiple_methods(path,method_names,referencePath,outpath = path)

####GSE115978_melanoma####
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
path <- "GSE115978_melanoma"
outPath <- paste(out_d,path,sep='/')
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","italk_deg","icellnet","nichenet")
path <- outPath

referencePath <-"I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
multiple_methods(path,method_names,referencePath,outpath = path)

####GSE123813_BCC_Cancer-skin(太多，删点)####
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
path <- "GSE123813_BCC_Cancer-skin"
out <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"

meta <- read.table(paste(d,path,"GSE123813_Basal Cell Carcinoma_meta_filter.txt",sep="/"),sep='\t',header = T,as.is = T)
head(meta)
table(meta$cluster)
dim(meta)
library(dplyr)
library(plyr)
meta2 <- meta %>% select(cell=cell.id,cell_type = cluster) %>% 
  filter(!(cell_type %in% c("CAFs","Endothelial","Melanocytes","Myofibroblasts")))
dim(meta2)

group <- read.table(paste(out,path,"GSE123813_BCC_Cancer_ComparisonGroup.txt",sep="/"),sep='\t',header = T,as.is = T)

meta3 <- merge(meta2,group,by="cell") 
table(meta3$cell_type,meta3$replicate)
table(meta3$replicate,meta3$compare_group)
#去掉没有肿瘤细胞的病人,且没有两组的样本
meta3_1 <- meta3[meta3$replicate %in% c("su004","su005","su006","su008"),]

meta3_sub <- ddply(meta3_1,.(cell_type,compare_group,replicate),function(x) {
  return(x[1:round(nrow(x)*0.5),])
})
table(meta3_sub$cell_type,meta3_sub$replicate)
#每个重复具有相同细胞类型
celllabel <- meta3_sub %>% select(cell,cell_type)
write.table(celllabel,paste(out,path,"cellLabel.txt",sep='/'),sep='\t',quote = F,row.names = F)
compare_group <- meta3_sub %>% select(cell,compare_group,replicate)
write.table(compare_group,paste(out,path,"compare_group.txt",sep='/'),sep='\t',quote = F,row.names = F)

count <- readRDS(paste(d,path,"bcc_scRNA_counts.rds",sep="/"))
count2 <- count[,celllabel$cell]
saveRDS(count2,paste(out,path,"ExpressionData.rds",sep='/'))
####run####
d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
path <- "GSE123813_BCC_Cancer-skin"

dataFile <-paste(d2,path,"ExpressionData.rds",sep='/')
metaFile <- paste(d2,path,"cellLabel.txt",sep='/')
groupFile <- paste(d2,path,"compare_group.txt",sep='/')

referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
dir.create(paste(out_d,path,sep='/'))
outPath <- paste(out_d,path,sep='/')

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_top_genes.R")
start <- Sys.time()
iTALK_method_top_genes(dataFile,metaFile,referencePath,gene_names = "symbol",top_genes = 50,
                       organism = "human",stats = "mean",outPath,
                       databaseFile = "union_LR_database")
end <- Sys.time()
end-start
#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/iTALK_method_DEG.R")
start <- Sys.time()
iTALK_method_DEG(dataFile,metaFile,groupFile,referencePath,gene_names = "symbol",organism = "human",
                 min_valid_cells = 0,min_gene_expressed = 0,
                 DEG_method = "Wilcox",q_cutoff = 0.05,databaseFile= "union_LR_database",outPath)
end <- Sys.time()
end-start
#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/celltalker_method.R")
start <- Sys.time()
celltalker_method(dataFile,metaFile,group_replicateFile = groupFile,referencePath,gene_names = "symbol",
                  organism = "human",cells.reqd = 0,
                  freq.pos.reqd = 0,freq.group.in.cluster = 0,
                  databaseFile = "union_LR_database",outPath)

end <- Sys.time()
end-start

#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/nichenet_method.R")
start <- Sys.time()
NicheNet_method(dataFile,metaFile,defined_receiver_cell = "false",genesetFile="empty",
                groupFile,referencePath,
                gene_names = "symbol",organism = "human",q_cutoff=0.05,
                databaseFile = "union_LR_database",outPath)

end <- Sys.time()
end-start
#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICELLNET_method.R")
start <- Sys.time()
ICELLNET_method(dataFile,metaFile,referencePath,gene_names = "symbol",organism = "human",top = 4,
                direction = "out",databaseFile = "union_LR_database",
                outPath)
end <- Sys.time()
end-start
#
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","italk_deg","celltalker","icellnet","nichenet")
path <- outPath

multiple_methods(path,method_names,outpath = path)

####GSE131907_lung-LUAD####
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
path <- "GSE131907_lung-LUAD"
outPath <- paste(out_d,path,sep='/')
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","italk_deg","celltalker","icellnet","nichenet")
path <- outPath

referencePath <-"I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
multiple_methods(path,method_names,referencePath,outpath = path)

####GSE146771_Colon_Cancer####
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_result"
path <- "GSE146771_Colon_Cancer"
outPath <- paste(out_d,path,sep='/')
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/multiple_methods.R")
#c("italk_top","italk_deg","celltalker","icellnet","nichenet")
method_names <- c("italk_top","icellnet")
path <- outPath

referencePath <-"I:/F/tumor_immune_interaction/computation/A_webserver/database/reference"
multiple_methods(path,method_names,referencePath,outpath = path)