#
#===================GSE72056_melanoma==============================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE72056_melanoma"
path <- file.path(d,file)
data <- read.table(paste0(path,"/GSE72056_melanoma_single_cell_revised_v2.txt"),sep='\t',header = T,as.is = T)

meta <- cbind.data.frame(colnames(data),t(data[1:3,]),stringsAsFactors =F)
meta2 <- meta[-1,]
colnames(meta2) <- as.matrix(meta[1,])
colnames(meta2)[3:4] <- c("malignant","non_malignant")
library(dplyr)
library(plyr)
#malignant(1=no,2=yes,0=unresolved)
#non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)

meta_filter <- meta2 %>% select(Cell,malignant,non_malignant) %>%
  filter(!as.numeric(malignant)==0) %>%
  filter(!(as.numeric(malignant)==1&as.numeric(non_malignant)==0)) %>%
  filter(!(as.numeric(malignant)==2&as.numeric(non_malignant)%in%c(1,3)))

meta_filter$celltype <- NA
meta_filter[as.numeric(meta_filter$malignant)==2&as.numeric(meta_filter$non_malignant)==0,]$celltype <- "Malignant"
meta_filter[as.numeric(meta_filter$malignant)==1&as.numeric(meta_filter$non_malignant)==1,]$celltype <- "T cell"
meta_filter[as.numeric(meta_filter$malignant)==1&as.numeric(meta_filter$non_malignant)==2,]$celltype <- "B cell"
meta_filter[as.numeric(meta_filter$malignant)==1&as.numeric(meta_filter$non_malignant)==3,]$celltype <- "Macrophage"
meta_filter[as.numeric(meta_filter$malignant)==1&as.numeric(meta_filter$non_malignant)==4,]$celltype <- "Endothelial"
meta_filter[as.numeric(meta_filter$malignant)==1&as.numeric(meta_filter$non_malignant)==5,]$celltype <- "CAF"
meta_filter[as.numeric(meta_filter$malignant)==1&as.numeric(meta_filter$non_malignant)==6,]$celltype <- "NK"

write.table(meta_filter,paste0(path,"/GSE72056_melanoma_meta_filter.txt"),sep='\t',quote=F,row.names=F)
meta_filter <- read.table(paste0(path,"/GSE72056_melanoma_meta_filter.txt"),sep='\t',header = T,as.is = T)        
table(meta_filter$celltype)
    
count <- data[-c(1:3),]
frenq <- table(as.matrix(count[,1]))
max(frenq)

uni_count <- count[which(as.matrix(count[,1]) %in% names(frenq[frenq!=max(frenq)])),]
uni_count2 <- uni_count[,-1]
rownames(uni_count2) <- as.matrix(uni_count[,1])

more_count <- count[which(as.matrix(count[,1]) %in% names(frenq[frenq==max(frenq)])),]
more_count2 <- apply(more_count[,-1],2,function(x){
  tapply(x, factor(more_count[,1]), function(x) mean(as.numeric(x)))
})

count2 <- rbind(uni_count2,more_count2)
count2 <- count2[,meta_filter$Cell]
identical(meta_filter$Cell,colnames(count2))
#TRUE
saveRDS(count2,paste0(path,"/GSE72056_melanoma_count.rds"))
count2 <- readRDS(paste0(path,"/GSE72056_melanoma_count.rds"))
#
#============================GSE146771_colon_cancer=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE146771_colon_cancer"
path <- file.path(d,file)

count <- read.table(paste0(path,"/CRC.Leukocyte.Smart-seq2.TPM.txt"),sep=" ",header = T,as.is = T)
dim(count)
#15179 10468
meta <- read.table(paste0(path,"/CRC.Leukocyte.Smart-seq2.Metadata.txt"),sep='\t',header = T,as.is = T)
head(meta)

library(dplyr)
meta_filter <- meta %>% select(cell = CellName,Sample,Tissue,cell_type = Global_Cluster) %>%
  filter(!(cell_type == "Malignant cell" & !(Tissue %in% "T")))

table(meta_filter$cell_type)

write.table(meta_filter,paste0(path,"/GSE146771_colon_cancer_meta_filter.txt"),sep='\t',quote=F,row.names=F)

count2 <- count[,meta_filter$cell]
saveRDS(count2,paste0(path,"/GSE146771_colon_cancer_count.rds"))
#
#============================GSE75688_breast cancer=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE75688_breast cancer"
path <- file.path(d,file)

samples <- read.table(paste0(path,"/GSE75688_final_sample_information.txt"),sep='\t',header = T,as.is = T)
head(samples)

library(dplyr)
samples_filter <- samples %>% select(sample:index3) %>% filter(!type=="Bulk") 
dim(samples_filter)
#515   5
celltype <- readxl::read_xlsx(paste0(path,"/GSE75688_umap.xlsx"))

meta <- samples_filter %>% select(sample,index3)
colnames(meta) <- c("cell","cell_type")
write.table(meta,paste0(path,"/GSE75688_GEO_processed_Breast_Cancer_meta_filter.txt"),sep='\t',quote = F,row.names = F)

meta <- read.table(paste0(path,"/GSE75688_GEO_processed_Breast_Cancer_meta_filter.txt"),sep='\t',header = T,as.is = T)
table(meta$cell_type)

data <- read.table(paste0(path,"/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"),sep='\t',header = T,as.is = T)
dim(data)
count <- data[,c("gene_name",meta$sample)]

frenq <- table(count$gene_name)
max(frenq)

uni_count <- count[count$gene_name %in% names(frenq[frenq==1]),]
uni_count2 <- uni_count[,-1]
rownames(uni_count2) <- as.matrix(uni_count[,1])

more_count <- count[count$gene_name %in% names(frenq[frenq > 1]),]
more_count2 <- apply(more_count[,-1],2,function(x){
  tapply(x,factor(more_count[,1]),function(x) mean(as.numeric(x)))
})

count2 <- rbind(uni_count2,more_count2)
saveRDS(count2,paste0(path,"/GSE75688_GEO_processed_Breast_Cancer_count.rds"))
#
#============================GSE84465_GBM=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE84465_GBM"
path <- file.path(d,file)

count <- read.table(paste0(path,"/GSE84465_GBM_All_data.csv"),sep=' ',header = T,as.is = T,check.names = F)
meta <- readxl::read_xlsx(paste0(path,"/GSE84465_GBM_umap.xlsx"))
table(meta$cell_type)
library(dplyr)
saveRDS(count,paste0(path,"/GSE84465_GBM_All_data.rds"))
#
#============================GSE103322_HNSCC=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE103322_HNSCC"
path <- file.path(d,file)

data <- read.table(paste0(path,"/GSE103322_HNSCC_all_data.txt"),sep='\t',header = T,as.is = T)

meta <- read.table(paste0(path,"/GSE103322_HNSCC_metafile.txt"),sep='\t',header = T,as.is = T)

library(dplyr)
meta_filter <- meta %>% select(-c(processed.by.Maxima.enzyme,Lymph.node)) %>%
  filter(!(as.numeric(classified..as.cancer.cell)==1&as.numeric(classified.as.non.cancer.cells)==1)) %>%
  filter(!(as.numeric(classified..as.cancer.cell)==0&as.numeric(classified.as.non.cancer.cells)==0))

meta_filter$cell_type <- meta_filter$non.cancer.cell.type
meta_filter[meta_filter$non.cancer.cell.type=="-Fibroblast",]$cell_type <- "Fibroblast"
meta_filter[meta_filter$non.cancer.cell.type==0,]$cell_type <- "Malignant"

write.table(meta_filter,paste0(path,"/GSE103322_HNSCC_meta_filter.txt"),sep='\t',quote = F,row.names = F)

meta_filter <- read.table(paste0(path,"/GSE103322_HNSCC_meta_filter.txt"),sep='\t',header = T,as.is = T)
table(meta_filter$cell_type)
count <- data[-c(1:5),]
count2 <- count[,-1]
rownames(count2) <- as.matrix(count[,1])

count2 <- count2[,meta_filter$cell]
saveRDS(count2,paste0(path,"/GSE103322_HNSCC_count.rds"))
#
#============================GSE115978_melanoma=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE115978_melanoma"
path <- file.path(d,file)

data <- read.table(paste0(path,"/GSE115978_counts.csv"),sep=',',header = T,as.is = T)
dim(data)
#[1] 23686  7187
count <- data[,-1]
rownames(count) <- as.matrix(data[,1])
dim(count)
#23686  7186
saveRDS(count,paste0(path,"/GSE115978_counts.rds"))

meta <- read.table(paste0(path,"/GSE115978_cell.annotations.csv"),sep=',',header = T,as.is = T)
dim(meta)
 
library(dplyr)
meta_filter <- meta %>% select(cells,cell.types) %>% 
  filter(!cell.types == "?")
dim(meta_filter)
colnames(meta_filter) <- c("cell","cell_type")

write.table(meta_filter,paste0(path,"/GSE115978_melanoma_meta_filter.txt"),sep='\t',quote = F,row.names = F)
#
#============================GSE117988_Merkel cell carcinoma=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE117988_Merkel cell carcinoma"
path <- file.path(d,file)
meta <- readxl::read_xlsx(paste0(path,"/GSE117988_umap.xlsx"))
head(meta)
table(meta$cell_type)
table(meta$cell_sample)
dim(meta)
#18396     5

library(dplyr)
index <- lapply(meta$cell_id,function(x) unlist(strsplit(x,"_"))[1]) %>% unlist

PBMC <- read.table(paste0(path,"/GSE117988_raw.expMatrix_PBMC.csv"),sep=',',header = T,as.is = T,check.names = F)
dim(PBMC)
colnames(PBMC)[1] <- "gene"
#17712     12875
tumor  <- read.table(paste0(path,"/GSE117988_raw.expMatrix_Tumor.csv"),sep=',',header = T,as.is = T,check.names = F)
dim(tumor)
colnames(tumor)[1] <- "gene"
#21861  7432

count <- merge(PBMC,tumor,by = "gene",all.y = TRUE)
count[is.na(count)] <- 0
count2 <- count[,-1]
rownames(count2) <- as.matrix(count[,1])

inter_sam <- intersect(index,colnames(count2))
meta_filter <- meta[which(index %in% inter_sam),]
index2 <- index[index %in% inter_sam]
meta_filter$cell <- index2
meta_filter2 <- meta_filter %>% select(cell,cell_type,cell_sample) %>%
  filter(!(grepl("Tumor",cell_sample)&!(cell_type %in% c("Tumor"))))
write.table(meta_filter2,paste0(path,"/GSE117988_Merkel cell carcinoma_meta_filter.txt"),sep='\t',quote = F,row.names = F)

meta_filter2 <- read.table(paste0(path,"/GSE117988_Merkel cell carcinoma_meta_filter.txt"),sep='\t',header = T,as.is = T)

table(meta_filter2$cell_sample)
table(meta_filter2$cell_type)

count3 <- count2[,meta_filter2$cell]
dim(count3)
#21861 18064

saveRDS(count3,paste0(path,"/GSE117988_Merkel cell carcinoma_count.rds"))
#
#============================GSE118056_Merkel cell carcinoma=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE118056_Merkel cell carcinoma"
path <- file.path(d,file)

meta <- readxl::read_xlsx(paste(path,"GSE118056_umap.xlsx",sep='/'))
dim(meta)
#11071     5
table(meta$cell_type)
table(meta$cell_sample)

count <- read.table(paste(path,"GSE118056_raw.expMatrix.csv",sep='/'),sep=',',header = T,as.is = T)
dim(count)
#25066 11072
count2 <- count[,-1]
rownames(count2) <- as.matrix(count[,1])
saveRDS(count2,paste(path,"GSE118056_raw.expMatrix.rds",sep='/'))
#
#============================GSE123813_Basal Cell Carcinoma=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

####scRNA-seq
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE123813_Basal Cell Carcinoma"
path <- file.path(d,file)
meta <- read.table(paste0(path,"/bcc_all_metadata.txt"),sep='\t',header = T,as.is = T,check.names = F,fileEncoding="UTF16")
head(meta)
table(meta$cluster)
dim(meta)

meta <- na.omit(meta)
dim(meta)

table(meta$sort)
library(dplyr)

meta_filter <- meta %>% select(cell.id,sort,cluster) %>%
  filter(!(grepl("CD45\\+",sort)&cluster %in% c("CAFs","Endothelial","Melanocytes","Myofibroblasts",
                                            "Tumor_1","Tumor_2"))) %>%
  filter(!(sort == "CD45- CD3-"&cluster %in% c("B_cells_1","B_cells_2","CD4_T_cells","CD8_act_T_cells",
                                              "CD8_ex_T_cells","CD8_mem_T_cells","DCs","Macrophages",
                                              "NK_cells","pDCs","Plasma_cells","Tcell_prolif","Tregs")))
meta_filter[meta_filter$cluster %in% c("B_cells_1","B_cells_2"),]$cluster <- "B_cells"
meta_filter[meta_filter$cluster %in% c("Tumor_1","Tumor_2"),]$cluster <- "Tumor"
write.table(meta_filter,paste0(path,"/meta_filter.txt"),sep='\t',quote = F,row.names = F)

meta_filter <- read.table(paste0(path,"/",file,"_meta_filter.txt"),sep='\t',header = T,as.is = T)

count <- read.table(paste0(path,"/bcc_scRNA_counts.txt"),sep='\t',header = T,as.is = T)
saveRDS(count,paste0(path,"/bcc_scRNA_counts.rds"))
#
#============================GSE125449_liver cancer=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE125449_liver cancer"
path <- file.path(d,file)

meta1 <- read.table(paste(path,"GSE125449_Set1_samples.txt",sep='/'),sep='\t',header = T,as.is = T,check.names = F)
head(meta1)

meta2 <- read.table(paste(path,"GSE125449_Set2_samples.txt",sep='/'),sep='\t',header = T,as.is = T,check.names = F)
head(meta2)

library(dplyr)
meta <- rbind(meta1,meta2) %>% select(cell = `Cell Barcode`,cell_type = Type)
dim(meta)
write.table(meta,paste(path,"GSE125449_liver_cancer_cellLabel.txt",sep='/'),sep='\t',quote = F,row.names = F)
library(Matrix)
####set1####
barcode.path1 <- paste0(path,"/GSE125449_Set1_barcodes.tsv")
features.path1 <- paste0(path,"/GSE125449_Set1_genes.tsv")
matrix.path1 <- paste0(path,"/GSE125449_Set1_matrix.mtx")

mat1 <- readMM(file = matrix.path1)
feature.names1 = read.delim(features.path1,header = FALSE, stringsAsFactors = FALSE)
barcode.names1 = read.delim(barcode.path1,header=FALSE, stringsAsFactors = FALSE)
colnames(mat1) = barcode.names1$V1
rownames(mat1) = feature.names1$V2
dim(mat1)
#20124  5115
mat1 <- cbind.data.frame(rownames(mat1),as.matrix(mat1))
colnames(mat1)[1] <- "gene"
####set2####
barcode.path2 <- paste0(path,"/GSE125449_Set2_barcodes.tsv")
features.path2 <- paste0(path,"/GSE125449_Set2_genes.tsv")
matrix.path2 <- paste0(path,"/GSE125449_Set2_matrix.mtx")

mat2 <- readMM(file = matrix.path2)
feature.names2 = read.delim(features.path2,header = FALSE, stringsAsFactors = FALSE)
barcode.names2 = read.delim(barcode.path2,header=FALSE, stringsAsFactors = FALSE)
colnames(mat2) = barcode.names2$V1
rownames(mat2) = feature.names2$V2
dim(mat2)
#19572  4831
mat2 <- cbind.data.frame(rownames(mat2),as.matrix(mat2))
colnames(mat2)[1] <- "gene"

count <- merge(mat1,mat2,by="gene")
dim(count)
count2 <- gene_mean(count)
dim(count2)
#18367  9946
saveRDS(count2,paste(path,"GSE125449_liver_cancer_count.rds",sep='/'))
####多个基因取均值，function####
gene_mean <- function(count){
  frenq <- table(count$gene)

  uni_count <- count[count$gene %in% names(frenq[frenq==1]),]
  uni_count2 <- uni_count[,-1]
  rownames(uni_count2) <- as.matrix(uni_count[,1])
  
  more_count <- count[count$gene %in% names(frenq[frenq > 1]),]
  more_count2 <- apply(more_count[,-1],2,function(x){
    tapply(x,factor(more_count[,1]),function(x) mean(as.numeric(x)))
  })
  
  count2 <- rbind(uni_count2,more_count2)
  return(count2)
}
#
#============================GSE134520_Early Gastric Cancer=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE134520_Early Gastric Cancer"
path <- file.path(d,file)

meta <- read.table(paste0(path,"/cellAnnotations.txt"),sep='\t',header = F,as.is = T)
dim(meta)
colnames(meta) <- c("cell","cell_type")
#2861    2
table(meta$cell_type)
write.table(meta,paste0(path,"/GSE134520_Early_Gastric_Cancer_cellLabels.txt"),sep="\t",quote = F,row.names = F)

count <- read.table(paste0(path,"/expdata_GSE134520_Early Gastric Cancer.txt"),sep='\t',header = T,as.is = T)
dim(count)
# 11962  2862
saveRDS(count,paste0(path,"/GSE134520_Early_Gastric_Cancer_count.rds"))
#
#============================GSE141982_glioma=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")
lncRNA <- read.table("GRCh38_lncRNA_name.txt",sep='\t',header = F,as.is = T)
lncRNA <- as.matrix(lncRNA)
dim(lncRNA)

mRNA <- read.table("GRCh38_mRNA_name.txt",sep='\t',header = F,as.is = T)
mRNA <- as.matrix(mRNA)

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE141982_glioma"
path <- file.path(d,file)

count <- read.table(paste0(path, "/expdata_GSE141982-glioma.txt"),sep='\t',header = T,as.is = T)
dim(count)
colnames(count) <- gsub("\\.","-",colnames(count))
#16628  6947
saveRDS(count,paste0(path, "/GSE141982_glioma_count.rds"))

meta <- read.table(paste0(path,"/cellAnnotations.txt"),sep='\t',header = F,as.is = T)
head(meta)
colnames(meta) <- c("cell","cell_type")
write.table(meta,paste0(path,"/GSE141982_glioma_metafile.txt"),sep='\t',quote = F,row.names = F)
#
#============================GSE145137_bladder_cancer=====================
setwd("I:/F/project_3/original_scRNA_lncRNA")

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_original"
file <- "GSE145137_bladder_cancer"
path <- file.path(d,file)

count <- read.table(paste(path,"GSM4307111_GEO_processed_BC159-T_3_log2TPM_matrix_final.txt",sep='/'),sep='\t',header = T,as.is = T)
dim(count)
#14233  2076
meta <- readxl::read_xlsx(paste0(path,"/GSE145137_BC159-T_3_All_SC_final_QC_Celltype_Information.xlsx"))
table(meta$cell_type)
library(dplyr)

meta_filter <- meta %>% select(cell=Samples,cell_type, cell_index) %>%
  filter(!(grepl("Basal tumor cells",cell_type)&cell_index != "Tumor")) %>%
  filter(!cell_type %in% c("Unknown1","Unknown2"))

write.table(meta_filter,paste0(path,"/GSE145137_BC159-T_3_meta_filter.txt"),sep='\t',quote = F,row.names = F)

count2 <- count[,c("gene",meta_filter$cell)]
count3 <- count2[,-1]
rownames(count3) <- as.matrix(count2[,1])
dim(count3)
#14233  1996
saveRDS(count3,paste0(path,"/GSE145137_BC159-T_3_count.rds"))
