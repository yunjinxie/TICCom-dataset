#==================================browse3癌症归类============================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search3_new/search3_update/")
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/download/cancer_catagory"

out <- "data_update"

cancer_type <- readxl::read_excel(paste0(d,"/TCGA+ICGC+EMBL+scRNA3_3.xlsx"))

d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search2_new"
union_LR <- read.table(paste0(d2,"/union_source_curated.txt"),sep='\t',header = T,as.is = T)
union_pairs <- paste(union_LR$ligand,union_LR$receptor,sep='_')
####TCGA####
#verify
verify <- read.table("search3_TCGA_verify.txt",sep='\t',header = T,as.is = T)

cancer_type_TCGA <- cancer_type[grepl("TCGA",cancer_type$Source2),]
cancers <- unique(verify$Cancer)
for(i in 1:length(cancers)){
  cc <- cancer_type_TCGA[cancer_type_TCGA$`Short name` == cancers[i],]$Rename2
  verify[verify$Cancer == cancers[i],]$Cancer <- cc
}
write.table(verify,paste0(out,"/search3_TCGA_verify.txt"),sep='\t',quote = F,row.names = F)
#LR
LR <- read.table("search3_TCGA_LR.txt",sep='\t',header = T,as.is = T)
pairs <- paste(LR$ligand,LR$receptor,sep='_')
LR$curated <- union_LR[match(pairs,union_pairs),]$Evidence
colnames(LR)[colnames(LR)=="curated"] <- "Evidence"

cancers <- unique(LR$cancer)
for(i in 1:length(cancers)){
  cc <- cancer_type_TCGA[cancer_type_TCGA$`Short name` == cancers[i],]$Rename2
  LR[LR$cancer == cancers[i],]$cancer <- cc
}
LR2 <- LR[,c(1:5,9,6:8)]
write.table(LR2,paste0(out,"/search3_TCGA_LR.txt"),sep='\t',quote = F,row.names = F)
####ICGC####
cancer_type_ICGC <- cancer_type[grepl("ICGC",cancer_type$Source2),]
#verify
verify2 <- read.table("search3_ICGC_verify.txt",sep='\t',header = T,as.is = T)

cancers <- unique(verify2$cancer_type)
for(i in 1:length(cancers)){
  cc <- cancer_type_ICGC[cancer_type_ICGC$`Short name` == cancers[i],]$Rename2
  verify2[verify2$cancer_type == cancers[i],]$cancer_type <- cc
}
verify2_out <- verify2[,c(1:6,10,8,9)]
colnames(verify2_out)[7] <- "Cancer" 
write.table(verify2_out,paste0(out,"/search3_ICGC_verify.txt"),sep='\t',quote = F,row.names = F)
#LR
LR2 <- read.table("search3_ICGC_LR.txt",sep='\t',header = T,as.is = T)
pairs <- paste(LR2$ligand,LR2$receptor,sep='_')
LR2$curated <- union_LR[match(pairs,union_pairs),]$Evidence
colnames(LR2)[colnames(LR2)=="curated"] <- "Evidence"

cancers <- unique(LR2$cancer_type)
for(i in 1:length(cancers)){
  cc <- cancer_type_ICGC[cancer_type_ICGC$`Short name` == cancers[i],]$Rename2
  LR2[LR2$cancer_type == cancers[i],]$cancer_type <- cc
}
LR2_out <- LR2[,c(1:5,10,9,7,8)]
colnames(LR2_out)[7] <- "cancer" 
write.table(LR2_out,paste0(out,"/search3_ICGC_LR.txt"),sep='\t',quote = F,row.names = F)
####EMBL####
cancer_type_EMBL <- cancer_type[grepl("EMBL",cancer_type$Source),]
#verify
verify3 <- read.table("search3_EMBL_verify.txt",sep='\t',header = T,as.is = T)
table(verify3$Cancer)
verify3$cancer_type <- "-"
verify3[verify3$Cancer == "fibrosarcoma",]$cancer_type <- "FS"
verify3[verify3$Cancer == "nasopharyngeal carcinoma",]$cancer_type <- "NPC"
verify3[verify3$Cancer == "non-small-cell lung cancer",]$cancer_type <- "NSCLC"
verify3[verify3$Cancer == "ovarian serous adenocarcinoma",]$cancer_type <- "OVSA"
verify3[verify3$Cancer == "pancreatic ductal carcinoma(PDAC)",]$cancer_type <- "PDAC"
verify3[verify3$Cancer == "papillary thyroid carcinoma",]$cancer_type <- "PTC"
verify3[verify3$Cancer == "Small Cell Lung Cancer",]$cancer_type <- "SCLC"
verify3[verify3$Cancer == "uterine leiomyosarcoma",]$cancer_type <- "ULMS"

cancers <- unique(verify3$cancer_type)
for(i in 1:length(cancers)){
  cc <- cancer_type_EMBL[cancer_type_EMBL$`Short name` == cancers[i],]$Rename2
  verify3[verify3$cancer_type == cancers[i],]$cancer_type <- cc
}
verify3_out <- verify3[,c(1:6,10,8,9)]
colnames(verify3_out)[7] <- "Cancer" 
write.table(verify3_out,paste0(out,"/search3_EMBL_verify.txt"),sep='\t',quote = F,row.names = F)
#LR
LR3 <- read.table("search3_EMBL_LR.txt",sep='\t',header = T,as.is = T)
pairs <- paste(LR3$ligand,LR3$receptor,sep='_')
LR3$curated <- union_LR[match(pairs,union_pairs),]$Evidence
colnames(LR3)[colnames(LR3)=="curated"] <- "Evidence"

table(LR3$cancer)
LR3$cancer_type <- "-"
LR3[LR3$cancer == "fibrosarcoma",]$cancer_type <- "FS"
LR3[LR3$cancer == "nasopharyngeal carcinoma",]$cancer_type <- "NPC"
LR3[LR3$cancer == "non-small-cell lung cancer",]$cancer_type <- "NSCLC"
LR3[LR3$cancer == "ovarian serous adenocarcinoma",]$cancer_type <- "OVSA"
LR3[LR3$cancer == "pancreatic ductal carcinoma(PDAC)",]$cancer_type <- "PDAC"
LR3[LR3$cancer == "papillary thyroid carcinoma",]$cancer_type <- "PTC"
LR3[LR3$cancer == "Small Cell Lung Cancer",]$cancer_type <- "SCLC"
LR3[LR3$cancer == "uterine leiomyosarcoma",]$cancer_type <- "ULMS"
LR3[LR3$cancer == "neuroblastoma",]$cancer_type <- "NBS"

cancers <- unique(LR3$cancer_type)
for(i in 1:length(cancers)){
  cc <- cancer_type_EMBL[cancer_type_EMBL$`Short name` == cancers[i],]$Rename2
  LR3[LR3$cancer_type == cancers[i],]$cancer_type <- cc
}
LR3_out <- LR3[,c(1:5,9,10,7,8)]
colnames(LR3_out)[7] <- "cancer" 
write.table(LR3_out,paste0(out,"/search3_EMBL_LR.txt"),sep='\t',quote = F,row.names = F)
####scRNA####
cancer_type_scRNA <- cancer_type[!grepl("EMBL|TCGA|ICGC",cancer_type$Source),]
#LR
LR_scRNA <- read.table("search3_predicted_scRNA_new.txt",sep='\t',header = T,as.is = T)
pairs <- paste(LR_scRNA$ligand,LR_scRNA$receptor,sep='_')
LR_scRNA$curated <- union_LR[match(pairs,union_pairs),]$Evidence
colnames(LR_scRNA)[colnames(LR_scRNA)=="curated"] <- "Evidence"

cancers <- unique(LR_scRNA$cancer)
for(i in 1:length(cancers)){
  cc <- cancer_type_scRNA[cancer_type_scRNA$`Short name` == cancers[i],]$Rename2
  cc2 <- unique(cc)
  LR_scRNA[LR_scRNA$cancer == cancers[i],]$cancer <- cc2
}
write.table(LR_scRNA,paste0(out,"/search3_predicted_scRNA_new.txt"),sep='\t',quote = F,row.names = F)

LR_scRNA <- read.table(paste0(out,"/search3_predicted_scRNA_new.txt"),sep='\t',header = T,as.is = T)
scRNA_cancer_type <- readxl::read_excel(paste0(out,"/scRNA.xlsx"))
scRNA_cancer_type2 <- unique(scRNA_cancer_type[,c(1,5)])
LR_scRNA$short <- "-"
cancer <- unique(LR_scRNA$cancer)
for(i in 1:length(cancer)){
  short <- scRNA_cancer_type2[scRNA_cancer_type2$`name in table` %in% cancer[i],2]
  LR_scRNA[LR_scRNA$cancer %in% cancer[i],]$short <- as.character(short)
}
colnames(LR_scRNA)[6] <- "iTALK_top"
colnames(LR_scRNA)[11] <- "curated"
head(LR_scRNA)
write.table(LR_scRNA,paste0(out,"/search3_predicted_scRNA_table.txt"),sep="\t",quote = F,row.names = F)
