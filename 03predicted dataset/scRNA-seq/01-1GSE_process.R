setwd("I:/F/tumor_immune_interaction/computation/A_webserver/3cancer_data/GSE103322_HNSCC/celltalker")

cell <- read.table("celltalker_Cell_Lables_File.txt",sep='\t',header=T,as.is = T)
group <- read.table("celltalker_Comparison_Sample_Replicate_File.txt",sep='\t',header=T,as.is = T)

table(cell$cell_type)
cell2 <- cell[cell$cell_type %in% c("B cell","cancer_cell","Dendritic","Macrophage","Mast","T cell"),]

cell_group <- merge(cell2,group,by="cell")
frenq <- table(cell_group$cell_type,cell_group$compare_group)
frenq

####处理GSE103322_HNSCC####
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/GSE103322_HNSCC")
exp <- read.table("GSE103322_HNSCC_all_data.txt",sep='\t',header = T,as.is = T)
sample <- t(exp[1:5,])
sample <- cbind(rownames(sample),sample)
sample2 <- sample[-1,]
colnames(sample2) <- sample[1,]
colnames(sample2)[1] <- "cell"

write.table(sample2,"GSE103322_HNSCC_metafile.txt",sep='\t',quote = F,row.names = F)

sample3 <- sample2[!(as.numeric(sample2[,4])==0&as.numeric(sample2[,5])==0),]

table(sample3[,6])
sample4 <- sample3[!sample3[,6] %in% c("-Fibroblast","Endothelial","Fibroblast","myocyte"),]

a <- sample4[sample4[,4]==1,]
table(a[,6])

sample4[sample4[,4]==1,6] <- "cancer cell"
colnames(sample4)[3] <- "compare_group"
colnames(sample4)[6] <- "cell_type"

sample4[as.numeric(sample4[,3])==1,3] <- "Lymph node"
sample4[as.numeric(sample4[,3])==0,3] <- "Primary"

frenq <- table(sample4[,"cell_type"],sample4[,"compare_group"])
frenq
write.table(sample4,"GSE103322_HNSCC_CellLabel.txt",sep='\t',quote = F,row.names = F)

####exp####
exp2 <- exp[-c(1:5),c("X",sample4[,1])]
exp3 <- exp2[,-1]
rownames(exp3) <- exp2[,1]
write.table(exp3,"GSE103322_HNSCC_ExpressionData.txt",sep='\t',quote = F,row.names = F)
#===================================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/3cancer_data/GSE115978_melanoma/celltalker")

cell <- read.table("celltalker_Cell_Lables_File.txt",sep='\t',header=T,as.is = T)
group <- read.table("celltalker_Comparison_Sample_Replicate_File.txt",sep='\t',header=T,as.is = T)

table(cell$cell_type)
cell2 <- cell[cell$cell_type %in% c("B.cell","Mal","NK","Macrophage","T.CD4","T.CD8","T.cell"),]

cell_group <- merge(cell2,group,by="cell")
frenq <- table(cell_group$cell_type,cell_group$compare_group)
frenq

#===================================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/GSE131907_Lung_Cancer")

cell <- read.table("GSE131907_Lung_Cancer_cell_annotation.txt",sep='\t',header=T,as.is = T)

table(cell$Cell_type)

cell2 <- cell[cell$Cell_subtype %in% "Malignant cells",]
cell3 <- cell[!cell$Cell_type %in% c("Endothelial cells","Epithelial cells","Fibroblasts","Oligodendrocytes","Undetermined"),]
cell4 <- rbind(cell2,cell3)
cell5 <- cell4[!cell4$Cell_subtype %in% "Undetermined",]
cell6 <- cell5[cell5$Sample_Origin != "PE",]

sample <- read.table("GSE131907_sample.txt",sep='\t',header=T,as.is = T)
sample2 <- sample[,c(1,10:11)]

cell_sample <- merge(cell6,sample2,by.x="Sample",by.y = "Sample_title")
cell_sample[cell_sample$Cell_type %in% "Epithelial cells","Cell_type"] <- "Malignant cells"

#在肿瘤内，转移和非转移,恶性细胞和NK细胞
cell_sample2 <- cell_sample[!cell_sample[,4] %in% c("nLN","nLung"),]
cell_sample3 <- cell_sample2[cell_sample2[,4] %in% c("mLN","tL/B")&cell_sample2[,5] %in% c("Malignant cells","NK cells"),]
frenq <- table(cell_sample3[,5],cell_sample3[,4],cell_sample3[,8])
#保留大于50个细胞的patient
cell_sample4 <- cell_sample3[cell_sample3[,8] %in% c("patient id: P1006","patient id: P1015","patient id: P1019",
                                                     "patient id: P1049","patient id: P1051","patient id: P1058"),]

memory.limit(32634*3)
exp <- readRDS("GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds")
exp2 <- exp[,cell_sample4[,2]]
dim(exp2)
#29634  4985
rm(exp)
write.table(exp2,"GSE131907_Lung_Cancer_ExpressionData.txt",sep = "\t",quote = F)

cell_label <- cell_sample4[,c(2,5)]
colnames(cell_label) <- c("cell","cell_type")
write.table(cell_label,"GSE131907_Lung_Cancer_cellLabel.txt",sep='\t',quote=F,row.names=F)

group <- cell_sample4[,c(2,4,8)]
colnames(group) <- c("cell","compare_group","replicate")
group[,3] <- gsub("patient id: ","",group[,3])
write.table(group,"GSE131907_Lung_Cancer_ComparisonGroup.txt",sep='\t',quote=F,row.names=F)


frenq <- table(cell_sample4[,5],cell_sample4[,4],cell_sample4[,8])
frenq


