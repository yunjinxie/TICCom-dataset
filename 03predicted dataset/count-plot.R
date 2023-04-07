####===================================download画图=======================================
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search3_new/search3_update/data_update"
####====================================================bulk====================================
####预测证实####
verified <- read.table(paste0(d,"/search3_bulk_verify_final.txt"),sep='\t',header = T,as.is = T)
verified_cell1 <- verified[,c(1,4,3,11)]
colnames(verified_cell1) <- c("Gene1","Gene2","cell_type","cancer")
verified_cell2 <- verified[,c(1,4,6,11)]
colnames(verified_cell2) <- c("Gene1","Gene2","cell_type","cancer")
verified2 <- unique(rbind(verified_cell1,verified_cell2))
verified2[verified2$cell_type == "APC",]$cell_type <- "DC cells"
library(plyr)
count_v <- ddply(verified2,.(cell_type,cancer),nrow)
# #堆积条形图
# cell <- unique(count_v$cell_type)
# cell_order <- c("T cells","B cells","DC cells","macrophage","MDSCs","neutrophils","mast cells","NK cells",
#                 "leukocytes","cancer cell")
# count_v$cell_type <- factor(count_v$cell_type,levels=cell_order)
# cols <- brewer.pal(12,"Set3")[c(1:8,10,11)]
# names(cols) <- cell_order
# pdf(paste0(d,"search3_Verified.pdf"))
# ggplot(count_v,aes(x= cancer,y=V1,fill=cell_type))+geom_bar(stat="identity",position = "stack")+coord_flip()+
#   labs(y="Number of tumor-immune cell interactions",x="")+
#   theme(axis.text.x = element_text(colour = "black",size=12,hjust = 1),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.title.y = element_text(size=12,colour = "black"))+
#   scale_fill_manual(values=cols)
# 
# dev.off()
####预测LR####
library(ggplot2)
library(ggforce)
library(RColorBrewer)
display.brewer.all()

LR <- read.table(paste0(d,"/search3_bulk_LR_final.txt"),sep="\t",header = T,as.is = T)
cancer <- unique(LR$cancer_type)
cancer_v_diff <- setdiff(cancer,unique(count_v$cancer))
# length(cancer_v_diff)
# #0
# cc <- matrix(c("A","A",cancer_v_diff,0,0),nrow=2,byrow = F)
# colnames(cc) <- c("cell_type","cancer","V1")
# count_v2 <- rbind(count_v,cc)
count_v2 <- count_v
count_cancer <- ddply(count_v2,.(cancer),function(x){sum(as.numeric(x[,3]))})
count_cancer_sort <- count_cancer[order(count_cancer$V1,decreasing = F),]
count_v2$cancer <- factor(count_v2$cancer,levels = count_cancer_sort$cancer)

cell_order <- c("T cells","B cells","DC cells","macrophage","MDSCs","neutrophils","mast cells","NK cells","leukocytes","lymphocytes","cancer cell")
count_v2$cell_type <- factor(count_v2$cell_type,levels = cell_order)
cols <- brewer.pal(12,"Set3")[c(1:8,10,11,9)]
names(cols) <- cell_order
#证实的
pdf(paste0(d,"/search3_bulk_verify.pdf"))
ggplot(count_v2,aes(x= cancer,y=as.numeric(V1),fill=cell_type))+geom_bar(stat="identity",position = "stack")+coord_flip()+
  theme(axis.text.x = element_text(colour = "black",size=12,hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12,colour = "black"))+
  labs(y="Number of tumor-immune cell interactions")+
  scale_fill_manual(values=cols)+
  ylim(60,0)
dev.off()
#LR
LR <- read.table(paste0(d,"/search3_bulk_LR_final.txt"),sep="\t",header = T,as.is = T)
count_LR <- ddply(LR,.(cancer_type,Evidence),nrow)
count_LR$cancer <- factor(count_LR$cancer,levels = count_cancer_sort$cancer)
pdf(paste0(d,"/search3_bulk_LR.pdf"))
ggplot(count_LR,aes(x= cancer,y=as.numeric(V1),fill=Evidence))+geom_bar(stat="identity",position = "stack")+coord_flip()+
  theme(axis.text.x = element_text(colour = "black",size=12,hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12,colour = "black"))+
  labs(y="Number of ligand-receptor interactions")+
  ylim(0,3000)+
  scale_fill_manual(values = c("manually curated" = "#FFDEAD","predicted" = "#FFF8DC"))
dev.off()
####每个癌症样本数
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search3_new/search3_update/data_update")
ct <- readxl::read_excel("TCGA+ICGC+EMBL+scRNA3_3.xlsx")

#TCGA
count_t <- read.table("sample_count/TCGA_count_info.txt",sep='\t',header = T,as.is = T)
ct_t <- ct[grepl("TCGA",ct$Source),]
count_t$cancer_type <- ct_t[match(count_t$cancer,ct_t$`Short name`),]$`Cancer Type`
count_t2 <- ddply(count_t,.(cancer_type),function(x) sum(x[,3]))
count_t2$source <- "TCGA"
colnames(count_t2)[2] <- "sample_num"
#ICGC
count_i <- read.table("sample_count/ICGC_count_info_30_log_05.txt",sep='\t',header = T,as.is = T)
count_i$cancer <- gsub("-.*","",count_i$cancer)
table(count_i$cancer)
count_i2 <- ddply(count_i,.(cancer),function(x) sum(x[,3]))
ct_i <- ct[grepl("ICGC",ct$Source),]
count_i2$cancer_type <- ct_i[match(count_i2$cancer,ct_i$`Short name`),]$`Cancer Type`
count_i3 <- ddply(count_i2,.(cancer_type),function(x) sum(x[,2]))
count_i3$source <- "ICGC"
colnames(count_i3)[2] <- "sample_num"

#EMBL
count_e1 <- read.table("sample_count/EMBL_cancer_count_info.txt",sep='\t',header = T,as.is = T)
count_e2 <- read.table("sample_count/EMBL_cancer-normal_count_info.txt",sep='\t',header = T,as.is = T)
count_e <- rbind(count_e1[,c(1,3)],count_e2[,c(1,3)])
count_e$short <- "-"
count_e[count_e$cancer == "fibrosarcoma",]$short <- "FS"
count_e[count_e$cancer == "nasopharyngeal carcinoma",]$short <- "NPC"
count_e[count_e$cancer == "neuroblastoma",]$short <- "NBS"
count_e[count_e$cancer == "non-small-cell lung cancer",]$short <- "NSCLC"
count_e[count_e$cancer == "ovarian serous adenocarcinoma",]$short <- "OVSA"
count_e[count_e$cancer == "pancreatic ductal carcinoma(PDAC)",]$short <- "PDAC"
count_e[count_e$cancer == "papillary thyroid carcinoma",]$short <- "PTC"
count_e[count_e$cancer == "Small Cell Lung Cancer",]$short <- "SCLC"
count_e[count_e$cancer == "uterine leiomyosarcoma",]$short <- "ULMS"

ct_e <- ct[grepl("EMBL",ct$Source),]
count_e$source <- ct_e[match(count_e$short,ct_e$`Short name`),]$Source2
count_e$cancer_type <- ct_e[match(count_e$short,ct_e$`Short name`),]$`Cancer Type`
count_e_2 <- count_e[,c(5,2,4)]

#all
all_count <- rbind(count_t2,count_i3,count_e_2)
write.table(all_count,"sample_count/all_source_count.txt",sep='\t',quote = F,row.names = F)
all_count <- read.table("sample_count/all_source_count.txt",sep='\t',header = T,as.is = T)

library(ggplot2)
cols <- brewer.pal(12,"Paired")[1:11]
names(cols) <- unique(all_count$source)

pdf("sample_count/all_source_count.pdf")
ggplot(all_count,aes(x=cancer_type,y=sample_num,fill= source))+geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,hjust=1),
        panel.grid = element_blank())+
  labs(x="",y="Sample number")+
  scale_fill_manual(values= cols)
dev.off()
#===============================================scRNA========================================
#=======================单细胞数据细胞个数=======================
# d <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed"
# file <- dir(d)
# file2 <- file[!grepl("TISCH|\\..*",file)]
# 
# nrows <- c()
# for(i in 1:length(file2)){
#   data <- read.table(paste0(d,"/",file2[i],"/cellLabel.txt"),sep='\t',header = T,as.is = T)
#   nrows <- rbind(nrows,cbind(file2[i],nrow(data)))
# }
# 
# d2 <- "I:/F/tumor_immune_interaction/computation/A_webserver/single_cell_RNAseq_processed/TISCH"
# file_2 <- dir(d2)
# nrows2 <- c()
# for(i in 1:length(file_2)){
#   data <- read.table(paste0(d2,"/",file_2[i],"/cellLabel.txt"),sep='\t',header = T,as.is = T)
#   nrows2 <- rbind(nrows2,cbind(file_2[i],nrow(data)))
# }
# 
# nrows_f <- rbind(nrows,nrows2)
# write.table(nrows_f,paste0(d,"/count_cell_number.txt"),sep='\t',quote = F,row.names = F)
# 
# library(ggplot2)
# d3 <- "I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search3_new/search3_update"
# count <- read.table(paste0(d3,"/scRNA_cell_number.txt"),sep='\t',header = T,as.is = T)
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search3_new/search3_update/data_update"
count <- read.table(paste0(d,"/sample_count/scRNA_cell_number.txt"),sep='\t',header = T,as.is = T)
ct_sc <- ct[!grepl("TCGA|ICGC|EMBL",ct$Source),]
count2 <- ddply(count,.(Short.name),function(x) sum(x[,3]))
count2$cancer_type <- ct_sc[match(count2$Short.name,ct_sc$`Short name`),]$`Cancer Type`
count2_2 <- ddply(count2,.(cancer_type),function(x) sum(x[,2]))

cancer_sort <- count2_2[order(count2_2$V1,decreasing = F),"cancer_type"]
pdf(paste0(d,"/scRNA_cell_number.pdf"))
ggplot(count2_2,aes(x=factor(cancer_type,levels = cancer_sort),y=V1))+geom_bar(stat = "identity",fill= "#7EC2C6")+coord_flip()+
  labs(x="")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylim(70000,0)
dev.off()
####癌症-免疫细胞####
LR <- read.table(paste0(d,"/search3_predicted_scRNA_table.txt"),sep='\t',header = T,as.is = T)
LR$cancer_type <- ct_sc[match(LR$dataset,ct_sc$Source),]$`Cancer Type`
table(LR$cancer_type)

LR_cell1 <- LR[,c(1,3,2,15)]
colnames(LR_cell1) <- c("ligand","receptor","cell_type","cancer")
LR_cell2 <- LR[,c(1,3,4,15)]
colnames(LR_cell2) <- c("ligand","receptor","cell_type","cancer")
LR2 <- unique(rbind(LR_cell1,LR_cell2))
head(LR2)
cell_type <- unique(LR2$cell_type)
sort(cell_type)
cell_order <- c("T cells","B cells","DC cells","ILC","NK cells","Myeloid cells","cancer cells","other cells")

LR2$cell_type2 <- "-"
LR2[LR2$cell_type %in% c("AC-like Malignant","Basal tumor cells","cancer cell","Mal","Malignant",
                         "Malignant cell","NB-like Malignant","Neoplatic_cells_3","OC-like Malignant",
                         "OPC-like Malignant","Tumor"),]$cell_type2 <- "cancer cells"

LR2[LR2$cell_type %in% c("Acinar","Astrocyte","CAF","Ductal","earlyEry","Endothelial","Endothelial cells",
                         "Endothelial_cells","Epithelial","Epithelial_cells","Fibroblasts","GMP","HPC-like",
                         "HSC","lateEry","Muscle cells","Myofibroblasts","Neuron","Oligodendrocyte",
                         "Oligodendrocytes","Others","Pericytes","Plasma","Plasma_cells","Prog","Stromal",
                         "Urothelial cells"),]$cell_type2 <- "other cells"

LR2[LR2$cell_type %in% c("B","B-memory","B-naive","B cell","B.cell","B_cells","Bcell","ProB"),]$cell_type2 <- "B cells"

LR2[LR2$cell_type %in% c("CD4-CCR7-TCF7","CD4-IL7R-MAL","CD4 T cell","CD4_T_cells","CD4Tconv",
                         "CD8-GNLY","CD8-SLC4A10-MAIT","CD8-transition","CD8 T cell","CD8_act_T_cells",
                         "CD8_ex_T_cells","CD8_mem_T_cells","CD8T","CD8Tex","CTL","Cycling-T","T","T cell",
                         "T cells","T.CD4","T.CD8","T.cell","Tcell","Tcell_prolif","TEC","Tprolif","Treg",
                         "Tregs"),]$cell_type2 <- "T cells"

LR2[LR2$cell_type %in% c("cDC","DC","DC-CD1C-cDC2","DC-CLEC9A-cDC1","DC-LAMP3","DCs","Dendritic","pDC",
                         "pDCs"),]$cell_type2 <- "DC cells"

LR2[LR2$cell_type %in% c("Macrophage","Macrophage-AREG","Macrophage-C1QC-PLTP","Macrophage-MARCO","Macrophage-SPP1-ACP5",
                         "Macrophages","Mast","Mono","Mono/Macro","Monocyte-FCGR3A-nonClassic","Monocyte-VCAN",
                         "Myeloid","Myeloid cell","Neutrophils","ProMono","TAM"),]$cell_type2 <- "Myeloid cells"

LR2[LR2$cell_type %in% c("NK","NK-FCGR3A","NK-KLRC1","NK cells","NK_cells"),]$cell_type2 <- "NK cells"
LR2[LR2$cell_type %in% "ILC",]$cell_type2 <- "ILC"
#Innate lymphoid cells (ILCs) are a heterogeneous population of non-B non-T lymphocytes that 
#originate from the common lymphoid progenitor but lack antigen-specific receptors
table(LR2$cell_type2)
#Innate lymphoid cells (ILCs)
library(plyr)
count <- ddply(LR2,.(cell_type2,cancer), nrow)

pdf(paste0(d,"/search3_predicted_scRNA_new.pdf"))
library(ggplot2)
ggplot(count,aes(x=cell_type2,y=factor(cancer,levels = cancer_sort),color=V1,size=V1))+geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=60,size=12,colour = "black",hjust=1),
        axis.text.y = element_text(size=12,colour = "black"),
        panel.grid = element_blank())+
  labs(x="",y="")+
  scale_color_gradientn(colours = colorRampPalette(c("#1597E5","#3CCF4E","#FFD124"))(20))
dev.off()
