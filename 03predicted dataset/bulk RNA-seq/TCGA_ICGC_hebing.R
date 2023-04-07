#===================TCGA(LGG+GBM)合并，批次效应========================
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/TCGA_test_data"

LGG_exp <- read.table(paste0(d,"/LGG_cancer_genes_expression_30_05_log2.txt"),sep = "\t",header = T,as.is = T)
GBM_exp <- read.table(paste0(d,"/GBM_cancer_genes_expression_30_05_log2.txt"),sep = "\t",header = T,as.is = T)

LGG_exp[1:3,1:3]
LGG_exp2 <- cbind(rownames(LGG_exp),LGG_exp)
colnames(LGG_exp2)[1] <- "gene"

GBM_exp[1:3,1:3]
GBM_exp2 <- cbind(rownames(GBM_exp),GBM_exp)
colnames(GBM_exp2)[1] <- "gene"

combined <- merge(LGG_exp2,GBM_exp2,by="gene")
combined[1:3,1:3]

####批次效应
library(ggplot2)
library(ggbiplot)
library(devtools)
install_github("vqv/ggbiplot")

group = ifelse(colnames(combined[,-1]) %in% colnames(LGG_exp),'LGG','GBM')
barcode=colnames(combined[,-1])

design = data.frame('Barcode'=barcode,'Group'=group)

exp_m_t=as.data.frame(t(combined[,-1]))
pca_result <- prcomp(exp_m_t,scale=T)

pdf("I:/F/tumor_immune_interaction/computation/A_webserver/pca_TCGA_LGG_GBM.pdf")
ggbiplot(pca_result, 
         var.axes=F,            # 是否为变量画箭头
         obs.scale = 1,         # 横纵比例 
         groups = design$Group, # 添加分组信息，将按指定的分组信息上色
         ellipse = T,           # 是否围绕分组画椭圆
         circle = F)
dev.off()

#combine
combined2 <- combined[,-1]
rownames(combined2) <- combined[,1]
save(combined2,file = "I:/F/tumor_immune_interaction/computation/A_webserver/Rfunction/function_relay/GBMLGG_cancer_genes_expression_30_05_log2.RData")

##============================ICGC liver cancer合并===============================
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM/30_log_05/reference"
LICA_expl <- load(paste0(d,"/LICA-FR_tumour_30_log_05.RData"))
LICA_exp <- eval(parse(text = LICA_expl))

LIRI_expl <- load(paste0(d,"/LIRI-JP_tumour_30_log_05.RData"))
LIRI_exp <- eval(parse(text = LIRI_expl))

LICA_exp[1:3,1:3]
LICA_exp2 <- cbind(rownames(LICA_exp),LICA_exp)
colnames(LICA_exp2)[1] <- "gene"

LIRI_exp[1:3,1:3]
LIRI_exp2 <- cbind(rownames(LIRI_exp),LIRI_exp)
colnames(LIRI_exp2)[1] <- "gene"

combined <- merge(LICA_exp2,LIRI_exp2,by="gene")
combined[1:3,1:3]

group = ifelse(colnames(combined[,-1]) %in% colnames(LICA_exp),'LICA','LIRI')
barcode=colnames(combined[,-1])

design = data.frame('Barcode'=barcode,'Group'=group)

exp_m_t=as.data.frame(t(combined[,-1]))
pca_result <- prcomp(exp_m_t,scale=T)

pdf("I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM/30_log_05/reference/pca_ICGC_LICA_LIRI.pdf")
ggbiplot(pca_result, 
         var.axes=F,            # 是否为变量画箭头
         obs.scale = 1,         # 横纵比例 
         groups = design$Group, # 添加分组信息，将按指定的分组信息上色
         ellipse = T,           # 是否围绕分组画椭圆
         circle = F)
dev.off()

combined2 <- combined[,-1]
rownames(combined2) <- combined[,1]
save(combined2,file = "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM/30_log_05/reference/LICA_LIRI_tumour_30_log_05.RData")
