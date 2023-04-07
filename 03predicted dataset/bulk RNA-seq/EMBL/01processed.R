#===========================EMBL_expression_atlas===============================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/EMBL_expression_atlas/")
####FPKM(疾病不好合并)###
####E-GEOD-26284####
id <- "E-GEOD-26284"
count <- read.table(paste0(id,"-query-results.fpkms.tsv"),sep='\t',header = T,as.is = T,check.names = F)
count[1:5,1:5]

samples <- read.table(paste0(id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T)
table(samples$Sample.Characteristic.disease.)
####chronic myelogenous leukemia (CML)####
CML_samples <- samples[samples$Sample.Characteristic.disease.== "chronic myelogenous leukemia (CML)",
                       c("Sample.Characteristic.disease.","Sample.Characteristic.cellular.component.","Factor.Value.RNA.","Sample.Characteristic.cell.line.")]
CML_samples2 <- paste(CML_samples[,2],CML_samples[,3],CML_samples[,4],sep=', ')
inter <- intersect(CML_samples2,colnames(count))
CML_count <- count[,colnames(count) %in% inter]
rownames(CML_count) <- count[,1]
CML_count[is.na(CML_count)] <- 0
CML_count <- as.data.frame(CML_count)
save(CML_count,file = "chronic myelogenous leukemia_count.RData")
write.table(CML_samples,"chronic myelogenous leukemia_samples.txt",sep='\t',quote = F,row.names = F)

####neuroblastoma####
neuroblastoma_samples <- samples[samples$Sample.Characteristic.disease.== "neuroblastoma",
                       c("Sample.Characteristic.disease.","Sample.Characteristic.cellular.component.","Factor.Value.RNA.","Sample.Characteristic.cell.line.")]
neuroblastoma_samples2 <- paste(neuroblastoma_samples[,2],neuroblastoma_samples[,3],neuroblastoma_samples[,4],sep=', ')
inter <- intersect(neuroblastoma_samples2,colnames(count))
neuroblastoma_count <- count[,colnames(count) %in% inter]
rownames(neuroblastoma_count) <- count[,1]
neuroblastoma_count[is.na(neuroblastoma_count)] <- 0
neuroblastoma_count <- as.data.frame(neuroblastoma_count)
save(neuroblastoma_count,file = "neuroblastoma_count.RData")
write.table(neuroblastoma_samples,"neuroblastoma_samples.txt",sep='\t',quote = F,row.names = F)
####normal####
normal_samples <- samples[samples$Sample.Characteristic.disease.== "normal",
                                 c("Sample.Characteristic.disease.","Sample.Characteristic.cellular.component.","Factor.Value.RNA.","Sample.Characteristic.cell.line.")]

#HMEC cell line : mammary
normal_samples_m <- normal_samples[normal_samples$Sample.Characteristic.cell.line.== "HMEC cell line",,drop=F]
normal_samples_m2 <- paste(normal_samples_m[,2],normal_samples_m[,3],normal_samples_m[,4],sep=', ')
inter <- intersect(normal_samples_m2,colnames(count))
normal_count_m <- count[,colnames(count) %in% inter]
rownames(normal_count_m) <- count[,1]
normal_count_m[is.na(normal_count_m)] <- 0
normal_count_m <- as.data.frame(normal_count_m)
save(normal_count_m,file = "normal_celllines_mammary_HMEC_count.RData")
write.table(normal_samples_m,"normal_celllines_mammary_HMEC_samples.txt",sep='\t',quote = F,row.names = F)

#IMR-90 : lung + NHLF cell line : lung
normal_samples_l <- normal_samples[normal_samples$Sample.Characteristic.cell.line. %in% c("IMR-90","NHLF cell line"),,drop=F]
normal_samples_l2 <- paste(normal_samples_l[,2],normal_samples_l[,3],normal_samples_l[,4],sep=', ')
inter <- intersect(normal_samples_l2,colnames(count))
normal_count_l <- count[,colnames(count) %in% inter]
rownames(normal_count_l) <- count[,1]
normal_count_l[is.na(normal_count_l)] <- 0
normal_count_l <- as.data.frame(normal_count_l)
save(normal_count_l,file = "normal_celllines_lung_IMR-90_NHLF_count.RData")
write.table(normal_samples_l,"normal_celllines_lung_IMR-90_NHLF_samples.txt",sep='\t',quote = F,row.names = F)

#NHEK cell line : 	normal human epidermal keratinocyte (skin)
normal_samples_s <- normal_samples[normal_samples$Sample.Characteristic.cell.line. =="NHEK cell line",,drop=F]
normal_samples_s2 <- paste(normal_samples_s[,2],normal_samples_s[,3],normal_samples_s[,4],sep=', ')
inter <- intersect(normal_samples_s2,colnames(count))
normal_count_s <- count[,colnames(count) %in% inter]
rownames(normal_count_s) <- count[,1]
normal_count_s[is.na(normal_count_s)] <- 0
normal_count_s <- as.data.frame(normal_count_s)
save(normal_count_s,file = "normal_celllines_skin_NHEK_count.RData")
write.table(normal_samples_s,"normal_celllines_skin_NHEK_samples.txt",sep='\t',quote = F,row.names = F)

####E-MTAB-2706####
id <- "E-MTAB-2706"
count <- readLines(paste0(id,"-query-results.fpkms.tsv"))

count <- read.table(paste0(id,"-query-results.fpkms.tsv"),sep='\t',header = F,as.is = T,check.names = F,skip = 4,fill=T)

count[1:5,1:5]

samples <- read.table(paste0(id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)

####E-MTAB-2770####
id <- "E-MTAB-2770"
count <- readLines(paste0(id,"-query-results.fpkms.tsv"))

samples <- read.table(paste0(id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)

####fibrosarcoma-E-MTAB-6823####
id <- "E-MTAB-6823"
samples <- read.table(paste0("fibrosarcoma/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)

count <- read.table(paste0("fibrosarcoma/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
####nasopharyngeal carcinoma-E-MTAB-7841####
id <- "E-MTAB-7841"
samples <- read.table(paste0("nasopharyngeal carcinoma/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)

count <- read.table(paste0("nasopharyngeal carcinoma/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
####neuroblastoma-E-MTAB-7025####
id <- "E-MTAB-7025"
samples <- read.table(paste0("neuroblastoma/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)

count <- read.table(paste0("neuroblastoma/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
####non-small-cell lung cancer-E-GEOD-81089####
id <- "E-GEOD-81089"
samples <- read.table(paste0("non-small-cell lung cancer/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)

count <- read.table(paste0("non-small-cell lung cancer/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
####ovarian serous adenocarcinoma-E-MTAB-7284####
id <- "E-MTAB-7284"
samples <- read.table(paste0("ovarian serous adenocarcinoma/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)

count <- read.table(paste0("ovarian serous adenocarcinoma/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
####pancreatic ductal carcinoma(PDAC)-E-GEOD-63776####
id <- "E-GEOD-63776"
samples <- read.table(paste0("pancreatic ductal carcinoma(PDAC)/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)

count <- read.table(paste0("pancreatic ductal carcinoma(PDAC)/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
####papillary thyroid carcinoma-E-GEOD-48850####
id <- "E-GEOD-48850"
samples <- read.table(paste0("papillary thyroid carcinoma/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)
cancer_samples <- samples[samples$Sample.Characteristic.disease. == "papillary thyroid carcinoma",]

count <- read.table(paste0("papillary thyroid carcinoma/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
cancer_count <- count[,colnames(count) %in% cancer_samples$Run]
cancer_count2 <- cbind(count[,c(1,2)],cancer_count)
write.table(cancer_count2,paste0("papillary thyroid carcinoma/",id,"-cancer_raw-counts.tsv"),sep='\t',quote = F,row.names = F)

####normal####
normal_samples <- samples[samples$Sample.Characteristic.disease. == "normal",]
normal_count <- count[,colnames(count) %in% normal_samples$Run]
normal_count2 <- cbind(count[,c(1,2)],normal_count)
write.table(normal_count2,paste0("papillary thyroid carcinoma/",id,"-normal_raw-counts.tsv"),sep='\t',quote = F,row.names = F)

####Small Cell Lung Cancer-E-GEOD-60052####
id <- "E-GEOD-60052"
samples <- read.table(paste0("Small Cell Lung Cancer/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)
cancer_samples <- samples[samples$Sample.Characteristic.disease. == "small cell lung carcinoma",]

count <- read.table(paste0("Small Cell Lung Cancer/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
cancer_count <- count[,colnames(count) %in% cancer_samples$Run]
cancer_count2 <- cbind(count[,c(1,2)],cancer_count)
write.table(cancer_count2,paste0("Small Cell Lung Cancer/",id,"-cancer_raw-counts.tsv"),sep='\t',quote = F,row.names = F)

####normal####
normal_samples <- samples[samples$Sample.Characteristic.disease. == "normal",]
normal_count <- count[,colnames(count) %in% normal_samples$Run]
normal_count2 <- cbind(count[,c(1,2)],normal_count)
write.table(normal_count2,paste0("Small Cell Lung Cancer/",id,"-normal_raw-counts.tsv"),sep='\t',quote = F,row.names = F)

####uterine leiomyosarcoma-E-MTAB-5762####
id <- "E-MTAB-5762"
samples <- read.table(paste0("uterine leiomyosarcoma/",id,"-experiment-design.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
table(samples$Sample.Characteristic.disease.)
cancer_samples <- samples[samples$Sample.Characteristic.disease. == "uterine leiomyosarcoma",]

count <- read.table(paste0("uterine leiomyosarcoma/",id,"-raw-counts.tsv"),sep='\t',header = T,as.is = T,fill=T,quote = "")
count[1:4,1:4]
cancer_count <- count[,colnames(count) %in% cancer_samples$Run]
cancer_count2 <- cbind(count[,c(1,2)],cancer_count)
write.table(cancer_count2,paste0("uterine leiomyosarcoma/",id,"-cancer_raw-counts.tsv"),sep='\t',quote = F,row.names = F)

####normal####
normal_samples <- samples[samples$Sample.Characteristic.disease. == "normal",]
normal_count <- count[,colnames(count) %in% normal_samples$Run]
normal_count2 <- cbind(count[,c(1,2)],normal_count)
write.table(normal_count2,paste0("uterine leiomyosarcoma/",id,"-normal_raw-counts.tsv"),sep='\t',quote = F,row.names = F)
