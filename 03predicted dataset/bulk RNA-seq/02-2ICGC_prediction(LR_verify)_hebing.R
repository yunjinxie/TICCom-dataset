#==================================整合ICGC=====================================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM/30_log_05/result")

d <- dir(getwd())
cancer0 <- d[!grepl(".txt",d)]
cancer0 <- cancer0[cancer0!="R"]
cancer <- setdiff(cancer0,c("LICA-FR","LIRI-JP")) 

####证实的####
file <- "TCGA_interaction_strength_gene_re.txt"

library(jsonlite)
library(dplyr)
write.table(cbind("Gene1Name","Gene2Name","cancer","Interaction_Strength","P_value"),
            paste0(getwd(),"/ICGC_interaction_strength_gene_result.txt"),sep='\t',quote = F,row.names = F,
            col.names = F,append = T)
for(i in 1:length(cancer)){
  filePath <- paste0(getwd(),"/",cancer[i],"/",file)
  
  data <- do.call(rbind,
                  lapply(paste(readLines(filePath,warn = F),
                                      collapse = ""),
                                jsonlite::fromJSON))
  data <- data[[1]]
  data$Cancer <- cancer[i]
  write.table(data,paste0(getwd(),"/ICGC_interaction_strength_gene_result.txt"),sep='\t',quote = F,row.names = F,
              col.names = F,append = T)
}

verify <- read.table("I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search1_new/search1_final/search_detail.txt",
                     sep='\t',header = T,as.is = T,fill=T)

verify2 <- verify %>% select(Gene1Name,Gene1Id,Cell,Gene2Name,Gene2Id,Cell_1 = Cell_1) %>% unique()
verify2$pair <- paste(verify2$Gene1Name,verify2$Gene2Name,sep='_')

data <- read.table(paste0(getwd(),"/ICGC_interaction_strength_gene_result.txt"),sep='\t',header = T,as.is = T,fill=T)
data$pair <- paste(data$Gene1Name,data$Gene2Name,sep='_')
data_merge <- merge(data,verify2,by = "pair")
data_merge2 <- data_merge %>% select(Gene1Name = Gene1Name.x,Gene1Id,Cell,Gene2Name = Gene2Name.x,Gene2Id,Cell_1,
                                     Cancer = cancer,Interaction_Strength,P_value) %>% 
  filter(as.numeric(P_value) < 0.05)
head(data_merge2)

data_merge2$cancer_type <- gsub("-.*","",data_merge2$Cancer)
write.table(data_merge2,paste0(getwd(),"/search3_ICGC_verify.txt"),sep='\t',quote = F,row.names = F)
####配体-受体####
file <- "ICCom_gene_re.txt"
library(jsonlite)
library(dplyr)
write.table(cbind("ligand","receptor","cancer","Interaction_Strength","P_value","Funcion"),
            paste0(getwd(),"/ICGC_ICCom_result.txt"),sep='\t',quote = F,row.names = F,
            col.names = F,append = T)
for(i in 1:length(cancer)){
  filePath <- paste0(getwd(),"/",cancer[i],"/",file)
  
  data <- do.call(rbind,
                  lapply(paste(readLines(filePath,warn = F),
                               collapse = ""),
                         jsonlite::fromJSON))
  data <- data[[1]]
  data$Cancer <- cancer[i]
  write.table(data,paste0(getwd(),"/ICGC_ICCom_result.txt"),sep='\t',quote = F,row.names = F,
              col.names = F,append = T)
}

####symbol2ENSG####
trans <- read.table("I:/F/tumor_immune_interaction/computation/A_webserver/symbolChangeData/ID_transfer_file.txt",
                    sep='\t',header = T,as.is = T)
data <- read.table(paste0(getwd(),"/ICGC_ICCom_result.txt"),sep='\t',header = T,as.is = T)
data_merge <- merge(data,trans,by.x = "ligand",by.y="Gene.name") %>%
  select(ligand,Ensembl_1 = Ensembl,receptor,Funcion,cancer,Interaction_Strength,P_value)
data_merge2 <- merge(data_merge,trans,by.x = "receptor",by.y = "Gene.name") %>%
  select(ligand,Ensembl_1,receptor,Ensembl_2 = Ensembl,Funcion,cancer,Interaction_Strength,P_value) %>%
  filter(as.numeric(P_value) < 0.05)

data_merge2$cancer_type <- gsub("-.*","",data_merge2$cancer)

write.table(data_merge2,paste0(getwd(),"/search3_ICGC_LR.txt"),sep='\t',quote = F,row.names = F)
