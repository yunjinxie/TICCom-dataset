#==============search3 bulk====================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search3_new/search3_update/data_update")
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/download/cancer_catagory"
cancer_type <- readxl::read_excel(paste0(d,"/TCGA+ICGC+EMBL+scRNA3_3.xlsx"))

write.table(cbind("Gene1Name","Gene1Id","Cell","Gene2Name","Gene2Id","Cell_1","cancer",
                  "Interaction_Strength","P_value","source"),"search3_bulk_verify.txt",sep='\t',quote = F,row.names = F,col.names = F,append = T)

write.table(cbind("ligand","Ensembl_1","receptor","Ensembl_2","Funcion","Evidence","cancer","Interaction_Strength",
                  "P_value","source"),"search3_bulk_LR.txt",sep='\t',quote = F,row.names = F,col.names = F,append = T)

####TCGA####
verify <- read.table("search3_TCGA_verify.txt",sep='\t',header = T,as.is = T)
head(verify)
verify$source <- "TCGA"
write.table(verify,"search3_bulk_verify.txt",sep='\t',quote = F,row.names = F,col.names = F,append = T)

LR <- read.table("search3_TCGA_LR.txt",sep='\t',header = T,as.is = T)
head(LR)
LR$source <- "TCGA"
write.table(LR,"search3_bulk_LR.txt",sep='\t',quote = F,row.names = F,col.names = F,append = T)

cancers <- union(verify$Cancer,LR$cancer)
cc1 <- cancer_type[match(cancers,cancer_type$Rename2),]$`Cancer Type`
asso1 <- cbind(cancers,cc1,"TCGA")

####ICGC####
verify <- read.table("search3_ICGC_verify.txt",sep='\t',header = T,as.is = T)
head(verify)
verify$source <- "ICGC"
write.table(verify,"search3_bulk_verify.txt",sep='\t',quote = F,row.names = F,col.names = F,append = T)

LR <- read.table("search3_ICGC_LR.txt",sep='\t',header = T,as.is = T)
head(LR)
LR$source <- "ICGC"
head(LR)
write.table(LR,"search3_bulk_LR.txt",sep='\t',quote = F,row.names = F,col.names = F,append = T)

cancers <- union(verify$Cancer,LR$cancer)
cc2 <- cancer_type[match(cancers,cancer_type$Rename2),]$`Cancer Type`
asso2 <- cbind(cancers,cc2,"ICGC")

####EMBL####
verify <- read.table("search3_EMBL_verify.txt",sep='\t',header = T,as.is = T)
head(verify)

LR <- read.table("search3_EMBL_LR.txt",sep='\t',header = T,as.is = T)
head(LR)
cancers <- union(verify$Cancer,LR$cancer)
cancer_type2 <- cancer_type[cancer_type$Source == "EMBL",]
cc3 <- cancer_type2[match(cancers,cancer_type2$Rename2),]
asso3 <- cbind(cancers,cc3$`Cancer Type`,cc3$Source2)

verify$source <- cc3[match(verify$Cancer,cc3$Rename2),]$Source2
write.table(verify,"search3_bulk_verify.txt",sep='\t',quote = F,row.names = F,col.names = F,append = T)

LR$source <- cc3[match(LR$cancer,cc3$Rename2),]$Source2
write.table(LR,"search3_bulk_LR.txt",sep='\t',quote = F,row.names = F,col.names = F,append = T)
####hebing####
search3_asso <- rbind(unname(asso1),unname(asso2),unname(asso3))
colnames(search3_asso) <- c("name_in_table","cancer","source")
write.table(search3_asso,"search3_bulk_association.txt",sep='\t',quote = F,row.names = F)

search3_asso2 <- as.data.frame(search3_asso)
search3_asso2$organ <- "-"
for(i in 1:nrow(search3_asso2)){
  index <- which(cancer_type$Rename2 %in% search3_asso2[i,1]&
                   grepl(search3_asso2[i,3],cancer_type$Source2))
  organ <- unique(cancer_type[index,,drop=F]$Organ)
  search3_asso2[i,]$organ <- organ
}
write.table(search3_asso2,"search3_bulk_organ.txt",sep='\t',quote = F,row.names = F)

####加两列####
search3_asso2 <- read.table("search3_bulk_organ.txt",sep='\t',header = T,as.is = T)
verify <- read.table("search3_bulk_verify.txt",sep='\t',header = T,as.is = T)
verify$cancer_type <- "-"
verify$tissue <- "-"
source <-unique(verify$source)
for(i in 1:length(source)){
  s <- source[i]
  v <- verify[verify$source == s,]
  cancer <- unique(v$cancer)
  for(j in 1:length(cancer)){
    ct <- search3_asso2[search3_asso2$source == s&search3_asso2$name_in_table == cancer[j],]$cancer
    tis <- search3_asso2[search3_asso2$source == s&search3_asso2$name_in_table == cancer[j],]$organ
    verify[verify$source == s&verify$cancer == cancer[j],]$cancer_type <- ct
    verify[verify$source == s&verify$cancer == cancer[j],]$tissue <- tis
  }
}
write.table(verify,"search3_bulk_verify_final.txt",sep='\t',quote = F,row.names = F)
#
LR2 <- read.table("search3_bulk_LR.txt",sep='\t',header = T,as.is = T)
LR2$cancer_type <- "-"
LR2$tissue <- "-"
source <-unique(LR2$source)
for(i in 1:length(source)){
  s <- source[i]
  v <- LR2[LR2$source == s,]
  cancer <- unique(v$cancer)
  for(j in 1:length(cancer)){
    ct <- search3_asso2[search3_asso2$source == s&search3_asso2$name_in_table == cancer[j],]$cancer
    tis <- search3_asso2[search3_asso2$source == s&search3_asso2$name_in_table == cancer[j],]$organ
    LR2[LR2$source == s&LR2$cancer == cancer[j],]$cancer_type <- ct
    LR2[LR2$source == s&LR2$cancer == cancer[j],]$tissue <- tis
  }
}
write.table(LR2,"search3_bulk_LR_final.txt",sep='\t',quote = F,row.names = F)

#####关联文件####
head(search3_asso2)

ss <- c()
cancer <- unique(search3_asso2$cancer)
cancer_sort <- sort(cancer)
for(i in 1:length(cancer_sort)){
  x <- search3_asso2[search3_asso2$cancer == cancer_sort[i],,drop=F]
  s <- paste(unique(x$source),collapse = ",")
  ss <- rbind(ss,paste('"',cancer_sort[i],'": ','"',s,'"',sep=''))
}
ss_output <- paste(ss,collapse = ",\n")
ss_output2 <- paste('{',ss_output,'}',sep='\n')
write.table(ss_output2,"bulk_parts2.txt",sep='\t',quote = F,row.names = F,col.names = F)
####name：value
nv <- c()
for(i in 1:length(cancer_sort)){
  nn <- paste('{name: ',"'",cancer_sort[i],"', ",'value: ',"'",cancer_sort[i],"'}",sep='')
  nv <- rbind(nv,nn)
}
nv2 <- paste(nv,collapse =',\n')
nv2_2 <- paste('[',nv2,']',sep='')
write.table(nv2_2,"bulk_name_value2.txt",sep='\t',quote = F,row.names = F,col.names = F)

#单细胞
scRNA <- readxl::read_excel("scRNA.xlsx")
ss <- c()
cancer <- unique(scRNA$`Cancer Type`)
for(i in 1:length(cancer)){
  x <- scRNA[scRNA$`Cancer Type` == cancer[i],,drop=F]
  s <- paste(unique(x$Source2),collapse = ",")
  ss <- rbind(ss,paste('"',cancer[i],'": ','"',s,'"',sep=''))
}
ss_output <- paste(ss,collapse = ",\n")
ss_output2 <- paste('{',ss_output,'}',sep='\n')
write.table(ss_output2,"singlecell_parts.txt",sep='\t',quote = F,row.names = F,col.names = F)

