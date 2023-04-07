#============根据ITALK和ICELLNET对剩下5个数据集配体功能注释==============================
####取并集####
setwd("H:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data")


file_name <- c("cellchatDB_2","celltalkDB_LR_database","nichenet_LR_database",
               "ramilowski_pairs","singlecellsignalR_LR_database",
               "icellnet_LR_database","iTALK_LR_database")

datasets <- list()
nrows <- c()
for(i in 1:length(file_name)){
  data <- read.table(paste(file_name[i],".txt",sep=''),sep='\t',header=T,as.is = T)
  head(data)
  data2 <- data[,c("ligand","receptor")]
  data2$pair <- paste(data2$ligand,data2$receptor,sep='_')
  head(data2)
  data2 <- unique(data2)
  nrows <- c(nrows,dim(data2)[1])
  
  datasets[[i]] <- data2
}
nrows2 <- cbind(file_name,nrows)
write.table(nrows,"classification_results/nrow_LR_datasets.txt",sep='\t',quote = F,row.names = F,col.names = F)
names(datasets) <- file_name

union <- Reduce(rbind,datasets)
dim(union)
union2 <- unique(union)
dim(union2)
#[1] 14338     3
write.table(union2,"union_LR_database.txt",sep='\t',quote = F,row.names = F)
union2 <- as.data.frame(union2)
####交集####
datasets2 <- list()
for(j in 1:length(datasets)){
  datasets2[[j]] <- datasets[[j]][,3]
}
inter_set <- Reduce(intersect,datasets2)
inter_set2 <- do.call(rbind,strsplit(inter_set,"_"))
inter_set3 <- cbind(inter_set2,inter_set)
colnames(inter_set3) <- c("ligand","receptor","pair")
write.table(inter_set3,"inter_LR_database.txt",sep='\t',quote = F,row.names=F)
####保留icellnet和italk中功能注释结果，交集取icellnet的####
icellnet <- read.table("icellnet_LR_database.txt",sep='\t',header = T,as.is = T)
icellnet$pair <- paste(icellnet$ligand,icellnet$receptor,sep='_')
head(icellnet)
icellnet2 <- unique(icellnet)
dim(icellnet2)
#456

italk <- read.table("iTALK_LR_database.txt",sep='\t',header=T,as.is = T)
italk$pair <- paste(italk$ligand,italk$receptor,sep='_')
head(italk)
dim(italk)
italk2 <- unique(italk)
dim(italk2)
#[1] 2575    4
table(italk2[,4])
####首字母大写####
#stringr::str_to_title()
italk2[,4] <- unlist(lapply(italk2[,4],function(x){
  paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
}))
write.table(italk2,"classification_results/iTALK_LR_database.txt",sep='\t',quote = F,row.names = F)
####italk与icellnet互作对合并####
icell_italk_union <- union(icellnet2[,3],italk2[,3])
length(icell_italk_union)
#2654
#交集取icellnet注释的结果
icell_italk_inter <- intersect(icellnet2[,3],italk2[,3])
length(icell_italk_inter)
#377
icellnet_class1 <- icellnet2[which(icellnet2$pair %in% icell_italk_inter),]

only_in_icellnet <- setdiff(icellnet2[,3],icell_italk_inter)
icellnet_class2 <- icellnet2[which(icellnet2$pair %in% only_in_icellnet),]

only_in_italk <- setdiff(italk2[,3],icell_italk_inter)
italk_class1 <- italk2[which(italk2$pair %in% only_in_italk),]

icell_italk_in_union <- rbind(icellnet_class1,icellnet_class2,italk_class1)
dim(icell_italk_in_union)
# 2654    4
write.table(icell_italk_in_union,"classification_results/icell_italk_in_union_gene_class.txt",sep='\t',quote = F,row.names = F)
####需要重新注释的互作对####
union_no_icell_italk <- union2[which(!union2$pair %in% icell_italk_in_union[,3]),] 
len <- nrow(union_no_icell_italk)
iclass_ligand<-matrix(NA,len,10)
iclass_receptor<-matrix(NA,len,10)

dir <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/data/go_term"
file <- list.files(dir)
go_term <- gsub(".txt","",file)
path <- file.path(dir,file)

union_gene_class<-cbind(union_no_icell_italk,iclass_ligand)
union_gene_class<-cbind(union_gene_class,iclass_receptor)

all_genes <- read.table("H:/F/tumor_immune_interaction/computation/A_webserver/database/data/c5.go.v7.2.symbols.gmt",sep='\t',header = F,as.is = T,fill = T)
all_genes <- as.matrix(all_genes)

go_term_genes <- lapply(path,function(x){
  go <- read.table(x,sep='\t',header = F,as.is = T)
  go <- as.matrix(go)
  gene <- c(all_genes[which(all_genes[,1] %in% go),-c(1,2)])
  gene <- gene[gene!=""]
})
names(go_term_genes)<- go_term
########ligand
col_index <- 4
for(i in 1:length(go_term_genes)){
  union_gene_class[which(union_gene_class$ligand %in% go_term_genes[[i]]),col_index]<-go_term[i]
  col_index <- col_index+1
}
########receptor
col_index <- 14
for(i in 1:length(go_term_genes)){
  union_gene_class[which(union_gene_class$receptor %in% go_term_genes[[i]]),col_index]<-go_term[i]
  col_index <- col_index+1
}
colnames(union_gene_class)[-c(1:3)] <- c(paste("ligand",1:10,sep=''),paste("receptor",1:10,sep=''))
write.table(union_gene_class,"classification_results/union_no_icellnet_italk_gene_class.txt",sep='\t',quote = F,row.names = F)
dim(union_gene_class)
#[1] 11684    23
####确定最终功能####
#####优先级：Notch signaling，Antigen binding，Neuropeptide，Hormone，cytokine（growth factor，interferon，interleukin，tumor necrosis factor，Chemokine）
len2 <- nrow(union_gene_class)
final_union_class<-union_gene_class[,1:3]
class<-matrix(NA,len2,1)
final_union_class<-cbind(final_union_class,class)

final_union_class[which(union_gene_class[,"ligand9"]=="Notch signaling"&union_gene_class[,"receptor9"]=="Notch signaling"),4]<-"Notch signaling"
final_union_class[which(is.na(final_union_class$class)),4]<-"Other"

final_union_class[which(final_union_class$class=="Other"&union_gene_class[,"ligand1"]=="Antigen binding"&union_gene_class[,"receptor1"]=="Antigen binding"),4]<-"Antigen binding"

final_union_class[which(final_union_class$class=="Other"&union_gene_class[,"ligand8"]=="Neuropeptide"&union_gene_class[,"receptor8"]=="Neuropeptide"),4]<-"Neuropeptide"

final_union_class[which(final_union_class$class=="Other"&union_gene_class[,"ligand5"]=="Hormone"&union_gene_class[,"receptor5"]=="Hormone"),4]<-"Hormone"
####筛选cytokine####
final_union_class[which(final_union_class$class=="Other"&union_gene_class[,"ligand2"]=="Chemokine"&union_gene_class[,"receptor2"]=="Chemokine"&is.na(union_gene_class[,"ligand4"])&
                          is.na(union_gene_class[,"receptor4"])&is.na(union_gene_class[,"ligand6"])&is.na(union_gene_class[,"receptor6"])&is.na(union_gene_class[,"ligand7"])&
                          is.na(union_gene_class[,"receptor7"])&is.na(union_gene_class[,"ligand10"])&is.na(union_gene_class[,"receptor10"])),4]<-"Chemokine"

final_union_class[which(final_union_class$class=="Other"&union_gene_class[,"ligand4"]=="Growth factor"&union_gene_class[,"receptor4"]=="Growth factor"&is.na(union_gene_class[,"ligand2"])&is.na(union_gene_class[,"receptor2"])&
                          is.na(union_gene_class[,"ligand6"])&is.na(union_gene_class[,"receptor6"])&is.na(union_gene_class[,"ligand7"])&
                          is.na(union_gene_class[,"receptor7"])&is.na(union_gene_class[,"ligand10"])&is.na(union_gene_class[,"receptor10"])),4]<-"Growth factor"

final_union_class[which(final_union_class$class=="Other"&union_gene_class[,"ligand6"]=="Interferon"&union_gene_class[,"receptor6"]=="Interferon"&is.na(union_gene_class[,"ligand4"])&is.na(union_gene_class[,"receptor4"])&
                          is.na(union_gene_class[,"ligand2"])&is.na(union_gene_class[,"receptor2"])&is.na(union_gene_class[,"ligand7"])&
                          is.na(union_gene_class[,"receptor7"])&is.na(union_gene_class[,"ligand10"])&is.na(union_gene_class[,"receptor10"])),4]<-"Interferon"

final_union_class[which(final_union_class$class=="Other"& union_gene_class[,"ligand7"]=="Interleukin"&union_gene_class[,"receptor7"]=="Interleukin"&is.na(union_gene_class[,"ligand2"])&is.na(union_gene_class[,"receptor2"])&
                          is.na(union_gene_class[,"ligand4"])&is.na(union_gene_class[,"receptor4"])&is.na(union_gene_class[,"ligand6"])&
                          is.na(union_gene_class[,"receptor6"])&is.na(union_gene_class[,"ligand10"])&is.na(union_gene_class[,"receptor10"])),4]<-"Interleukin"

final_union_class[which(final_union_class$class=="Other"&union_gene_class[,"ligand10"]=="Tumor necrosis factor"&union_gene_class[,"receptor10"]=="Tumor necrosis factor"&is.na(union_gene_class[,"ligand2"])&is.na(union_gene_class[,"receptor2"])&
                          is.na(union_gene_class[,"ligand4"])&is.na(union_gene_class[,"receptor4"])&is.na(union_gene_class[,"ligand6"])&
                          is.na(union_gene_class[,"receptor6"])&is.na(union_gene_class[,"ligand7"])&is.na(union_gene_class[,"receptor7"])),4]<-"Tumor necrosis factor"

final_union_class[which(final_union_class$class=="Other"&union_gene_class[,"ligand3"]=="Cytokine"&union_gene_class[,"receptor3"]=="Cytokine"),4]<- "Cytokine"

final_union_class[which(substr(final_union_class$ligand,1,3)=="TNF"&substr(final_union_class$receptor,1,3)=="TNF"),4]<-"Tumor necrosis factor"
colnames(final_union_class)[4] <- "classification"
write.table(final_union_class,"classification_results/union_no_icellnet_italk_gene_class_final.txt",sep='\t',quote = F,row.names = F)

####整合####
final_union_class2 <- rbind(final_union_class,icell_italk_in_union)
write.table(final_union_class2,"classification_results/union_LR_database.txt",sep='\t',quote = F,row.names = F)
#####提取单个数据集分类结果####
file_name <- c("cellchatDB_2","celltalkDB_LR_database","nichenet_LR_database",
               "ramilowski_pairs","singlecellsignalR_LR_database",
               "icellnet_LR_database","iTALK_LR_database")

cellchatDB_class <- final_union_class2[which(final_union_class2$pair %in% datasets[[1]][,3]),]
dim(cellchatDB_class)

celltalkDB_class <- final_union_class2[which(final_union_class2$pair %in% datasets[[2]][,3]),]
dim(celltalkDB_class)

nichenet2_class <- final_union_class2[which(final_union_class2$pair %in% datasets[[3]][,3]),]
dim(nichenet2_class)

ramilowski_class <- final_union_class2[which(final_union_class2$pair %in% datasets[[4]][,3]),]
dim(ramilowski_class)

singlecellsignalR_class <- final_union_class2[which(final_union_class2$pair %in% datasets[[5]][,3]),]
dim(singlecellsignalR_class)

write.table(cellchatDB_class,"classification_results/cellchatDB_LR_database.txt",sep='\t',quote = F,row.names = F)
write.table(celltalkDB_class,"classification_results/celltalkDB_LR_database.txt",sep='\t',quote = F,row.names = F)
write.table(nichenet2_class,"classification_results/nichenet_LR_database.txt",sep='\t',quote = F,row.names = F)
write.table(ramilowski_class,"classification_results/ramilowski_LR_database.txt",sep='\t',quote = F,row.names = F)
write.table(singlecellsignalR_class,"classification_results/singlecellsignalR_LR_database.txt",sep='\t',quote = F,row.names = F)
####交集注释结果####
inter_set3_class <- final_union_class2[which(final_union_class2$pair %in% inter_set3[,3]),]
dim(inter_set3_class)
#294   4
write.table(inter_set3_class,"classification_results/inter_LR_database.txt",sep='\t',quote = F,row.names = F)
