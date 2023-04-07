#============================对小鼠LR数据集功能分类===========================
####提取小鼠功能####
setwd("H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse")

gene2go_mouse <- read.table("mouse_gene2go",sep='\t',header=F,as.is = T,fill=T,strip.white = T,quote = "")
gene2go_mouse <- unique(gene2go_mouse)
go <- unique(gene2go_mouse[,2])

id_transfer <- read.table("mouse_ID_transfer_file.txt",sep='\t',header=T,as.is = T,check.names=F)

go_term_dir <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/分类结果/new/go_term"
file <- list.files(go_term_dir)
path <- file.path(go_term_dir,file)

for(i in 1:length(path)){
  go <- read.table(path[i],sep='\t',header=F,as.is = T)
  go <- as.matrix(go)
  go2 <- gsub("GO_","",go)
  go2 <- gsub("_"," ",go2)
  go3 <- tolower(go2)
  
  cat(i)
  cat("#####################################\n")
  for(j in 1:length(go3)){
    gene <- gene2go_mouse[grepl(paste("^",go3[j],"$",sep=''),gene2go_mouse[,2],ignore.case = T),1,drop=F]
    
    if(nrow(gene)!=0){
      colnames(gene) <- "Entrez ID"
      gene2symbol <- merge(gene,id_transfer,by="Entrez ID")
      symbol <- gene2symbol[,"Gene name"]
      output_s <- c(go[j],symbol)
      write.table(t(output_s),"mouse_go_symbol.txt",sep = "\t",quote = F,row.names = F,col.names = F,append = T)
    }
  }
}

####功能注释####
setwd("H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse")

output_dir <- "classification"
dir.create(output_dir)
####RNAMagnet####
RNAMagnet <- read.table("ensmbl/RNAMagnet_split.txt",sep='\t',header=T,as.is = T)
# RNAMagnet2 <- RNAMagnet[,c(1,3)]
# pair <- paste(RNAMagnet2[,1],RNAMagnet2[,2],sep='_')
# identical(RNAMagnet$pair,pair)
# #[1] TRUE
##############
go_term_dir <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/data/go_term"
file <- list.files(go_term_dir)
path <- file.path(go_term_dir,file)

go_term <- gsub(".txt","",file)

all_genes <- readLines("mouse_go_symbol.txt")
all_go <- lapply(all_genes,function(x){unlist(strsplit(x,"\t"))[1]})
all_go <- do.call(rbind,all_go)

go_term_genes <- lapply(path,function(x){
  go <- read.table(x,sep='\t',header = F,as.is = T)
  go <- as.matrix(go)
  
  index <- which(all_go %in% go)
  if(length(index)!=0){
    genes <- lapply(all_genes[index],function(x) unlist(strsplit(x,"\t"))[-1])
    genes <- unlist(genes)
  }
})
names(go_term_genes)<- go_term

len <- nrow(RNAMagnet)
RNAMagnet_ligand<-matrix(NA,len,10)
RNAMagnet_receptor<-matrix(NA,len,10)

RNAMagnet_class<-cbind(RNAMagnet,RNAMagnet_ligand)
RNAMagnet_class<-cbind(RNAMagnet_class,RNAMagnet_receptor)

########ligand
col_index <- 6
for(i in 1:length(go_term_genes)){
  RNAMagnet_class[which(RNAMagnet_class[,"ligand"] %in% go_term_genes[[i]]),col_index]<-go_term[i]
  col_index <- col_index+1
}
########receptor
col_index <- 16
for(i in 1:length(go_term_genes)){
  RNAMagnet_class[which(RNAMagnet_class[,"receptor"] %in% go_term_genes[[i]]),col_index]<-go_term[i]
  col_index <- col_index+1
}
colnames(RNAMagnet_class)[-c(1:5)] <- c(paste("ligand",1:10,sep=''),paste("receptor",1:10,sep=''))
write.table(RNAMagnet_class,"RNAMagnet_class.txt",sep='\t',quote = F,row.names = F)
dim(RNAMagnet_class)

#####优先级：Notch signaling，Antigen binding，Neuropeptide，Hormone，cytokine（growth factor，interferon，interleukin，tumor necrosis factor,Chemokine)
len2 <- nrow(RNAMagnet_class)
final_class<-RNAMagnet_class[,1:5]
class<-matrix(NA,len2,1)
final_class<-cbind(final_class,class)

final_class[which(RNAMagnet_class[,"ligand9"]=="Notch signaling"&RNAMagnet_class[,"receptor9"]=="Notch signaling"),6]<-"Notch signaling"
final_class[which(is.na(final_class$class)),6]<-"Other"

final_class[which(final_class$class=="Other"&RNAMagnet_class[,"ligand1"]=="Antigen binding"&RNAMagnet_class[,"receptor1"]=="Antigen binding"),6]<-"Antigen binding"

final_class[which(final_class$class=="Other"&RNAMagnet_class[,"ligand8"]=="Neuropeptide"&RNAMagnet_class[,"receptor8"]=="Neuropeptide"),6]<-"Neuropeptide"

final_class[which(final_class$class=="Other"&RNAMagnet_class[,"ligand5"]=="Hormone"&RNAMagnet_class[,"receptor5"]=="Hormone"),6]<-"Hormone"
##########筛选cytokine############
final_class[which(final_class$class=="Other"&RNAMagnet_class[,"ligand2"]=="Chemokine"&RNAMagnet_class[,"receptor2"]=="Chemokine"&is.na(RNAMagnet_class[,"ligand4"])&
                          is.na(RNAMagnet_class[,"receptor4"])&is.na(RNAMagnet_class[,"ligand6"])&is.na(RNAMagnet_class[,"receptor6"])&is.na(RNAMagnet_class[,"ligand7"])&
                          is.na(RNAMagnet_class[,"receptor7"])&is.na(RNAMagnet_class[,"ligand10"])&is.na(RNAMagnet_class[,"receptor10"])),6]<-"Chemokine"

final_class[which(final_class$class=="Other"&RNAMagnet_class[,"ligand4"]=="Growth factor"&RNAMagnet_class[,"receptor4"]=="Growth factor"&is.na(RNAMagnet_class[,"ligand2"])&is.na(RNAMagnet_class[,"receptor2"])&
                          is.na(RNAMagnet_class[,"ligand6"])&is.na(RNAMagnet_class[,"receptor6"])&is.na(RNAMagnet_class[,"ligand7"])&
                          is.na(RNAMagnet_class[,"receptor7"])&is.na(RNAMagnet_class[,"ligand10"])&is.na(RNAMagnet_class[,"receptor10"])),6]<-"Growth factor"

final_class[which(final_class$class=="Other"&RNAMagnet_class[,"ligand6"]=="Interferon"&RNAMagnet_class[,"receptor6"]=="Interferon"&is.na(RNAMagnet_class[,"ligand4"])&is.na(RNAMagnet_class[,"receptor4"])&
                    is.na(RNAMagnet_class[,"ligand2"])&is.na(RNAMagnet_class[,"receptor2"])&is.na(RNAMagnet_class[,"ligand7"])&
                    is.na(RNAMagnet_class[,"receptor7"])&is.na(RNAMagnet_class[,"ligand10"])&is.na(RNAMagnet_class[,"receptor10"])),6]<-"Interferon"

final_class[which(final_class$class=="Other"& RNAMagnet_class[,"ligand7"]=="Interleukin"&RNAMagnet_class[,"receptor7"]=="Interleukin"&is.na(RNAMagnet_class[,"ligand2"])&is.na(RNAMagnet_class[,"receptor2"])&
                    is.na(RNAMagnet_class[,"ligand4"])&is.na(RNAMagnet_class[,"receptor4"])&is.na(RNAMagnet_class[,"ligand6"])&
                    is.na(RNAMagnet_class[,"receptor6"])&is.na(RNAMagnet_class[,"ligand10"])&is.na(RNAMagnet_class[,"receptor10"])),6]<-"Interleukin"

final_class[which(final_class$class=="Other"&RNAMagnet_class[,"ligand10"]=="Tumor necrosis factor"&RNAMagnet_class[,"receptor10"]=="Tumor necrosis factor"&is.na(RNAMagnet_class[,"ligand2"])&is.na(RNAMagnet_class[,"receptor2"])&
                    is.na(RNAMagnet_class[,"ligand4"])&is.na(RNAMagnet_class[,"receptor4"])&is.na(RNAMagnet_class[,"ligand6"])&
                    is.na(RNAMagnet_class[,"receptor6"])&is.na(RNAMagnet_class[,"ligand7"])&is.na(RNAMagnet_class[,"receptor7"])),6]<-"Tumor necrosis factor"

final_class[which(final_class$class=="Other"&RNAMagnet_class[,"ligand3"]=="Cytokine"&RNAMagnet_class[,"receptor3"]=="Cytokine"),6]<- "Cytokine"

final_class[which(substr(final_class$ligand,1,3)=="Tnf"&substr(final_class$receptor,1,3)=="Tnf"),6]<-"Tumor necrosis factor"
colnames(final_class)[6] <- "classification"
write.table(final_class,"classification/RNAMagnetDB_LR_database.txt",sep='\t',quote = F,row.names = F)
####celltalkDB####
celltalkDB <- read.table("ensmbl/celltalkDB_mouse_LR_database0.txt",sep='\t',header=T,as.is = T)
# celltalkDB <- celltalkDB0[,c(1,3)]
# celltalkDB$pair <- paste(celltalkDB[,1],celltalkDB[,2],sep='_')
# identical(celltalkDB0$pair,pair)
# #[1] TRUE

len <- nrow(celltalkDB)
celltalkDB_ligand<-matrix(NA,len,10)
celltalkDB_receptor<-matrix(NA,len,10)

celltalkDB_class<-cbind(celltalkDB,celltalkDB_ligand)
celltalkDB_class<-cbind(celltalkDB_class,celltalkDB_receptor)

########ligand
col_index <- 6
for(i in 1:length(go_term_genes)){
  celltalkDB_class[which(celltalkDB_class[,"ligand"] %in% go_term_genes[[i]]),col_index]<-go_term[i]
  col_index <- col_index+1
}
########receptor
col_index <- 16
for(i in 1:length(go_term_genes)){
  celltalkDB_class[which(celltalkDB_class[,"receptor"] %in% go_term_genes[[i]]),col_index]<-go_term[i]
  col_index <- col_index+1
}
colnames(celltalkDB_class)[-c(1:5)] <- c(paste("ligand",1:10,sep=''),paste("receptor",1:10,sep=''))
write.table(celltalkDB_class,"celltalkDB_class.txt",sep='\t',quote = F,row.names = F)
dim(celltalkDB_class)

#####优先级：Notch signaling，Antigen binding，Neuropeptide，Hormone，cytokine（growth factor，interferon，interleukin，tumor necrosis factor，Chemokine）
len2 <- nrow(celltalkDB_class)
final_class<-celltalkDB_class[,1:5]
class<-matrix(NA,len2,1)
final_class<-cbind(final_class,class)


final_class[which(celltalkDB_class[,"ligand9"]=="Notch signaling"&celltalkDB_class[,"receptor9"]=="Notch signaling"),6]<-"Notch signaling"
final_class[which(is.na(final_class$class)),6]<-"Other"

final_class[which(final_class$class=="Other"&celltalkDB_class[,"ligand1"]=="Antigen binding"&celltalkDB_class[,"receptor1"]=="Antigen binding"),6]<-"Antigen binding"

final_class[which(final_class$class=="Other"&celltalkDB_class[,"ligand8"]=="Neuropeptide"&celltalkDB_class[,"receptor8"]=="Neuropeptide"),6]<-"Neuropeptide"

final_class[which(final_class$class=="Other"&celltalkDB_class[,"ligand5"]=="Hormone"&celltalkDB_class[,"receptor5"]=="Hormone"),6]<-"Hormone"
##########筛选cytokine############
final_class[which(final_class$class=="Other"&celltalkDB_class[,"ligand2"]=="Chemokine"&celltalkDB_class[,"receptor2"]=="Chemokine"&is.na(celltalkDB_class[,"ligand4"])&
                    is.na(celltalkDB_class[,"receptor4"])&is.na(celltalkDB_class[,"ligand6"])&is.na(celltalkDB_class[,"receptor6"])&is.na(celltalkDB_class[,"ligand7"])&
                    is.na(celltalkDB_class[,"receptor7"])&is.na(celltalkDB_class[,"ligand10"])&is.na(celltalkDB_class[,"receptor10"])),6]<-"Chemokine"

final_class[which(final_class$class=="Other"&celltalkDB_class[,"ligand4"]=="Growth factor"&celltalkDB_class[,"receptor4"]=="Growth factor"&is.na(celltalkDB_class[,"ligand2"])&is.na(celltalkDB_class[,"receptor2"])&
                    is.na(celltalkDB_class[,"ligand6"])&is.na(celltalkDB_class[,"receptor6"])&is.na(celltalkDB_class[,"ligand7"])&
                    is.na(celltalkDB_class[,"receptor7"])&is.na(celltalkDB_class[,"ligand10"])&is.na(celltalkDB_class[,"receptor10"])),6]<-"Growth factor"

final_class[which(final_class$class=="Other"&celltalkDB_class[,"ligand6"]=="Interferon"&celltalkDB_class[,"receptor6"]=="Interferon"&is.na(celltalkDB_class[,"ligand4"])&is.na(celltalkDB_class[,"receptor4"])&
                    is.na(celltalkDB_class[,"ligand2"])&is.na(celltalkDB_class[,"receptor2"])&is.na(celltalkDB_class[,"ligand7"])&
                    is.na(celltalkDB_class[,"receptor7"])&is.na(celltalkDB_class[,"ligand10"])&is.na(celltalkDB_class[,"receptor10"])),6]<-"Interferon"

final_class[which(final_class$class=="Other"& celltalkDB_class[,"ligand7"]=="Interleukin"&celltalkDB_class[,"receptor7"]=="Interleukin"&is.na(celltalkDB_class[,"ligand2"])&is.na(celltalkDB_class[,"receptor2"])&
                    is.na(celltalkDB_class[,"ligand4"])&is.na(celltalkDB_class[,"receptor4"])&is.na(celltalkDB_class[,"ligand6"])&
                    is.na(celltalkDB_class[,"receptor6"])&is.na(celltalkDB_class[,"ligand10"])&is.na(celltalkDB_class[,"receptor10"])),6]<-"Interleukin"

final_class[which(final_class$class=="Other"&celltalkDB_class[,"ligand10"]=="Tumor necrosis factor"&celltalkDB_class[,"receptor10"]=="Tumor necrosis factor"&is.na(celltalkDB_class[,"ligand2"])&is.na(celltalkDB_class[,"receptor2"])&
                    is.na(celltalkDB_class[,"ligand4"])&is.na(celltalkDB_class[,"receptor4"])&is.na(celltalkDB_class[,"ligand6"])&
                    is.na(celltalkDB_class[,"receptor6"])&is.na(celltalkDB_class[,"ligand7"])&is.na(celltalkDB_class[,"receptor7"])),6]<-"Tumor necrosis factor"

final_class[which(final_class$class=="Other"&celltalkDB_class[,"ligand3"]=="Cytokine"&celltalkDB_class[,"receptor3"]=="Cytokine"),6]<- "Cytokine"

final_class[which(substr(final_class$ligand,1,3)=="Tnf"&substr(final_class$receptor,1,3)=="Tnf"),6]<-"Tumor necrosis factor"
colnames(final_class)[6] <- "classification"
write.table(final_class,"classification/celltalkDB_mouse_LR_database.txt",sep='\t',quote = F,row.names = F)
####union&inter####
celltalkDB_classs <- read.table("classification/celltalkDB_mouse_LR_database.txt",sep='\t',header = T,as.is = T)
RNAMagnetDB_class <- read.table("classification/RNAMagnetDB_LR_database.txt",sep='\t',header = T,as.is = T)

union_class <- unique(rbind(celltalkDB_classs,RNAMagnetDB_class))
dim(union_class)
write.table(union_class,"classification/unionDB_mouse_LR_database.txt",sep='\t',quote = F,row.names = F)

inter_class <- na.omit(RNAMagnetDB_class[match(celltalkDB_classs$pair,RNAMagnetDB_class$pair),])##1157
dim(unique(inter_class))
write.table(inter_class,"classification/interDB_mouse_LR_database.txt",sep='\t',quote = F,row.names = F)
