#=================RNAMagnet LR database拆分==========================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/database/mouse")

RNAMagnet0 <- read.table("RNAMagnet_mouse_LR_database0.txt",sep='\t',header=T,as.is = T)
RNAMagnet1 <- RNAMagnet0[,c(2,3)]
dim(RNAMagnet0)
RNAMagnet[21,]

spe_index <- which(grepl("&|\\|",RNAMagnet1[,1])|grepl("&|\\|",RNAMagnet1[,2]))
spe_pairs <- RNAMagnet1[spe_index,]

pairs <- c()
for(i in 1:nrow(spe_pairs)){
  gene1 <- unlist(strsplit(spe_pairs[i,1],"&|\\|"))
  gene2 <- unlist(strsplit(spe_pairs[i,2],"&|\\|"))
  
  for(n in 1:length(gene1)){
    for(m in 1:length(gene2)){
      pairs <-rbind(pairs,cbind(gene1[n],gene2[m]))
    }
  }
}

RNAMagnet <- rbind(as.matrix(unname(RNAMagnet1[-spe_index,])),pairs)
RNAMagnet <- as.data.frame(RNAMagnet)
RNAMagnet$pair <- paste(RNAMagnet[,1],RNAMagnet[,2],sep='_')
colnames(RNAMagnet)[c(1,2)] <- c("ligand","receptor")
write.table(RNAMagnet,"RNAMagnet_split.txt",sep='\t',row.names = F,quote = F)
#=======================转化成ensmbl===============================
setwd("H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse")

file_name <- c("celltalkDB_mouse_LR_database0","RNAMagnet_split")
file_name2 <- paste(file_name,".txt",sep='')

dataset_dir <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse"
gene_listPath <- file.path(dataset_dir,file_name2)

referencePath <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/data/data2"
# ID_transfer_data <- read.table(paste(referencePath,"mouse_ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
# ID_transfer_data <- as.matrix(ID_transfer_data)
# 
# for(i in 1:length(gene_listPath)){
#   gene_list <- read.table(gene_listPath[i],sep='\t',header = T,as.is = T,fill = T,strip.white = T)
#   gene_list2 <- as.matrix(gene_list)
#   
#   diff <- setdiff(union(gene_list2[,1],gene_list2[,2]),ID_transfer_data[,3])
#   if(length(diff)!=0){
#     write.table(diff,paste("H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse/mouse_diff_",i,".txt",sep=""),sep='\t',quote = F,row.names = F,col.names = F)
#   }
# }

output_dir <- "ensmbl"
dir.create(output_dir)

for(i in 1:length(gene_listPath)){
  gene_list <- read.table(gene_listPath[i],sep='\t',header = T,as.is = T,fill = T,strip.white = T)
  gene_list2 <- as.matrix(gene_list)
  ############
  gene1_trans <- ID_transfer(gene_list = gene_list2[,1],ID_transfer_Path = referencePath)
  if(nrow(gene1_trans)==0){
    cat("#######1:")
    cat(i)
    cat("######\n")
  }
  
  gene2 <- gene1_trans[match(unique(gene1_trans[,1]),gene1_trans[,1]),]

  cat(max(table(gene2[,1])))
  cat(max(table(gene2[,2])))
  
  colnames(gene2) <- c("ligand","Gene1Id")
  gene_list3 <- merge(gene_list2,gene2,by='ligand')
  
  gene2_trans <- ID_transfer(gene_list = gene_list2[,2],ID_transfer_Path = referencePath)
  if(nrow(gene2_trans)==0){
    cat("#######2:")
    cat(i)
    cat("######\n")
  }
  
  gene2_2 <- gene2_trans[match(unique(gene2_trans[,1]),gene2_trans[,1]),]
  
  cat(max(table(gene2_2[,1])))
  cat(max(table(gene2_2[,2])))
  
  colnames(gene2_2) <- c("receptor","Gene2Id")
  gene_list4 <- merge(gene_list3,gene2_2,by='receptor')
  
  pairs2 <- gene_list4[,c("ligand","Gene1Id","receptor","Gene2Id","pair"),drop=F]
  
  dim(pairs2)
  dim(unique(pairs2))
  
  t <- table(pairs2$pair)
  cat(max(t))
  
  pairs2[pairs2$pair %in% names(t[t==2]),]
  cat(file_name2[i])
  write.table(pairs2,paste("ensmbl/",file_name2[i],sep=''),sep='\t',quote = F,row.names = F)
}
##################################
ID_transfer <- function(gene_list,ID_transfer_Path=NULL){
  
  ID_transfer_data <- read.table(paste(ID_transfer_Path,"mouse_ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
  
  gene <- gene_list[1]
  
  if(grepl("ENSG",gene[1])){
    gene_names = "ensembl"
  }else if(!grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])){
    gene_names = "entrez"
  }else{
    gene_names = "symbol"
  }
  
  if(gene_names == "ensembl"){
    res<- unique(ID_transfer_data[ID_transfer_data[,2] %in% gene_list,c(2,2)])
  }else if(gene_names == "entrez"){
    res <- unique(ID_transfer_data[ID_transfer_data[,1] %in% gene_list,c(1,2)])
  }else{
    res <- unique(ID_transfer_data[ID_transfer_data[,3] %in% gene_list,c(3,2)])
  }
  return(res)
}
