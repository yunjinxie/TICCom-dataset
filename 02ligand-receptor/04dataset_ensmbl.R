setwd("H:/F/tumor_immune_interaction/computation/A_webserver/database/data/data2")

file_name <- c("cellchatDB_LR_database","celltalkDB_LR_database","nichenet_LR_database",
               "ramilowski_LR_database","singlecellsignalR_LR_database",
               "icellnet_LR_database","iTALK_LR_database","inter_LR_database",
               "union_LR_database")
file_name2 <- paste(file_name,".txt",sep='')

dataset_dir <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data/classification_results"
gene_listPath <- file.path(dataset_dir,file_name2)

referencePath <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/data/data2"

for(i in 1:length(gene_listPath)){
  gene_list <- read.table(gene_listPath[i],sep='\t',header = T,as.is = T,fill = T,strip.white = T)
  gene_list2 <- as.matrix(gene_list)
  ############
  gene1_trans <- ID_transfer(gene_list = gene_list2[,1],ID_transfer_Path = referencePath)
  if(nrow(gene1)==0){
    cat("#######1:")
    cat(i)
    cat("######\n")
  }

  gene2 <- gene1_trans[match(unique(gene1_trans[,1]),gene1_trans[,1]),]
  gene3 <- gene2[match(unique(gene2[,2]),gene2[,2]),,drop=F]
  
  cat(max(table(gene3[,1])))
  cat(max(table(gene3[,2])))
  
  colnames(gene3) <- c("ligand","Gene1Id")
  gene_list3 <- merge(gene_list2,gene3,by='ligand')
  
  gene2_trans <- ID_transfer(gene_list = gene_list2[,2],ID_transfer_Path = referencePath)
  if(nrow(gene1)==0){
    cat("#######2:")
    cat(i)
    cat("######\n")
  }
  
  gene2_2 <- gene2_trans[match(unique(gene2_trans[,1]),gene2_trans[,1]),]
  gene3_2 <- gene2_2[match(unique(gene2_2[,2]),gene2_2[,2]),,drop=F]
  
  cat(max(table(gene3_2[,1])))
  cat(max(table(gene3_2[,2])))
  
  colnames(gene3_2) <- c("receptor","Gene2Id")
  gene_list4 <- merge(gene_list3,gene3_2,by='receptor')
  
  pairs2 <- gene_list4[,c("ligand","Gene1Id","receptor","Gene2Id","pair","classification"),drop=F]
  
  dim(pairs2)
  dim(unique(pairs2))
  
  t <- table(pairs2$pair)
  cat(max(t))
  
  pairs2[pairs2$pair %in% names(t[t==2]),]
  cat(file_name2[i])
  write.table(pairs2,file_name2[i],sep='\t',quote = F,row.names = F)
}
##################################
ID_transfer <- function(gene_list,ID_transfer_Path=NULL){
  
  ID_transfer_data <- read.table(paste(ID_transfer_Path,"ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
  
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
