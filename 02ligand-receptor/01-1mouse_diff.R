setwd("H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse")

file_name <- c("celltalkDB_mouse_LR_database0","RNAMagnet_split")
file_name2 <- paste(file_name,".txt",sep='')

dataset_dir <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse"
gene_listPath <- file.path(dataset_dir,file_name2)

referencePath <- "H:/F/tumor_immune_interaction/computation/A_webserver/database/data/data2"
ID_transfer_data <- read.table(paste(referencePath,"mouse_ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
ID_transfer_data <- as.matrix(ID_transfer_data)

for(i in 1:length(gene_listPath)){
  gene_list <- read.table(gene_listPath[i],sep='\t',header = T,as.is = T,fill = T,strip.white = T)
  gene_list2 <- as.matrix(gene_list)

  diff <- setdiff(union(gene_list2[,1],gene_list2[,2]),ID_transfer_data[,3])
  if(length(diff)!=0){
    write.table(diff,paste("H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse/mouse_diff_",i,".txt",sep=""),sep='\t',quote = F,row.names = F,col.names = F)
  }
}

RM <- read.table("H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse/RNAMagnet_split.txt",sep='\t',header=T,as.is = T)
index1 <- which(RM[,1]=="Cd244")
length(index1)
RM[index1,1] <- gsub("Cd244","Cd244a",RM[index1,1])

index2 <- which(RM[,2]=="Cd244")
length(index2)
RM[index2,2] <- gsub("Cd244","Cd244a",RM[index2,2])

index3 <- grep("_Cd244$|^Cd244_",RM[,3])
length(index3)
RM[index3,3] <- gsub("Cd244","Cd244a",RM[index3,3])

RM <- unique(RM)
write.table(RM,"H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse/RNAMagnet_split.txt",sep='\t',quote = F,row.names = F)

#==================================
cell <- read.table("H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse/celltalkDB_mouse_LR_database0.txt",sep='\t',header=T,as.is = T)
index1 <- which(cell[,1]=="Il1f5")
length(index1)
cell[index1,1] <- gsub("Il1f5","Il36rn",cell[index1,1])

index2 <- which(cell[,2]=="Il1f5")
length(index2)
cell[index2,2] <- gsub("Il1f5","Il36rn",cell[index2,2])

index3 <- grep("_Il1f5$|^Il1f5_",cell[,3])
length(index3)
cell[index3,3] <- gsub("Il1f5","Il36rn",cell[index3,3])

cell <- unique(cell)
write.table(cell,"H:/F/tumor_immune_interaction/computation/A_webserver/database/mouse/celltalkDB_mouse_LR_database0.txt",sep='\t',quote = F,row.names = F)
