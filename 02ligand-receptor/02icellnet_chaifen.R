setwd("I:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data")

data <- read.table("icellnet_LR_database4.txt",sep='\t',header = T,as.is = T) 
data <- as.matrix(data)

head(data)

pair_all <- c()

for(i in 1:nrow(data)){
  pair <- data[i,c("Ligand.1","Receptor.1","Family"),drop=F]
  pair_all <- rbind(pair_all,unname(pair))
  
  if(!is.na(data[i,"Ligand.2"])){
    
    pair2 <- data[i,c("Ligand.2","Receptor.1","Family"),drop=F]
    pair_all <- rbind(pair_all,unname(pair2))
    
    if(!is.na(data[i,"Receptor.2"])){
      pair3 <- data[i,c("Ligand.2","Receptor.2","Family"),drop=F]
      pair4 <- data[i,c("Ligand.1","Receptor.2","Family"),drop=F]
      pair_all <- rbind(pair_all,unname(pair3),unname(pair4))
    }
    
    if(!is.na(data[i,"Receptor.3"])){
      pair5 <- data[i,c("Ligand.1","Receptor.3","Family"),drop=F]
      pair6 <- data[i,c("Ligand.2","Receptor.3","Family"),drop=F]
      pair_all <- rbind(pair_all,unname(pair5),unname(pair6))
    }
    
  }else{
    
    if(!is.na(data[i,"Receptor.2"])){
      pair2 <- data[i,c("Ligand.1","Receptor.2","Family"),drop=F]
      pair_all <- rbind(pair_all,unname(pair2))
    }
    if(!is.na(data[i,"Receptor.3"])){
      pairs3 <- data[i,c("Ligand.1","Receptor.3","Family"),drop=F]
      pair_all <- rbind(pair_all,unname(pairs3))
    }
  }

}
pair <- paste(pair_all[,1],pair_all[,2],sep='_')
pair_all2 <- cbind(pair_all[,c(1,2)],pair,pair_all[,3,drop=F])
colnames(pair_all2) <- c("ligand","receptor","pair","classification")
dim(pair_all2)

pair_all2_final <- unique(pair_all2)
dim(pair_all2_final)

table(pair_all2_final[,4])
write.table(pair_all2_final,"I:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data/icellnet_LR_database.txt",sep='\t',quote = F,row.names = F)
