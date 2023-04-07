#===============================配体-受体数据集标记(manually curated/predicted)===============================
rm(list=ls())
gc()

setwd("I:/F/tumor_immune_interaction/computation/A_webserver/database/data/data2")

####分是否literature supported或者curated 来源数据库！！！
####union####
data <- read.table("I:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search2/union_source.txt",
                   sep='\t',header = T,as.is = T)
data_pairs <- paste(data[,1],data[,3],sep="_")
####human####
####cellchatDB_human####
data$cellchatDB_human_curated <- "-"
index <- which(data$cellchatDB_human == "Yes")
data[index,]$cellchatDB_human_curated <- "manually curated"

####celltalkDB_human####
data$celltalkDB_human_curated <- "-"
index <- which(data$celltalkDB_human == "Yes")
data[index,]$celltalkDB_human_curated <- "manually curated"

####icellnetDB_human####
data$icellnetDB_human_curated <- "-"
index <- which(data$icellnetDB_human == "Yes")
data[index,]$icellnetDB_human_curated <- "manually curated"

####ramilowskiDB_human####
LR_orig <- readxl::read_excel("I:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data/ramilowski.xlsx",
                              sheet = "All.Pairs")
table(LR_orig$Pair.Evidence)
LR_orig_curated <- LR_orig[LR_orig$Pair.Evidence == "literature supported",]
index <- which(data$ramilowskiDB_human == "Yes")
inter <- intersect(LR_orig_curated$Pair.Name,data_pairs[index])

data$ramilowskiDB_human_curated <- "-"
data[data_pairs %in% inter & data$ramilowskiDB_human == "Yes",]$ramilowskiDB_human_curated <- "manually curated"

putative_pairs <- LR_orig[LR_orig$Pair.Evidence %in% c("EXCLUDED",
                                                      "EXCLUDED not ligand","EXCLUDED not receptor","putative"),]$Pair.Name
inter_putative <- intersect(data_pairs[index],putative_pairs)
data[data_pairs %in% inter_putative & data$ramilowskiDB_human == "Yes",]$ramilowskiDB_human_curated <- "predicted"
####singlecellsignalRDB_human####
data$singlecellsignalRDB_human_curated <- "-"
index <- which(data$singlecellsignalRDB_human == "Yes")

load("I:/F/tumor_immune_interaction/computation/SingleCellSignalR/LRdb.rda")
PMID <- LRdb[LRdb$PMIDs != "",]
PMID_pairs <- paste(PMID$ligand,PMID$receptor,sep='_')
inter <- intersect(data_pairs[index],PMID_pairs)
data[data_pairs %in% inter & data$singlecellsignalRDB_human == "Yes",]$singlecellsignalRDB_human_curated <- "manually curated"

####前4个取个并集
index <- which(data$cellchatDB_human_curated == "manually curated"|data$celltalkDB_human_curated == "manually curated"|
                 data$icellnetDB_human_curated == "manually curated"|data$ramilowskiDB_human_curated == "manually curated")
curated_pairs_4 <- paste(data[index,]$ligand,data[index,]$receptor,sep='_')
########
LRdb_other <- LRdb[LRdb$PMIDs == "",]
other_pairs <- paste(LRdb_other$ligand,LRdb_other$receptor,sep='_')
inter_mc <- intersect(other_pairs,curated_pairs_4)
data[data_pairs %in% inter_mc & data$singlecellsignalRDB_human == "Yes",]$singlecellsignalRDB_human_curated <- "manually curated"

other_diff <- setdiff(other_pairs,curated_pairs_4)
data[data_pairs %in% other_diff & data$singlecellsignalRDB_human == "Yes",]$singlecellsignalRDB_human_curated <- "predicted"

####nichenetDB_human####
LR_orig <- readRDS("I:/F/tumor_immune_interaction/computation/NicheNet_network/lr_network.rds")
table(LR_orig$database)
LR_orig_curated <- LR_orig[LR_orig$database %in% c("guide2pharmacology","kegg","ramilowski"),]
LR_orig_curated <- as.data.frame(LR_orig_curated)
LR_orig_curated_pairs <- paste(LR_orig_curated$from,LR_orig_curated$to,sep="_")

index <- which(data$nichenetDB_human == "Yes")
inter <- intersect(LR_orig_curated_pairs,data_pairs[index])

data$nichenetDB_human_curated <- "-"
data[data_pairs %in% inter & data$nichenetDB_human == "Yes",]$nichenetDB_human_curated <- "manually curated"

others <- LR_orig[LR_orig$database %in% c("ppi_prediction","ppi_prediction_go"),]
others <- as.data.frame(others)
others_pairs <- paste(others$from,others$to,sep="_")

#前5个curated并集
index <- which(data$cellchatDB_human_curated == "manually curated"|data$celltalkDB_human_curated == "manually curated"|
                 data$icellnetDB_human_curated == "manually curated"|data$ramilowskiDB_human_curated == "manually curated"|
                 data$singlecellsignalRDB_human_curated == "manually curated")
curated_pairs_5 <- paste(data[index,]$ligand,data[index,]$receptor,sep='_')
inter_others <- intersect(others_pairs,curated_pairs_5)
data[data_pairs %in% inter_others & data$nichenetDB_human == "Yes",]$nichenetDB_human_curated <- "manually curated"

diff_others <- setdiff(others_pairs,curated_pairs_5)
data[data_pairs %in% diff_others & data$nichenetDB_human == "Yes",]$nichenetDB_human_curated <- "predicted"

####iTALKDB_human####
data$iTALKDB_human_curated <- "-"
index <- which(data$iTALKDB_human == "Yes")
LR_l <- load("I:/F/tumor_immune_interaction/computation/iTALK-master/data/LR_database.rda")
LR <- eval(parse(text = LR_l))
pairs <- paste(LR$Ligand.ApprovedSymbol,LR$Receptor.ApprovedSymbol,sep='_')

#前6个并集
index <- which(data$cellchatDB_human_curated == "manually curated"|data$celltalkDB_human_curated == "manually curated"|
                 data$icellnetDB_human_curated == "manually curated"|data$ramilowskiDB_human_curated == "manually curated"|
                 data$singlecellsignalRDB_human_curated == "manually curated"|
                 data$nichenetDB_human_curated == "manually curated")
curated_pairs_6 <- paste(data[index,]$ligand,data[index,]$receptor,sep='_')

inter <- intersect(pairs,curated_pairs_6)
data[data_pairs %in% inter & data$iTALKDB_human == "Yes",]$iTALKDB_human_curated <- "manually curated"

pairs_diff <- setdiff(pairs,curated_pairs_6)
data[data_pairs %in% pairs_diff & data$iTALKDB_human == "Yes",]$iTALKDB_human_curated <- "predicted"

####house mouse####
####celltalkDB_mouse####
data$celltalkDB_mouse_curated <- "-"
index <- which(data$celltalkDB_mouse == "Yes")
data[index,]$celltalkDB_mouse_curated <- "manually curated"

####RNAMagnetDB_mouse####
# load("I:/F/tumor_immune_interaction/computation/RNA_Magnet/ligandsReceptors_3.0.0.rda")
# RNAMagnet <- ligandsReceptors_3.0.0
# 
# spe_index <- which(grepl("&|\\|",RNAMagnet[,2])|grepl("&|\\|",RNAMagnet[,3]))
# spe_pairs <- RNAMagnet[spe_index,]
# 
# pairs <- c()
# for(i in 1:nrow(spe_pairs)){
#   gene1 <- unlist(strsplit(spe_pairs[i,2],"&|\\|"))
#   gene2 <- unlist(strsplit(spe_pairs[i,3],"&|\\|"))
#   ManualAnnotation <- spe_pairs[i,"ManualAnnotation"]
#   Source <- spe_pairs[i,"Source"]
#   
#   for(n in 1:length(gene1)){
#     for(m in 1:length(gene2)){
#       pairs <-rbind(pairs,cbind(gene1[n],gene2[m],Source,ManualAnnotation))
#     }
#   }
# }
# 
# RNAMagnet2 <- rbind(as.matrix(unname(RNAMagnet[-spe_index,c(2,3,4,5)])),pairs)
# RNAMagnet2 <- as.data.frame(RNAMagnet2)
# RNAMagnet2$pair <- paste(RNAMagnet2[,1],RNAMagnet2[,2],sep='_')
# colnames(RNAMagnet2)[c(1,2)] <- c("ligand","receptor")
# write.table(RNAMagnet2,"I:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data/RNAMagnet_split.txt",sep='\t',row.names = F,quote = F)
# 
# RNAMagnet2 <- read.table("I:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data/RNAMagnet_split.txt",sep='\t',header = T,as.is = T)

RNAMagnet2 <- read.table("I:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data/RNAMagnet_split2.txt",sep='\t',header = T,as.is = T)
pairs <- paste(RNAMagnet2$ligand,RNAMagnet2$receptor,sep='_')
#celltalkDB_mouse
data$RNAMagnetDB_mouse_curated <- "-"
RNAMagnet2_curated <- RNAMagnet2[RNAMagnet2$Source %in% c("Baccin","cellphoneDB","Ramilowski"),]
curated_pairs <- paste(RNAMagnet2_curated$ligand,RNAMagnet2_curated$receptor,sep='_')
data[data_pairs %in% curated_pairs & data$RNAMagnetDB_mouse == "Yes",]$RNAMagnetDB_mouse_curated <- "manually curated"

other <- RNAMagnet2[RNAMagnet2$Source %in% "Julia Schnell",]
other_pairs <- paste(other$ligand,other$receptor,sep="_")

celltalkDB_mouse <- data[data$celltalkDB_mouse == "Yes",]
celltalkDB_mouse_pair <- paste(celltalkDB_mouse$ligand,celltalkDB_mouse$receptor,sep='_')
inter <- intersect(other_pairs,celltalkDB_mouse_pair)
data[data_pairs %in% inter & data$RNAMagnetDB_mouse == "Yes",]$RNAMagnetDB_mouse_curated <- "manually curated"

other_diff <- setdiff(other_pairs,celltalkDB_mouse_pair)
data[data_pairs %in% other_diff & data$RNAMagnetDB_mouse == "Yes",]$RNAMagnetDB_mouse_curated <- "predicted"

write.table(data,"curate_confirm/union_source_curated_respectively.txt",sep='\t',quote = F,row.names = F)
########
data <- read.table("curate_confirm/union_source_curated_respectively.txt",sep='\t',header = T,as.is = T)
library(dplyr)
data2 <- data %>% select(ligand:intersectDB_mouse)
data2$Evidence <- "-"

len <- apply(data,1,function(x){
  length(which(x=="manually curated"))
})

dd <- data[len==0,]
data2[len!=0,]$Evidence <- "manually curated"
data2[len==0,]$Evidence <- "predicted"
write.table(data2,"curate_confirm/union_source_curated.txt",sep='\t',quote = F,row.names = F)

####union####
data2$pair <- paste(data2$ligand,data2$receptor,sep='_')
union_LR <- data2[data2$unionDB_human == "Yes",] %>% select(ligand:Gene2Id,
                                                            pair,classification,Evidence)
save(union_LR,file = "curate_confirm/union_LR_database.RData")

inter_LR <- data2[data2$intersectDB_human == "Yes",] %>% select(ligand:Gene2Id,
                                                            pair,classification,Evidence)
save(inter_LR,file = "curate_confirm/inter_LR_database.RData")
####小鼠####
union_LR_m <- data2[data2$unionDB_mouse == "Yes",] %>% select(ligand:Gene2Id,
                                                            pair,classification,Evidence)
save(union_LR_m,file = "curate_confirm/unionDB_mouse_LR_database.RData")

inter_LR_m <- data2[data2$intersectDB_mouse == "Yes",] %>% select(ligand:Gene2Id,
                                                                pair,classification,Evidence)
save(inter_LR_m,file = "curate_confirm/interDB_mouse_LR_database.RData")

LR_datasets_h <- c("cellchatDB","celltalkDB","iTALKDB","nichenetDB",
                 "ramilowskiDB","singlecellsignalRDB","icellnetDB")
data$pair <- paste(data$ligand,data$receptor,sep='_')
for(i in LR_datasets_h){
  h <- paste0(i,"_human")
  hc <- paste0(i,"_human_curated")
  LR <- data[data[,h]=="Yes",] %>% select(ligand:Gene2Id,
                                   pair,classification,
                                   Evidence = all_of(hc))
  out <- paste0(i,"_LR_database")
  save(LR,file = paste0("curate_confirm/",out,".RData"))
}
LR_datasets_m <- c("celltalkDB","RNAMagnetDB")
for(i in LR_datasets_m){
  h <- paste0(i,"_mouse")
  hc <- paste0(i,"_mouse_curated")
  LR <- data[data[,h]=="Yes",] %>% select(ligand:Gene2Id,
                                          pair,classification,
                                          Evidence = all_of(hc))
  out <- paste0(i,"_mouse_LR_database")
  save(LR,file = paste0("curate_confirm/",out,".RData"))
}

#####ICELLNET####
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/data/data2"
LR <- read.table(paste0(d,"/icellnet_LR_database4.txt"),sep='\t',header = T,as.is = T,check.names = F)
colnames(LR)
LR$Evidence <- "manually curated"
setwd(d)
save(LR,file = paste0("curate_confirm/icellnet_LR_database4.RData"))
