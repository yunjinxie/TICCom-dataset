#############配体-受体数据集合并
setwd("H:/F/tumor_immune_interaction/computation/A_webserver/database/data/data2")

human_union <- read.table("union_LR_database.txt",sep='\t',header=T,as.is = T)
human_union <- as.matrix(human_union)
dim(human_union)
#14190
mouse_union <- read.table("unionDB_mouse_LR_database.txt",sep='\t',header=T,as.is = T)
mouse_union <- as.matrix(mouse_union)

dim(mouse_union)#3650
dim(unique(mouse_union))

union <- unique(rbind(human_union,mouse_union))
dim(union)
union_ID_pair <- paste(union[,2],union[,4],sep='_')
#intersect(human_union[,"pair"],mouse_union[,"pair"])

source <- matrix("-",nrow = nrow(union),ncol=13)
union_source <- cbind(union,source)
colnames(union_source)[7:19] <- c("cellchatDB_human","celltalkDB_human","icellnetDB_human","iTALKDB_human","nichenetDB_human",
                                  "ramilowskiDB_human","singlecellsignalRDB_human","unionDB_human","intersectDB_human","celltalkDB_mouse","RNAMagnetDB_mouse","unionDB_mouse","intersectDB_mouse")

human_file_names <- c("cellchatDB_LR_database.txt","celltalkDB_LR_database.txt",
                "icellnet_LR_database.txt","iTALK_LR_database.txt","nichenet_LR_database.txt",
                "ramilowski_LR_database.txt","singlecellsignalR_LR_database.txt")

mouse_file_names <- c("celltalkDB_mouse_LR_database.txt","RNAMagnetDB_LR_database.txt")

for(i in 1:length(human_file_names)){
  data <- read.table(human_file_names[i],sep='\t',header=T,as.is = T)
  data <- as.matrix(data)
  
  ID_pair <- paste(data[,2],data[,4],sep='_')
  
  union_source[union_ID_pair %in% ID_pair,(i+6)] <- "Yes"
  cat(i)
}

union_source[union_source[,7]=="Yes"|union_source[,8]=="Yes"|union_source[,9]=="Yes"|union_source[,10]=="Yes"|union_source[,11]=="Yes"|
                              union_source[,12]=="Yes"|union_source[,13]=="Yes","unionDB_human"] <- "Yes"

union_source[union_source[,7]=="Yes"&union_source[,8]=="Yes"&union_source[,9]=="Yes"&union_source[,10]=="Yes"&union_source[,11]=="Yes"&
                                  union_source[,12]=="Yes"&union_source[,13]=="Yes","intersectDB_human"] <- "Yes"

for(i in 1:length(mouse_file_names)){
  data <- read.table(mouse_file_names[i],sep='\t',header=T,as.is = T)
  data <- as.matrix(data)
  
  ID_pair_m <- paste(data[,2],data[,4],sep='_')
  
  union_source[union_ID_pair %in% ID_pair_m,(i+15)] <- "Yes"
  cat(i)
}

union_source[union_source[,16]=="Yes"|union_source[,17]=="Yes","unionDB_mouse"] <- "Yes"

union_source[union_source[,16]=="Yes"&union_source[,17]=="Yes","intersectDB_mouse"] <- "Yes"

# aa <- union_source[union_source[,"unionDB_human"]=="Yes",]
# bb <- union_source[union_source[,"intersectDB_human"]=="Yes",]
# inter_LR <- bb[,c(1,2,3,4,5,6)]
# dim(inter_LR)

# write.table(inter_LR,"H:/F/tumor_immune_interaction/computation/A_webserver/database/data/data2/inter_LR_database.txt",sep='\t',quote = F,row.names = F)

# aa2 <- union_source[union_source[,"unionDB_mouse"]=="Yes",]
# 
# setdiff(aa2[,"pair"],mouse_union[,"pair"])
# 
# bb2 <- union_source[union_source[,"intersectDB_mouse"]=="Yes",]
# dim(bb2)
union_source_final <- subset(union_source,select = -pair)

write.table(union_source_final,"H:/F/tumor_immune_interaction/computation/A_webserver/download/TICCom_data/search2/union_source.txt",sep='\t',quote = F,row.names = F)

