##R调用python
library(reticulate)
py_config()#安装的python版本环境查看，显示anaconda和numpy的详细信息。
py_available()#[1] TRUE   #检查您的系统是否安装过Python
py_module_available("pandas")#检查“pandas”是否安装
py_module_available("umap")

library(yulab.utils)
install_zip("I:/F/tumor_immune_interaction/computation/umapr-master.zip") 

library(yulab.utils)
install_zip("I:/F/tumor_immune_interaction/computation/mindr-1.2.zip") 

hsa<-CellChatDB.human
head(hsa)
interact <- CellChatDB.human$interaction
head(interact)
com<-CellChatDB.human$complex
head(com)
cofactor<-CellChatDB.human$cofactor
head(cofactor)
geneinfo<-CellChatDB.human$geneInfo
head(geneinfo)


write.table(interact,"I:/F/tumor_immune_interaction/computation/A_webserver/database/data/cellchatDB.txt",sep='\t',quote = F,row.names = F)
#====================================处理cellchat中多亚基结构======================
interact <- read.table("H:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data/cellchatDB.txt",sep='\t',header=T,as.is = T)
####ligand_receptor与interact（第一列）一致####
pair <- paste(interact[,3],interact[,4],sep='_')
ident_interact <- interact[interact[,1]==pair,c(1,3,4,2,10)]
####受体是多个，配体是一个####
ligand_index <- which(grepl('_',ident_interact[,2]))
ligand_1 <- ident_interact[-ligand_index,]
receptor_ident <- lapply(ligand_1[,3],function(x){unlist(strsplit(x,"_"))})

res_1 <- c()
for(i in 1:nrow(ligand_1)){
  res0 <- cbind(ligand_1[i,2],receptor_ident[[i]])
  res01 <- cbind(res0,ligand_1[i,4],ligand_1[i,5])
  res_1 <- rbind(res_1,res01)
}
res_1 <- unique(res_1)
#7420
####配体是多个,受体是一个####
ligand_2 <- ident_interact[ligand_index,]
ligand_ident <- lapply(ligand_2[,2],function(x){unlist(strsplit(x,"_"))})
res_2 <- c()
for(i in 1:nrow(ligand_2)){
  res0 <- cbind(ligand_ident[[i]],ligand_2[i,3])
  res01 <- cbind(res0,ligand_1[i,4],ligand_1[i,5])
  res_2 <- rbind(res_2,res01)
}
res_2 <- unique(res_2)
res <- rbind(res_1,res_2)
#7436
####ligand_receptor与interact（第一列）不一致####
diff_interact <- interact[interact[,1]!=pair,c(1,3,4,2,10)]
#99
d1<-diff_interact[grep("TGFbR",diff_interact[,3]),]
#12
receptor <- t(apply(d1,1,function(x){unlist(strsplit(x[1],"_"))[c(2,3)]}))
res13 <- c()
for(i in 1:nrow(d1)){
    res11 <-rbind(c(d1[i,2],receptor[i,1]),
                   c(d1[i,2],receptor[i,2]))
    res12 <- cbind(res11,d1[i,4],d1[i,5])
    res13 <- rbind(res13,res12)
}
##########
d2<-diff_interact[!grepl("TGFbR",diff_interact[,3])&grepl("complex|receptor",diff_interact[,3]),]
#3
receptor_d2 <- apply(d2,1,function(x){unlist(strsplit(x[1],"_"))[-1]})
res2 <- c()
for(i in 1:nrow(d2)){
  res21 <- cbind(d2[i,2],receptor_d2[[i]])
  res22 <- cbind(res21,d2[i,4],d2[i,5])
  res2 <- rbind(res22,res2)
}
#######
d3<-diff_interact[!grepl("TGFbR",diff_interact[,3])&!grepl("complex|receptor",diff_interact[,3]),]
#84
###
index <- apply(d3,1,function(x){grepl(paste("^",x[2],"_",sep=''),x[1])})
d3_ligand <- d3[index,]
ligand_i_d3_ligand <- which(grepl("_",d3_ligand[,2]))
#empty
receptor_d3_ligand <- lapply(d3_ligand[,3],function(x){unlist(strsplit(x,"_"))})
res3 <- c()
for(i in 1:nrow(d3_ligand)){
  res31 <- cbind(d3_ligand[i,2],receptor_d3_ligand[[i]])
  res32 <- cbind(res31,d3_ligand[i,4],d3_ligand[i,5])
  res3 <- rbind(res3,res32)
}
res3 <- unique(res3)
#216
d3_ligand_2 <- d3[!index,]
d3_ligand_2[d3_ligand_2[,2]=="Activin AB",2] <- "INHBA"
d3_ligand_2[d3_ligand_2[,2]=="Inhibin A",2] <- "INHA"
d3_ligand_2[d3_ligand_2[,2]=="Inhibin B",2] <- "INHBB"
d3_ligand_2[d3_ligand_2[,2]=="IL12AB",2] <- "IL12B"
d3_ligand_2[d3_ligand_2[,2]=="LTa1b2",2] <- "LTA_LTB"
d3_ligand_3 <- d3_ligand_2[!grepl("complex",d3_ligand_2[,2]),]

ligand_2_index <- which(grepl("_",d3_ligand_3[,2]))
d3_ligand_31 <- d3_ligand_3[-ligand_2_index,]
receptor_d3_ligand31 <- lapply(d3_ligand_31[,3],function(x){unlist(strsplit(x,"_"))})
res4_1 <- c()
for(i in 1:nrow(d3_ligand_31)){
  res41 <- cbind(d3_ligand_31[i,2],receptor_d3_ligand31[[i]])
  res42 <- cbind(res41,d3_ligand_31[i,4],d3_ligand_31[i,5])
  res4_1 <- rbind(res4_1,res42)
}
res4_1 <- unique(res4_1)

d3_ligand_32 <- d3_ligand_3[ligand_2_index,]
ligand_d3_ligand32 <- strsplit(as.character(d3_ligand_32[2]),"_")
res4_2 <- cbind(ligand_d3_ligand32[[1]],d3_ligand_32[,3],
                 d3_ligand_32[,4],d3_ligand_32[,5])

d3_ligand_4 <- d3_ligand_2[grepl("complex",d3_ligand_2[,2]),]
d3_ligand_4_ligand <- lapply(d3_ligand_4[,1],function(x){unlist(strsplit(x,"_"))[1]})
d3_ligand_4_receptor <- lapply(d3_ligand_4[,1],function(x){unlist(strsplit(x,"_"))[-1]})

res4_3 <- c()
for(i in 1:nrow(d3_ligand_4)){
  res4_31 <- cbind(d3_ligand_4_ligand[[i]],d3_ligand_4_receptor[[i]])
  res4_32 <- cbind(res4_31,d3_ligand_4[i,4],d3_ligand_4[i,5])
  res4_3 <- rbind(res4_3,res4_32)
}
res4 <- rbind(res4_1,res4_2,res4_3)
#######
res_final <- unique(rbind(res,res13,res2,res3,res4))
colnames(res_final)<- c("ligand","receptor","pathway_name","annotation")
write.table(res_final,"H:/F/tumor_immune_interaction/computation/A_webserver/database/data/original data/cellchatDB_2.txt",sep='\t',quote=F,row.names=F)
######
