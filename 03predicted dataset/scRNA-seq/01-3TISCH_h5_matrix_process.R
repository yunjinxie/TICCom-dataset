#############
library(Seurat)
library(hdf5r)
setwd("I:/F/tumor_immune_interaction/computation/A_webserver")
data <- Read10X_h5("ALL_GSE132509/ALL_GSE132509_expression.h5")

data[1:4,1:4]

data@i[1:10]
str(data)

mm <- as_matrix(data)
dim(mm)
write.table(mm,"ALL_GSE132509/ALL_GSE132509_Expression_Data.txt",sep='\t',quote = F)

meta <- read.table("ALL_GSE132509/ALL_GSE132509_CellMetainfo_table.tsv",sep='\t',header=T,as.is = T,fill=T,quote = )
head(meta)

Cell_Lables_File <- meta[,c("Cell","Celltype..major.lineage.")]
colnames(Cell_Lables_File) <- c("cell","cell_type")
write.table(Cell_Lables_File,"ALL_GSE132509/ALL_GSE132509_Cell_Lables_File.txt",sep='\t',quote = F,row.names=F)

Comparison_Sample_Replicate_File <- meta[,c("Cell","Source","Sample")]
colnames(Comparison_Sample_Replicate_File) <- c("cell","compare_group","replicate")
write.table(Comparison_Sample_Replicate_File,"ALL_GSE132509/ALL_GSE132509_Comparison_Sample_Replicate_File.txt",sep='\t',quote = F,row.names=F)

#===============================================
data <- Read10X_h5("AML_GSE116256/AML_GSE116256_expression.h5")

data[1:4,1:4]

data@i[1:10]
str(data)

mm <- as_matrix(data)
dim(mm)
#[1] 19288 38348
write.table(mm,"AML_GSE116256/ALL_GSE132509_Expression_Data.txt",sep='\t',quote = F)

meta <- read.table("AML_GSE116256/AML_GSE116256_CellMetainfo_table.tsv",sep='\t',header=T,as.is = T,fill=T,quote = )
head(meta)

Cell_Lables_File <- meta[,c("Cell","Celltype..major.lineage.")]
colnames(Cell_Lables_File) <- c("cell","cell_type")
write.table(Cell_Lables_File,"AML_GSE116256/ALL_GSE132509_Cell_Lables_File.txt",sep='\t',quote = F,row.names=F)

Comparison_Sample_Replicate_File <- meta[,c("Cell","Source","Sample")]
colnames(Comparison_Sample_Replicate_File) <- c("cell","compare_group","replicate")
write.table(Comparison_Sample_Replicate_File,"AML_GSE116256/ALL_GSE132509_Comparison_Sample_Replicate_File.txt",sep='\t',quote = F,row.names=F)

#===============================================
data <- Read10X_h5("CRC_GSE139555/CRC_GSE139555_expression.h5")

data[1:4,1:4]

data@i[1:10]
str(data)

mm <- as_matrix(data)
dim(mm)
#[1] 19288 38348
write.table(mm,"CRC_GSE139555/CRC_GSE139555_Expression_Data.txt",sep='\t',quote = F)

meta <- read.table("CRC_GSE139555/CRC_GSE139555_CellMetainfo_table.tsv",sep='\t',header=T,as.is = T,fill=T,quote = )
head(meta)

Cell_Lables_File <- meta[,c("Cell","Celltype..major.lineage.")]
colnames(Cell_Lables_File) <- c("cell","cell_type")
write.table(Cell_Lables_File,"CRC_GSE139555/CRC_GSE139555_Cell_Lables_File.txt",sep='\t',quote = F,row.names=F)

Comparison_Sample_Replicate_File <- meta[,c("Cell","Source","Sample")]
colnames(Comparison_Sample_Replicate_File) <- c("cell","compare_group","replicate")
write.table(Comparison_Sample_Replicate_File,"CRC_GSE139555/CRC_GSE139555_Comparison_Sample_Replicate_File.txt",sep='\t',quote = F,row.names=F)

#========================================
file_name <- c("CRC_GSE139555","CRC_GSE146771_10X","CRC_GSE146771_Smartseq2","KIRC_GSE145281_aPDL1")

for(i in 1:length(file_name)){
  data <- Read10X_h5(paste(file_name[i],"/",file_name[i],"_expression.h5",sep=''))
  
  data[1:4,1:4]
  
  data@i[1:10]
  str(data)
  
  mm <- as_matrix(data)
  dim(mm)
  #[1] 19288 38348
  write.table(mm,paste(file_name[i],"/",file_name[i],"_Expression_Data.txt",sep=""),sep='\t',quote = F)
  
  meta <- read.table(paste(file_name[i],"/",file_name[i],"_CellMetainfo_table.tsv",sep=""),sep='\t',header=T,as.is = T,fill=T,quote = )
  head(meta)
  
  Cell_Lables_File <- meta[,c("Cell","Celltype..major.lineage.")]
  colnames(Cell_Lables_File) <- c("cell","cell_type")
  write.table(Cell_Lables_File,paste(file_name[i],"/",file_name[i],"_Cell_Lables_File.txt",sep=""),sep='\t',quote = F,row.names=F)
  
  Comparison_Sample_Replicate_File <- meta[,c("Cell","Source","Sample")]
  colnames(Comparison_Sample_Replicate_File) <- c("cell","compare_group","replicate")
  write.table(Comparison_Sample_Replicate_File,paste(file_name[i],"/",file_name[i],"_Comparison_Sample_Replicate_File.txt",sep=""),sep='\t',quote = F,row.names=F)
  
  cat(i)
  cat("\n##############################\n")
}
#========================================================
file_name2 <- c("HNSC_GSE103322","PAAD_CRA001160")

for(i in 1:length(file_name2)){
  data <- Read10X_h5(paste(file_name2[i],"/",file_name2[i],"_expression.h5",sep=''))
  
  data[1:4,1:4]
  
  data@i[1:10]
  str(data)
  
  mm <- as_matrix(data)
  dim(mm)
  #[1] 19288 38348
  write.table(mm,paste(file_name2[i],"/",file_name2[i],"_Expression_Data.txt",sep=""),sep='\t',quote = F)
  
  meta <- read.table(paste(file_name2[i],"/",file_name2[i],"_CellMetainfo_table.tsv",sep=""),sep='\t',header=T,as.is = T,fill=T,quote = )
  head(meta)
  
  Cell_Lables_File <- meta[,c("Cell","Celltype..major.lineage.")]
  colnames(Cell_Lables_File) <- c("cell","cell_type")
  write.table(Cell_Lables_File,paste(file_name2[i],"/",file_name2[i],"_Cell_Lables_File.txt",sep=""),sep='\t',quote = F,row.names=F)
  
  Comparison_Sample_Replicate_File <- meta[,c("Cell","Source","Patient")]
  colnames(Comparison_Sample_Replicate_File) <- c("cell","compare_group","replicate")
  write.table(Comparison_Sample_Replicate_File,paste(file_name2[i],"/",file_name2[i],"_Comparison_Sample_Replicate_File.txt",sep=""),sep='\t',quote = F,row.names=F)
  
  cat(i)
  cat("\n##############################\n")
}
#################稀疏矩阵转化为普通矩阵
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
