library(dplyr)
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC"

file <- dir(d)
file2 <- file[grepl("exp_seq.*tsv$",file)]
file2

out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count"

info <- c()
for(i in 1:length(file2)){

    exp <- read.table(paste0(d,"/",file2[i]),sep='\t',header = T,as.is = T)
    colnames(exp)
    head(exp[,c("gene_id","icgc_sample_id","raw_read_count","gene_model")])

     sam <- unique(exp$icgc_sample_id)
     gene <- unique(exp$gene_id)
     count <- matrix(0,ncol=length(sam),nrow=length(gene))
     rownames(count) <- gene
     colnames(count) <- sam
    
     for(j in 1:length(sam)){
       cc <- exp[exp$icgc_sample_id == sam[j],c("gene_id","raw_read_count")]
       count[cc[,1],sam[j]] <- cc[,2]
     }

    saveRDS(count,paste0(out_d,"/",file2[i],"_count.rds"))
    count <- readRDS(paste0(out_d,"/",file2[i],"_count.rds"))
    
    info <- rbind(info,cbind(file2[i],unique(exp$assembly_version),unique(exp$gene_model),nrow(count),ncol(count)))
    
    cat(paste0(file2[i],"\n"))
 
}
write.table(info,paste0(out_d,"/count_info.txt"),sep='\t',quote = F,row.names = F)

####symbol2ENSG####
gene_gtf <- read.table("J:/xieyunjin/lncRNA_synergisty_project/gencode.v19.annotation.gtf/all_gene_symbol.txt",sep='\t',header = F,as.is = T)
colnames(gene_gtf) <- c("ENSG","symbol")

out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count"
info <- read.table(paste0(out_d,"/count_info.txt"),sep='\t',header = T,as.is = T)
file <- info[info[,3]!="Ensembl",1]

for(i in 1:length(file)){
  dataFile <- paste0(out_d,"/",file[i],"_count.rds")
  data <- readRDS(dataFile)
  data2 <- cbind(rownames(data),data)
  colnames(data2)[1] <- "symbol"
  
  data_merge <- merge(gene_gtf,data2,by="symbol")
  frenq <- table(data_merge[,"ENSG"])
  
  if(max(frenq)==1){
    data_merge2 <- data_merge[,-c(1,2)]
    rownames(data_merge2) <- as.matrix(data_merge[,2])
  }else{
    uni_data <- data_merge[data_merge[,2] %in% names(frenq[frenq==1]),]
    uni_data2 <- uni_data[,-c(1,2)]
    rownames(uni_data2) <- as.matrix(uni_data[,2])
    
    more_data <- data_merge[data_merge[,2] %in% names(frenq[frenq!=1]),,drop=F]
    more_data2 <- apply(more_data[,-c(1,2)],2,function(x){
      tapply(x,factor(more_data[,2]),function(x) mean(as.numeric(x)))
    })
    
    data_merge2 <- rbind(uni_data2,more_data2)
  }
  
  saveRDS(data_merge2,paste0(out_d,"/",file[i],"_count_ensmbl.rds"))
  cat(paste0(file[i],"\n"))
}

####exp_seq.PRAD-CA.tsv####
file <- "exp_seq.PRAD-CA.tsv_count.rds"
out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count"

dataFile <- paste0(out_d,"/",file)
data <- readRDS(dataFile)
data2 <- cbind(rownames(data),data)
colnames(data2)[1] <- "symbol"

data_merge <- merge(gene_gtf,data2,by="symbol")
frenq <- table(data_merge[,"ENSG"])

if(max(frenq)==1){
  data_merge2 <- data_merge[,-c(1,2)]
  rownames(data_merge2) <- as.matrix(data_merge[,2])
}

saveRDS(data_merge2,paste0(out_d,"/",file,"_ensmbl.rds"))


####ENST_2_ENSG####
####提取GRCh37中的转录本####
gtf <- readLines("J:/xieyunjin/lncRNA_synergisty_project/gencode.v19.annotation.gtf/gencode.v19.annotation1.gtf")
transcript <- gtf[grepl("\ttranscript?\t",gtf)]
transcript2 <- strsplit(transcript,"\t")
transcript2_mat <- do.call(rbind,transcript2)
transcript2_mat_9 <- strsplit(transcript2_mat[,9],";")
transcript2_mat_9_mat <- do.call(rbind,
                                 lapply(lapply(transcript2_mat_9,unlist),`length<- `,
                                        max(lengths(transcript2_mat_9))))

transcript2_mat_9_mat <- do.call(rbind,
                            lapply(lapply(transcript2_mat_9,unlist),`length<-`,
                                   max(lengths(transcript2_mat_9))))

transcript3 <- cbind(transcript2_mat[,1:8],transcript2_mat_9_mat)
write.table(transcript3,"J:/xieyunjin/lncRNA_synergisty_project/gencode.v19.annotation.gtf/transcript_gtf_mat.txt",sep='\t',
            quote = F,row.names = F,col.names = F)

transcript_2_gene <- transcript3[,c(9,10)]
transcript_2_gene[,1] <- stringr::str_match(transcript_2_gene[,1],'gene_id.*"(.*?)\\"')[,2]
transcript_2_gene[,2] <- stringr::str_match(transcript_2_gene[,2],'transcript_id.*"(.*?)\\"')[,2]
write.table(transcript_2_gene,"J:/xieyunjin/lncRNA_synergisty_project/gencode.v19.annotation.gtf/transcript_2_gene.txt",
            sep='\t',quote = F,row.names = F,col.names = F)

####转化####
trans_gene <- read.table("J:/xieyunjin/lncRNA_synergisty_project/gencode.v19.annotation.gtf/transcript_2_gene.txt",
                          sep='\t',header = F,as.is = T)
colnames(trans_gene) <- c("gene","transcript")
trans_gene[,1] <- gsub("\\..*","",trans_gene[,1])
trans_gene[,2] <- gsub("\\..*","",trans_gene[,2])

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count"
file <- "exp_seq.PBCA-US.tsv_count.rds"

data <- readRDS(paste(d,file,sep='/'))
data2 <- cbind(rownames(data),data)
colnames(data2)[1] <- "transcript"

data_merge <- merge(trans_gene,data2,by= "transcript") 
frenq <- table(data_merge$gene)

if(max(frenq)==1){
  data_merge2 <- data_merge[,-c(1,2)]
  rownames(data_merge2) <- as.matrix(data_merge[,2])
}else{
  uni_data <- data_merge[data_merge[,2] %in% names(frenq[frenq==1]),]
  uni_data2 <- uni_data[,-c(1,2)]
  rownames(uni_data2) <- as.matrix(uni_data[,2])
  
  more_data <- data_merge[data_merge[,2] %in% names(frenq[frenq!=1]),,drop=F]
  more_data2 <- apply(more_data[,-c(1,2)],2,function(x){
    tapply(x,factor(more_data[,2]),function(x) mean(as.numeric(x)))
  })
  
  data_merge2 <- rbind(uni_data2,more_data2)
}
saveRDS(data_merge2,paste0(d,"/",file,"_ensmbl.rds"))

