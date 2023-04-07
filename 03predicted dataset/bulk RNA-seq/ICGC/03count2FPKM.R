#==================================count计算FPKM====================================
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("J:/xieyunjin/lncRNA_synergisty_project/gencode.v19.annotation.gtf/gencode.v19.annotation1.gtf",format="auto")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
n=t(as.data.frame(exons_gene_lens))

out_d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count"

write.table(n,paste0(out_d,'/gene_length.txt'),col.names=F,row.names=T,quote=F,sep='\t')
length=read.table(paste0(out_d,'/gene_length.txt'),header=F,sep='\t',check.names=F)
names(length)<-c('gene','length')
length$gene <- gsub("\\..*","",length$gene)

file <- dir(out_d)
file2 <- file[grepl(".rds",file)]
for(i in 4:length(file2)){
  
  tryCatch({
    count <- readRDS(paste(out_d,file2[i],sep='/'))
    count2 <- matrix(as.numeric(as.matrix(count)),nrow=nrow(count))
    rownames(count2) <- rownames(count)
    colnames(count2) <- colnames(count)
    
    count2 <- as.data.frame(count2)
    count2 <- cbind(rownames(count2),count2)
    colnames(count2)[1] <- "gene"
    count2$gene <- gsub("\\..*","",count2$gene)
    #匹配gene_count与gene_length
    merge<-merge(count2,length,by = 'gene') 
    samples <- setdiff(colnames(merge),c("gene","length"))
    #计算每个样本的mapped reads数
    mapped_reads <- colSums(merge[,samples])
    #计算FPKM值
    FPKM <- merge[,samples]/(10^-9*matrix(as.numeric(merge[,"length"]),ncol=1) %*% matrix(mapped_reads,nrow=1))
    frenq <- table(merge$gene)
    if(max(frenq)==1){
      rownames(FPKM) <- merge$gene
    }else{
      uni_data <- FPKM[merge$gene %in% names(frenq[frenq==1]),]
      rownames(uni_data) <- merge[merge$gene %in% names(frenq[frenq==1]),,drop=F]$gene
      
      more_data <- FPKM[merge$gene %in% names(frenq[frenq!=1]),,drop=F]
      more_gene <- merge[merge$gene %in% names(frenq[frenq!=1]),,drop=F]$gene
      more_data2 <- apply(more_data,2,function(x){
        tapply(x,factor(more_gene),function(x) mean(as.numeric(x)))
      })
      
      FPKM <- rbind(uni_data,more_data2)
    }

    FPKM_out <- paste0(out_d,"/FPKM") 
    dir.create(FPKM_out)
    
    saveRDS(FPKM,paste0(FPKM_out,"/",file2[i],'_FPKM.rds'))
    cat(paste0(file2[i],"\n"))
  },error=function(e){
    cat(paste0(file2[i],": wrong\n"))
  })
  
}
