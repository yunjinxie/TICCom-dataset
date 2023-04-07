##################human 构建GENE ID 转化文件################
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/symbolChangeData")
#NCBI ftp下载的
gene2ensembl <- read.table("gene2ensembl",sep='\t',header=F,as.is = T)
human_gene2ensembl <- gene2ensembl[gene2ensembl[,1]=="9606",]
human_gene2ensembl2 <- unique(human_gene2ensembl[,c(1,2,3)])
colnames(human_gene2ensembl2) <- c("Species","Entrez ID","Ensembl")

# frenq1 <- table(human_gene2ensembl2[,2])
# entrez <- names(frenq1)[frenq1>1]
# entrez_morethan1 <- human_gene2ensembl2[human_gene2ensembl2[,2] %in% entrez,]

# frenq2 <- table(human_gene2ensembl2[,3])
# ENSG <- names(frenq2)[frenq2>1]
# ENSG_morethan1 <- human_gene2ensembl2[human_gene2ensembl2[,3] %in% ENSG,]


gene2accession <- read.table("human_symbol2ID2.txt",sep='\t',header=T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
frenq <- table(gene2accession[,2])
aa <- gene2accession[gene2accession[,2] %in% names(frenq[frenq>1]),]

# frenq3 <- table(gene2accession[,1])
# names(frenq3[frenq3>1])


# gene2accession_mouse <- read.table("mouse_symbol2ID2.txt",sep='\t',header=F,as.is = T,fill=T,strip.white = T,quote = "")
# frenq3 <- table(gene2accession_mouse[,3])
# names(frenq3[frenq3>1])

# aa2 <- gene2accession_mouse[gene2accession_mouse[,2] %in% names(frenq2[frenq2>1]),]

aa3 <- merge(human_gene2ensembl2,gene2accession,by="Entrez ID") 
aa4 <- aa3[,c(1,3,4)]
write.table(aa4,"ID_transfer_file.txt",sep='\t',quote = F,row.names = F)

ff <- table(aa4$Ensembl)
aa5 <- aa4[aa4[,2] %in% names(ff[ff>1]),]
aa5[order(aa5[,2]),]

##################mouse 构建GENE ID 转化文件################
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/symbolChangeData")

gene2ensembl <- read.table("gene2ensembl",sep='\t',header=F,as.is = T)
mouse_gene2ensembl <- gene2ensembl[gene2ensembl[,1]=="10090",]
mouse_gene2ensembl2 <- unique(mouse_gene2ensembl[,c(1,2,3)])
colnames(mouse_gene2ensembl2) <- c("Species","Entrez ID","Ensembl")

# frenq1 <- table(human_gene2ensembl2[,1])
# entrez <- names(frenq1)[frenq1>1]
# entrez_morethan1 <- human_gene2ensembl2[human_gene2ensembl2[,1] %in% entrez,]

gene2accession <- read.table("mouse_symbol2ID2.txt",sep='\t',header=T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
frenq <- table(gene2accession[,2])
aa <- gene2accession[gene2accession[,2] %in% names(frenq[frenq>1]),]

# gene2accession_mouse <- read.table("mouse_symbol2ID2.txt",sep='\t',header=F,as.is = T,fill=T,strip.white = T,quote = "")
# frenq2 <- table(gene2accession_mouse[,2])
# aa2 <- gene2accession_mouse[gene2accession_mouse[,2] %in% names(frenq2[frenq2>1]),]

aa3 <- merge(mouse_gene2ensembl2,gene2accession,by="Entrez ID") 
aa4 <- aa3[,c(1,3,4)]
write.table(aa4,"mouse_ID_transfer_file.txt",sep='\t',quote = F,row.names = F)

ff <- table(aa4$Ensembl)
aa5 <- aa4[aa4[,2] %in% names(ff[ff>1]),]
aa5[order(aa5[,2]),]
