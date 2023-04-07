#=============移除FPKM在超过30%（包括30%）样本为0的基因，log2(FPKM+0.05)转化====================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM")

file <- dir(getwd())
file2 <- file[grepl(".rds",file)]
out_d <- "30_log_05"
dir.create(paste0(getwd(),"/",out_d))

nrow <- c()
for(i in 1:length(file2)){
  tryCatch({
    data <- readRDS(paste0(getwd(),"/",file2[i]))
    cat(paste(class(data),"\n"))
    len <- apply(data,1,function(x) length(which(x==0)))
    data2 <- data[len < round(ncol(data)*0.3),,drop=F]
    
    out_file <- gsub("\\.rds","",file2[i])
  
    data3 <- log2(data2+0.05)
    nrow <- rbind(nrow,cbind(out_file,nrow(data2),ncol(data2)))
    
    cat(paste0(rownames(data3)[1],"\n"))
    cat(paste0(colnames(data3)[1],"\n"))
    
    saveRDS(data3,paste0(out_d,"/",out_file,"_30_log_05.rds"))
    cat(paste0(out_file,"\n"))
  },error=function(e){
    cat(paste0(out_file,": wrong\n"))
  })
    
}
write.table(nrow,paste0(out_d,"/nrow.txt"),sep = "\t",quote = F,row.names = F)

#####样本编号：icgc_specimen_id-icgc_donor_id-specimen_type(N or P)####
#改！！！！！！！！！！！
#N： Normal，Primary tumour

setwd("I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM/30_log_05")

file <- dir(getwd())
file2 <- file[grepl(".rds",file)]
cancer <- stringr::str_match(file2,"exp_seq\\.(.*?)\\..*")[,2]
d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM/30_log_05/clinic"
# specimen <- dir(d)
# specimen2 <- specimen[grepl(".gz",specimen)]
# library(R.utils)
# for(j in 1:length(specimen2)){
#   gunzip(paste0(d,"/",specimen2[j]))
# }
library(dplyr)
library(plyr)
nrows <- c()
for(i in 1:length(file2)){
  
  specimen <- read.table(paste0(d,"/","specimen.",cancer[i],".tsv"),sep='\t',header = T,as.is = T,fill=T,na.strings = "")
  samples <- read.table(paste0(d,"/","sample.",cancer[i],".tsv"),sep='\t',header = T,as.is = T,fill=T,na.strings = "")
  
  #type <- as.data.frame(table(specimen$specimen_type))
  #write.table(cbind(cancer[i],type),"specimen_type.txt",sep='\t',quote = F,row.names = F,append = T,col.names = F)
  tumour <- specimen[grepl("tumour",specimen$specimen_type,ignore.case = T),
                    c("icgc_specimen_id","icgc_donor_id","specimen_type")]
  samples <- samples[,c("icgc_sample_id","icgc_specimen_id")]

  tumour_merge <- merge(tumour,samples,by= "icgc_specimen_id")

  count <- readRDS(file2[i])
  ####tumour####
  count_sam <- tumour_merge[match(colnames(count),tumour_merge$icgc_sample_id),]
  count_sam <- na.omit(count_sam)
  count_sam$submitter_id <- paste(count_sam$icgc_donor_id,count_sam$icgc_specimen_id,sep='-')
  
  frenq <- table(count_sam$submitter_id)
  if(max(frenq)==1){
    tumour_count <- count[,count_sam$icgc_sample_id]
    colnames(tumour_count) <- count_sam$submitter_id
    
  }else{
    uni_sam <- count_sam[count_sam$submitter_id %in% names(frenq[frenq==1]),,drop=F]
    uni_count <- count[,uni_sam$icgc_sample_id]
    colnames(uni_count) <- uni_sam$submitter_id
    
    more_sam <- count_sam[count_sam$submitter_id %in% names(frenq[frenq!=1]),,drop=F]
    more_count <- rbind(more_sam$submitter_id,count[,more_sam$icgc_sample_id,drop=F])
    
    more_count2 <- apply(more_count[-1,],1, function(x){
      tapply(x,factor(more_count[1,]),function(x) mean(as.numeric(x)))
    })
    
    if(is.null(nrow(more_count2))){
      more_count2 <- as.data.frame(more_count2)
      colnames(more_count2) <- levels(factor(more_count[1,]))
    }
    
    tumour_count <- cbind(uni_count,more_count2)
  }

  save(tumour_count,file=paste0(getwd(),"/reference/",cancer[i],"_tumour_30_log_05.RData"))
  ####clinic####
  clinic <- read.table(paste0(d,"/","donor.",cancer[i],".tsv"),sep='\t',header = T,as.is = T,fill=T,na.strings = "")
  clinic <- clinic[,c("icgc_donor_id","donor_vital_status","donor_survival_time")]
  clinic2 <- merge(count_sam,clinic,by="icgc_donor_id")
  clinic2 <- na.omit(clinic2)
  clinic2 <- clinic2 %>% 
    select(submitter_id,vital_status = donor_vital_status,days = donor_survival_time)
  clinic3 <- clinic2[match(colnames(tumour_count),clinic2$submitter_id),] 
  clinic3 <- na.omit(clinic3)
  tumour_count_clinic <- tumour_count[,clinic3$submitter_id,drop=F]
  save(tumour_count_clinic,file = paste0(getwd(),"/reference/",cancer[i],"_clinic_count.RData"))
  
  status <- as.data.frame(table(clinic3$vital_status))
  if(nrow(status)==2){
    clinic3[clinic3$vital_status == "deceased",]$vital_status <- 1
    clinic3[clinic3$vital_status == "alive",]$vital_status <- 0
  }else{
    if(status[,1] == "alive"){
      clinic3[clinic3$vital_status == "alive",]$vital_status <- 0
    }else{
      clinic3[clinic3$vital_status == "deceased",]$vital_status <- 1
    }
  }
  
  save(clinic3,file = paste0(getwd(),"/reference/",cancer[i],"_clinic.RData"))
  
  ####normal####
  normal <- specimen[grepl("normal",specimen$specimen_type,ignore.case = T),
                     c("icgc_specimen_id","icgc_donor_id","specimen_type")]
  normal_merge <- merge(normal,samples,by = "icgc_specimen_id")
  
  normal_sam <- normal_merge[match(colnames(count),normal_merge$icgc_sample_id),]
  normal_sam <- na.omit(normal_sam)
  
  if(nrow(normal_sam)!=0){
    normal_sam$submitter_id <- paste(normal_sam$icgc_donor_id,normal_sam$icgc_specimen_id,sep='-')
    
    frenq <- table(normal_sam$submitter_id)
    if(max(frenq)==1){
      normal_count <- count[,normal_sam$icgc_sample_id]
      colnames(normal_count) <- normal_sam$submitter_id
      
    }else{
      uni_sam <- normal_sam[normal_sam$submitter_id %in% names(frenq[frenq==1]),,drop=F]
      uni_count <- count[,uni_sam$icgc_sample_id]
      colnames(uni_count) <- uni_sam$submitter_id
      
      more_sam <- normal_sam[normal_sam$submitter_id %in% names(frenq[frenq!=1]),,drop=F]
      more_count <- rbind(more_sam$submitter_id,count[,more_sam$icgc_sample_id,drop=F])
      
      more_count2 <- apply(more_count[-1,],1, function(x){
        tapply(x,factor(more_count[1,]),function(x) mean(as.numeric(x)))
      })
      
      if(is.null(nrow(more_count2))){
        more_count2 <- as.data.frame(more_count2)
        colnames(more_count2) <- levels(factor(more_count[1,]))
      }
      normal_count <- cbind(uni_count,more_count2)
    }
    save(normal_count,file = paste0(getwd(),"/reference/",cancer[i],"_normal_30_log_05.RData"))
    nrows <- rbind(nrows,cbind(cancer[i],nrow(tumour_count),ncol(tumour_count),nrow(clinic3),ncol(normal_count)))
  }else{
    nrows <- rbind(nrows,cbind(cancer[i],nrow(tumour_count),ncol(tumour_count),nrow(clinic3),0))
  }
  cat(paste0(cancer[i],"\n"))
  
}

colnames(nrows) <- c("cancer","gene","tumor","clinic","normal")
write.table(nrows,"count_info_30_log_05.txt",sep='\t',quote = F,row.names = F)
