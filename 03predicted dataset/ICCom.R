##################
ICCom <- function(expPath = NULL,Path = NULL,compare_expPath = NULL,p_value = 0.1,referencePath,databaseFile = "iTALK_LR_database",
                  gene_names = "symbol",organism = "human",outpath){
  
  suppressMessages({library(jsonlite)})
  
  stopMessage <- NULL
  
  UserUploadData <- tail(unlist(strsplit(expPath,"/")),n=1)
  #webPath <- gsub(paste0("/",UserUploadData),"",expPath)
  #ID_new <- paste(tail(unlist(strsplit(webPath,"/")),n=1),collapse = "/")
  ID_new <- UserUploadData
  ## updata Job_ID.txt, append a new job and state is run
  jobFile <- paste0(referencePath,"/Job_ID.txt")
  if(file.info(jobFile)$size==0){
    Task_ID_new <- cbind(ID_new,"run")
    write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
  }else{
    Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
    Task_ID_new <- rbind(Task_ID,c(ID_new,"run"))
    write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
  }
  
  tryCatch({
    set.seed(4)
    
    if(p_value != "empty"){
      show_level = "gene" 
    }else if(compare_expPath !="empty"){
      show_level = "patient"
    }
    
    #database <- read.table(paste(referencePath,"/",databaseFile,".txt",sep=''),sep='\t',header=T,as.is=T)
    l <- load(paste(referencePath,"/",databaseFile,".RData",sep=''))
    database <- eval(parse(text = l))
    
    if(file.info(expPath)$size==0){
      stopMessage <- "wrong data."
    }else{
      #exp0 <- read.table(expPath,sep='\t',header = T,as.is = T,fill = T,strip.white = T,check.names=F)
      #exp0 <- readRDS(expPath)
      expl0 <- load(expPath)
      exp0 <- eval(parse(text = expl0))
      if((class(exp0[,1]) == "character") | (class(exp0[,1]) == "factor")){
        rownames(exp0) <- as.matrix(exp0[,1])
        exp0 <- exp0[,-1]
      }
      exp <- ensmbl_entrez_to_symbol(data = exp0,gene_names,organism,ID_transfer_Path = referencePath)
      #rank
      genes <- rownames(exp)
      sam <- colnames(exp)
      exp_rank <- lapply(1:ncol(exp),function(x) find_exp_rank(sam_exp = exp[,x],genes = genes,sam = sam[x]))
      exp_rank2 <- Reduce(function(x,y) merge(x,y,by="GeneID"),exp_rank)
      exp_rank3 <- exp_rank2[,-1]
      rownames(exp_rank3) <- as.matrix(exp_rank2[,1])
      
      if(show_level == "gene"){
      
        score_result <- Interaction_strength_gene(exp_rank = exp_rank3,pairs_immune = database)
        plot_data <- score_result[[1]]
        detail <- score_result[[2]]
        
        if(is.null(plot_data)){
          stopMessage <- "no results."
        }else{
          plot_data_sig <- plot_data[as.numeric(plot_data[,"P value"]) < as.numeric(p_value),,drop=F]
          if(nrow(plot_data_sig)==0){
            stopMessage <- "no results."
          }else{
            #===========================================output table
            rownames(plot_data_sig) <- NULL
            plot_data_sig[,"Interaction Strength"] <- round(as.numeric(plot_data_sig[,"Interaction Strength"]),4)
            
            plot_data_df <- as.data.frame(plot_data_sig)
            output <- list(data = plot_data_df)
            jsoncars <- toJSON(output, pretty=TRUE)
            cat(jsoncars, file = paste(outpath,'ICCom_gene_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
            #==============================================pie plot
            row_name <- paste(plot_data_sig[,"ligand"],plot_data_sig[,"receptor"],sep='_')
            plot_data_sig2 <- cbind(plot_data_sig,row_name)
            
            group <- cut(as.numeric(plot_data_sig2[,4]),5)
            group_name <- levels(group)
            plot_data2_group <- cbind(plot_data_sig2,group)
            
            xuhao <- unique(plot_data2_group[,"group"])
            for(i in xuhao){
              plot_data2_group[plot_data2_group[,"group"]==i,"group"] <- group_name[as.numeric(i)]
            }
            
            group_count <- aggregate(plot_data2_group[,"row_name"],list(plot_data2_group[,"group"]),length)   
            pie_output <- t(apply(group_count,1,function(x) paste("{value: ",x[2],",name: '",x[1],"'}",sep = '')))
            pie_output2 <- paste(pie_output,collapse = ",")
            write.table(pie_output2,paste(outpath,"ICCom_r_re_pie.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
            #================================================riverplot
            plot_data_sig2_sort <- plot_data_sig2[order(plot_data_sig2[,"Interaction Strength"],decreasing = T),,drop=F]
            
            if(nrow(plot_data_sig2_sort)>=50){
              plot_data_sig2_sort <- plot_data_sig2_sort[1:50,]
            }else{
              plot_data_sig2_sort <- plot_data_sig2_sort
            }
            jsoncars <- plot_riverplot(plot_data = plot_data_sig2_sort)
            cat(jsoncars, file = paste(outpath,'ICCom_r_re_riverplot.json',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
            write.table("successful",paste(outpath,"ICCom_r_erro.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
          }
        }
        if(is.null(stopMessage)){
          download_table="ICCom_gene_re.txt"
          Res_pie_file="ICCom_r_re_pie.txt"
          Res_river_file="ICCom_r_re_riverplot.json"
          resut_merge=paste0("{",
                             '"download_table" :','"',download_table,'",',
                             '"Res_pie_file" :','"',Res_pie_file,'",',
                             '"Res_river_file" :','"',Res_river_file,'",',
                             '"error_attention" :','"no',
                             '"}')
          write.table(resut_merge,paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
          
          Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
          Task_ID_new <- Task_ID
          Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "ICCom_gene+success"
          write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
          
          result_message <- paste0("{",'"error_attention" :','"success"',"}")
          return(result_message)
        }
        ####====================================================patient  
      }else if(show_level == "patient"){
        if(compare_expPath!= "empty"){
          
          if(file.info(compare_expPath)$size==0){
            stopMessage <- "wrong data."
          }else{
            compare_exp0 <- read.table(compare_expPath,sep='\t',header=T,as.is = T,fill=T,strip.white = T,check.names=F)
            if((class(compare_exp0[,1]) == "character") | (class(compare_exp0[,1]) == "factor")){
              rownames(compare_exp0) <- as.matrix(compare_exp0[,1])
              compare_exp0 <- compare_exp0[,-1]
            }
            compare_exp <- ensmbl_entrez_to_symbol(compare_exp0,gene_names,organism,ID_transfer_Path = referencePath)
            #rank
            genes <- rownames(compare_exp)
            sam <- colnames(compare_exp)
            compare_exp_rank <- lapply(1:ncol(compare_exp),function(x) find_exp_rank(sam_exp = compare_exp[,x],genes = genes,sam = sam[x]))
            compare_exp_rank2 <- Reduce(function(x,y) merge(x,y,by="GeneID"),compare_exp_rank)
            compare_exp_rank3 <- compare_exp_rank2[,-1]
            rownames(compare_exp_rank3) <- as.matrix(compare_exp_rank2[,1])
          }
        }else{
          compare_exp_rank3 = NULL
        }
        
        if(is.null(stopMessage)){
          score_result <- Interaction_strength_patient(exp_rank = exp_rank3,compare_exp_rank = compare_exp_rank3,pairs_immune = database)
          if(is.null(score_result)){
            stopMessage <- "no results."
          }else{
            ####=====================================================boxplot
            box_res <- plot_box_data(all = score_result)
            write.table(box_res,paste(outpath,"/ICCom_patient_boxplot.txt",sep=''),sep='\t',quote = F,row.names = F,col.names = F)
            write.table("successful",paste(outpath,"ICCom_r_erro.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
          }
        }
        if(is.null(stopMessage)){
          boxplot_table="ICCom_patient_boxplot.txt"
          resut_merge=paste0("{",
                             '"boxplot_table" :','"',boxplot_table,'",',
                             '"error_attention" :','"no',
                             '"}')
          write.table(resut_merge,paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
          
          Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
          Task_ID_new <- Task_ID
          Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "ICCom_patient+success"
          write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
          
          result_message <- paste0("{",'"error_attention" :','"success"',"}")
          return(result_message)
        }
      }else{
        stopMessage <- "wrong choices."
      }
    }
    write.table(paste0("{",'"error_attention" :','"',stopMessage,'"}'),paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
    
    Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
    Task_ID_new <- Task_ID
    Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "dead"
    write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
    
    result_message <- paste0("{",'"error_attention" :','"',stopMessage,'"}')
    return(result_message)
  },error = function(e){
    stopMessage <- "task quit."
    write.table(paste0("{",'"error_attention" :','"',stopMessage,'"}'),paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
    
    Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
    Task_ID_new <- Task_ID
    Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "dead"
    write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
    
    result_message <- paste0("{",'"error_attention" :','"',stopMessage,'"}')
    return(result_message)
  })
}
######################function
ensmbl_entrez_to_symbol <- function(data,gene_names,organism,ID_transfer_Path=NULL){
  gene <- rownames(data)
  if(organism == "human"){
    if(gene_names == "symbol"){
      if(grepl("ENSG",gene[1])){
        gene_names = "ensembl"
      }else if(!grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])){
        gene_names = "entrez"
      }
    }else if(gene_names == "ensembl"){
      if(!grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])){
        gene_names = "entrez"
      }else if(grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])& (!grepl("ENSG",gene[1]))){
        gene_names = "symbol"
      }
    }else if(gene_names == "entrez"){
      if(grepl("ENSG",gene[1])){
        gene_names = "ensembl"
      }else if(grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])& (!grepl("ENSG",gene[1]))){
        gene_names = "symbol"
      }
    }
  }else if(organism == "mouse"){
    if(gene_names == "symbol"){
      if(grepl("ENSMUSG",gene[1])){
        gene_names = "ensembl"
      }else if(!grepl(paste("[",paste(letters,collapse ='|'),"]",sep=''),gene[1],ignore.case=TRUE)){
        gene_names = "entrez"
      }
    }else if(gene_names == "ensembl"){
      if(!grepl(paste("[",paste(letters,collapse ='|'),"]",sep=''),gene[1],ignore.case=TRUE)){
        gene_names = "entrez"
      }else if(grepl(paste("[",paste(letters,collapse ='|'),"]",sep=''),gene[1])){
        gene_names = "symbol"
      }
    }else if(gene_names == "entrez"){
      if(grepl("ENSMUSG",gene[1])){
        gene_names = "ensembl"
      }else if(grepl(paste("[",paste(letters,collapse ='|'),"]",sep=''),gene[1])){
        gene_names = "symbol"
      }
    }
  }
  
  if(organism == "human"){
    l2 <- load(paste(ID_transfer_Path,"ID_transfer_file.RData",sep='/'))
    ID_transfer_data <- eval(parse(text = l2))
    #ID_transfer_data <- read.table(paste(ID_transfer_Path,"ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
  }else{
    l2 <- load(paste(ID_transfer_Path,"mouse_ID_transfer_file.RData",sep='/'))
    ID_transfer_data <- eval(parse(text = l2))
    #ID_transfer_data <- read.table(paste(ID_transfer_Path,"mouse_ID_transfer_file.txt",sep='/'),sep='\t',header = T,as.is = T,fill=T,strip.white = T,quote = "",check.names = F)
  }
  
  if(gene_names == "ensembl"){
    if(grepl("\\.",gene[1])){
      gene <- gsub("\\..*","",gene)
    }
    res <- unique(ID_transfer_data[ID_transfer_data[,2] %in% gene,c(2,3)])
    res2 <- res[match(rownames(data),res[,1]),]
    res2 <- na.omit(res2)
    frenq <- table(res2[,2])
    
    if(max(frenq)>1){
      exp_gene <- data[res2[,1],,drop=F]
      exp_gene2 <- cbind(res2,exp_gene)
      
      uni_exp <- exp_gene2[match(names(frenq[frenq==1]),exp_gene2[,2]),-c(1,2),drop=F]
      rownames(uni_exp) <- names(frenq[frenq==1])
      
      more_exp <- exp_gene2[exp_gene2[,2] %in% names(frenq[frenq>1]),,drop=F]
      more_exp2 <- apply(more_exp[,-c(1,2)],2,function(x) {
        tapply(x, factor(more_exp[,2]), function(x) mean(as.numeric(x))) 
      })
      data2 <- rbind(uni_exp,more_exp2)
    }else{
      data2 <- data[rownames(data) %in% res2[,1],,drop=F]
      rownames(data2) <- res2[match(rownames(data2),res2[,1]),2]
    }
  }else if(gene_names == "entrez"){
    res <- unique(ID_transfer_data[ID_transfer_data[,1] %in% gene,c(1,3)])
    res2 <- res[match(rownames(data),res[,1]),]
    res2 <- na.omit(res2)
    data2 <- data[as.character(res2[,1]),]
    rownames(data2) <- res2[,2]
    
  }else if(gene_names == "symbol"){
    data2 <- data
  }
  return(data2)
}
###########
find_exp_rank <- function(sam_exp,genes,sam){
  sam_exp2 <- cbind(genes,sam_exp)
  
  sam_exp_sort <- sam_exp2[order(as.numeric(sam_exp2[,2])),c(1,2)]
  zhici <- as.data.frame(cut(seq(as.numeric(sam_exp_sort[,2])), 10, labels = FALSE))
  sam_exp_sort_zhici <- cbind(sam_exp_sort,zhici)
  
  colnames(sam_exp_sort_zhici) <- c("GeneID","Exp",sam)
  
  return(sam_exp_sort_zhici[,c(1,3)])
}
######################patient
Interaction_strength_patient <- function(exp_rank,compare_exp_rank = NULL,pairs_immune){
  
  
  gene1 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune[,"ligand"]]
  gene2 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune[,"receptor"]]
  ###############
  index1 <- which(pairs_immune[,"ligand"] %in% gene1)
  index2 <- which(pairs_immune[,"receptor"]%in% gene2)
  index <- intersect(index1,index2)
  
  if(length(index)==0){
    return(plot_data = NULL)
  }
  
  pairs <- unique(pairs_immune[index,])
  nlen <- nrow(pairs)
  
  score_result <- lapply(1:nlen,function(x){
    zhici1 <- as.numeric(exp_rank[pairs[x,"ligand"],])
    zhici2 <- as.numeric(exp_rank[pairs[x,"receptor"],])
    score <- (zhici1+zhici2) * (10-abs(zhici2-zhici1))
    return(score)
  })
  score_result <- do.call(rbind,score_result)
  colnames(score_result) <- colnames(exp_rank)
  rownames(score_result) <- pairs[,"pair"]
  
  ####
  sample_level <- colMeans(score_result)
  
  if(!is.null(compare_exp_rank)==T){
    gene1 <- rownames(compare_exp_rank)[rownames(compare_exp_rank) %in% pairs_immune[,"ligand"]]
    gene2 <- rownames(compare_exp_rank)[rownames(compare_exp_rank) %in% pairs_immune[,"receptor"]]
    
    index1 <- which(pairs_immune[,"ligand"] %in% gene1)
    index2 <- which(pairs_immune[,"receptor"] %in% gene2)
    index <- intersect(index1,index2)
    
    if(length(index)!=0){
      normal_pairs <- unique(pairs_immune[index,])
      normal_nlen <- nrow(normal_pairs)
      
      normal_score_result <- lapply(1:normal_nlen,function(x){
        zhici1 <- as.numeric(compare_exp_rank[normal_pairs[x,"ligand"],])
        zhici2 <- as.numeric(compare_exp_rank[normal_pairs[x,"receptor"],])
        score <- (zhici1+zhici2) * (10-abs(zhici2-zhici1))
        return(score)
      })
      normal_score_result <- do.call(rbind,normal_score_result)
      colnames(normal_score_result) <- colnames(compare_exp_rank)
      rownames(normal_score_result) <- normal_pairs[,"pair"]
      normal_sample_level <- colMeans(normal_score_result)
      #############boxplot
    }else{
      normal_sample_level <- 0
    }
  }else{
    normal_sample_level <- 0
  }

  plot_data <- list(cancer = "User cancer",cancer_box = sample_level,normal_box = normal_sample_level)

  return(plot_data = plot_data)
}
####################
plot_box_data <- function(all = plot_data){
  
    if(length(all[[3]]) == 1){
      p <- 1
    }else{
      p <- t.test(as.numeric(all[[2]]),as.numeric(all[[3]]),alternative = "two.sided")$p.value
    }
    xx <- paste("User cancer(p=",round(p,4),")",sep='')
    cancer_s <-paste("'",xx,"'",sep='')
    
    real_box_s0 <- paste(round(as.numeric(all[[2]]),4),collapse=',')
    real_box_s <- paste("[[",real_box_s0,"]]",sep='')
    
    random_box_s0 <- paste(round(as.numeric(all[[3]]),4),collapse=',')
    random_box_s <- paste("[[",random_box_s0,"]]",sep='')
    return(t(c(cancer_s,real_box_s,random_box_s)))
}
######################gene
Interaction_strength_gene <- function(exp_rank,pairs_immune){

  gene1 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune[,"ligand"]]
  gene2 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune[,"receptor"]]
  
  index1 <- which(pairs_immune[,"ligand"] %in% gene1)
  index2 <- which(pairs_immune[,"receptor"] %in% gene2)
  index <- intersect(index1,index2)
  
  if(length(index)==0){
    return(list(plot_data = NULL,detail = NULL))
  }
  
  pairs <- unique(pairs_immune[index,])
  nlen <- nrow(pairs)
  
  score_result <- lapply(1:nlen,function(x){
    zhici1 <- as.numeric(exp_rank[pairs[x,"ligand"],])
    zhici2 <- as.numeric(exp_rank[pairs[x,"receptor"],])
    score <- (zhici1+zhici2) * (10-abs(zhici2-zhici1))
    return(score)
  })
  score_result <- do.call(rbind,score_result)
  colnames(score_result) <- colnames(exp_rank)
  rownames(score_result) <- pairs[,"pair"]
  
  ####gene interaction level
  gene_pairs_level <-  rowMeans(score_result)
  
  # ####random 10000
   r_score_result <- lapply(1:10000,function(x){
      r1 <- sample(1:nrow(exp_rank),1)
      r2 <- sample(1:nrow(exp_rank),1)

      zhici1 <- as.numeric(exp_rank[r1,])
      zhici2 <- as.numeric(exp_rank[r2,])
      score <- (zhici1+zhici2) * (10-abs(zhici2-zhici1))

      r_mean <- mean(score)
      return(r_mean)
    })
    r_score <- unlist(r_score_result)
  
  #########p value
  p_value <- lapply(1:length(gene_pairs_level),function(x){
    length(which(r_score > gene_pairs_level[x]))/10000
  })
  p_value <- unlist(p_value)
  gene_pairs_pvalue <- cbind(gene_pairs_level,p_value)
  
  gene_level_result <- cbind(pairs[,"ligand"],pairs[,"receptor"],"-",gene_pairs_pvalue,pairs[,"classification"])
  colnames(gene_level_result) <- c("ligand","receptor","Cancer","Interaction Strength","P value","Funcion")

  return(list(plot_data = gene_level_result,detail = score_result))
}
##############
plot_riverplot <- function(plot_data){
  plot_data <- as.matrix(plot_data)

  plot_data[,"ligand"] <- paste(plot_data[,"ligand"],"(ligand)",sep='')
  plot_data[,"receptor"] <- paste(plot_data[,"receptor"],"(receptor)",sep='')
  
  node <- matrix(unique(c(plot_data[,"ligand"],plot_data[,"receptor"],plot_data[,c("Funcion")])),ncol=1)
  colnames(node) <- "name"
  node <- as.data.frame(node)
  
  index <- rep(1,times = nrow(plot_data))
  plot_data <- cbind(index,plot_data)
  
  data_ligand <- aggregate(plot_data[,"index"],list(plot_data[,"ligand"],plot_data[,c("Funcion")]),length)
  colnames(data_ligand) <- c("source","target","value")
  
  data_Funcion <- aggregate(plot_data[,"index"],list(plot_data[,c("Funcion")],plot_data[,"receptor"]),length)
  colnames(data_Funcion) <- c("source","target","value")
  
  link <- Reduce(rbind,list(data_ligand,data_Funcion))
  
  output <- list(nodes=node,links=link)
  jsoncars <- toJSON(output, pretty=TRUE)
  return(jsoncars)
}

