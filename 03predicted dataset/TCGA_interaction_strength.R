#show_level : Patient (sample),Gene Interaction,
###########interaction strength
TCGA_interaction_strength <- function(IDPath=NULL,expPath=NULL,compare_expPath = NULL,gene_listPath=NULL,cancer=NULL,referencePath=NULL,show_level="gene",outpath){
 
   suppressMessages({
    library(jsonlite)
    library(ggplot2)
  })
  
  stopMessage <- NULL
  #ID
  ID_new <- tail(unlist(strsplit(IDPath,"/")),n=1)
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
    
    if(gene_listPath!="empty"){
      gene_list <- read.table(gene_listPath,sep='\t',header = F,as.is = T,fill = T,strip.white = T)
      gene_list <- as.matrix(gene_list)
      ############ABL1-CSF1;LCAM-SDCBP
      gene_list2 <- Reduce(rbind,strsplit(gene_list,"_|:|&"))
      
      if(is.null(ncol(gene_list2))){
        if(length(gene_list2)==2){
          gene_list2 <- matrix(gene_list2,nrow=1,ncol=2)
          colnames(gene_list2) <- c("UserID1","UserID2")
          gene1 <- ID_transfer(gene_list = gene_list2[,1],ID_transfer_Path = referencePath)
          gene2 <- ID_transfer(gene_list = gene_list2[,2],ID_transfer_Path = referencePath)
          if(nrow(gene1)==0|nrow(gene2)==0){
            stopMessage <- "no matched gene symbols."
          }else{
            colnames(gene1) <- c("UserID1","Gene1Name","Gene1Id")
            gene_list3 <- merge(gene_list2,gene1,by='UserID1')
            colnames(gene2) <- c("UserID2","Gene2Name","Gene2Id")
            gene_list4 <- merge(gene_list3,gene2,by='UserID2')
            if(nrow(gene_list4)==0){
              stopMessage <- "no matched gene symbols."
            }else{
              pairs2 <- gene_list4[,c("Gene1Name","Gene1Id","Gene2Name","Gene2Id"),drop=F]
            }
          }
        }else{
          stopMessage <- "too few gene interactions input."
        }
      }else{
        colnames(gene_list2) <- c("UserID1","UserID2")
        gene1 <- ID_transfer(gene_list = gene_list2[,1],ID_transfer_Path = referencePath)
        gene2 <- ID_transfer(gene_list = gene_list2[,2],ID_transfer_Path = referencePath)
        if(nrow(gene1)==0|nrow(gene2)==0){
          stopMessage <- "no matched gene symbols."
        }else{
          colnames(gene1) <- c("UserID1","Gene1Name","Gene1Id")
          gene_list3 <- merge(gene_list2,gene1,by='UserID1')
          colnames(gene2) <- c("UserID2","Gene2Name","Gene2Id")
          gene_list4 <- merge(gene_list3,gene2,by='UserID2')
          if(nrow(gene_list4)==0){
            stopMessage <- "no matched gene symbols."
          }else{
            pairs2 <- gene_list4[,c("Gene1Name","Gene1Id","Gene2Name","Gene2Id"),drop=F]
          }
        }
      }
    }else{
      l <- load(paste(referencePath,"immune_cancer_genePairs.RData",sep='/'))
      pairs <- eval(parse(text = l))
      pairs2 <- pairs
    }
    
    if(is.null(stopMessage)){
      #=====================================================================
      show_level <- unlist(strsplit(show_level,","))
      
      if(sum(cancer %in% "empty")==0){
        cancer <- unlist(strsplit(cancer,","))
        TCGARankPath <- paste(referencePath,"/",cancer,"_sort.RData",sep='')
        exp_rank <- lapply(TCGARankPath,function(x){
          expl <- load(x)
          e_rank <- eval(parse(text = expl))
          return(e_rank)
        })
        ###########################
        if(sum(show_level %in% "patient")==1){
          
          normal_cancer <- c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC","LUAD",
                             "LUSC","PAAD","PCPG","PRAD","READ","SARC","STAD","THCA","THYM","UCEC")
          
          compare_exp_rank <- lapply(cancer,function(x) {
            inter_cancer <- normal_cancer[normal_cancer %in% x]
            if(length(inter_cancer)==0){
              return(NULL)
            }else{
              compare_TCGARankPath <- paste(referencePath,"/",inter_cancer,"_normal_sort.RData",sep='')
              com_expl <- load(compare_TCGARankPath)
              data <- eval(parse(text = com_expl))
              return(data)
            }
          })
          
          score_result <- lapply(1:length(cancer),function(x) Interaction_strength_patient(exp_rank = exp_rank[[x]],compare_exp_rank = compare_exp_rank[[x]],pairs_immune = pairs2,cancer = cancer[x]))
          names(score_result) <- cancer
          
          is_null <- unlist(lapply(score_result,function(x) is.null(x[[1]])))
          cancer2 <- cancer[!is_null]
          
          if(length(cancer2)==0){
            stopMessage <- "no results."
          }else{
            plot_data <- lapply(cancer2,function(x) score_result[[x]][[1]])
            names(plot_data) <- cancer2
            box_res <- plot_box_data(all = plot_data,choose_cancer = TRUE)
            write.table(box_res,paste(outpath,"/TCGA_interaction_strength_patient_boxplot.txt",sep=''),sep='\t',quote = F,row.names = F,col.names = F)
          }
        }else if(sum(show_level %in% "gene")==1){
          
          random_strength_path <- paste(referencePath,"/",cancer,"_random_strength.RData",sep='')
          random_strength <- lapply(random_strength_path,function(x){
            randl <- load(x)
            e_rand <- eval(parse(text = randl))
            return(e_rand)
          })
          score_result <- lapply(1:length(cancer),function(x) Interaction_strength_gene(exp_rank = exp_rank[[x]],pairs_immune = pairs2,cancer = cancer[x],r_score = random_strength[[x]],outpath))
          names(score_result) <- cancer
          is_null <- unlist(lapply(score_result,function(x) is.null(x[[1]])))
          cancer2 <- cancer[!is_null]
          
          if(length(cancer2)==0){
            stopMessage <- "no results."
          }else{
            plot_data <- lapply(cancer2,function(x) score_result[[x]][[1]])
            plot_data <- do.call(rbind,plot_data)
            detail <- lapply(cancer2,function(x) score_result[[x]][[2]])
            names(detail) <- cancer2
            ####=========================================output table
            rownames(plot_data) <- NULL
            plot_data[,"Interaction Strength"] <- round(as.numeric(plot_data[,"Interaction Strength"]),4)
            plot_data[,"P value"] <- signif(as.numeric(plot_data[,"P value"]),4)
            
            plot_data_df <- as.data.frame(plot_data)
            output <- list(data = plot_data_df)
            jsoncars <- toJSON(output, pretty=TRUE)
            cat(jsoncars, file = paste(outpath,'TCGA_interaction_strength_gene_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
            ####=============================================heatmap
            row_name <- paste(plot_data[,"Gene1Name"],plot_data[,"Gene2Name"],sep='_')
            plot_data2 <- cbind(plot_data,row_name)
            write.table(plot_data2[,c(1:5)],paste(outpath,"/TCGA_interaction_strength_gene_heatmap_download.txt",sep=''),sep='\t',row.names = F,quote = F)
            
            if(length(unique(as.matrix(plot_data2[,"Cancer"])))==1){
              plot_data2 <- plot_data2[order(as.numeric(plot_data2[,"Interaction Strength"]),decreasing = T),,drop=F]
              
              if(nrow(plot_data2) >= 25){
                plot_data_final <- plot_data2[1:25,]
              }else{
                plot_data_final <- plot_data2
              }
            }else{
              sig <- matrix(0,nrow = nrow(plot_data2))
              sig[as.numeric(plot_data2[,"P value"]) < 0.05,] <- 1
              plot_data_sig <- cbind(plot_data2,sig)
              colnames(plot_data_sig)[7] <- "sig"
              sig_count <- aggregate(plot_data_sig[,"sig"],list(plot_data_sig[,"row_name"]),function(x) length(which(x==1)))
              sig_count_sort <- sig_count[order(sig_count[,2],decreasing = T),,drop=F]
              
              if(nrow(sig_count_sort)>=25){
                sig_count_sort2 <- sig_count_sort[1:25,]
              }else{
                sig_count_sort2 <- sig_count_sort
              }
              sig_count_sort2 <- as.matrix(sig_count_sort2)
              plot_data_final <- plot_data2[plot_data2[,"row_name"] %in% sig_count_sort2[,1],]
            }
            colnames(plot_data_final)[4] <- "Value"
            
            heatmap_res <- plot_all_heatmap(plot_score = plot_data_final,choose_cancer = TRUE)
            write.table(heatmap_res,paste(outpath,"/TCGA_interaction_strength_gene_heatmap.txt",sep=''),sep='\t',row.names = F,col.names = F,quote = F)
            
            colnames(plot_data2)[4] <- "Value"
            heatmap_res_all <- plot_all_heatmap(plot_score = plot_data2,choose_cancer = TRUE)
            write.table(heatmap_res_all,paste(outpath,"/TCGA_interaction_strength_gene_heatmap_all.txt",sep=''),sep='\t',row.names = F,col.names = F,quote = F)
            ####======================================================violin
            p <- lapply(cancer2,function(x){
              plot_cancer <- plot_data2[plot_data2[,"Cancer"]==x,,drop=F]
              plot_cancer <- as.matrix(plot_cancer)
              plot_violin(plot_matrix = plot_cancer,cancer = x,cancer_strength = detail[[x]],outpath = outpath,level_name="gene")
            })
          }
          ##########################
        }else if(length(which(c("APC","B cells","DC cells","leukocyte","lymphocytes","macrophage","mast cells",
                                "MDSCs","myeloid cells","neutrophils","NK cells","T cells") %in% show_level))!=0){
          
          if(gene_listPath!="empty"){
            stopMessage <- "no immune cells information."
          }else{
            score_result <- lapply(1:length(cancer),function(x) Interaction_strength_immune_cell(exp_rank = exp_rank[[x]],pairs_immune = pairs2,cancer = cancer[x],show_level = show_level))
            
            names(score_result) <- cancer
            is_null <- unlist(lapply(score_result,function(x) is.null(x[[1]])))
            cancer2 <- cancer[!is_null]
            
            if(length(cancer2)==0){
              stopMessage <- "no results."
            }else{
              plot_data <- lapply(cancer2,function(x) score_result[[x]][[1]])
              plot_data <- do.call(rbind,plot_data)
              detail <- lapply(cancer2,function(x) score_result[[x]][[2]])
              names(detail) <- cancer2
              ####========================================================output table
              rownames(plot_data) <- NULL
              plot_data[,"Interaction Strength"] <- round(as.numeric(plot_data[,"Interaction Strength"]),4)
              plot_data_df <- as.data.frame(plot_data)
              output <- list(data = plot_data_df)
              jsoncars <- toJSON(output, pretty=TRUE)
              cat(jsoncars, file = paste(outpath,'TCGA_interaction_strength_immcell_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
              ####========================================================violin
              p <- lapply(cancer2,function(x){
                plot_cancer <- as.matrix(plot_data[plot_data[,"Cancer"]==x,,drop=F])
                row_name <- paste(plot_cancer[,"Gene1Name"],plot_cancer[,"Gene2Name"],sep='_')
                plot_cancer2 <- cbind(plot_cancer,row_name)
                plot_violin(plot_matrix = plot_cancer2,cancer = x,cancer_strength = detail[[x]],outpath = outpath,level_name="immune cell")
              })
              ####=========================================================heatmap
              cancer_immune <- aggregate(plot_data[,"Interaction Strength"],list(plot_data[,"Cancer"],plot_data[,"Immune cell"]),function(x) mean(as.numeric(x)))
              colnames(cancer_immune) <- c("Cancer","row_name","Value")
              cancer_immune <- as.matrix(cancer_immune)
              cancer_immune[,"Value"] <- round(as.numeric(cancer_immune[,"Value"]),4)
              
              heatmap_res <- plot_all_heatmap(plot_score = cancer_immune,choose_cancer = TRUE)   
              write.table(heatmap_res,paste(outpath,"/TCGA_interaction_strength_immcell_heatmap.txt",sep=''),sep='\t',row.names = F,col.names = F,quote = F)
              write.table(plot_data,paste(outpath,"/TCGA_interaction_strength_immcell_heatmap_download.txt",sep=''),sep='\t',row.names = F,quote = F)
            }
          }
        }
      }else if(expPath!="empty"){
        #exp0 <- read.table(expPath,sep='\t',header=T,as.is = T,fill=T,strip.white = T,check.names=F)
        #exp0 <- readRDS(expPath)
        expl0 <- load(expPath)
        exp0 <- eval(parse(text = expl0))
        if((class(exp0[,1]) == "character") | (class(exp0[,1]) == "factor")){
          rownames(exp0) <- as.matrix(exp0[,1])
          exp0 <- exp0[,-1]
        }
        gene <- rownames(exp0)
        gene2 <- ID_transfer(gene_list = gene,ID_transfer_Path = referencePath)
        if(nrow(gene2)==0){
          stopMessage <- "no matched gene ENSGs."
        }else{
          gene3 <- gene2[match(rownames(exp0),gene2[,1]),]
          gene3 <- na.omit(gene3)
          frenq <- table(gene3[,3])
          
          if(max(frenq)>1){
            exp_gene <- exp0[gene3[,1],,drop=F]
            exp_gene2 <- cbind(gene3,exp_gene)
            
            uni_exp <- exp_gene2[match(names(frenq[frenq==1]),exp_gene2[,3]),-c(1,2,3),drop=F]
            rownames(uni_exp) <- names(frenq[frenq==1])
            
            more_exp <- exp_gene2[exp_gene2[,3] %in% names(frenq[frenq>1]),,drop=F]
            more_exp2 <- apply(more_exp[,-c(1,2,3)],2,function(x) {
              tapply(x, factor(more_exp[,3]), function(x) mean(as.numeric(x))) 
            })
            exp <- rbind(uni_exp,more_exp2)
          }else{
            exp <- exp0[rownames(exp0) %in% gene3[,1],,drop=F]
            rownames(exp) <- gene3[match(rownames(exp),gene3[,1]),3]
          }
          #rank
          genes <- rownames(exp)
          sam <- colnames(exp)
          exp_rank <- lapply(1:ncol(exp),function(x) find_exp_rank(sam_exp = exp[,x],genes = genes,sam = sam[x]))
          exp_rank2 <- Reduce(function(x,y) merge(x,y,by="GeneID"),exp_rank)
          exp_rank3 <- exp_rank2[,-1]
          rownames(exp_rank3) <- as.matrix(exp_rank2[,1])
          
          save(exp_rank3,file = paste0(outpath,"/sort.RData"))
          ########################################
          if(sum(show_level %in% "patient")==1){
            
            if(compare_expPath!= "empty"){
              compare_exp0 <- read.table(compare_expPath,sep='\t',header=T,as.is = T,fill=T,strip.white = T,check.names=F)
              if((class(compare_exp0[,1]) == "character") | (class(compare_exp0[,1]) == "factor")){
                rownames(compare_exp0) <- as.matrix(compare_exp0[,1])
                compare_exp0 <- compare_exp0[,-1]
              }
              gene <- rownames(compare_exp0)
              gene2 <- ID_transfer(gene_list = gene,ID_transfer_Path = referencePath)
              if(nrow(gene2)==0){
                stopMessage <- "no matched gene ENSGs."
              }
              gene3 <- gene2[match(rownames(compare_exp0),gene2[,1]),]
              gene3 <- na.omit(gene3)
              frenq <- table(gene3[,3])
              
              if(max(frenq)>1){
                exp_gene <- compare_exp0[gene3[,1],,drop=F]
                exp_gene2 <- cbind(gene3,exp_gene)
                
                uni_exp <- exp_gene2[match(names(frenq[frenq==1]),exp_gene2[,3]),-c(1,2,3),drop=F]
                rownames(uni_exp) <- names(frenq[frenq==1])
                
                more_exp <- exp_gene2[exp_gene2[,3] %in% names(frenq[frenq>1]),,drop=F]
                more_exp2 <- apply(more_exp[,-c(1,2,3)],2,function(x) {
                  tapply(x, factor(more_exp[,3]), function(x) mean(as.numeric(x))) 
                })
                compare_exp <- rbind(uni_exp,more_exp2)
              }else{
                compare_exp <- compare_exp0[rownames(compare_exp0) %in% gene3[,1],,drop=F]
                rownames(compare_exp) <- gene3[match(rownames(compare_exp),gene3[,1]),3]
              }
              #rank
              genes <- rownames(compare_exp)
              sam <- colnames(compare_exp)
              compare_exp_rank <- lapply(1:ncol(compare_exp),function(x) find_exp_rank(sam_exp = compare_exp[,x],genes = genes,sam = sam[x]))
              compare_exp_rank2 <- Reduce(function(x,y) merge(x,y,by="GeneID"),compare_exp_rank)
              compare_exp_rank3 <- compare_exp_rank2[,-1]
              rownames(compare_exp_rank3) <- as.matrix(compare_exp_rank2[,1])
            }else{
              compare_exp_rank3 = NULL
            }
            
            score_result <- Interaction_strength_patient(exp_rank = exp_rank3,compare_exp_rank = compare_exp_rank3,pairs_immune = pairs2)
            plot_data <- score_result[[1]]
            
            if(is.null(plot_data)){
              stopMessage <- "no results."
            }else{
              ####================================================boxplot
              box_res <- plot_box_data(all = plot_data,choose_cancer = FALSE)  
              write.table(box_res,paste(outpath,"/TCGA_interaction_strength_patient_boxplot.txt",sep=''),sep='\t',quote = F,row.names = F,col.names = F)
            }
            ##############################
          }else if(sum(show_level %in% "gene")==1){
            
            score_result <- Interaction_strength_gene(exp_rank = exp_rank3,pairs_immune = pairs2,outpath = outpath)
            plot_data <- score_result[[1]]
            detail <- score_result[[2]]

            if(is.null(plot_data)){
              stopMessage <- "no results."
            }else{
              ####==============================================output table
              rownames(plot_data) <- NULL
              plot_data[,"Interaction Strength"] <- round(as.numeric(plot_data[,"Interaction Strength"]),4)
              plot_data[,"P value"] <- signif(as.numeric(plot_data[,"P value"]),4)
              
              plot_data_df <- as.data.frame(plot_data)
              output <- list(data = plot_data_df)
              jsoncars <- toJSON(output, pretty=TRUE)
              cat(jsoncars, file = paste(outpath,'TCGA_interaction_strength_gene_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
              ####==================================================heatmap
              row_name <- paste(plot_data[,"Gene1Name"],plot_data[,"Gene2Name"],sep='_')
              plot_data2 <- cbind(plot_data,row_name)
              write.table(plot_data2[,c(1:5),drop=F],paste(outpath,"/TCGA_interaction_strength_gene_heatmap_download.txt",sep=''),sep='\t',row.names = F,quote = F)
              
              plot_data2 <- plot_data2[order(as.numeric(plot_data2[,"Interaction Strength"]),decreasing = T),,drop=F]
              
              if(nrow(plot_data2) >= 25){
                plot_data_final <- plot_data2[1:25,]
              }else{
                plot_data_final <- plot_data2
              }
              colnames(plot_data_final)[4] <- "Value"
              
              heatmap_res <- plot_all_heatmap(plot_score = plot_data_final,choose_cancer = FALSE)
              write.table(heatmap_res,paste(outpath,"/TCGA_interaction_strength_gene_heatmap.txt",sep=''),sep='\t',row.names = F,col.names = F,quote = F)
              
              colnames(plot_data2)[4] <- "Value"
              heatmap_res_all <- plot_all_heatmap(plot_score = plot_data2,choose_cancer = FALSE)
              write.table(heatmap_res_all,paste(outpath,"/TCGA_interaction_strength_gene_heatmap_all.txt",sep=''),sep='\t',row.names = F,col.names = F,quote = F)
              ####=====================================================violin plot
              plot_data2 <- as.matrix(plot_data2)
              plot_violin(plot_matrix = plot_data2,cancer_strength = detail,outpath = outpath,level_name ="gene")
            }
          }else if(length(which(c("APC","B cells","DC cells","leukocyte","lymphocytes","macrophage","mast cells",
                                  "MDSCs","myeloid cells","neutrophils","NK cells","T cells") %in% show_level))!=0){
            
            if(gene_listPath!="empty"){
              stopMessage <- "no immune cells information."
            }else{
              score_result <- Interaction_strength_immune_cell(exp_rank = exp_rank3,pairs_immune = pairs2,show_level = show_level)
              
              plot_data <- score_result[[1]]
              detail <- score_result[[2]]
              
              if(is.null(plot_data)){
                stopMessage <- "no results."
              }else{
                ####====================================================output table
                rownames(plot_data) <- NULL
                plot_data[,"Interaction Strength"] <- round(as.numeric(plot_data[,"Interaction Strength"]),4)
                
                plot_data_df <- as.data.frame(plot_data)
                output <- list(data = plot_data_df)
                jsoncars <- toJSON(output, pretty=TRUE)
                cat(jsoncars, file = paste(outpath,'TCGA_interaction_strength_immcell_re.txt',sep='/'), fill = FALSE, labels = NULL, append = FALSE)
                ####=======================================================violin
                row_name <- paste(plot_data[,"Gene1Name"],plot_data[,"Gene2Name"],sep='_')
                plot_data2 <- cbind(plot_data,row_name)
                plot_data2 <- as.matrix(plot_data2)
                
                plot_violin(plot_matrix = plot_data2,cancer_strength = detail,outpath = outpath,level_name="immune cell")
                ####=======================================================heatmap
                cancer_immune <- aggregate(plot_data[,"Interaction Strength"],list(plot_data[,"Cancer"],plot_data[,"Immune cell"]),function(x) mean(as.numeric(x)))
                colnames(cancer_immune) <- c("Cancer","row_name","Value")
                cancer_immune <- as.matrix(cancer_immune)
                cancer_immune[,"Value"] <- round(as.numeric(cancer_immune[,"Value"]),4)
                
                heatmap_res <- plot_all_heatmap(plot_score = cancer_immune,choose_cancer = FALSE)   
                write.table(heatmap_res,paste(outpath,"/TCGA_interaction_strength_immcell_heatmap.txt",sep=''),sep='\t',row.names = F,col.names = F,quote = F)
                write.table(plot_data,paste(outpath,"/TCGA_interaction_strength_immcell_heatmap_download.txt",sep=''),sep='\t',row.names = F,quote = F)
              }
            }
          }
        }
      }
    }
    if(is.null(stopMessage)){
      
      if(sum(show_level %in% "gene")==1){
        download_table="TCGA_interaction_strength_gene_re.txt"
        heat_download_table="TCGA_interaction_strength_gene_heatmap_download.txt"
        Res_heat_all_file="TCGA_interaction_strength_gene_heatmap_all.txt"
        Res_heat_file="TCGA_interaction_strength_gene_heatmap.txt"
        resut_merge=paste0("{",
                           '"download_table" :','"',download_table,'",',
                           '"heat_download_table" :','"',heat_download_table,'",',
                           '"Res_heat_all_file" :','"',Res_heat_all_file,'",',
                           '"Res_heat_file" :','"',Res_heat_file,'",',
                           '"error_attention" :','"no',
                           '"}')
        
        success <- "TCGA_interaction_strength_gene+success"
        
      }else if(sum(show_level %in% "patient")==1){

        Res_boxplot_file="TCGA_interaction_strength_patient_boxplot.txt"
        resut_merge=paste0("{",
                           '"Res_boxplot_file" :','"',Res_boxplot_file,'",',
                           '"error_attention" :','"no',
                           '"}')
        success <- "TCGA_interaction_strength_patient+success"
        
      }else if(length(which(c("APC","B cells","DC cells","leukocyte","lymphocytes","macrophage","mast cells",
                              "MDSCs","myeloid cells","neutrophils","NK cells","T cells") %in% show_level))!=0){
        
        download_table="TCGA_interaction_strength_immcell_re.txt"
        heat_download_table="TCGA_interaction_strength_immcell_heatmap_download.txt"
        Res_heat_file="TCGA_interaction_strength_immcell_heatmap.txt"
        resut_merge=paste0("{",
                           '"download_table" :','"',download_table,'",',
                           '"heat_download_table" :','"',heat_download_table,'",',
                           '"Res_heat_file" :','"',Res_heat_file,'",',
                           '"error_attention" :','"no',
                           '"}')
        success <- "TCGA_interaction_strength_imm+success"
        
      }
      write.table(resut_merge,paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- success
      write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
      
      result_message <- paste0("{",'"error_attention" :','"success"',"}")
      return(result_message)
    }else{
      write.table(paste0("{",'"error_attention" :','"',stopMessage,'"}'),paste(outpath,"job_re.txt",sep='/'),sep='\t',quote = F,row.names = F,col.names = F)
      
      Task_ID <- as.matrix(read.table(paste0(referencePath,"/Job_ID.txt"),header = F,sep = "\t",check.names = F,fill = T))
      Task_ID_new <- Task_ID
      Task_ID_new[which(Task_ID_new[,1]==ID_new),2] <- "dead"
      write.table(Task_ID_new,paste0(referencePath,"/Job_ID.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
      
      result_message <- paste0("{",'"error_attention" :','"',stopMessage,'"}')
      return(result_message)
    }
  },error=function(e){
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
######################patient
Interaction_strength_patient <- function(exp_rank,compare_exp_rank = NULL,pairs_immune,cancer = NULL){
  pairs_immune$pairs2 <- paste(pairs_immune$Gene1Name,pairs_immune$Gene2Name,sep='_')
  
  gene1 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune$Gene1Id]
  gene2 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune$Gene2Id]
  
  index1 <- which(pairs_immune$Gene1Id %in% gene1)
  index2 <- which(pairs_immune$Gene2Id %in% gene2)
  index <- intersect(index1,index2)
  
  if(length(index)==0){
    return(list(plot_data = NULL,detail = NULL))
  }
  
  pairs <- unique(pairs_immune[index,c(1:4),drop=F])
  nlen <- nrow(pairs)
  
  score_result <- lapply(1:nlen,function(x){
    zhici1 <- as.numeric(exp_rank[pairs[x,,drop=F]$Gene1Id,])
    zhici2 <- as.numeric(exp_rank[pairs[x,,drop=F]$Gene2Id,])
    score <- (zhici1+zhici2) * (10-abs(zhici2-zhici1))
    return(score)
  })
  score_result <- do.call(rbind,score_result)
  colnames(score_result) <- colnames(exp_rank)
  rownames(score_result) <- paste(pairs[,,drop=F]$Gene1Name,pairs[,,drop=F]$Gene2Name,sep='_')
  
  ####
  sample_level <- colMeans(score_result)
  
  if(!is.null(compare_exp_rank)==T){
    gene1 <- rownames(compare_exp_rank)[rownames(compare_exp_rank) %in% pairs_immune$Gene1Id]
    gene2 <- rownames(compare_exp_rank)[rownames(compare_exp_rank) %in% pairs_immune$Gene2Id]
    
    index1 <- which(pairs_immune$Gene1Id %in% gene1)
    index2 <- which(pairs_immune$Gene2Id %in% gene2)
    index <- intersect(index1,index2)
    
    if(length(index)!=0){
      normal_pairs <- unique(pairs_immune[index,c(1:4),drop=F])
      normal_nlen <- nrow(normal_pairs)
      
      normal_score_result <- lapply(1:normal_nlen,function(x){
        zhici1 <- as.numeric(compare_exp_rank[normal_pairs[x,,drop=F]$Gene1Id,])
        zhici2 <- as.numeric(compare_exp_rank[normal_pairs[x,,drop=F]$Gene2Id,])
        score <- (zhici1+zhici2) * (10-abs(zhici2-zhici1))
        return(score)
      })
      normal_score_result <- do.call(rbind,normal_score_result)
      colnames(normal_score_result) <- colnames(compare_exp_rank)
      rownames(normal_score_result) <- paste(normal_pairs[,,drop=F]$Gene1Name,normal_pairs[,,drop=F]$Gene2Name,sep='_')
      normal_sample_level <- colMeans(normal_score_result)
      #############boxplot
    }else{
      normal_sample_level <- 0
    }
  }else{
    normal_sample_level <- 0
  }
  if(!is.null(cancer)==T){
    plot_data <- list(cancer = cancer,cancer_box = sample_level,normal_box = normal_sample_level)
  }else{
    plot_data <- list(cancer = "User cancer",cancer_box = sample_level,normal_box = normal_sample_level)
  }
  return(list(plot_data = plot_data,detail = score_result))
}
######################gene
Interaction_strength_gene <- function(exp_rank,pairs_immune,cancer = NULL,r_score = NULL,outpath){
  pairs_immune$pairs2 <- paste(pairs_immune$Gene1Name,pairs_immune$Gene2Name,sep='_')
  
  gene1 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune$Gene1Id]
  gene2 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune$Gene2Id]
  
  index1 <- which(pairs_immune$Gene1Id %in% gene1)
  index2 <- which(pairs_immune$Gene2Id %in% gene2)
  index <- intersect(index1,index2)
  
  if(length(index)==0){
    return(list(plot_data = NULL,detail = NULL))
  }
  
  pairs <- unique(pairs_immune[index,c(1:4),drop=F])
  nlen <- nrow(pairs)
  
  score_result <- lapply(1:nlen,function(x){
    zhici1 <- as.numeric(exp_rank[pairs[x,,drop=F]$Gene1Id,])
    zhici2 <- as.numeric(exp_rank[pairs[x,,drop=F]$Gene2Id,])
    score <- (zhici1+zhici2) * (10-abs(zhici2-zhici1))
    return(score)
  })
  score_result <- do.call(rbind,score_result)
  colnames(score_result) <- colnames(exp_rank)
  rownames(score_result) <- paste(pairs[,,drop=F]$Gene1Name,pairs[,,drop=F]$Gene2Name,sep='_')
  
  ####gene interaction level
  gene_pairs_level <-  rowMeans(score_result)
  # ####random 10000
  if(!is.null(cancer)==T){
    r_score <- as.numeric(as.matrix(r_score))
  }else{
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
    
  }
  save(r_score,file = paste0(outpath,"/random_strength.RData"))
  #########p value
  p_value <- lapply(1:length(gene_pairs_level),function(x){
    length(which(r_score > gene_pairs_level[x]))/10000
  })
  p_value <- unlist(p_value)
  gene_pairs_pvalue <- cbind(gene_pairs_level,p_value)
  
  if(!is.null(cancer)==T){
    gene_level_result <- cbind(pairs[,"Gene1Name"],pairs[,"Gene2Name"],cancer,gene_pairs_pvalue)
    colnames(gene_level_result) <- c("Gene1Name","Gene2Name","Cancer","Interaction Strength","P value")
  }else{
    gene_level_result <- cbind(pairs[,"Gene1Name"],pairs[,"Gene2Name"],"-",gene_pairs_pvalue)
    colnames(gene_level_result) <- c("Gene1Name","Gene2Name","Cancer","Interaction Strength","P value")
  }
  
  return(list(plot_data = gene_level_result,detail = score_result))
}
######################immune cell
Interaction_strength_immune_cell <- function(exp_rank,pairs_immune,cancer = NULL,show_level = show_level){
  pairs_immune$pairs2 <- paste(pairs_immune$Gene1Name,pairs_immune$Gene2Name,sep='_')
  
  gene1 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune$Gene1Id]
  gene2 <- rownames(exp_rank)[rownames(exp_rank) %in% pairs_immune$Gene2Id]
  
  index1 <- which(pairs_immune$Gene1Id %in% gene1)
  index2 <- which(pairs_immune$Gene2Id %in% gene2)
  index <- intersect(index1,index2)
  
  if(length(index)==0){
    return(list(plot_data = NULL,detail = NULL))
  }
  
  pairs <- unique(pairs_immune[index,c(1:4),drop=F])
  nlen <- nrow(pairs)
  
  score_result <- lapply(1:nlen,function(x){
    zhici1 <- as.numeric(exp_rank[pairs[x,,drop=F]$Gene1Id,])
    zhici2 <- as.numeric(exp_rank[pairs[x,,drop=F]$Gene2Id,])
    score <- (zhici1+zhici2) * (10-abs(zhici2-zhici1))
    return(score)
  })
  score_result <- do.call(rbind,score_result)
  colnames(score_result) <- colnames(exp_rank)
  rownames(score_result) <- paste(pairs[,,drop=F]$Gene1Name,pairs[,,drop=F]$Gene2Name,sep='_')
  
  immune_cells <- show_level

  score_result2 <- cbind(rownames(score_result),score_result)
  colnames(score_result2)[1] <- "pairs"
  score_result2_immune <- merge(pairs_immune[,c("pairs2","Immucell_type")],score_result2,by.x = "pairs2",by.y = "pairs")
  
  score_result2_immune_choose <- score_result2_immune[score_result2_immune[,"Immucell_type"] %in% immune_cells,,drop=F]
  
  if(nrow(score_result2_immune_choose)==0){
    return(list(plot_data = NULL,detail = NULL))
  }
  score_result2_immune_choose2 <- matrix(as.numeric(as.matrix(score_result2_immune_choose[,-c(1,2)])),nrow=nrow(score_result2_immune_choose))
  immune_level <- rowMeans(score_result2_immune_choose2)
  
  score_result2_immune_level <- cbind(score_result2_immune_choose[,c(1,2)],immune_level)
  g12 <- do.call(rbind,lapply(score_result2_immune_level[,"pairs2"],function(x) unlist(strsplit(x,"_"))))
  
  if(!is.null(cancer)==T){
    immune_cell_level_result <- cbind(g12,cancer,score_result2_immune_level[,c(2,3)])
  }else{
    immune_cell_level_result <- cbind(g12,"-",score_result2_immune_level[,c(2,3)])
  }
  colnames(immune_cell_level_result) <- c("Gene1Name","Gene2Name","Cancer","Immune cell","Interaction Strength")
  
  return(list(plot_data = immune_cell_level_result,detail = score_result))
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
#######################################
ID_transfer <- function(gene_list,ID_transfer_Path=NULL){
  
  l2 <- load(paste(ID_transfer_Path,"ID_transfer_file.RData",sep='/'))
  ID_transfer_data <- eval(parse(text = l2))
  
  gene <- gene_list[1]
  
  if(grepl("ENSG",gene[1])){
    gene_names = "ensembl"
  }else if(!grepl(paste("[",paste(toupper(letters),collapse ='|'),"]",sep=''),gene[1])){
    gene_names = "entrez"
  }else{
    gene_names = "symbol"
  }
  
  if(gene_names == "ensembl"){
    res<- unique(ID_transfer_data[ID_transfer_data[,2] %in% gene_list,c(2,3,2)])
  }else if(gene_names == "entrez"){
    res <- unique(ID_transfer_data[ID_transfer_data[,1] %in% gene_list,c(1,3,2)])
  }else{
    res <- unique(ID_transfer_data[ID_transfer_data[,3] %in% gene_list,c(3,3,2)])
  }
  return(res)
}
########################plot heatmap
plot_all_heatmap <- function(plot_score,choose_cancer = TRUE){
  usp <- unique(plot_score[,"row_name"])
  uspy <- cbind(usp,c(0:(length(usp)-1)))
  
  if(choose_cancer == TRUE){
    score_cancers <- unique(plot_score[,"Cancer"])
    scx <- cbind(score_cancers,c(0:(length(score_cancers)-1)))
  }else{
    score_cancers <- "User cancer"
    scx <- cbind("User cancer",c(0:(length(score_cancers)-1)))
  }
  
  plot_score2 <- cbind(plot_score,x=0,y=0)
  ####
  for(i in 1:nrow(uspy)){
    plot_score2[which(plot_score2[,"row_name"] %in% uspy[i,1]),"y"] <- uspy[i,2]
  }
  
  for(i in 1:nrow(scx)){
    plot_score2[which(plot_score2[,"Cancer"] %in% scx[i,1]),"x"] <- scx[i,2]
  }
  ####[y,x,value]
  xyAxis0 <- paste("[",plot_score2[,"y"],",",plot_score2[,"x"],",",plot_score2[,"Value"],"]",sep='')
  xyAxis1 <- paste(xyAxis0,collapse =',')
  xyAxis <- paste("[",xyAxis1,"]",sep = "")
  y0<- paste(usp,collapse="','")
  y <- paste("['",y0,"']",sep='')
  x0<- paste(score_cancers,collapse="','")
  x <- paste("['",x0,"']",sep='')
  
  value_min <- min(as.numeric(plot_score2[,"Value"]))
  value_max <- max(as.numeric(plot_score2[,"Value"]))
  
  return(t(c(x,y,xyAxis,value_min,value_max)))
}
###########boxplot
plot_box_data <- function(all = plot_data,choose_cancer = TRUE){

  if(choose_cancer == TRUE){
    all_cancer <- names(all)
    
    all_test <- lapply(all_cancer,function(x){
      if(length(all[[x]][[2]])==1|length(all[[x]][[3]]) == 1){
        p <- 1
      }else if(length(unique(all[[x]][[2]]))==1&length(unique(all[[x]][[3]]))==1){
        p <- 1
      }else{
        p <- t.test(as.numeric(all[[x]][[2]]),as.numeric(all[[x]][[3]]),alternative = "two.sided")$p.value
      }
      
      cancer_box0 <- paste(round(as.numeric(all[[x]][[2]]),4),collapse=',')
      cancer_box <- paste("[",cancer_box0,"]",sep = '')
      
      normal_box0 <- paste(round(as.numeric(all[[x]][[3]]),4),collapse=',')
      normal_box <- paste("[",normal_box0,"]",sep = '')
    
      return(list(p_value = p,cancer_box = cancer_box,normal_box = normal_box))
    })
    
    p_value <- lapply(all_test,function(x) x$p_value)
    p_value <- do.call(rbind,p_value)
    
    xx <- paste(all_cancer,"(p=",round(p_value,4),")",sep='')
    all_cancer0 <- paste(xx,collapse = "','")
    all_cancer_s <- paste("'",all_cancer0,"'",sep='')
   
    cancer_box <- lapply(all_test,function(x) x$cancer_box)
    cancer_box <- do.call(rbind,cancer_box)
    cancer_box_s0 <- paste(cancer_box,collapse=',')
    cancer_box_s <- paste("[",cancer_box_s0,"]",sep='')
   
    normal_box <- lapply(all_test,function(x) x$normal_box)
    normal_box <- do.call(rbind,normal_box)
    normal_box_s0 <- paste(normal_box,collapse=',')
    normal_box_s <- paste("[",normal_box_s0,"]",sep='')
    return(t(c(all_cancer_s,cancer_box_s,normal_box_s)))
  }else{
    
    if(length(all[[2]])==1|length(all[[3]]) == 1){
      p <- 1
    }else if(length(unique(all[[2]]))==1&length(unique(all[[3]]))==1){
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
}
###################violin
plot_violin <- function(plot_matrix,cancer = NULL,cancer_strength,outpath,level_name){
  
  if(!is.null(cancer)==T){
    lapply(1:nrow(plot_matrix),function(x){
      
      pair_name <- plot_matrix[x,"row_name"]
      pair_exp <- as.numeric(cancer_strength[pair_name,])
      
      data <- data.frame(gene = unname(pair_name),exp = pair_exp)
      
      file_name <- paste(cancer,plot_matrix[x,1],plot_matrix[x,2],sep="_")
      if(level_name == "immune cell"){
        file_name2 <- paste(outpath,"/interaction_strength_imm_",file_name,".png",sep='')
      }else{
        file_name2 <- paste(outpath,"/interaction_strength_",file_name,".png",sep='')
      }
      
      
      #png(file = file_name2,units= "in",width = 5,height = 5,res = 300)
      
      plot_title <- paste("The interaction strength of ",pair_name," in ",cancer,sep='')
      
       p <- ggplot(data,aes(x=gene,y=exp))+geom_violin(fill = "pink2")+
        theme(axis.title.x = element_blank())+
        labs(y="Interaction strength",title = plot_title)+
        theme(plot.title = element_text(hjust = 0.5))  
       
       #print(p)

       ggsave(file_name2,p)

    })
  }else{
    lapply(1:nrow(plot_matrix),function(x){
      
      pair_name <- plot_matrix[x,"row_name"]
      pair_exp <- as.numeric(cancer_strength[pair_name,])
      
      data <- data.frame(gene = unname(pair_name),exp = pair_exp)
      
      file_name <- paste("-",plot_matrix[x,1],plot_matrix[x,2],sep="_")
      
      if(level_name == "immune cell"){
        file_name2 <- paste(outpath,"/interaction_strength_imm_",file_name,".png",sep='')
      }else{
        file_name2 <- paste(outpath,"/interaction_strength_",file_name,".png",sep='')
      }
      
      #png(file = file_name2,units= "in",width = 5,height = 5,res = 300)
      
      plot_title <- paste("The interaction strength of ",pair_name,sep='')
      
      p <- ggplot(data,aes(x=gene,y=exp))+geom_violin(fill = "pink2")+
        theme(axis.title.x = element_blank())+
        labs(y="Interaction strength",title = plot_title)+
        theme(plot.title = element_text(hjust = 0.5))  
      
      #print(p)

      ggsave(file_name2,p)

    })
  }
}
