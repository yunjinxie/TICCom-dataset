#=======================ICGC run(重跑)=============================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM/30_log_05")

#cs <- c("CLLE-ES","LIRI-JP","MALY-DE","PACA-AU","PACA-CA","PBCA-US","RECA-EU")

cs <- "LICA_LIRI"

for(i in 1:length(cs)){
  d <- "I:/F/tumor_immune_interaction/computation/A_webserver/ICGC/count/FPKM/30_log_05/reference"
  IDPath <- expPath <- paste0(getwd(),"/reference/",cs[i],"_tumour_30_log_05.RData")
  compare_expPath <- gene_listPath <- cancer<- "empty"
  referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference/reference"
  
  cc <- cs[i]
  outpath <- paste0(getwd(),"/result/",cc)
  dir.create(outpath)
  ###
  source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/TCGA_interaction_strength.R")
  start <- Sys.time()
  TCGA_interaction_strength(IDPath = IDPath,expPath = expPath,compare_expPath = compare_expPath,gene_listPath = gene_listPath,
                            cancer = cancer,referencePath = referencePath,show_level="gene",outpath = outpath)
  end <- Sys.time()
  end-start
  ####
  source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICCom.R")
  
  start <- Sys.time()
  ICCom(expPath = expPath,compare_expPath =compare_expPath,p_value = 0.1,referencePath = referencePath,
        databaseFile = "union_LR_database",gene_names = "symbol",organism = "human",
        outpath =outpath)
  end <- Sys.time()
  end-start
  
}

####TCGA####
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/sjy_TCGA/TCGA_analysis/result/interaction_strength/33cancer_gene")

cs <- "GBMLGG"

d <- "I:/F/tumor_immune_interaction/computation/A_webserver/Rfunction/function_relay"
IDPath <- expPath <- paste0(d,"/",cs,"_cancer_genes_expression_30_05_log2.RData")
compare_expPath <- gene_listPath <- cancer<- "empty"
referencePath <- "I:/F/tumor_immune_interaction/computation/A_webserver/database/reference/reference"

cc <- cs
outpath <- paste0(getwd(),"/GBMLGG")
dir.create(outpath)
###
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/TCGA_interaction_strength.R")
start <- Sys.time()
TCGA_interaction_strength(IDPath = IDPath,expPath = expPath,compare_expPath = compare_expPath,gene_listPath = gene_listPath,
                          cancer = cancer,referencePath = referencePath,show_level="gene",outpath = outpath)
end <- Sys.time()
end-start
####
source("I:/F/tumor_immune_interaction/computation/A_webserver/RfunctionNew/TICCom_new2/run/ICCom.R")

start <- Sys.time()
ICCom(expPath = expPath,compare_expPath =compare_expPath,p_value = 0.1,referencePath = referencePath,
      databaseFile = "union_LR_database",gene_names = "symbol",organism = "human",
      outpath =outpath)
end <- Sys.time()
end-start




