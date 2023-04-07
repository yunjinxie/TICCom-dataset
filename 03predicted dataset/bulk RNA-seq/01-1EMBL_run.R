#=======================EMBL run=============================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/EMBL_expression_atlas/30_log_05")

cs <- c("fibrosarcoma","nasopharyngeal carcinoma","neuroblastoma","non-small-cell lung cancer",
                  "ovarian serous adenocarcinoma","pancreatic ductal carcinoma(PDAC)",
                  "papillary thyroid carcinoma","Small Cell Lung Cancer","uterine leiomyosarcoma")

for(i in 2:length(cs)){

  IDPath <- expPath <- paste0(getwd(),"/",cs[i],"_cancer_FPKM_30_log_05.RData")
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
