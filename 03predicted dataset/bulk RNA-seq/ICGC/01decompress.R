#=======================================ICGC==========================
setwd("I:/F/tumor_immune_interaction/computation/A_webserver/ICGC")

d <- dir(getwd())
d2 <- d[grepl("exp",d)]

array <- d2[grepl("array",d2)]
seq <- d2[grepl("seq",d2)]

array_cancer <- stringr::str_match(array,"exp_array\\.(.*)-.*")[,2]
seq_cancer <- stringr::str_match(seq,"exp_seq\\.(.*)-.*")[,2]

library(R.utils)
for(i in 1:length(d2)){
  gunzip(d2[i])
}
