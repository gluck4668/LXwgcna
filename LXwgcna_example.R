
if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXwgcna")

library(LXwgcna)

??LXwgcan

#-----------------------
data(meta_data_example)
data(experiment_info_example)

#-----------------------
rm(list=ls())

data_file <- "didymin-model-old orgin.xlsx"
expriment_info <- "data-test1.xlsx"
key_phenotype <- "TG"

LXwgcna (data_file,expriment_info,key_phenotype)
