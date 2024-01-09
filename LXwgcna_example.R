
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

data_file <- "test_data.xlsx"
expriment_info <- "test_trait.xlsx"

key_phenotype <- "BIS"

abnormal_sample_height = NA #如果有异常值，可以输入一个合适的heitht以去掉；如果没有异常值，填NA。

LXwgcna (data_file,expriment_info,key_phenotype)
