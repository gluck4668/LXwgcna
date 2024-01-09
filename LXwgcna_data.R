
library(openxlsx)

meta_data_example <- read.xlsx("meta_data_example.xlsx")
experiment_info_example <- read.xlsx("experiment_info_example.xlsx")

usethis::use_data(meta_data_example,overwrite = T)
usethis::use_data(experiment_info_example,overwrite = T)

rm(list=ls())

data(meta_data_example)
data(experiment_info_example)

