\name{LXwgcna}
\alias{LXwgcna}
\title{LXwgcna::LXwgcna}
\usage{
LXwgcna (data_file,expriment_ifo,key_phenotype)
}
\description{
LXwgcna can be used to assess the phenotype-related hub metabolites. The codes can be seen at
"https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/"
}
\examples{

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

data_file <- "meta_data_example.xlsx"
expriment_ifo <- "experiment_info_example.xlsx"
key_phenotype <- "col_expression"

LXwgcna (data_file,expriment_ifo,key_phenotype)


}
