###################################################################################################################
#Running tests
###################################################################################################################
#automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doParallel",
  "dplyr",
  "tidyr",
  "tidyverse",
  "data.table",
  "glmnet",
  "stringr",
  "MASS", 
  "stats",
  "bigsnpr",
  "RcppArmadillo",
  "Rcpp",
  "Matrix"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
      )
    )
}

source(".../TransferTWAS_func.R")
sourceCpp(file='.../TransferTWAS_cpp.cpp')

TransferTWAS(fold = 5, 
          n_tune = 5, 
          n_cores=48,
          gene_name="ENSG00000169174.10",
          dir_geno=".../ENSG00000169174.10.rds",
          dir_exp=".../ENSG00000169174.10",
          dir_output=".../output",
          tar_tis_name="Liver")
