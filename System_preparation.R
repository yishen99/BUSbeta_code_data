##############################################################
################## System preparation ########################
##############################################################
rm(list=ls())

### Install basic packages, if necessary
if (!require("aricode", quietly = TRUE)) 
  install.packages("aricode", version="1.0.3")
if (!require("circlize", quietly = TRUE))
  install.packages("circlize", version="0.4.16")
if (!require("coda", quietly = TRUE))
  install.packages("coda", version="0.19.4.1")
if (!require("ggbreak", quietly = TRUE)) 
  install.packages("ggbreak", version="0.1.4")
if (!require("ggplot2", quietly = TRUE)) 
  install.packages("ggplot2", version="3.5.1")
if (!require("invgamma", quietly = TRUE)) 
  install.packages("invgamma", version="1.1")
if (!require("matrixStats", quietly = TRUE)) 
  install.packages("matrixStats", version="1.5.0")
if (!require("remotes", quietly = TRUE))
  install.packages("remotes", version="2.5.0")
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr", version="1.3.1")
if (!require("truncnorm", quietly = TRUE)) 
  install.packages("truncnorm", version="1.0.9")
if (!require("umap", quietly = TRUE))
  install.packages("umap", version="0.2.10.0")
if (!require("vroom", quietly = TRUE)) 
  install.packages("vroom", version="1.6.5")


### install the Bioconductor package manager, if necessary
if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")


### BUS
if (!require("BUScorrect", quietly = TRUE))
  BiocManager::install("BUScorrect")

### ComBat
if (!require("sva", quietly = TRUE))
  BiocManager::install("sva")

### BEclear
if (!require("BEclear", quietly = TRUE))
  BiocManager::install("BEclear")

### GMM
if (!require("mclust", quietly = TRUE)) 
  install.packages("mclust", version="6.1.1")

### MetaSparseKmeans
if (!require("MetaSparseKmeans", quietly = TRUE))
  remotes::install_github("Caleb-Huo/MetaSparseKmeans")


### Attention!
# For Windows users, please note that the version of Rtools 
# needs to be compatible with the version of R!

