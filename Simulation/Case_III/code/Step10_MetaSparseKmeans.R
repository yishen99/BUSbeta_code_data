##############################################################
#################     Simulation Case III   ##################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

if (!require("remotes", quietly = TRUE))
  install.packages("remotes")
if (!require("MetaSparseKmeans", quietly = TRUE))
  remotes::install_github("Caleb-Huo/MetaSparseKmeans")

library(MetaSparseKmeans)
library(aricode)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y
load("../input_data/Y.RData")
load("../input_data/Z_real.RData")
B <- length(Y)
J <- ncol(Y[[1]])
K <- 5 ### From BIC analysis
n_vec <- sapply(Y, nrow)


## Set the same seed as BUSbeta in simulation case III
seed = 84
set.seed(seed)


## MetaSparseKmeans clustering label
Z_MSkms = MetaSparseKmeans(Y, K = 5, wbounds=10, lambda=0.5) $ Cs


## Compute ARI
ARI_MSkms <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_MSkms[[1]]), 3), 
                        Batch2 = round(ARI(Z_real[[2]], Z_MSkms[[2]]), 3), 
                        Batch3 = round(ARI(Z_real[[3]], Z_MSkms[[3]]), 3),
                        Overall = round(ARI(unlist(Z_real), unlist(Z_MSkms)), 3))


## Save ARI
write.csv(ARI_MSkms, row.names = "MetaSparseKmeans ARI", file = "../result_data/ARI_MetaSparseKmeans.csv")


