##############################################################
#################     Simulation Case I    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

if (!require("mclust", quietly = TRUE))
  install.packages("mclust")

library(mclust)
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


## Set the same seed as BUSbeta in simulation case I
seed = 2002
set.seed(seed)


## GMM clustering label
Z_gmm <- list()
for(b in 1:B){
  Z_gmm[[b]] <- Mclust(Y[[b]]) $ classification
}


## Compute ARI
ARI_GMM <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_gmm[[1]]), 3), 
                      Batch2 = round(ARI(Z_real[[2]], Z_gmm[[2]]), 3), 
                      Batch3 = round(ARI(Z_real[[3]], Z_gmm[[3]]), 3),
                      Overall = "-")


## Save ARI
write.csv(ARI_GMM, row.names = "GMM ARI", file = "../result_data/ARI_GMM.csv")

