##############################################################
################    Down syndrome data    ####################
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


## Read data Y_preprocessed
Y_preprocessed <- read.table("../input_data/Y_preprocessed.txt", sep = "\t", header = T)
Y_preprocessed <- as.data.frame(Y_preprocessed)

B <- length(unique(Y_preprocessed $ batch))
J <- ncol(Y_preprocessed) - 2
K <- 2 ### From BIC analysis
n_vec <- as.vector(table(Y_preprocessed $ batch))

## Data Y list
Y <- list()
for(b in 1:B){
  ind = which(Y_preprocessed $ batch == b)
  Y[[b]] = as.matrix(Y_preprocessed[ind, c(-J-1, -J-2)])
}


## Read true subgroups Z_real
load("../input_data/Z_real.RData")


## Set seed for BUSbeta in application to the Down syndrome dataset
seed = 8
set.seed(seed)


## GMM clustering label
Z_gmm <- list()
for(b in 1:B){
  Z_gmm[[b]] <- Mclust(Y[[b]]) $ classification
}


## Compute ARI
ARI_GMM <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_gmm[[1]]), 3), 
                      Batch2 = round(ARI(Z_real[[2]], Z_gmm[[2]]), 3), 
                      Overall = "-")


## Save ARI
write.csv(ARI_GMM, row.names = "GMM ARI", file = "../result_data/ARI_GMM.csv")

