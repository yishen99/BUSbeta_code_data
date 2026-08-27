##############################################################
#################     Simulation Case II    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

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


## Set the same seed as BUSbeta in simulation case II
seed = 64
set.seed(seed)


## Kmeans clustering label
Z_kms <- list()
for(b in 1:B){
  Z_kms[[b]] <- kmeans(Y[[b]], centers = K) $ cluster
}

## Compute ARI
ARI_Kmeans <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_kms[[1]]), 3), 
                         Batch2 = round(ARI(Z_real[[2]], Z_kms[[2]]), 3), 
                         Batch3 = round(ARI(Z_real[[3]], Z_kms[[3]]), 3),
                         Overall = "-")


## Save ARI
write.csv(ARI_Kmeans, row.names = "Kmeans ARI", file = "../result_data/ARI_Kmeans.csv")




