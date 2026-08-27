##############################################################
################    Kabuki syndrome data    ##################
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


## Set seed for BUSbeta in application to the Kabuki syndrome dataset
seed = 11
set.seed(seed)


## Kmeans clustering label
Z_kms <- list()
for(b in 1:B){
  Z_kms[[b]] <- kmeans(Y[[b]], centers = K) $ cluster
}

## Compute ARI
ARI_Kmeans <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_kms[[1]]), 3), 
                         Batch2 = round(ARI(Z_real[[2]], Z_kms[[2]]), 3), 
                         Overall = "-")


## Save ARI
write.csv(ARI_Kmeans, row.names = "Kmeans ARI", file = "../result_data/ARI_Kmeans.csv")




