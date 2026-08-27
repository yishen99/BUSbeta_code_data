##############################################################
################    Kabuki syndrome data    ##################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("BUScorrect", quietly = TRUE))
  BiocManager::install("BUScorrect")

library(BUScorrect)
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


## BUS clustering label
Y_trs = list(t(Y[[1]]), t(Y[[2]]))
BUSfits = BUSgibbs(Y_trs, n.subtypes = K, n.iterations = 1000, hyperparameters = c(2, 0.1, 0.5, 2, 2, 1, 2, 1, 1, 1, 1), showIteration = FALSE)
Z_BUS = BUSfits $ Subtypes


## Compute ARI
ARI_BUS <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_BUS[[1]]), 3), 
                      Batch2 = round(ARI(Z_real[[2]], Z_BUS[[2]]), 3), 
                      Overall = round(ARI(unlist(Z_real), unlist(Z_BUS)), 3))


## Save ARI
write.csv(ARI_BUS, row.names = "BUS ARI", file = "../result_data/ARI_BUS.csv")


###################### batch correction #########################

alpha_BUS = BUSfits $ alpha ### J
mu_BUS = BUSfits $ mu ### J * K
gamma_BUS = BUSfits $ gamma ### J * B
sigma_BUS = BUSfits $ sigma_sq ### J * B

Y_BUS <- list()
for(b in 1:B){
  Y_BUS[[b]] <- matrix(NA, nrow = n_vec[b], ncol = J)
  for(i in 1:n_vec[b]){
    Y_BUS[[b]][i, ] = alpha_BUS + mu_BUS[, Z_BUS[[b]][i]] + (Y[[b]][i, ]-alpha_BUS-mu_BUS[, Z_BUS[[b]][i]]-gamma_BUS[, b]) / (sigma_BUS[, b]/sigma_BUS[, 1])
  }
}


## Save Corrected data and clustering label
save(Y_BUS, file = "../result_data/Y_BUS.RData")
save(Z_BUS, file = "../result_data/Z_BUS.RData")



