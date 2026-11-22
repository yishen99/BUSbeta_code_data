##############################################################
#################     Simulation Case III    #################
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


## BUS clustering label
Y_trs = list(t(Y[[1]]), t(Y[[2]]), t(Y[[3]]))
BUSfits = BUSgibbs(Y_trs, n.subtypes = 5, n.iterations = 1000, hyperparameters = c(2, 0.1, 0.5, 2, 2, 1, 2, 1, 1, 1, 1), showIteration = FALSE)
Z_BUS = BUSfits $ Subtypes


## Compute ARI
ARI_BUS <- data.frame(Batch1 = round(ARI(Z_real[[1]], Z_BUS[[1]]), 3), 
                      Batch2 = round(ARI(Z_real[[2]], Z_BUS[[2]]), 3), 
                      Batch3 = round(ARI(Z_real[[3]], Z_BUS[[3]]), 3),
                      Overall = round(ARI(unlist(Z_real), unlist(Z_BUS)), 3))



## Save ARI
write.csv(ARI_BUS, row.names = "BUS ARI", file = "../result_data/ARI_BUS.csv")


###################### batch correction #########################

alpha_BUS = BUSfits $ alpha ### J
mu_BUS = BUSfits $ mu ### J * K
gamma_BUS = BUSfits $ gamma ### J * B
sigma_BUS = BUSfits $ sigma_sq ### J * B

Y_BUS <- list()
Y_BUS[[1]] = Y[[1]]
for(b in 2:B){
  Y_BUS[[b]] <- matrix(NA, nrow = n_vec[b], ncol = J)
  for(i in 1:n_vec[b]){
    Y_BUS[[b]][i, ] = alpha_BUS + mu_BUS[, Z_BUS[[b]][i]] + (Y[[b]][i, ]-alpha_BUS-mu_BUS[, Z_BUS[[b]][i]]-gamma_BUS[, b]) / (sigma_BUS[, b]/sigma_BUS[, 1])
  }
}



## Save Corrected data and clustering label
save(Y_BUS, file = "../result_data/Y_BUS.RData")
save(Z_BUS, file = "../result_data/Z_BUS.RData")



