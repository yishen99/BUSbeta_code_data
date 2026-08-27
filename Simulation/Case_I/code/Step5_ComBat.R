##############################################################
#################     Simulation Case I    ###################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("sva", quietly = TRUE))
  BiocManager::install("sva")

library(sva)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y
load("../input_data/Y.RData")
B <- length(Y)
J <- ncol(Y[[1]])
K <- 5 ### From BIC analysis
n_vec <- sapply(Y, nrow)


## Set the same seed as BUSbeta in simulation case I
seed = 2002
set.seed(seed)


###################### batch correction #########################

Y_combined = t(rbind(Y[[1]], Y[[2]], Y[[3]]))
Y_batch = c(rep(1, n_vec[1]), rep(2, n_vec[2]), rep(3, n_vec[3])) 

Y_ComBat_corrected = ComBat(dat = Y_combined, batch = Y_batch)

Y_ComBat = list()
Y_ComBat[[1]] = t(Y_ComBat_corrected[, 1:n_vec[1]])
Y_ComBat[[2]] = t(Y_ComBat_corrected[, (n_vec[1]+1):(n_vec[1]+n_vec[2])])
Y_ComBat[[3]] = t(Y_ComBat_corrected[, (n_vec[1]+n_vec[2]+1):(n_vec[1]+n_vec[2]+n_vec[3])])


## Save Corrected data
save(Y_ComBat, file = "../result_data/Y_ComBat.RData")



