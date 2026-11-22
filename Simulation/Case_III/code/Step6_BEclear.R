##############################################################
#################     Simulation Case III    #################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("BEclear", quietly = TRUE))
  BiocManager::install("BEclear")

library(BEclear)


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


## Set the same seed as BUSbeta in simulation case III
seed = 84
set.seed(seed)


###################### batch correction #########################

Y_combined <- t(rbind(Y[[1]], Y[[2]], Y[[3]]))
rownames(Y_combined) = paste0("j", 1:nrow(Y_combined))
colnames(Y_combined) = paste0("i", 1:ncol(Y_combined))
Y_BEclear_samples = data.frame(sample_id = colnames(Y_combined), batch_id = c(rep(1, n_vec[1]), rep(2, n_vec[2]), rep(3, n_vec[3])))

Y_BEclear_corrected = correctBatchEffect(Y_combined, Y_BEclear_samples, adjusted=T, method="fdr") $ correctedPredictedData

Y_BEclear <- list()
Y_BEclear[[1]] = t(Y_BEclear_corrected[, 1:n_vec[1]])
Y_BEclear[[2]] = t(Y_BEclear_corrected[, (n_vec[1]+1):(n_vec[1]+n_vec[2])])
Y_BEclear[[3]] = t(Y_BEclear_corrected[, (n_vec[1]+n_vec[2]+1):(n_vec[1]+n_vec[2]+n_vec[3])])


## Save Corrected data
save(Y_BEclear, file = "../result_data/Y_BEclear.RData")

