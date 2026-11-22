##############################################################
################    Kabuki syndrome data    ##################
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



## Set seed for BUSbeta in application to the Kabuki syndrome dataset
seed = 11
set.seed(seed)


###################### batch correction #########################

Y_combined = t(rbind(Y[[1]], Y[[2]]))
Y_batch = c(rep(1, nrow(Y[[1]])), rep(2, nrow(Y[[2]]))) 

Y_ComBat_corrected = ComBat(dat = Y_combined, batch = Y_batch)

Y_ComBat = list()
Y_ComBat[[1]] = t(Y_ComBat_corrected[, 1:n_vec[1]])
Y_ComBat[[2]] = t(Y_ComBat_corrected[, (n_vec[1]+1):(n_vec[1]+n_vec[2])])


## Save Corrected data
save(Y_ComBat, file = "../result_data/Y_ComBat.RData")



