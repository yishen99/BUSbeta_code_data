##############################################################
################    Kabuki syndrome data    ##################
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


## Read data Y_preprocessed
Y_preprocessed <- read.table("../input_data/Y_preprocessed.txt", sep = "\t", header = T)
Y_preprocessed <- as.data.frame(Y_preprocessed)

B <- length(unique(Y_preprocessed $ batch))
J <- ncol(Y_preprocessed) - 2
K <- 2 ### From BIC analysis
n_vec <- as.vector(table(Y_preprocessed $ batch))


## Set seed for BUSbeta in application to the Kabuki syndrome dataset
seed = 11
set.seed(seed)


###################### batch correction #########################

Y_BEclear_data <- t(Y_preprocessed[,c(-J-1, -J-2)])
Y_BEclear_samples = data.frame(sample_id = rownames(Y_preprocessed), batch_id = Y_preprocessed[, J+1])

Y_BEclear_corrected = correctBatchEffect(Y_BEclear_data, Y_BEclear_samples, adjusted=TRUE, method="fdr") $ correctedPredictedData

Y_BEclear <- list()
for(b in 1:B){
  ind = which(Y_preprocessed$batch==b)
  Y_BEclear[[b]] = as.matrix(t(Y_BEclear_corrected[, ind]))
}

## Save Corrected data
save(Y_BEclear, file = "../result_data/Y_BEclear.RData")



