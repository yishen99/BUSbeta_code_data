##############################################################
################    Kabuki syndrome data    ##################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

########################## Note ##############################
## Please first install required R packages.
## R: System_preparation.R

## Please download the two data sources from GEO with accession
## numbers GSE218186 and GSE97362 to the "input_data" folder. 
## GSE218186_series_matrix.txt:
## https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218186/matrix/GSE218186_series_matrix.txt.gz
## GSE218186_RAW.txt:
## https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE218186&format=file
## GSE97362_series_matrix.txt:
## https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97362/matrix/GSE97362_series_matrix.txt.gz
##############################################################

if (!require("vroom", quietly = TRUE))
  install.packages("vroom")

library(vroom)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


######################## Batch 1 ########################

Y1 <- vroom("../input_data/GSE218186_RAW.txt", delim = "\t", skip = 1, col_names = FALSE)
Y1 = as.data.frame(Y1)[, c(1, 1:62 * 2)]

Y1_first <- vroom("../input_data/GSE218186_series_matrix.txt", delim = "\t", skip = 36, n_max = 33, col_names = FALSE)
Y1_first = as.data.frame(Y1_first)

## Raw data in batch 1
Y1_proc = Y1[, 2:ncol(Y1)]
rownames(Y1_proc) = Y1[, 1]
colnames(Y1_proc) = Y1_first[2, -1]

## True subgroup in batch 1
Y1_subject = c(NA)
Y1_subject[which(Y1_first[11, -1]=="disease state: Kabuki syndrome")] = 1 # 7 KS
Y1_subject[which(Y1_first[11, -1]=="disease state: Control")] = 0 # 55 normal

######################## Batch 2 ########################


Y2 <- vroom("../input_data/GSE97362_series_matrix.txt", delim = "\t", skip = 66, col_names = FALSE)
Y2 = as.data.frame(Y2)
Y2 = Y2[-nrow(Y2), ]
cpg_batch2_NULL = which(is.na(Y2), arr.ind = TRUE)
Y2 = Y2[-cpg_batch2_NULL[, 1], ]

Y2_first <- vroom("../input_data/GSE97362_series_matrix.txt", delim = "\t", skip = 30, n_max = 35, col_names = FALSE)
Y2_first = as.data.frame(Y2_first)

## Raw data in batch 2
Y2_proc = Y2[2:nrow(Y2), 2:ncol(Y2)]
rownames(Y2_proc) = Y2[2:nrow(Y2), 1]
colnames(Y2_proc) = Y2[1, -1]
KABUKI_ind = c(which(Y2_first[12, -1]=="disease state: Kabuki"),94,97,98,101,102,104,105,108)
control_ind = which(Y2_first[12, -1]=="disease state: Control")
Y2_proc = Y2_proc[, c(KABUKI_ind, control_ind)]

## True subgroup in batch 2
Y2_subject = c(NA)
Y2_subject[1:length(KABUKI_ind)] = rep(1, times=length(KABUKI_ind)) # 19 KS
Y2_subject[(length(KABUKI_ind)+1):(length(KABUKI_ind)+length(control_ind))] = rep(0, times=length(KABUKI_ind)+length(control_ind)) # 125 normal





## Select informative CpG sites
cpg_Kabuki <- read.table("../input_data/informative_cpg_Kabuki.txt", header = T)[, 3] # 148 CpG sites
cpg_batch1_ind = na.omit(match(cpg_Kabuki, rownames(Y1_proc)))
cpg_batch2_ind = na.omit(match(cpg_Kabuki, rownames(Y2_proc)))

Y1_proc = Y1_proc[cpg_batch1_ind, ]
Y2_proc = Y2_proc[cpg_batch2_ind, ]
cpg_batch1_0 = which(apply(Y1_proc, MARGIN=1, min)==0)
cpg_batch2_0 = which(apply(Y2_proc, MARGIN=1, min)==0)
if(length(cpg_batch1_0)>0) Y1_proc = Y1_proc[-cpg_batch1_0, ] 
if(length(cpg_batch2_0)>0) Y2_proc = Y2_proc[-cpg_batch2_0, ] 


## Preprocessed data: Y_preprocessed
Y1_proc = t(Y1_proc)
Y2_proc = t(Y2_proc)
Y_preprocessed = rbind(Y1_proc, Y2_proc)
Y_preprocessed = cbind(Y_preprocessed, batch = c(rep(1, times = nrow(Y1_proc)), rep(2, times = nrow(Y2_proc))))
Y_preprocessed = cbind(Y_preprocessed, subject = c(Y1_subject, Y2_subject)) 


## Save the prepreocessed data
write.table(Y_preprocessed, file = "../input_data/Y_preprocessed.txt", sep = "\t", col.names = T, row.names = T)

## Save true subgroups Z_real
Z_real = list(Y1_subject, Y2_subject)
save(Z_real, file = "../input_data/Z_real.RData")




