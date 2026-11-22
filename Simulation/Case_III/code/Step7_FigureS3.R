##############################################################
#################     Simulation Case III    #################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
library(ggplot2)

font.size = 20


################# Visualization function #####################

visualize_by_subgroup_no_batch <- function(Data, feature_ind,
                                           title_name,
                                           color_key_range = seq(0, 1, 0.1), 
                                           is.Y = NULL,
                                           subgroup_list = NULL,
                                           batch_colors = NULL) {
  
  K <- length(Data)
  n_vec <- sapply(Data, nrow)
  total_n <- sum(n_vec)
  
  df <- do.call(cbind, lapply(Data, t))
  df <- df[feature_ind, ]
  
  if(is.Y==TRUE) subgroup_labels <- rep(paste0("A", 1:K), times = n_vec)
  else subgroup_labels <- rep(paste0("C", 1:K), times = n_vec)
  
  if (!is.null(subgroup_list)) {
    batch_vec <- unlist(subgroup_list)
    if (length(batch_vec) != total_n) stop("Subgroups don't match sample number.")
  } else {
    batch_vec <- rep("NA", total_n)
  }
  
  anno_df <- data.frame(
    Batch = batch_vec,
    Subgroup = subgroup_labels
  )
  
  if(is.Y==TRUE) subgroup_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0", "#F1937B", "#8DB0FE")[1:K], paste0("A", 1:K))
  else subgroup_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0", "#F1937B", "#8DB0FE")[1:K], paste0("C", 1:K))
  if (is.null(batch_colors)) {
    batch_levels <- unique(batch_vec)
    n_batch <- length(batch_levels)
    gray_shades <- gray.colors(n_batch, start = 0.75, end = 0.2)
    batch_colors <- setNames(gray_shades, batch_levels)
  }
  
  ha <- HeatmapAnnotation(
    Batch = anno_df$Batch,
    Subgroup = anno_df$Subgroup,
    col = list(
      Batch = batch_colors,
      Subgroup = subgroup_colors
    ),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      Subgroup = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      ),
      Batch = list(
        legend = FALSE
      )
    )
  )
  
  ht <- Heatmap(df,
                name = "Value",
                col = colorRamp2(color_key_range, colorRampPalette(c("#3B4CC0", "white", "#B40426"))(length(color_key_range))),
                top_annotation = ha,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                column_title = title_name,
                heatmap_legend_param = list(
                  title = "Value",
                  title_gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = font.size),
                  labels_gp = gpar(fontsize = 16),
                  legend_height = unit(4, "cm")
                )
  )
  draw(ht, heatmap_legend_side = "left")
  decorate_heatmap_body("Value", {
    grid.rect(gp = gpar(col = "black", lwd = 1, fill = NA))
  })
}



visualize_by_batch <- function(Data, feature_ind,
                               title_name,
                               color_key_range = seq(0, 1, 0.1),
                               is.Y = NULL,
                               subgroup_list = NULL,
                               subgroup_colors = NULL) {
  
  B <- length(Data)
  K <- length(unique(unlist(subgroup_list)))
  n_vec <- sapply(Data, nrow)
  total_n <- sum(n_vec)
  
  df <- do.call(cbind, lapply(Data, t))
  df <- df[feature_ind, ]
  
  batch_labels <- rep(1:B, times = n_vec)
  if(is.Y==TRUE) subgroup_labels <- paste0("A", unlist(subgroup_list))
  else subgroup_labels <- paste0("C", unlist(subgroup_list))
  
  anno_df <- data.frame(
    Batch = batch_labels,
    Subgroup = subgroup_labels
  )
  
  if (is.null(subgroup_colors)) {
    subgroup_levels <- unique(unlist(subgroup_list))
    n_subgroup <- length(subgroup_levels)
    if(is.Y==TRUE) subgroup_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0", "#F1937B", "#8DB0FE")[1:n_subgroup], paste0("A", 1:n_subgroup))
    else subgroup_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0", "#F1937B", "#8DB0FE")[1:n_subgroup], paste0("C", 1:n_subgroup))
  }
  gray_shades <- gray.colors(B, start = 0.75, end = 0.2)
  batch_colors <- setNames(gray_shades, 1:B)
  
  ha <- HeatmapAnnotation(
    df = anno_df,
    col = list(Batch = batch_colors, Subgroup = subgroup_colors),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      Subgroup = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      ),
      Batch = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      )
    )
  )
  
  ht <- Heatmap(df,
                name = "Value",
                col = colorRamp2(color_key_range, colorRampPalette(c("#3B4CC0", "white", "#B40426"))(length(color_key_range))),
                top_annotation = ha,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                column_title = title_name,
                heatmap_legend_param = list(
                  title = "Value",
                  title_gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = font.size),
                  labels_gp = gpar(fontsize = 16),
                  legend_height = unit(4, "cm")
                )
  )
  draw(ht, heatmap_legend_side = "left")
  decorate_heatmap_body("Value", {
    grid.rect(gp = gpar(col = "black", lwd = 1, fill = NA))
  })
}



visualize_by_subgroup <- function(Data, feature_ind,
                                  title_name,
                                  color_key_range = seq(0, 1, 0.1), 
                                  is.Y = NULL,
                                  batch_list = NULL,
                                  batch_colors = NULL) {
  
  K <- length(Data)
  n_vec <- sapply(Data, nrow)
  total_n <- sum(n_vec)
  
  df <- do.call(cbind, lapply(Data, t))
  df <- df[feature_ind, ]
  
  if(is.Y==TRUE) subgroup_labels <- rep(paste0("A", 1:K), times = n_vec)
  else subgroup_labels <- rep(paste0("C", 1:K), times = n_vec)
  
  if (!is.null(batch_list)) {
    batch_vec <- unlist(batch_list)
    if (length(batch_vec) != total_n) stop("Subgroups don't match sample number.")
  } else {
    batch_vec <- rep("NA", total_n)
  }
  
  anno_df <- data.frame(
    Batch = batch_vec,
    Subgroup = subgroup_labels
  )
  
  if(is.Y==TRUE) subgroup_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0", "#F1937B", "#8DB0FE")[1:K], paste0("A", 1:K))
  else subgroup_colors <- setNames(c("#E5E5E5", "#B40426", "#3B4CC0", "#F1937B", "#8DB0FE")[1:K], paste0("C", 1:K))
  if (is.null(batch_colors)) {
    batch_levels <- sort(unique(batch_vec))
    n_batch <- length(batch_levels)
    gray_shades <- gray.colors(n_batch, start = 0.75, end = 0.2)
    batch_colors <- setNames(gray_shades, batch_levels)
  }
  
  ha <- HeatmapAnnotation(
    Batch = anno_df$Batch,
    Subgroup = anno_df$Subgroup,
    col = list(
      Batch = batch_colors,
      Subgroup = subgroup_colors
    ),
    annotation_name_gp = gpar(fontsize = 0),
    annotation_legend_param = list(
      Subgroup = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      ),
      Batch = list(
        title_gp = gpar(fontsize = font.size, fontface = "bold"),
        labels_gp = gpar(fontsize = font.size),
        legend_direction = "vertical",
        grid_height = unit(0.8, "cm"),
        grid_width = unit(0.5, "cm"),     
        gap = unit(3, "mm")               
      )
    )
  )
  
  ht <- Heatmap(df,
                name = "Value",
                col = colorRamp2(color_key_range, colorRampPalette(c("#3B4CC0", "white", "#B40426"))(length(color_key_range))),
                top_annotation = ha,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_column_names = FALSE,
                show_row_names = FALSE,
                column_title = title_name,
                heatmap_legend_param = list(
                  title = "Value",
                  title_gp = gpar(fontfamily = "sans", fontface = "bold", fontsize = font.size),
                  labels_gp = gpar(fontsize = 16),
                  legend_height = unit(4, "cm")
                )
  )
  draw(ht, heatmap_legend_side = "left")
  decorate_heatmap_body("Value", {
    grid.rect(gp = gpar(col = "black", lwd = 1, fill = NA))
  })
}



batch_label_by_subgroup  <- function(Z_list) {
  batch_label = list()
  for(k in 1:K){
    temp = c()
    for(b in 1:B){
      if(length(which(Z_list[[b]]==k))>0) temp = c(temp, rep(b, times = length(which(Z_list[[b]]==k))))
    }
    batch_label[[k]] = temp
  }
  return(batch_label)
}


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

## Read True subgroup mean
load("../input_data/truemean_by_subgroup.RData")

## Read data Y_BUSbeta
load("../result_data/Y_BUSbeta.RData")
load("../result_data/Z_BUSbeta.RData")
# Reorder the clustering labels (without affecting the ARI result)
map <- c(1, 3, 5, 4, 2)
for(b in 1:B){
  Z_BUSbeta[[b]] <- map[Z_BUSbeta[[b]]]
}

## Read data Y_BUS
load("../result_data/Y_BUS.RData")
load("../result_data/Z_BUS.RData")
# Reorder the clustering labels (without affecting the ARI result)
map <- c(2, 3, 1, 4, 5)
for(b in 1:B){
  Z_BUS[[b]] <- map[Z_BUS[[b]]]
}

## Read data Y_ComBat
load("../result_data/Y_ComBat.RData")

## Read data Y_BEclear
load("../result_data/Y_BEclear.RData")



######## Raw data and corrected data organized by subgroups #########

Y_by_subgroup <- list()
for(k in 1:K){
  temp = Y[[1]][which(Z_real[[1]]==k), ]
  if(B>=2){
    for(b in 2:B) temp = rbind(temp, Y[[b]][which(Z_real[[b]]==k), ])
  }
  Y_by_subgroup[[k]] = temp
}

Yhat_by_subgroup <- list()
for(k in 1:K){
  temp = Y_BUSbeta[[1]][which(Z_BUSbeta[[1]]==k), ]
  if(B>=2){
    for(b in 2:B) temp = rbind(temp, Y_BUSbeta[[b]][which(Z_BUSbeta[[b]]==k), ])
  }
  Yhat_by_subgroup[[k]] = temp
}

YBUS_by_subgroup <- list()
for(k in 1:K){
  temp = Y_BUS[[1]][which(Z_BUS[[1]]==k), ]
  if(B>=2){
    for(b in 2:B) temp = rbind(temp, Y_BUS[[b]][which(Z_BUS[[b]]==k), ])
  }
  YBUS_by_subgroup[[k]] = temp
}

YComBat_by_subgroup <- list()
for(k in 1:K){
  temp = Y_ComBat[[1]][which(Z_real[[1]]==k), ]
  if(B>=2){
    for(b in 2:B) temp = rbind(temp, Y_ComBat[[b]][which(Z_real[[b]]==k), ])
  }
  YComBat_by_subgroup[[k]] = temp
}

YBEclear_by_subgroup <- list()
for(k in 1:K){
  temp = Y_BEclear[[1]][which(Z_real[[1]]==k), ]
  if(B>=2){
    for(b in 2:B) temp = rbind(temp, Y_BEclear[[b]][which(Z_real[[b]]==k), ])
  }
  YBEclear_by_subgroup[[k]] = temp
}


################# Visualization #####################


## True subgroup mean
png(filename="../figures/FigureS3(a).png", width = 500, height = 300)
visualize_by_subgroup_no_batch(Data = truemean_by_subgroup, feature_ind = 1:J, title_name = "True subgroup means", is.Y = TRUE)
dev.off()


### Y by subgroup
png(filename="../figures/FigureS3(b).png", width = 500, height = 300)
visualize_by_subgroup(Data = Y_by_subgroup, feature_ind = 1:J, title_name = "Raw simulated data", is.Y = TRUE, batch_list = batch_label_by_subgroup(Z_real))
dev.off()


### Y by batch
png(filename="../figures/FigureS3(c).png", width = 500, height = 300)
visualize_by_batch(Data = Y, feature_ind = 1:J, title_name = "Raw simulated data", is.Y = TRUE, subgroup_list = Z_real)
dev.off()


### Y_BUSbeta
png(filename="../figures/FigureS3(e).png", width = 500, height = 300)
visualize_by_subgroup(Data = Yhat_by_subgroup, feature_ind = 1:J, title_name = "BUSbeta", is.Y = FALSE, batch_list = batch_label_by_subgroup(Z_BUSbeta))
dev.off()


### Y_BUS
png(filename="../figures/FigureS3(f).png", width = 500, height = 300)
visualize_by_subgroup(Data = YBUS_by_subgroup, feature_ind = 1:J, title_name = "BUS", is.Y = FALSE, batch_list = batch_label_by_subgroup(Z_BUS))
dev.off()


### Y_ComBat
png(filename="../figures/FigureS3(g).png", width = 500, height = 300)
visualize_by_subgroup(Data = YComBat_by_subgroup, feature_ind = 1:J, title_name = "ComBat", is.Y = TRUE, batch_list = batch_label_by_subgroup(Z_real))
dev.off()


### Y_BEclear
png(filename="../figures/FigureS3(h).png", width = 500, height = 300)
visualize_by_subgroup(Data = YBEclear_by_subgroup, feature_ind = 1:J, title_name = "BEclear", is.Y = TRUE, batch_list = batch_label_by_subgroup(Z_real))
dev.off()



