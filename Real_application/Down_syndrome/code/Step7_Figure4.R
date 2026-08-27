##############################################################
################    Down syndrome data    ####################
##############################################################

########################## Note ##############################
## Please set the working directory to the source file 
## location.
##############################################################

if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")
if (!require("umap", quietly = TRUE))
  install.packages("umap")
if (!require("ggbreak", quietly = TRUE))
  install.packages("ggbreak")

library(ggplot2)
library(tidyr)
library(umap)
library(ggbreak)

font.size = 28


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read data Y_preprocessed
Y_preprocessed <- read.table("../input_data/Y_preprocessed.txt", sep = "\t", header = T)
Y_preprocessed <- as.data.frame(Y_preprocessed)

## Read data Y_BUSbeta
load("../result_data/Y_BUSbeta.RData")

## Read data Y_BUS
load("../result_data/Y_BUS.RData")

## Read data Y_ComBat
load("../result_data/Y_ComBat.RData")

## Read data Y_BEclear
load("../result_data/Y_BEclear.RData")


B <- length(unique(Y_preprocessed $ batch))
J <- ncol(Y_preprocessed) - 2
K <- 2 ### From BIC analysis
n_vec <- as.vector(table(Y_preprocessed $ batch))


## Set seed for BUSbeta in application to the Down syndrome dataset
seed = 8
set.seed(seed)


#################  Preprocessed data  ####################

umap_Y = umap(Y_preprocessed[, c(-J-1, -J-2)], min_dist = 0.5, spread = 1) $ layout
umap_Y = data.frame(umap_Y)
plot_Y = cbind(umap_Y, Y_preprocessed[, c(J+1, J+2)])


plot_Y$group <- interaction(plot_Y$batch, plot_Y$subject)
plot_Y$group_label <- factor(plot_Y$group,
                             levels = c("1.0", "2.0", "1.1", "2.1"),
                             labels = c("Subgroup A1 in batch 1",
                                        "Subgroup A1 in batch 2",
                                        "Subgroup A2 in batch 1",
                                        "Subgroup A2 in batch 2")
)

color_values <- c("#AA4466", "#DDCC77", "#44AA99", "#332288")
color_labels <- levels(plot_Y$group_label)

umap_plot_Y <- ggplot(plot_Y, aes(x = X1, y = X2, colour = group_label)) +
  geom_point(size = 3) +
  annotate("segment", x = -8, xend = 7, y = -14, yend = -14, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  annotate("segment", x = -8, xend = -8, y = -14, yend = 12, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  scale_color_manual(
    values = color_values,
    labels = color_labels,
    drop = FALSE,
    guide = guide_legend(override.aes = list(size = 7))
  ) +
  scale_y_break(c(-10.5, -1.5), scales = "fixed", space = 0) +
  scale_y_break(c(2, 7), scales = "fixed", space = 0) + 
  labs(x = NULL, y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = font.size),
    legend.key.height = unit(50, "pt") 
  )

ggsave("../figures/Figure4(a).png", umap_plot_Y, width=9, height=6)


#################  BUSbeta  ####################

YBUSbeta_proc = rbind(Y_BUSbeta[[1]], Y_BUSbeta[[2]])
YBUSbeta_proc = cbind(YBUSbeta_proc, Y_preprocessed[, c(J+1, J+2)])

umap_YBUSbeta = umap(YBUSbeta_proc[, c(-J-1, -J-2)], min_dist = 0.5, spread = 1) $ layout
umap_YBUSbeta = data.frame(umap_YBUSbeta)
plot_YBUSbeta = cbind(umap_YBUSbeta, Y_preprocessed[, c(J+1, J+2)])

plot_YBUSbeta$group <- interaction(plot_YBUSbeta$batch, plot_YBUSbeta$subject)
plot_YBUSbeta$group_label <- factor(plot_YBUSbeta$group,
                                    levels = c("1.0", "2.0", "1.1", "2.1"),
                                    labels = c("Subgroup C1 in batch 1",
                                               "Subgroup C1 in batch 2",
                                               "Subgroup C2 in batch 1",
                                               "Subgroup C2 in batch 2")
)

color_values <- c("#AA4466", "#DDCC77", "#44AA99", "#332288")
color_labels <- levels(plot_YBUSbeta$group_label)

umap_plot_YBUSbeta <- ggplot(plot_YBUSbeta, aes(x = X1, y = X2, colour = group_label)) +
  geom_point(size = 3) +
  annotate("segment", x = -11, xend = 13, y = -7, yend = -7, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  annotate("segment", x = -11, xend = -11, y = -7, yend = 9, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  scale_color_manual(
    values = color_values,
    labels = color_labels,
    drop = FALSE,
    guide = guide_legend(override.aes = list(size = 7))
  ) +
  scale_x_break(c(-6, 8), scales = "fixed", space = 0) +
  labs(x = NULL, y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = font.size),
    legend.key.height = unit(50, "pt") 
  )

ggsave("../figures/Figure4(c).png", umap_plot_YBUSbeta, width=9, height=6)


#################  BUS  ####################

YBUS_proc = rbind(Y_BUS[[1]], Y_BUS[[2]])
YBUS_proc = cbind(YBUS_proc, Y_preprocessed[, c(J+1, J+2)])

umap_YBUS = umap(YBUS_proc[, c(-J-1, -J-2)], min_dist = 0.5, spread = 1) $ layout
umap_YBUS = data.frame(umap_YBUS)
plot_YBUS = cbind(umap_YBUS, Y_preprocessed[, c(J+1, J+2)])


plot_YBUS$group <- interaction(plot_YBUS$batch, plot_YBUS$subject)
plot_YBUS$group_label <- factor(plot_YBUS$group,
                                levels = c("1.0", "2.0", "1.1", "2.1"),
                                labels = c("Subgroup C1 in batch 1",
                                           "Subgroup C1 in batch 2",
                                           "Subgroup C2 in batch 1",
                                           "Subgroup C2 in batch 2")
)

color_values <- c("#AA4466", "#DDCC77", "#44AA99", "#332288")
color_labels <- levels(plot_YBUS$group_label)

umap_plot_YBUS <- ggplot(plot_YBUS, aes(x = X1, y = X2, colour = group_label)) +
  geom_point(size = 3) +
  annotate("segment", x = -9, xend = 8, y = -9, yend = -9, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  annotate("segment", x = -9, xend = -9, y = -9, yend = 10, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  scale_color_manual(
    values = color_values,
    labels = color_labels,
    drop = FALSE, 
    guide = guide_legend(override.aes = list(size = 7))
  ) +
  scale_x_break(c(-5.5, 3.5), scales = "fixed", space = 0) + 
  labs(x = NULL, y = NULL, color = NULL) + 
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = font.size),
    legend.key.height = unit(50, "pt"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave("../figures/Figure4(d).png", umap_plot_YBUS, width=9, height=6)


#################  ComBat  ####################

YCombat_proc = rbind(Y_ComBat[[1]], Y_ComBat[[2]])
YCombat_proc = cbind(YCombat_proc, Y_preprocessed[, c(J+1, J+2)])

umap_YCombat = umap(YCombat_proc[, c(-J-1, -J-2)], min_dist = 0.5, spread = 1) $ layout
umap_YCombat = data.frame(umap_YCombat)
plot_YCombat = cbind(umap_YCombat, Y_preprocessed[, c(J+1, J+2)])

plot_YCombat$group <- interaction(plot_YCombat$batch, plot_YCombat$subject)
plot_YCombat$group_label <- factor(plot_YCombat$group,
                                   levels = c("1.0", "2.0", "1.1", "2.1"),
                                   labels = c("Subgroup A1 in batch 1",
                                              "Subgroup A1 in batch 2",
                                              "Subgroup A2 in batch 1",
                                              "Subgroup A2 in batch 2")
)

color_values <- c("#AA4466", "#DDCC77", "#44AA99", "#332288")
color_labels <- levels(plot_YCombat$group_label)

umap_plot_YCombat <- ggplot(plot_YCombat, aes(x = X1, y = X2, colour = group_label)) +
  geom_point(size = 3) +
  annotate("segment", x = -13, xend = 15, y = -2.5, yend = -2.5, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  annotate("segment", x = -13, xend = -13, y = -2.5, yend = 2.5, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  scale_color_manual(
    values = color_values,
    labels = color_labels,
    drop = FALSE,
    guide = guide_legend(override.aes = list(size = 7))
  ) +
  scale_x_break(c(-8, 10), scales = "fixed", space = 0) + 
  labs(x = NULL, y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = font.size),
    legend.key.height = unit(50, "pt")
  )

ggsave("../figures/Figure4(e).png", umap_plot_YCombat, width=9, height=6)


#################  BEclear  ####################

YBEclear_proc = rbind(Y_BEclear[[1]], Y_BEclear[[2]])
YBEclear_proc = cbind(YBEclear_proc, Y_preprocessed[, c(J+1, J+2)])

umap_YBEclear = umap(YBEclear_proc[, c(-J-1, -J-2)], min_dist = 0.5, spread = 1) $ layout
umap_YBEclear = data.frame(umap_YBEclear)
plot_YBEclear = cbind(umap_YBEclear, Y_preprocessed[, c(J+1, J+2)])

plot_YBEclear$group <- interaction(plot_YBEclear$batch, plot_YBEclear$subject)
plot_YBEclear$group_label <- factor(plot_YBEclear$group,
                                    levels = c("1.0", "2.0", "1.1", "2.1"),
                                    labels = c("Subgroup A1 in batch 1",
                                               "Subgroup A1 in batch 2",
                                               "Subgroup A2 in batch 1",
                                               "Subgroup A2 in batch 2")
)

color_values <- c("#AA4466", "#DDCC77", "#44AA99", "#332288")
color_labels <- levels(plot_YBEclear$group_label)

umap_plot_YBEclear <- ggplot(plot_YBEclear, aes(x = X1, y = X2, colour = group_label)) +
  geom_point(size = 3) +
  annotate("segment", x = -20, xend = 23, y = -7, yend = -7, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  annotate("segment", x = -20, xend = -20, y = -7, yend = 8, 
           color = "white", 
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  scale_color_manual(
    values = color_values,
    labels = color_labels,
    drop = FALSE,
    guide = guide_legend(override.aes = list(size = 7))
  ) +
  scale_x_break(c(-13.5, 17.5), scales = "fixed", space = 0) +
  labs(x = NULL, y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = font.size),
    legend.key.height = unit(50, "pt") 
  )



ggsave("../figures/Figure4(f).png", umap_plot_YBEclear, width=9, height=6)




