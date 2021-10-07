#
# Author: Edward B. Irvine
#
# Description:
# Makes PCA score and loading plot for delta-sigH NHP vaccination cohort antibody profiling data. 
#
# Created: 5/21/20
#
# Modified: 10/6/21
#

###################
###### Housekeeping -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###################

# Load required libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(factoextra)
library(ggforce)

# Read in data
deepak_data <- read.csv("SigH Titers.csv")










#######################
###### Pre-process data -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######################

deepak_data$Color <- rep("#66CCFE", nrow(deepak_data))
deepak_data$Color[deepak_data$Group == "SigH"] <- "#FC9A13"

# Drop variables
deepak_data <- deepak_data %>%
  dplyr::select(!matches("EBOV"))

# Make Wk7 subset
wk7_keep <- c("Group", "Color", "Wk7")
wk7_data <- deepak_data %>%
  dplyr::select(matches(wk7_keep))

# Get rid of missing data
wk7_data <- na.omit(wk7_data)

# Prep aesthetics for plots
colnames(wk7_data) <- str_replace(colnames(wk7_data), "Wk7.", "")
wk7_data$Color <- as.factor(wk7_data$Color)










##########
###### PCA -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##########

# Run PCA
wk7_data_pca <- prcomp(as.matrix(wk7_data[3:(ncol(wk7_data))]))

# PCA score plot
tiff("sigH_PCA_score.tiff", units = "in", width = 7, height = 7, res = 300)
fviz_pca_ind(wk7_data_pca,
  label = "none",
  habillage = wk7_data$Group,
  repel = TRUE,
  geom = "point",
  pointsize = 8,
  alpha = 0.8,
  invisible = "quali",
  pch = 19) + 
  xlim(-4.35, 4.35) +
  ylim(-4.35, 4.35) +
  theme_linedraw() +
  scale_color_manual(values = c("#FFADD6", "#940003")) +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        aspect.ratio = 7/7) 
dev.off()
  
# PCA loading plot
tiff("sigH_PCA_loading.tiff", units = "in", width = 6.97, height = 6.58, res = 300)
fviz_pca_var(wk7_data_pca,
  geom = c("text", "arrow"),
  repel = TRUE,
  col.var="contrib",
  labelsize = 5) +
  scale_color_gradient2(name = "Contribution", low = "white", mid = "lightgrey", high = "black") +
  theme_linedraw() +
  xlim(-1, 1) +
  ylim(-1, 1) +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        aspect.ratio = 7/7) 
dev.off()

  
  

