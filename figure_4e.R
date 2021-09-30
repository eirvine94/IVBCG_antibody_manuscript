#
# Author: Edward B. Irvine
#
# Description: 
# LASSO PLS-DA analysis from humoral profiling from the BCG route study.
# Includes LASSO feature selection, PLS-DA classification, and an assessment of model significance.  
#
# Modified: 11/5/20
#
# Last updated: 9/30/21
#

###################
###### Housekeeping -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###################

# Load required libraries
library(plyr)
library(dplyr)
library(glmnet)
library(mixOmics)
library(tibble) 

# Read in data
fold_data <- read.csv("fold_z_data.csv")










#######################
###### Pre-process data -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######################

# Set the data frame row names to the associated monkey ID (stored in column 1)
rownames(fold_data) <- fold_data$ID
fold_data["ID"] <- list(NULL)

# Make binary outcome variable based on lung Mtb burden. "Protected" = 0; "Not protected" = 1
fold_data$outcome <- rep(0, nrow(fold_data))
fold_data$outcome[fold_data$Mtb.burden > 1000] <- 1
fold_data$outcome <- as.factor(fold_data$outcome)










#######################
###### Define functions -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######################

# LASSO feature selection 
lasso <- function(dataset, nrep) {
  # Description:
  # Implements lasso by leveraging the glmnet package. 
  # To mitigate overfitting, bootstrap samples are iteratively taken from each group to fit the lasso model. 
  # Bootstrapping the lasso estimator of the regression parameters allows the generation of an inclusion probability for each antibody feature.
  # Defined as the percentage of times, or the probability the lasso regression coefficient corresponding to a particular antibody feature is non-zero.
  #
  # Input:
  # dataset -- data frame in which the first column represents the class identity
  # nrep -- number of bootstrap datasets generated
  #
  # Returns: 
  # Named numeric vector containing the inclusion probability for each feature
  #
  
  # Initialize variables
  colnames(dataset)[1] <- "Group"
  lambda_grid <- 10^seq(10, -2, length = 100)
  count <- 1
  coeff_frame <- data.frame()
  coeff_frame_final <- data.frame(matrix(NA, nrow = ncol(dataset), ncol = 1))
  
  for (i in 1:nrep) {
    
    ### Create bootstrap dataset
    set.seed(i)
    boot_sample <- data.frame()
    sample_index <- sample(rownames(dataset), replace = TRUE)
    sample_index <- unlist(sample_index, use.names = FALSE)
    sample_index <- sort(sample_index)
    sample_index
    for (ind in sample_index) {
      boot_sample <- rbind(boot_sample, dataset[ind, ])
    }
    
    ### Reformat data for LASSO model
    x_sample <- model.matrix(boot_sample$Group ~ ., boot_sample)[,-1]
    y_sample <- boot_sample$Group
    
    # Fit a model, tuning over a grid of lambda values for cross validation
    cv.fit <- cv.glmnet(x_sample, y_sample,  alpha = 1, nfolds = nrow(dataset), family = "binomial", lambda = lambda_grid, grouped = FALSE)
    
    # Identify optimal lambda value selected by CV
    best_lam <- cv.fit$lambda.min
    
    # Fit new LASSO model using the optimal value for lambda
    final.mod <- glmnet(x_sample, y_sample, alpha = 1, family = "binomial", lambda = best_lam)
    
    # Extract coefficients from best model and bind to growing data frame
    coeff_frame <- data.frame(coef(final.mod)[, 1])
    coeff_frame_final <- cbind(coeff_frame_final, coeff_frame)
    
  }
  
  # Clean final data frame
  coeff_frame_final[1] <- list(NULL)
  coeff_frame_final <- data.frame(t(coeff_frame_final))
  coeff_frame_final[1] <- list(NULL)
  
  # Find variable inclusion probabilities
  inclusion_probs <- (colSums(coeff_frame_final != 0) / nrep)
  
  return(inclusion_probs)
}

# PLS-DA learning algorithm
tune_plsda <- function(dataset, feature_prob, cutoffs, nrep) {
  # Description:
  # Generate and compute the cross-validation accuracy for different PLS-DA models, each of which uses a distinct variable inclusion probability threshold. 
  # Cross-validation accuracy is balanced by group size.
  #
  # Inputs:
  # dataset -- data frame in which the first column represents the class identity 
  # feature_prob -- named numeric vector containing the inclusion probability for each feature
  # cutoffs -- numeric vector containing the different thresholds for variable inclusion probability you would like to test
  # nrep -- number of times to repeat cross-validation
  #
  # Returns: 
  # A data frame with the cross-validation accuracy and optimal number of components for each variable inclusion probability threshold.
  #
  
  # Initialize variables
  colnames(dataset)[1] <- "Group"
  accuracy_frame <- data.frame(matrix(nrow = length(cutoffs), ncol = 3))
  colnames(accuracy_frame) <- c("cutoff", "accuracy", "ncomp")
  accuracy_frame$cutoff <- cutoffs
  count <- 1

  for (cutoff in cutoffs) {
    
    # Create subset of data containing only the selected features
    features <- c(feature_prob >= cutoff)
    
    if (sum(features) <= 1) {
      break
    }
    
    feature_data <- dataset[ , colnames(dataset)[which(features) + 1]] 
    feature_data <- data.frame(cbind(dataset$Group, feature_data))
    colnames(feature_data)[1] <- "Group"
    
    # Prep data for PLS-DA
    x_plsda <- as.matrix(feature_data[ , 2:ncol(feature_data)]) 
    y_plsda <- as.factor(feature_data$Group) 
    
    ### Fit PLS-DA models
    # Select number of models to fit
    ncomp <- min(nrow(feature_data), ncol(feature_data)) - 1
    if (ncomp == 1) {
      ncomp <- 2
    } 
    if (ncomp > 20) {
      ncomp <- 20
    }
    
    # Fit models
    plsda.fit <- mixOmics::plsda(x_plsda, y_plsda, ncomp = ncomp, scale = FALSE)
    
    # Perform 5-fold repeated CV on each of the models (with different numbers of components)
    perf.plsda <- perf(plsda.fit, dist = "centroids.dist", validation = "Mfold", folds = 5, nrepeat = nrep, auc = FALSE, progressBar = TRUE)
    
    # Identify model with optimal number of components. If 1 component is optimal, convert to the 2 component model.
    best_size <- which.min(perf.plsda$error.rate$BER)
    
    if (best_size == 1) {
      best_size <- 2
    } 
    
    accuracy_frame$ncomp[count] <- best_size
    
    # Identify balanced error rate of best model. This is balanced for differences in group size. 
    model_error <- perf.plsda$error.rate$BER[paste("comp", best_size, sep = ""), "centroids.dist"]
    
    # Identify (balanced) accuracy of best model
    model_accuracy <- 1 - model_error
    accuracy_frame$accuracy[count] <- model_accuracy
    
    count <- count + 1
    
  }
  
  return(accuracy_frame)
}

# Generate permuted model
permute_model <- function(dataset, x_matrix, ncomp, nrep) {
  # Description:
  # Iteratively generate permuted PLS-DA models, where the class labels are randomly suffled during each iteration. 
  #
  # Inputs:
  # dataset -- data frame in which the first column represents the class identity 
  # x_matrix -- a matrix of your dataset of interest, excluding the column containing the class labels (this should be the same as the first parameter in the "mixOmics::plsda" function)
  # ncomp -- number of components used to generate your real model
  # nrep -- number of iterations
  #
  # Returns: 
  # Data frame comprised of the model accuracy across all the iterations.
  #
  
  colnames(dataset)[1] <- "Group"
  
  # Initialize vector to store balanced error rate from each permuted model (balanced by class size)
  permute_error_vec <- c()
  
  # Each time through the loop, the error rate for a permuted model is generated
  for (i in 1:nrep) {
    
    # Randomly shuffle the group labels
    set.seed(i)
    y_permute <- as.factor(sample(dataset$Group))
    
    # Generate PLS-DA models with with ncomp components
    plsda.permute <- mixOmics::plsda(x_matrix, y_permute, ncomp = ncomp, scale = FALSE)
    
    # Evaluate performance of models, and identify the optimal number of components in the model via 5-fold CV. If 1 component is optimal, convert to the 2 component model.
    perf.plsda.permute <- perf(plsda.permute, dist = "centroids.dist", validation = "Mfold", folds = 5, auc = FALSE, progressBar = TRUE)
    
    # Identify error rate of optimal model, and store in vector
    permute_error_vec[i] <- perf.plsda.permute$error.rate$BER[paste("comp", ncomp, sep = ""), "centroids.dist"]
  }
  
  # Calculate and return a data frame with the accuracy of the permuted model
  permute_accuracy_vec <- 1 - permute_error_vec
  permute_accuracy_vec <- data.frame(permute_accuracy_vec)
  colnames(permute_accuracy_vec)[1] <- "accuracy"
  
  return(data.frame(permute_accuracy_vec))
} 










############################
##### Implement LASSO PLS-DA -------------------------------------------------------------------------------------------------------------------------------------------------------------------
############################

nrep <- 100

# Make all subset
all_keep <- c("Plasma", "BAL", "outcome")
all_fold <- fold_data %>% 
  dplyr::select(matches(all_keep)) %>%
  dplyr::select(outcome, everything())



###
### LASSO
###
inclusion_probs_all <- lasso(all_fold, nrep)
inclusion_probs_all 



###
### PLS-DA
###

cutoffs_all_lasso <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)
accuracy_frame_all_lasso <- tune_plsda(all_fold, inclusion_probs_all, cutoffs_all_lasso, nrep)
accuracy_frame_all_lasso 

best_cutoff_lasso <- 0.45
best_ncomp_lasso <- 3

# Create subset of data containing only the selected features
features_all_lasso <- c(inclusion_probs_all >= best_cutoff_lasso)
feature_data_all_lasso <- all_fold[ , colnames(all_fold)[which(features_all_lasso) + 1]]
feature_data_all_lasso <- cbind(all_fold$outcome, feature_data_all_lasso)
colnames(feature_data_all_lasso)[1] <- "Group"

# Prep data for PLS-DA
x_plsda_all_lasso <- as.matrix(feature_data_all_lasso[ , 2:length(feature_data_all_lasso)])
y_plsda_all_lasso <- as.factor(feature_data_all_lasso$Group) 

# Fit PLS-DA
plsda.fit.final_all_lasso <- mixOmics::plsda(x_plsda_all_lasso, y_plsda_all_lasso, ncomp = best_ncomp_lasso, scale = FALSE)

# Plot PLS-DA
scores <- data.frame(cbind(plsda.fit.final_all_lasso$Y, plsda.fit.final_all_lasso$variates$X))
scores$V1 <- as.factor(scores$V1)
tiff("plsda_scores_plot.tiff", units = "in", width = 7, height = 7, res = 300)
scores_plot <- scores %>%
  ggplot(aes(x = comp1, y = comp2, color = V1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 6, alpha = 0.8) +
  theme_linedraw() +
  xlab("Scores on LV1") +
  ylab("Scores on LV2") +
  scale_color_manual("Group", labels = c("Protected", "Susceptible"), values = c("black", "lightgrey")) +
  stat_ellipse(level = 0.95) +
  theme(plot.title = element_text(size = 25, face = "bold"), 
        axis.text = element_text(size = 20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        axis.title = element_text(size = 24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=20),
        aspect.ratio = 7/7)
scores_plot
dev.off()



###
### VIP analysis
###
vip_all_fold_lasso <- data.frame(vip(plsda.fit.final_all_lasso))
vip_all_fold_LV1_lasso <- vip_all_fold_lasso %>% 
  rownames_to_column("feature") %>%
  dplyr::select(feature, comp1) %>%
  arrange(desc(comp1)) %>%
  column_to_rownames("feature")
vip_all_fold_LV1_lasso










###############################
##### Assess model significance -----------------------------------------------------------------------------------------------------------------------------------------------------------------
###############################

# Generate permuted model accuracy distribution
lasso_permuted_all <- permute_model(all_fold, x_plsda_all_lasso, ncomp = 3, nrep = 100)
mean(lasso_permuted_all$accuracy)

# Pull out CV accuracy from real model
set.seed(10)
perf.plsda.lasso.all <- perf(plsda.fit.final_all_lasso, dist = "centroids.dist", validation = "Mfold", folds = 5, nrepeat = 100, auc = TRUE, progressBar = TRUE)
lasso_all_error <- t(data.frame(perf.plsda.lasso.all$error.rate.all$BER))
lasso_all_error <- data.frame(lasso_all_error[ , "comp3"])
lasso_all_accuracy <- 1 - lasso_all_error
colnames(lasso_all_accuracy)[1] <- "accuracy"
mean(lasso_all_accuracy$accuracy)

# Perform Mann-Whitney U test to assess model significance, comparing your real model to the permuted model
lasso_all_mann_whit <- wilcox.test(lasso_all_accuracy$accuracy, lasso_permuted_all$accuracy)
lasso_all_mann_whit




