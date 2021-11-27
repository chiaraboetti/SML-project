library(ggplot2)
library(scatterplot3d)
library(dplyr)
library(knitr)
library(smotefamily)
#####################################################

blood.train = read.csv("dataset/blood_training.csv")
blood.test = read.csv("dataset/blood_test.csv")

# To color differently in the plot
index = which(blood.train$label == 1)
col.ind = rep("steelblue", dim(blood.train)[1])
col.ind[index] = "red"


############################################################
# ??? # Rough PCA
############################################################

# Let us drop the two categorical variables; moreover, since
# our data have the same scale, we do not specify anything else.
blood.pca = blood.train %>%
  select(-c(DepMap_ID, label)) %>%
  prcomp()

# As expected, we obtained 817 PCs. This is ok, since we are
# doing selection variables: we hope that the first 20 PCs are
# enough to obtain a proper classification

var_explained = blood.pca$sdev^2/sum(blood.pca$sdev^2)
round(cumsum(var_explained[1:312]), 4)
# The first 312 PCs provides the 75% of explained variance

screeplot(blood.pca, type = "lines", main = "Screeplot")
# Little elbow with 3 PCs

as.data.frame(blood.pca$rotation[,1:3]) %>%
  filter(abs(PC1) > 0.07 | abs(PC2) > 0.07 | abs(PC3) > 0.05) %>%
  kable()

as.data.frame(blood.pca$sdev[1:3]) 

as.data.frame(blood.pca$x) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(color = col.ind, lwd = 2) +
  ggtitle("PC1 vs. PC2") + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100,1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100,1), "%"))

# There is not a clear cut btw cancer obs and non-carcer obs.
# However, it can be notice that cancer obs tend to stay in 
# the left-up part of the 2-dim PCs plot.
# In particular, blood cancer cells have:
#   - less than 5 Pc1 values (except for 1 obs) 
#   - positive PC2 values (or near to zero)

scatterplot3d(blood.pca$x[,1:3], color = col.ind, pch = 16,
              main = "PC1 vs. PC2 vs. PC3")
legend("bottomright", legend = c("Blood Cancer", "Other"),
       col =  c("red", "steelblue"), pch = 16) 


# Similar notes also for the 3-dim PCs plot.
# This is an excellent result: we have gained a higlhy reduction 
# in dimensionality with just PCA.

# Thus, clouds of points, but hopefully first PCs can provide
# enough information

boxplot(blood.pca$x[,1], main = "PC1", outpch = 16, outcol = col.ind)
boxplot(blood.pca$x[,2], main = "PC2", outpch = 16, outcol = col.ind)
boxplot(blood.pca$x[,3], main = "PC3", outpch = 16, outcol = col.ind)
# outliers are all non-blood-cancer cells


# Let us try with the svd (= singular value decomposition) function
# and let us see the condition number
blood.S = blood.train %>%
  select(-c(DepMap_ID, label)) %>%
  svd()

# REM: the condition number k is the largest singular value 
#      divided by the smallest
# This provides a general thumb rule to detecting multicollinearity

which(blood.S$d < 0.1)
# None of the square roots of the eigenvalues of the X'X are close to zero

max(blood.S$d)/min(blood.S$d)
# k is 213.6117 >> 30
# --> big multicollinearity problem (as intuitively one expects)


############################################################
# ??? # PCA with adjustment for minority class 
############################################################
# First, let us balance the minority class using the SMOTE function 

length(which(blood.train$label == 1))/nrow(blood.train)
# 11% of obs are blood-cancer labelled

genData = SMOTE(blood.train[,-1], blood.train$label, dup_size = 1)
blood.train_bal = genData$data

length(which(blood.train_bal$label == 1))/nrow(blood.train_bal)
# about 20%

# write.csv(blood.train_bal,"C:/Users/user/Desktop/SML/Project/blood_training_balnced.csv",
#           row.names = F)

colnames(blood.train_bal)[c(1, 17394, 17395)]
which(blood.train_bal$label != blood.train_bal$class)

# To color differently in the plot
index_bal = which(blood.train_bal$label == 1)
col.ind_bal = rep("steelblue", dim(blood.train_bal)[1])
col.ind_bal[index_bal] = "red"

blood.pca_bal = blood.train_bal %>%
  select(-c(label, class)) %>%
  prcomp()

var_explained = blood.pca_bal$sdev^2/sum(blood.pca_bal$sdev^2)
round(cumsum(var_explained[1:294]), 4)
# The first 294 PCs provides the 75% of explained variance
# --> better than before
# Now the first PC explains more than the one of imbalanced dataset

as.data.frame(blood.pca_bal$x) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(color = col.ind_bal, lwd = 2) +
  ggtitle("PC1 vs. PC2") + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100,1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100,1), "%"))


scatterplot3d(blood.pca_bal$x[,1:3], color = col.ind_bal, pch = 16,
              main = "PC1 vs. PC2 vs. PC3")
legend("bottomright", legend = c("Blood Cancer", "Other"),
       col =  c("red", "steelblue"), pch = 16)

# Still a cloud of points, but now blood-cancer obs are explained
# differently by the PCs

boxplot(blood.pca$x[,1], main = "PC1", outpch = 16, outcol = col.ind_bal)
boxplot(blood.pca$x[,2], main = "PC2", outpch = 16, outcol = col.ind_bal)
boxplot(blood.pca$x[,3], main = "PC3", outpch = 16, outcol = col.ind_bal)

# Now outliers are blood-cancer cells