library(ggplot2)
library(scatterplot3d)
library(dplyr)
library(knitr)
library(smotefamily)
#####################################################

blood = read.csv("dataset/blood_dataset.csv")

# To color differently in the plot
index = which(blood$label == 1)
col.ind = rep("steelblue", dim(blood)[1])
col.ind[index] = "red"

colnames(blood)[1:4]
# Need to remove X and DepMap_ID

############################################################
#     Rough PCA                                            #
############################################################

# Let us drop the two categorical variables; moreover, since
# our data have the same scale, we do not specify anything else.
blood.pca = blood %>%
  select(-c(X, DepMap_ID, label)) %>%
  prcomp()

# As expected, we obtained 1022 PCs. This is ok, since we are
# doing selection variables: we hope that the first 50 or 100 PCs
# are enough to obtain a proper classification

var_explained = blood.pca$sdev^2/sum(blood.pca$sdev^2)
round(cumsum(var_explained[1:370]), 4)
# The first 369 PCs provides the 75% of explained variance

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
#   - more than -5 PC1 values (except for 1 obs) 
#   - negative PC2 values (or near to zero)

scatterplot3d(blood.pca$x[,1:3], color = col.ind, pch = 16,
              main = "PC1 vs. PC2 vs. PC3")
legend("bottomright", legend = c("Blood Cancer", "Other"),
       col =  c("red", "steelblue"), pch = 16) 


# Similar notes also for the 3-dim PCs plot.
# This is an excellent result: we have gained a highly reduction 
# in dimensionality with just PCA.

# Thus, clouds of points, but hopefully first PCs can provide
# enough information

boxplot(blood.pca$x[,1], main = "PC1", outpch = 16, outcol = col.ind)
boxplot(blood.pca$x[,2], main = "PC2", outpch = 16, outcol = col.ind)
boxplot(blood.pca$x[,3], main = "PC3", outpch = 16, outcol = col.ind)
# outliers are all non-blood-cancer cells


# Let us try with the svd (= singular value decomposition) function
# and let us see the condition number
blood.S = blood %>%
  select(-c(X,DepMap_ID, label)) %>%
  svd()

# REM: the condition number k is the largest singular value 
#      divided by the smallest
# This provides a general thumb rule to detecting multicollinearity

which(blood.S$d < 0.1)
# None of the square roots of the eigenvalues of the X'X are close to zero

max(blood.S$d)/min(blood.S$d)
# k is 253.2206 >> 30
# --> big multicollinearity problem (as intuitively one expects)


############################################################
#       PCA with adjustment for minority class             #
############################################################
# First, let us balance the minority class using the SMOTE function 

length(which(blood$label == 1))/nrow(blood)
# 11% of obs are blood-cancer labelled

genData = SMOTE(blood[,-c(1,2)], blood$label, dup_size = 3)
blood_bal = genData$data

length(which(blood_bal$label == 1))/nrow(blood_bal)
# about 30%

colnames(blood_bal)[c(1, 17394, 17395)]
which(blood_bal$label != blood_bal$class)

# To color differently in the plot
index_bal = which(blood_bal$label == 1)
col.ind_bal = rep("steelblue", dim(blood_bal)[1])
col.ind_bal[index_bal] = "red"

blood.pca_bal = blood_bal %>%
  select(-c(label, class)) %>%
  prcomp()

var_explained = blood.pca_bal$sdev^2/sum(blood.pca_bal$sdev^2)
round(cumsum(var_explained[1:305]), 4)
# The first 203 PCs provides the 75% of explained variance
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

boxplot(blood.pca_bal$x[,1], main = "PC1", outpch = 16, outcol = col.ind_bal)
boxplot(blood.pca_bal$x[,2], main = "PC2", outpch = 16, outcol = col.ind_bal)
boxplot(blood.pca_bal$x[,3], main = "PC3", outpch = 16, outcol = col.ind_bal)
# Now outliers are blood-cancer cells


############################################################
#           PCA for the filtered dataset                   #
############################################################
blood_train = read.csv("dataset/blood_filt_train.csv")
blood_test = read.csv("dataset/blood_filt_test.csv")

blood_fil = rbind(blood_train, blood_test)

# To color differently in the plot
index_fil = which(blood_fil$label == 1)
col.ind_fil = rep("steelblue", dim(blood_fil)[1])
col.ind_fil[index_fil] = "red"

colnames(blood_fil)

blood.pca_fil = blood_fil %>%
  select(-c(X, DepMap_ID, label)) %>%
  prcomp()

var_explained = blood.pca_fil$sdev^2/sum(blood.pca_fil$sdev^2)
round(cumsum(var_explained[1:22]), 4)
# The first 22 PCs provides the 75% of explained variance
# --> very good result!

as.data.frame(blood.pca_fil$x) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(color = col.ind_fil, lwd = 2) +
  ggtitle("PC1 vs. PC2") + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100,1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100,1), "%"))


scatterplot3d(blood.pca_fil$x[,1:3], color = col.ind_fil, pch = 16,
              main = "PC1 vs. PC2 vs. PC3")
legend("bottomright", legend = c("Blood Cancer", "Other"),
       col =  c("red", "steelblue"), pch = 16)

# Now there is a clear distinction btw non-cancer and cancer obs

boxplot(blood.pca_fil$x[,1], main = "PC1", outpch = 16, outcol = col.ind_bal)
boxplot(blood.pca_fil$x[,2], main = "PC2", outpch = 16, outcol = col.ind_bal)
boxplot(blood.pca_fil$x[,3], main = "PC3", outpch = 16, outcol = col.ind_bal)

as.data.frame(blood.pca_fil$rotation[,1:3]) %>%
  filter(abs(PC1) > 0.15 | abs(PC2) > 0.15 | abs(PC3) > 0.15) %>%
  kable()
