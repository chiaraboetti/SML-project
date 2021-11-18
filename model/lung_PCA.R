library(ggplot2)
library(scatterplot3d)
library(dplyr)
library(knitr)
#####################################################

lung.train = read.csv("dataset/lung_training.csv")
lung.test = read.csv("dataset/lung_test.csv")

# To color differently in the plot
index = which(lung.train$label == 1)
col.ind = rep("steelblue", dim(lung.train)[1])
col.ind[index] = "orange"


############################################################
# ♫ # Rough PCA
############################################################

# A full-rank dispersion matrix S [variance-covariance] cannot be 
# estimated using a number of observations n smaller than or equal 
# to the number of descriptors p.
# When n<=p, since there are n???1 DF in total, the rank of the 
# resulting S matrix of order p is (n-1).
# In such a case, the eigen-decomposition of S produces (n-1) real
# and p*(n-1) null eigenvalues. 
# Positioning n objects while respecting their distances requires
# (n-1) dimensions only.
# A PCA where n<=p produces (n-1) eigenvalues larger than 0 and
# the (n-1) corresponding eigenvectors and principle components.


# Because PCA works best with numerical data, we exclude the 
# two categorical variables
# Since our data have the same scale, we do not specify anything
lung.pca = lung.train %>%
  select(-c(DepMap_ID, label)) %>%
  prcomp()
# As expected, we obtained 817 PCs. This is ok, since we are
# doing selection variables: we hope that the first 20 PCs are
# enough to obtain a proper classification

var_explained = lung.pca$sdev^2/sum(lung.pca$sdev^2)
round(cumsum(var_explained[1:312]), 4)
# The first 312 PCs provides the 75% of explained variance

screeplot(lung.pca, type = "lines", main = "Screeplot")
# Little elbow with 3 PCs

as.data.frame(lung.pca$rotation[,1:3]) %>%
  filter(abs(PC1) > 0.07 | abs(PC2) > 0.07 | abs(PC3) > 0.05) %>%
  kable()

as.data.frame(lung.pca$sdev[1:3]) 

as.data.frame(lung.pca$x) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(color = col.ind, lwd = 2) +
  ggtitle("PC1 vs. PC2") + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100,1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100,1), "%"))


scatterplot3d(lung.pca$x[,1:3], color = col.ind, pch = 16,
              main = "PC1 vs. PC2 vs. PC3")
legend("bottomright", legend = c("Lung Cancer", "Other"),
       col =  c("orange", "steelblue"), pch = 16) 

# Both two and three PCs show a cloud of points, which means they
# do not provide enough information

boxplot(lung.pca$x[,1], main = "PC1")
boxplot(lung.pca$x[,2], main = "PC2")
boxplot(lung.pca$x[,3], main = "PC3")


# Let us try with the svd (= singular value decomposition) function
# and let us see the condition number
lung.S = lung.train %>%
  select(-c(DepMap_ID, label)) %>%
  svd()

# REM: the condition number k is the largest singular value 
#      divided by the smallest
# This provides a general thumb rule to detecting multicollinearity

which(lung.S$d < 0.1)
# None of the square roots of the eigenvalues of the X'X are close to zero
max(lung.S$d)/min(lung.S$d)
# k is 213.6117 >> 30
# --> big multicollinearity problem (as intuitively one expects)


############################################################
# ♫ # PCA with adjustment for minority class 
############################################################

lung.train_bal = read.csv("dataset/lung_training_balnced.csv")

# To color differently in the plot
index = which(lung.train_bal$label == 1)
col.ind = rep("steelblue", dim(lung.train_bal)[1])
col.ind[index] = "orange"

lung.pca_bal = lung.train_bal %>%
  select(-c(DepMap_ID, label)) %>%
  prcomp()

var_explained = lung.pca_bal$sdev^2/sum(lung.pca_bal$sdev^2)
round(cumsum(var_explained[1:312]), 4)
# The first 312 PCs provides the 75% of explained variance

as.data.frame(lung.pca_bal$x) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(color = col.ind, lwd = 2) +
  ggtitle("PC1 vs. PC2") + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100,1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100,1), "%"))


# Conclusion: there is no improvement
