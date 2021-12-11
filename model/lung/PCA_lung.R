# PCA

library(ggplot2)
library(scatterplot3d)
library(dplyr)
library(knitr)
library(smotefamily)
#####################################################

# load the datase
lung = read.csv("dataset/lung_dataset.csv")

index = which(lung$label == 1)
col.ind = rep("steelblue", dim(lung)[1])
col.ind[index] = "orange"


## Rough PCA 
lung.pca = lung %>%
  select(-c(X, DepMap_ID, label)) %>%
  prcomp()

var_explained = lung.pca$sdev^2/sum(lung.pca$sdev^2)
round(cumsum(var_explained[1:370]), 4)
# The first 369 PCs provides the 75% of explained variance


## Look at the first 3 PCs
screeplot(lung.pca, type = "lines", main = "Screeplot")

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
# do not provide enough information.

boxplot(lung.pca$x[,1], main = "PC1", outpch = 16, outcol = col.ind)
boxplot(lung.pca$x[,2], main = "PC2", outpch = 16, outcol = col.ind)
boxplot(lung.pca$x[,3], main = "PC3", outpch = 16, outcol = col.ind)

## Singular value decomposition
lung.S = lung %>%
  select(-c(X, DepMap_ID, label)) %>%
  svd()

which(lung.S$d < 0.1)

# Condition number k as general thumb rule to detecting multicollinearity
max(lung.S$d)/min(lung.S$d)
# k is 253.2206 >> 30
# --> big multicollinearity problem (as intuitively one expects)


## PCA with adjustment for minority class 
length(which(lung$label == 1))/nrow(lung)
# --> 12% of obs are lung cancer 

genData = SMOTE(lung[,-c(1,2)], lung$label, dup_size = 5)
lung_bal = genData$data
length(which(lung_bal$label == 1))/nrow(lung_bal)
# --> 45% of obs are lung cancer

which(lung_bal$label != lung_bal$class)

index_bal = which(lung_bal$label == 1)
col.ind_bal = rep("steelblue", dim(lung_bal)[1])
col.ind_bal[index_bal] = "orange"

lung.pca_bal = lung_bal %>%
  select(-c(label, class)) %>%
  prcomp()

var_explained = lung.pca_bal$sdev^2/sum(lung.pca_bal$sdev^2)
round(cumsum(var_explained[1:300]), 4)
# The first 300 PCs provides the 75% of explained variance

as.data.frame(lung.pca_bal$x) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(color = col.ind_bal, lwd = 2) +
  ggtitle("PC1 vs. PC2") + 
  labs(x = paste0("PC1: ", round(var_explained[1]*100,1), "%"),
       y = paste0("PC2: ", round(var_explained[2]*100,1), "%"))

scatterplot3d(lung.pca_bal$x[,1:3], color = col.ind_bal, pch = 16,
              main = "PC1 vs. PC2 vs. PC3")
legend("bottomright", legend = c("Lung Cancer", "Other"),
       col =  c("orange", "steelblue"), pch = 16)

# Still a cloud of points: PCs do not explain well Lung-cancer cells

boxplot(lung.pca_bal$x[,1], main = "PC1", outpch = 16, outcol = col.ind_bal)
boxplot(lung.pca_bal$x[,2], main = "PC2", outpch = 16, outcol = col.ind_bal)
boxplot(lung.pca_bal$x[,3], main = "PC3", outpch = 16, outcol = col.ind_bal)
