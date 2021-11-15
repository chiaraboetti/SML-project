library(ggplot2)
library(scatterplot3d)
library(dplyr)
library(knitr)
#####################################################

lung.train = read.csv("dataset/lung_training.csv")
lung.test = read.csv("dataset/lung_test.csv")
# data.frame(which(is.na(lung.train), arr.ind = TRUE))

# To color differently in the plot
index = which(lung.train$label == 1)
col.ind = rep("steelblue", dim(lung.train)[1])
col.ind[index] = "orange"


############################################################
# ??? # Rough PCA
############################################################

# A full-rank dispersion matrix S [variance-covariance] cannot be 
# estimated using a number of observations n smaller than or equal 
# to the number of descriptors p.
# When n<=p, since there are n???1 DF in total, the rank of the 
# resulting S matrix of order p is (n???1).
# In such a case, the eigen-decomposition of S produces (n???1) real
# and p???(n???1) null eigenvalues. 
# Positioning n objects while respecting their distances requires
# (n???1) dimensions only.
# A PCA where n<=p produces (n???1) eigenvalues larger than 0 and
# the (n???1) corresponding eigenvectors and principle components.


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
round(cumsum(var_explained[1:20]), 4)
# Very very bad situation

screeplot(lung.pca, type = "lines", main = "Screeplot")


boxplot(lung.pca$x[,1], main = "PC1")
boxplot(lung.pca$x[,2], main = "PC2")
boxplot(lung.pca$x[,3], main = "PC3")


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


# library(plotly)
# plot_ly(x = lung.pca$x[,1], y = lung.pca$x[,2], z =lung.pca$x[,3], 
#         type="scatter3d", mode = "markers", size = 2, color = col.ind) %>% 
#   layout(scene = list(xaxis = list(title = "PC1"),
#                       yaxis = list(title = "PC2"),
#                       zaxis = list(title = "PC3")))


############################################################
# ??? # PCA with adjustment for minority class 
############################################################
