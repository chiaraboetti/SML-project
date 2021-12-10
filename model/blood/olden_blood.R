library(tidyverse)
library(neuralnet)
library(nnet)
library(NeuralNetTools)
library(caret)

blood = read_csv("../../dataset/blood.csv")
blood = blood[, -c(1, 2)]

set.seed(42)
splitIndex = createDataPartition(blood$label, p=0.6, list=FALSE)
data1 = blood[splitIndex,] 
data2 = blood[-splitIndex,]

# train whole NN
set.seed(2311)
splitIndex = createDataPartition(data1$label, p=0.75, list=FALSE)
train = data1[splitIndex,] 
test = data1[-splitIndex,]

names(train) = c(names(train)[1:17393],"IsBlood")

# Set up formula
n = names(train)
f = as.formula(paste("IsBlood ~",
                     paste(n[!n %in% c("IsBlood")], collapse = " + ")))

importance = rep(0, 17393)
for (i in 1:8){
  
  nn_i = neuralnet(f,
                 data = train,
                 hidden = c(400, 300, 1),
                 act.fct = "sigmoid",
                 linear.output = FALSE,
                 lifesign = "minimal")
  
  d = olden(nn)$data
  importance = importance + d[order(as.numeric(row.names(d))), ]$importance
  
  
  
}
write.csv(importance,'olden_multiclass.csv')




