library(tidyverse)
library(neuralnet)
library(nnet)
library(NeuralNetTools)
library(caret)

blood = read_csv("../../dataset/blood_1.csv")
data = blood[, -c(1, 2, 3)]

# train whole NN
set.seed(42)
splitIndex = createDataPartition(data$label, p=0.8, list=FALSE)
train = data[splitIndex,] 
test = data[-splitIndex,]

names(train) = c(names(train)[1:17393],"IsBlood")

# Set up formula
n = names(train)
f = as.formula(paste("IsBlood ~",
                     paste(n[!n %in% c("IsBlood")], collapse = " + ")))

importance = rep(0, 17393)
blood_obs = train[train$IsBlood == 1, ]
for (i in 1:1){
  set.seed(i)
  no_blood_obs = sample_n(train[train$IsBlood == 0, ], 60)
  df = rbind(blood_obs, no_blood_obs)
  df = df[sample(nrow(df)), ]
  
  nn_i = neuralnet(f,
                 data = df,
                 hidden = c(400, 300, 1),
                 act.fct = "logistic",
                 linear.output = FALSE,
                 lifesign = "minimal")
  
  d = olden(nn_i)$data
  importance = importance + d[order(as.numeric(row.names(d))), ]$importance
}

write.csv(importance,'olden_blood.csv')

