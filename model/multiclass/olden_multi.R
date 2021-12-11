library(tidyverse)
library(neuralnet)
library(nnet)
library(NeuralNetTools)
library(caret)

multiclass = read_csv("../../dataset/multi_1.csv")
data = multiclass[, -c(1, 2, 3)]


# train whole NN
set.seed(42)
splitIndex = createDataPartition(data$label, p=0.8, list=FALSE)
train = data[splitIndex,] 
test = data[-splitIndex,]

train = cbind(train[,-c(17394)], class.ind(as.factor(train$label)))
names(train) = c(names(train)[1:17393],"Gastrointestinal","Genitals","Muscle_Bone",
                 "Neuro","Breast","Head_Neck","Blood", "Genitourinary", "Lung")


# Set up formula
n = names(train)
f = as.formula(paste("Gastrointestinal + Genitals + Muscle_Bone + Neuro + 
                     Breast + Head_Neck + Blood + Genitourinary + Lung ~",
                     paste(n[!n %in% c("Gastrointestinal","Genitals","Muscle_Bone","Neuro","Breast","Head_Neck","Blood",
                                       "Genitourinary", "Lung")], collapse = " + ")))

nn = neuralnet(f,
                data = train,
                hidden = c(400, 300, 9),
                act.fct = "logistic",
                linear.output = FALSE,
                lifesign = "minimal")

predict = function(data){
  prediction = data.frame(neuralnet::compute(nn,data.frame(data[,1:17393]))$net.result)
  labels = c("Gastrointestinal","Genitals","Muscle_Bone","Neuro","Breast","Head_Neck","Blood",
             "Genitourinary", "Lung")
  prediction_label = data.frame(max.col(prediction)) %>% 
    mutate(prediction=labels[max.col.prediction.]) %>% 
    select(2) %>% 
    unlist()
  
  table(data$label, prediction_label)
}

predict(test)

importance = olden(nn)$data
write.csv(importance,'olden_multiclass.csv')





