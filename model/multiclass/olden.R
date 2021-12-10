library(tidyverse)
library(neuralnet)
library(nnet)
library(NeuralNetTools)
library(caret)

multiclass = read_csv("../../dataset/multiclass_dataset_no_weird_obs.csv")
multiclass = multiclass[, -c(1, 2)]

set.seed(42)
splitIndex = createDataPartition(multiclass$label, p=0.6, list=FALSE)
data1 = multiclass[splitIndex,] 
data2 = multiclass[-splitIndex,]

# train whole NN
set.seed(2311)
splitIndex = createDataPartition(data1$label, p=0.75, list=FALSE)
train = data1[splitIndex,] 
test = data1[-splitIndex,]

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




