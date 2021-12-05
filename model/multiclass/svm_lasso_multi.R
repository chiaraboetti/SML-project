# LASSO-SVM 
library(sparseSVM)
library(caret)

# load the dataset
dataset = read.csv("../../dataset/multiclass_dataset_no_weird_obs.csv")
dataset = dataset[, -c(1, 2)]

set.seed(42)
split_dummy = sample(c(rep(0, 0.7 * nrow(dataset)), rep(1, 0.3 * nrow(dataset))))
train = dataset[split_dummy == 0, ]  
test = dataset[split_dummy == 1, ]  

# OVO classification function
multiSVM = function(data){
  
  indexes = as.vector(combn(0:8, 2))
  trained_models = c()
  
  for (n in 1:36){
      
      i = indexes[2*n - 1]
      j = indexes[2*n]
      data_ij = data[data$label %in% c(i, j), ]
      train_X = as.matrix(data_ij[, -c(17394)])
      train_y = data_ij$label
      
      cv.SVM = cv.sparseSVM(train_X, train_y, seed = 42)
      lasso.SVM = sparseSVM(train_X, train_y, lambda = cv.SVM$lambda.min)
      trained_models = c(trained_models, list(lasso.SVM))
    }
    
  return(trained_models)
}


multiSVMpredict_obs = function(observation){
  predictions = c()
  for (model_i in my_model){
    pred_y = predict(object = model_i, X=observation)
    predictions = c(predictions, pred_y)
  }
  return(as.factor(predictions))
}

multiSVM.single_pred = function(observation){
  global_pred = multiSVMpredict_obs(observation)
  tmp = as.data.frame(table(global_pred))
  return(which.max(tmp$Freq) - 1)
}

multiSVM.predict = function(dataset){
  pred = c()
  for (row in 1:nrow(dataset)){
    obs = dataset[row, ]
    pred_y = multiSVM.single_pred(obs)
    pred = c(pred, pred_y)
  }
  
  return(pred)
}


# testing
test_X = as.matrix(test[, -c(17394)])
test_y = test$label


my_model = multiSVM(train)
pred_y = multiSVM.predict(test_X)
confusionMatrix(factor(pred_y), factor(test_y))

