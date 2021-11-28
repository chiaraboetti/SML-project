# LASSO-SVM 
library(sparseSVM)
library(caret)
library(smotefamily)
library(tidyverse)
library(caret)
library(glmnet)

# load the dataset
train = read.csv("../../dataset/blood_training.csv")
test = read.csv("../../dataset/blood_test.csv")

# balancing the minority class using the SMOTE function 
length(which(train$label == 1))/nrow(train)
genData = SMOTE(train[,-1], train$label, dup_size = 6)
blood.training_balanced = genData$data
length(which(blood.training_balanced$label == 1))/nrow(blood.training_balanced)

train_X = as.matrix(train[, -c(1, 17395)])
train_y = train$label

test_X = as.matrix(test[, -c(1,17395)])
test_y = test$label


# Find the best lambda using cross-validation
set.seed(42) 
cv.lasso = cv.glmnet(train_X, train_y, alpha = 1, family = "binomial")
model =  glmnet(train_X, train_y, alpha = 1, family = "binomial",
                lambda = cv.lasso$lambda.min)

plot(cv.lasso)

# Make predictions on the test data
probabilities = model %>% predict(newx = test_X)
predicted.classes = ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(factor(pred_y), factor(test_y))


