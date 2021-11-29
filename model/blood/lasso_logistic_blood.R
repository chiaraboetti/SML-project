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
#length(which(train$label == 1))/nrow(train)
# genData = SMOTE(train[,-1], train$label, dup_size = 6)
# blood.training_balanced = genData$data
# length(which(blood.training_balanced$label == 1))/nrow(blood.training_balanced)

train_X = as.matrix(train[, -c(1, 17395)])
train_y = train$label

test_X = as.matrix(test[, -c(1,17395)])
test_y = test$label


# Find the best lambda using cross-validation
set.seed(42) 
cv.lasso = cv.glmnet(train_X, as.factor(train_y), alpha = 1, family = "binomial")
plot(cv.lasso)

which(cv.lasso$lambda == cv.lasso$lambda.min)
model =  glmnet(train_X, as.factor(train_y), alpha = 1, family = "binomial", 
                lambda = 0.008)



# Make predictions on the test data
probabilities = model %>% predict(newx = test_X)
predicted.classes = ifelse(probabilities > 0.5, 1, 0)
confusionMatrix(factor(predicted.classes), factor(test_y))
sum(model$beta != 0)

