# LASSO-SVM 
library(sparseSVM)
library(caret)
library(smotefamily)

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

# choose the best alpha (ElasticNet) and lambda (for shrikage) with 10-fold-cv
cv.SVM = cv.sparseSVM(train_X, train_y, seed = 42)
plot(cv.SVM)
log(cv.SVM$lambda.min)

# fit this SVM
lasso.SVM = sparseSVM(train_X, train_y, lambda = cv.SVM$lambda.min)
pred_y = predict(object = lasso.SVM, X = test_X)
confusionMatrix(factor(pred_y), factor(test_y))

sum(lasso.SVM$weights != 0)
# We have selected 114 variables
