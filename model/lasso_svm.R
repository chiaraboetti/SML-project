# LASSO-SVM 
library(sparseSVM)
library(caret)

# load the dataset
train = read.csv("../dataset/lung_training_balanced.csv")
test = read.csv("../dataset/lung_test.csv")

train_X = as.matrix(train[, -c(1, 17395)])
train_y = train$label

test_X = as.matrix(test[, -c(1,17395)])
test_y = test$label

# choose the best alpha (ElasticNet) and lambda (for shrikage) with 10-fold-cv
cv.SVM = cv.sparseSVM(train_X, train_y, seed = 42)
plot(cv.SVM)
log(cv.SVM$lambda.min)

# fit Lasso regularization to achieve sparsity
lasso.SVM = sparseSVM(train_X, train_y, alpha = 1, lambda = 0.03)
sum(lasso.SVM$weights > 0)
pred_y = predict(object = lasso.SVM, X = test_X)
confusionMatrix(data=pred_y, reference = test_y)

