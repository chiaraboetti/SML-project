# LASSO-SVM 
library(sparseSVM)
library(caret)

# load the dataset
dataset = read.csv("../../dataset/blood_dataset.csv")

# splitting into two datasets
set.seed(42)
split_dummy = sample(c(rep(0, 0.7 * nrow(dataset)), rep(1, 0.3 * nrow(dataset))))
train = dataset[split_dummy == 0, ]  
test = dataset[split_dummy == 1, ]  

train_X = as.matrix(train[, -c(1, 2, 17396)])
train_y = train$label

test_X = as.matrix(test[, -c(1, 2, 17396)])
test_y = test$label

# choose the best alpha (ElasticNet) and lambda (for shrinkage) with 10-fold-cv
cv.SVM = cv.sparseSVM(train_X, train_y, seed = 42)
plot(cv.SVM)
log(cv.SVM$lambda.min)

# fit this SVM
lasso.SVM = sparseSVM(train_X, train_y, lambda = cv.SVM$lambda.min)
pred_y = predict(object = lasso.SVM, X = test_X)
confusionMatrix(factor(pred_y), factor(test_y))

sum(lasso.SVM$weights != 0)
# We have selected 108 variables and achieved great accuracy!!


