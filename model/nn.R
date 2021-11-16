library(smotefamily)
library(keras)



lung.train = read.csv("../dataset/lung_training.csv")
lung.test = read.csv("../dataset/lung_test.csv")
genData = SMOTE(lung.train[,-1],lung.train[,17395], dup_size = 1)
lung.training_balanced = genData$data #192 positive obs over 917 (about 21%)


network = keras_model_sequential() %>%
  layer_dense(units = 512,
              activation = "relu",
              kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001),
              input_shape = ) %>%
  layer_dense(units = 1, activation = "softmax")

network %>% compile(
  optimizer = "rmsprop",
  loss = "mse",
  metrics = c("accuracy")
)

network %>% fit(train_images,
                train_labels,
                epochs = 5,
                batch_size = 128)

metrics = network %>% evaluate(test_images, test_labels)

metrics

