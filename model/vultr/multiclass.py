import tensorflow as tf
import keras
import pandas as pd
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import backend as K
from sklearn.model_selection import train_test_split
from keras.wrappers.scikit_learn import KerasClassifier
import eli5
from eli5.sklearn import PermutationImportance

multiclass = pd.read_csv('../../dataset/multiclass_dataset_no_weird_obs.csv',header= 0)
first_split, second_split = train_test_split (multiclass, test_size=403,
                                              train_size=604, random_state=42,
                                              shuffle=True, stratify=None,
)
Training, Test = train_test_split(first_split, train_size = 480,
                                   test_size = 124, random_state=42,
                                   shuffle=True, stratify=None,
)

X_test = np.array(Test.drop(['label','DepMap_ID','Unnamed: 0'], axis = 1))
y_test = pd.get_dummies(Test, columns = ['label']).iloc[:,17395:17404]
X_train = np.array(Training.drop(['label','Unnamed: 0','DepMap_ID'], axis = 1))
y_train = pd.get_dummies(Training, columns = ['label']).iloc[:,17395:17404]

X_train_reshaped = X_train.reshape(480,17393)
X_test_reshaped = X_test.reshape(124,17393)


def categorical_focal_loss(alpha, gamma=2.):
    """
    Softmax version of focal loss.
    When there is a skew between different categories/labels in your data set, you can try to apply this function as a
    loss.
           m
      FL = ∑  -alpha * (1 - p_o,c)^gamma * y_o,c * log(p_o,c)
          c=1
      where m = number of classes, c = class and o = observation
    Parameters:
      alpha -- the same as weighing factor in balanced cross entropy. Alpha is used to specify the weight of different
      categories/labels, the size of the array needs to be consistent with the number of classes.
      gamma -- focusing parameter for modulating factor (1-p)
    Default value:
      gamma -- 2.0 as mentioned in the paper
      alpha -- 0.25 as mentioned in the paper
    References:
        Official paper: https://arxiv.org/pdf/1708.02002.pdf
        https://www.tensorflow.org/api_docs/python/tf/keras/backend/categorical_crossentropy
    Usage:
     model.compile(loss=[categorical_focal_loss(alpha=[[.25, .25, .25]], gamma=2)], metrics=["accuracy"], optimizer=adam)
    """

    alpha = np.array(alpha, dtype=np.float32)

    def categorical_focal_loss_fixed(y_true, y_pred):
        """
        :param y_true: A tensor of the same shape as `y_pred`
        :param y_pred: A tensor resulting from a softmax
        :return: Output tensor.
        """
        y_true = tf.cast(y_true, tf.float32)
        # Clip the prediction value to prevent NaN's and Inf's
        epsilon = K.epsilon()
        y_pred = K.clip(y_pred, epsilon, 1. - epsilon)

        # Calculate Cross Entropy
        cross_entropy = -y_true * K.log(y_pred)

        # Calculate Focal Loss
        loss = alpha * K.pow(1 - y_pred, gamma) * cross_entropy

        # Compute mean loss in mini_batch
        return K.mean(K.sum(loss, axis=-1))

    return categorical_focal_loss_fixed

def nn_model():
    model = Sequential()

    #1°layer
    model.add(Dense(600, input_shape=(17393,)))
    model.add(Activation('relu'))
    #model.add(Dropout(0.25))

    #2°layer
    model.add(Dense(300))
    model.add(Activation('relu'))
    model.add(Dropout(0.2))

    #3°layer
    model.add(Dense(9))
    model.add(Activation('softmax'))

    #categorical_crossentropy loss per classificazione (penalizza molto le previsioni che sono sbagliate ma con alta probabilità)
    model.compile(optimizer = 'adam', metrics=['accuracy'], loss="categorical_crossentropy")

    return model

my_model = KerasClassifier(build_fn=nn_model)    
my_model.fit(X_train_reshaped,y_train)

perm = PermutationImportance(my_model, random_state=3).fit(X_test_reshaped, y_test)
avg_importance = pd.DataFrame(sorted(perm.feature_importances_, reverse=True))
np.save('nn_multi_importance', avg_importance) 