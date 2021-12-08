from imblearn.under_sampling import RandomUnderSampler
from tensorflow import keras
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder, StandardScaler
import eli5
from eli5.sklearn import PermutationImportance

iris = load_iris()
X = iris['data']
y = iris['target']
names = iris['target_names']
feature_names = iris['feature_names']

# One hot encoding
enc = OneHotEncoder()
Y = enc.fit_transform(y[:, np.newaxis]).toarray()

# Scale data to have mean 0 and variance 1 
# which is importance for convergence of the neural network
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Split the data set into training and testing
X_train, X_test, Y_train, Y_test = train_test_split(
    X_scaled, Y, test_size=0.5, random_state=2)

n_features = X.shape[1]
n_classes = Y.shape[1]

sum_weight = 0

def nn_model():
    model_i = keras.Sequential([
        keras.layers.Dense(8, activation="relu", input_shape=(4,)),
        keras.layers.Dense(4, activation="relu"),
        keras.layers.Dropout(0.3),
        keras.layers.Dense(3, activation="softmax"),
    ])

    model_i.compile(optimizer='adam', loss="categorical_crossentropy", metrics=['accuracy'])

    return model_i

for i in range(10):

    # define the model
    rus = RandomUnderSampler(random_state=i)
    X_res, y_res = rus.fit_resample(X_train, Y_train)

    my_model = KerasClassifier(build_fn=nn_model)    
    my_model.fit(X_res,y_res)

    perm = PermutationImportance(my_model, random_state=i).fit(X_test, Y_test)
    sum_weight += perm.feature_importances_

avg_importance = pd.DataFrame(sorted(sum_weight/10, reverse=True))
np.save('iris_importance', avg_importance) 