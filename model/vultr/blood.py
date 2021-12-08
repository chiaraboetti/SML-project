from imblearn.under_sampling import RandomUnderSampler
from tensorflow import keras
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import train_test_split
import eli5
from eli5.sklearn import PermutationImportance

dataset = pd.read_csv('../../dataset/blood_dataset.csv').iloc[:, 2:]
data1, data2 =  train_test_split(dataset, test_size=0.4, random_state=42)
train1, test1  = train_test_split(data1, test_size=0.2, random_state=46)
train2, test2  = train_test_split(data2, test_size=0.2, random_state=32)
train1_X = train1.iloc[:,:-1]
train1_y = train1.iloc[:,-1]
test1_X = test1.iloc[:,:-1]
test1_y = test1.iloc[:,-1]

def nn_model():
    model_i = keras.Sequential([
        keras.layers.Dense(400, activation="relu", input_shape=(17393,)),
        keras.layers.Dense(300, activation="relu"),
        keras.layers.Dropout(0.3),
        keras.layers.Dense(1, activation="sigmoid"),
    ])

    model_i.compile(optimizer='adam', loss="binary_crossentropy", metrics=['accuracy'])

    return model_i

sum_weight = 0

for i in range(1):

    # define the model
    rus = RandomUnderSampler(random_state=i)
    X_res, y_res = rus.fit_resample(train1_X, train1_y)

    my_model = KerasClassifier(build_fn=nn_model)    
    my_model.fit(X_res,y_res)

    perm = PermutationImportance(my_model, random_state=i).fit(test1_X, test1_y)
    sum_weight += perm.feature_importances_
    print("done: " + str(i))

avg_importance = pd.DataFrame(sorted(sum_weight/1, reverse=True))

np.save('avg_importance', avg_importance) 
imp_blood = np.load('avg_importance.npy')