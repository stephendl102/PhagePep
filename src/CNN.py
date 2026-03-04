import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import math
import tensorflow as tf
from scipy import stats
#from normalizelibrary import *
from normalizelibrary3 import *
from sortmatrix3_11_24 import *
from sklearn.model_selection import train_test_split

# Main PHASTpep file code
from tensorflow.keras.regularizers import l2
from tensorflow.keras.layers import Dropout
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense

def Forest(Output, Input):
    print('starting Forest')
    print(Input)

    X = pd.read_excel(Input)

    y = X.iloc[:,-2]
    print(y)
    X = X.iloc[:,2:len(X.columns)-2]
    print(X)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    X_train = np.array(X_train)
    X_test = np.array(X_test)
    y_train = np.array(y_train)
    y_test = np.array(y_test)

    # Reshape input matrix for CNN input
    X_train = X_train.reshape(X_train.shape[0], X_train.shape[1], 1)  # Reshape to (m, n, 1) for 1D convolution
    X_test = X_test.reshape(X_test.shape[0], X_test.shape[1], 1)  # Reshape to (m, n, 1) for 1D convolution

    # Define the CNN model
    model = Sequential([
        #Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(X_train.shape[1], 1), kernel_regularizer=l2(0.01)),
        Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(X_train.shape[1], 1)),
        MaxPooling1D(pool_size=2),
        Flatten(),
        #Dense(64, activation='relu', kernel_regularizer=l2(0.01)),
        Dense(64, activation='relu'),
        Dense(1)  # Output layer with one neuron for regression
    ])

    # Compile the model
    model.compile(optimizer='adam', loss='mean_squared_error')

    # Print model summary
    model.summary()

    # Train the model
    model.fit(X_train, y_train, epochs=10, batch_size=32, validation_split=0.2)

    # Make predictions
    y_pred = model.predict(X_test)

    OldNew = pd.DataFrame()
    OldNew['y_test'] = y_test
    OldNew['y_pred'] = y_pred
    print(OldNew)

    OldNew.to_excel(Output)

    # Example usage:
    # Forest("output.xlsx", "input.xlsx")

    OldNew.to_excel(Output)
    return


print(snakemake.input)
print(snakemake.output)
Forest(snakemake.output[0], snakemake.input[0])
