import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import math
from scipy import stats
#from normalizelibrary import *
from normalizelibrary3 import *
from sortmatrix3_11_24 import *

# Main PHASTpep file code
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

def Forest(Output,Input):
        print('starting Forest')
        print(Input)

        X = pd.read_excel(Input)
        print(X)

        y = X.iloc[:,-2]
        print(y)
        X = X.iloc[:,2:len(X.columns)-2]
        print(X)

        # Assuming X is your feature matrix of size m by n and y is your target vector of size m
        # Split the data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

        # Initialize the Random Forest Regressor model
        rf_model = RandomForestRegressor(n_estimators=100, random_state=42)

        # Train the model
        rf_model.fit(X_train, y_train)

        # Make predictions on the testing set
        y_pred = rf_model.predict(X_test)

        # Evaluate the model (optional)
        mse = mean_squared_error(y_test, y_pred)
        print("Mean Squared Error:", mse)
        print(y_test)
        OldNew = pd.DataFrame()
        OldNew['y_test'] = y_test
        OldNew['y_pred'] = y_pred
        print(OldNew)

        OldNew.to_excel(Output)
        return


print(snakemake.input)
print(snakemake.output)
Forest(snakemake.output[0], snakemake.input[0])
