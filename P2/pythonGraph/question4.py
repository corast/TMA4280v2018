#script to create graphs from the csv files.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def loadData(file):
    df = pd.read_csv(file, sep=";", header = 0)
    return df.values

#We now need to split this dataset into two parts.

print(df)

def graphs_n_time(data):
    #split the dataset into X and Y for easy plotting.
    X = data[:,0]
    Y = dat[:,1]
    print(X)
    print(Y)
    pass
