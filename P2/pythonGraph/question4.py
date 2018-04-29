#script to create graphs from the csv files.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def loadData(file):
    df = pd.read_csv(file, sep=";", header = 0)
    return df.values



def graph_speedup(data1, data2):
    #sp = T(1)/T(p) [time with processors p T(p)]
    fig, ax = plt.subplots()
    dim = data1.shape[0]#count rows of data.
    X = data1[:,0]
    Y = data1[:,1]
    T1 = Y[0]
    S = []
    for x in range(Y.shape[0]):
        S.append(T1/Y[x])
    #plot X with S.
    ax.plot(X,S, label="n=2^10")

    dim = data2.shape[0]#count rows of data.
    X2 = data2[:,0]
    Y2 = data2[:,1]
    T12 = Y2[0]
    S2 = []
    for x in range(Y2.shape[0]):
        S2.append(T12/Y2[x])
    #plot X2 with S2.
    ax.plot(X2,S2, label="n=2^12")

    O = list(range(0,dim+1)) #plot optimal line x/y
    ax.plot(O, O,'k--',label="ideal speedup")

    plt.xlabel("number of processes p")
    plt.ylabel("speedup (T1/Tp)")
    #plt.xscale('log')
    plt.ylim(ymax=dim/2, ymin=0)
    plt.xlim(xmax=dim, xmin=0)
    plt.grid()
    plt.title("speedup with increased processes")
    ax.legend(loc="upper left")
    plt.show()


def plot_complexity(data1, data4, data8):
    #plot the run time with p = 1 for differnt values of n, against the predicted complexity function O(n^2 log(n))
    fig, ax = plt.subplots()
    Y = data1[:,0]
    X = data1[:,1]
    ax.plot(X,Y,label="p = 1")
    Y1 = data4[:,0]
    X1 = data4[:,1]
    ax.plot(X1,Y1,label="p = 4")

    Y2 = data8[:,0]
    X2 = data8[:,1]
    ax.plot(X2,Y2, label="p = 8")


    ax.legend(loc="upper left")
    plt.ylabel("time (sec)")
    plt.xlabel("n")
    plt.title("Fixed p, different n")

    plt.show()

def plot_walltime(data1, data2):
    #plot the run time with p = 1 for differnt values of n, against the predicted complexity function O(n^2 log(n))
    fig, ax = plt.subplots()
    X = data1[:,0]
    Y = data1[:,1]
    ax.plot(X,Y, label = "n=2^10")
    
    X1 = data2[:,0]
    Y1 = data2[:,1]
    ax.plot(X1,Y1, label = "n=2^12")



    ax.legend(loc="upper center")
    plt.ylabel("time (sec)")
    plt.yscale('log')
    plt.xlabel("p")
    plt.title("walltime on p")
    plt.grid()
    plt.show()

def plot_complexity_function(data1):
    
    X = data1[:,1]
    
    # we need to plow the n^2 log(n) function from this.
    Y = [o(x) for x in X]
    plt.plot(X,Y)
    plt.xlabel("n")
    plt.title("growth of O(n^2 log(n))")
    plt.show()
    
def o(x):
    return x**2*np.log2(x)

def plot_efficiency(data1, data2):
    #sp = T(1)/T(p) [time with processors p T(p)]
    fig, ax = plt.subplots()
    dim = data1.shape[0]#count rows of data.
    X = data1[:,0]
    Y = data1[:,1]
    T1 = Y[0]
    S = []
    for x in range(Y.shape[0]):
        S.append(T1/(Y[x] * X[x]))
    #plot X with S.
    ax.plot(X,S, label="n=2^10")

    dim = data2.shape[0]#count rows of data.
    X2 = data2[:,0]
    Y2 = data2[:,1]
    T12 = Y2[0]
    S2 = []
    for x in range(Y2.shape[0]):
        S2.append(T12/(Y2[x] * X[x]))
    #plot X2 with S2.
    ax.plot(X2,S2, label="n=2^12")

    plt.xlabel("number of processes p")
    plt.ylabel("Efficiency %")
    plt.title("Efficiency as np increases")
    ax.legend(loc="upper center")
    plt.grid()
    plt.show()


data10 = loadData("../output/question4-p-10.csv");
data12 = loadData("../output/question4-p-12.csv");

datap1 = loadData("../output/question4-1-n-1.csv")
datap4 = loadData("../output/question4-4-n-1.csv")
datap8 = loadData("../output/question4-8-n-1.csv")

#plot_complexity(datap1,datap4,datap8)
#plot_walltime(data10, data12)
#graph_speedup(data10,data12)
#plot_complexity_function(datap1)
#plot_efficiency(data10, data12)