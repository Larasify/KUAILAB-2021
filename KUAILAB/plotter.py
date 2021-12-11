import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotqp(q,p):
    fig = plt.figure()
    ax = Axes3D(fig)
    q=np.transpose(q)
    p=np.transpose(p)
    ax.scatter(q[0],q[1],q[2],marker='o')
    ax.scatter(p[0],p[1],p[2],marker='^')
    plt.show()

def plotpoints(qlist,plist,brutepmatrix,improvedpmatrix):
    initqmatrix = np.array(qlist)
    initpmatrix = np.array(plist)
    print("Plotting initial matrices")
    plotqp(initqmatrix,initpmatrix) #INITIAL MATRICES PLOTTED
    print("Plotting matrices made by bruteforce algorithm")
    plotqp(initqmatrix, brutepmatrix) #BRUTE FORCE P MATRIX 
    print("Plotting matrices made by my implementation")
    plotqp(initqmatrix,improvedpmatrix) #MY ALGORITHMN P MATRIX