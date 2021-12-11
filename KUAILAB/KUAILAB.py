import numpy as np
from plotter import plotpoints
import math
import time
import os.path
from os import path
#CONFIGURATION
plotting = True
BruteforceREthreshold = 2
MyImplementationREthreshold = 1
Qfilepath="Q_data.txt"
Pfilepath="P_data.txt"

def brute_force():
    print("Starting brute_force algorithm...")
    start = time.time()
    f = open("bruteforce.txt","w")
    #CONVERT INPUTS TO  TO NUMPY MATRICES
    qmatrix = np.array(qlist)
    pmatrix = np.array(plist)
    #INITIALISE CALCULATION VARIABLES
    RE = float('inf')
    EK1 = 0
    accumulatedrot = np.array([0]);
    accumulatedtrans = np.array([0,0,0]);
    iteration = 0
    N = len(qmatrix) 

    while BruteforceREthreshold<RE:
        print("Iteration "+str(iteration+1),file=f)
        #ROTATION AND TRANSLATION 
        if EK1!=0:
            RT = np.column_stack((R,t))
            RT = np.row_stack((RT,[0,0,0,1]))
            pmatrix = RT.dot(np.row_stack((np.transpose(pmatrix),np.ones(N))))
            pmatrix = np.transpose(pmatrix)
            pmatrix = np.delete(pmatrix,3,1)


        #FIND Pi closest to each Qi
        tempqmatrix = np.array(qmatrix)
        temppmatrix = np.array(pmatrix)
        pairlist = []
        for qpoint in tempqmatrix:
            pairs = [[0,0,0],[0,0,0],[0]]
            distance = float('inf')
            tempppoint = [0,0,0]
            for x in range(len(qmatrix)):
                tempdis = math.sqrt((qpoint[0]-temppmatrix[x][0])**2 + (qpoint[1]-temppmatrix[x][1])**2 +  (qpoint[2]-temppmatrix[x][2])**2)
                if tempdis < distance:
                    distance = tempdis
                    tempppoint = temppmatrix[x]
            pairs[0] = qpoint
            pairs[1] = tempppoint
            pairs[2] = distance
            pairlist.append(pairs)
        #SORT PAIRS BY DISTANCE AND CHOOSE NEXT ITERATION PAIRS
        pairlist.sort(key = lambda x:x[2])

        comQ = pairlist[iteration][0]
        comP = pairlist[iteration][1]
        iteration +=1

        #MEAN SUBTRACTED MATRICES (q')
        meansubQlist = []
        meansubPlist = []
        for point in qmatrix:
            qchange = point - comQ
            meansubQlist.append(qchange)
        for point in pmatrix:
            pchange = point - comP
            meansubPlist.append(pchange)


        #CONVERT MEANSUB LIST TO NUMPY MATRIX
        meansubQmatrix = np.transpose(np.array(meansubQlist)) #WE NEED THE POINTS AS COLLUMNS SO WE WILL TRANSPOSE 
        meansubPmatrix = np.array(meansubPlist) #WE NEED COORDINATES FROM P AS ROWS SO WE WILL NOT TRANSPOSE

        #CROSS CORRELATION MATRIX AND SVD
        W = meansubQmatrix.dot(meansubPmatrix)
        U,S,VH = np.linalg.svd(W,full_matrices = True) #GET U,S and VT(VH) from SVD FUNCTION

        #ROTATION MATRIX AND TRANSLATION VECTOR
        R = U.dot(VH) #R=U*VT
        t = comQ - R.dot(comP.T)

        #ACCUMULATED ROTATIONMATRIX AND TRANSLATION VECTOR
        if EK1 == 0:
            accumulatedrot = R
            accumulatedtrans = t
        else:
            accumulatedrot = accumulatedrot.dot(R)
            accumulatedtrans += t

        #WRITE OUTPUTS
        print("Rotation matrix:\n",R,file=f)
        print("Translation vector:\n",t,file=f)
        print("Accumulated Rotation Matrix:\n",accumulatedrot,file=f)
        print("Accumulated Translation Vector:\n",accumulatedtrans,file=f)
        #CALCULATE EK
        suml2 = 0
        for point in qmatrix:
            qi = point
            pi = [0,0,0]
            shortest_distance = float('inf')
            #closest pi point
            for ppoint in pmatrix:
                distance = math.sqrt((qi[0]-ppoint[0])**2 + (qi[1]-ppoint[1])**2 +  (qi[2]-ppoint[2])**2)
                if distance<shortest_distance:
                    shortest_distance = distance
                    pi = ppoint
            a = qi-(R.dot(pi))-t
            l2 = 0
            for j in range(len(a)):
                l2 += a[j]*a[j]
            suml2 +=l2
        EK = suml2 / N

        #CALCULATE RE and EK-1
        if EK1!=0:
            RE = abs(100*(EK1-EK)/EK1)
            print("Relative Error Change:","%.2f" % RE+"%","\nDistance Error(E):","%.28f" % EK,file=f)
        EK1 = EK
        print("-------------------------------------------------------------------------\n",file=f)

    end = time.time()
    print("Time elapsed:",end-start,"seconds")
    f.close()
    return pmatrix



def my_imp():
    print("Starting my implementation...")
    start = time.time()
    f = open("myimplementation.txt","w")
    #CONVERT INPUTS TO  TO NUMPY MATRICES
    qmatrix = np.array(qlist)
    pmatrix = np.array(plist)
    #INITIALISE CALCULATION VARIABLES
    RE = float('inf')
    EK1 = 0
    accumulatedrot = np.array([0]);
    accumulatedtrans = np.array([0,0,0]);
    iteration = 0
    N = len(qmatrix) #LENGTH OF MATRIX
    while MyImplementationREthreshold<RE:
        print("Iteration "+str(iteration+1),file=f)
        #ROTATION AND TRANSLATION 
        if EK1!=0:
            RT = np.column_stack((R,t))
            RT = np.row_stack((RT,[0,0,0,1]))
            pmatrix = RT.dot(np.row_stack((np.transpose(pmatrix),np.ones(N))))
            pmatrix = np.transpose(pmatrix)
            pmatrix = np.delete(pmatrix,3,1)


        #IF ON THE CHOSEN ITERATION STEP CALCULATE CENTRE OF MASS AND GO ON AS IF IT IS THE OPTIMAL POINTS
        #IF IT IS NOT ON THE CHOSEN ITERATION GO ON WITH THE STANDARD ITERATION
        #THIS WAY WE SKIP THE WHOLE ITERATIVE STEP ONCE EVERY 3 ITERATIONS 
        #AND WHEN WE DO SKIP WE MAKE BIGGER MOVES THAN THE BRUTE FORCE ITERATIONS WHICH GETS US CLOSER 
        #AND LETS US SKIP BUNCH OF EARLY SMALL ITERATIONS
        sumQ = [0,0,0]
        sumP = [0,0,0]
        for point in qmatrix:
            sumQ += point
        for point in pmatrix:
            sumP += point

        tempqmatrix = np.array(qmatrix)
        temppmatrix = np.array(pmatrix)
        comQ = 0;
        comP = 0;
        if iteration%6 !=5:
            pairlist = []
            for qpoint in tempqmatrix:
                pairs = [[0,0,0],[0,0,0],[0]]
                distance = float('inf')
                tempppoint = [0,0,0]
                for x in range(len(qmatrix)):
                    tempdis = math.sqrt((qpoint[0]-temppmatrix[x][0])**2 + (qpoint[1]-temppmatrix[x][1])**2 +  (qpoint[2]-temppmatrix[x][2])**2)
                    if tempdis < distance:
                        distance = tempdis
                        tempppoint = temppmatrix[x]
                pairs[0] = qpoint
                pairs[1] = tempppoint
                pairs[2] = distance
                pairlist.append(pairs)
            pairlist.sort(key = lambda x:x[2])
            comQ = pairlist[iteration][0]
            comP = pairlist[iteration][1]
        elif iteration%6 == 5:
            comQ = sumQ/N
            comP = sumP/N
        iteration +=1

        #MEAN SUBTRACTED MATRICES (q')
        meansubQlist = []
        meansubPlist = []

        for point in qmatrix:
            qchange = point - comQ
            meansubQlist.append(qchange)
        for point in pmatrix:
            pchange = point - comP
            meansubPlist.append(pchange)


        #CONVERT MEANSUB LIST TO NUMPY MATRIX
        meansubQmatrix = np.transpose(np.array(meansubQlist)) #WE NEED THE POINTS AS COLLUMNS SO WE WILL TRANSPOSE 
        meansubPmatrix = np.array(meansubPlist) #WE NEED COORDINATES FROM P AS ROWS SO WE WILL NOT TRANSPOSE

        #CROSS CORRELATION MATRIX AND SVD
        W = meansubQmatrix.dot(meansubPmatrix)
        U,S,VH = np.linalg.svd(W,full_matrices = True)
        #ROTATION MATRIX AND TRANSLATION VECTOR
        R = U.dot(VH)
        t = comQ - R.dot(comP.T)

        #ACCUMULATED ROTATIONMATRIX AND TRANSLATION VECTOR
        if EK1 == 0:
            accumulatedrot = R
            accumulatedtrans = t
        else:
            accumulatedrot = accumulatedrot.dot(R)
            accumulatedtrans += t

        #WRITE OUTPUTS
        print("Rotation matrix:\n",R,file=f)
        print("Translation vector:\n",t,file=f)
        print("Accumulated Rotation Matrix:\n",accumulatedrot,file=f)
        print("Accumulated Translation Vector:\n",accumulatedtrans,file=f)
        #CALCULATE E
        suml2 = 0
        for point in qmatrix:
            qi = point
            pi = [0,0,0]
            shortest_distance = float('inf')
            #closest pi point
            for ppoint in pmatrix:
                distance = math.sqrt((qi[0]-ppoint[0])**2 + (qi[1]-ppoint[1])**2 +  (qi[2]-ppoint[2])**2)
                if distance<shortest_distance:
                    shortest_distance = distance
                    pi = ppoint
            a = qi-(R.dot(pi))-t
            l2 = 0
            for j in range(len(a)):
                l2 += a[j]*a[j]
            suml2 +=l2
        EK = suml2 / N
        #CALCULATE RE and EK-1
        if EK1!=0:
            RE = abs(100*(EK1-EK)/EK1)
            print("Relative Error Change:","%.2f" % RE+"%","\nDistance Error(E):","%.28f" % EK,file=f)
        EK1 = EK
        print("-------------------------------------------------------------------------\n",file=f)
        
    end = time.time()
    print("Time elapsed:",end-start,"seconds")
    f.close()
    return pmatrix



###MAIN START###
#READ THE INPUT
if not path.exists(Qfilepath):
    Qfilepath = input("Please Enter Q data path: ")
if not path.exists(Pfilepath):
    Pfilepath = input("Please Enter P data path: ")

qlist = []
plist = []
for line in open(Qfilepath,"r").readlines():
    if line != "\n":
        temp = line.rstrip("\n").split("\t")
        qlist.append(list(map(float,temp)))
for line in open(Pfilepath,"r").readlines():
    if line != "\n":
        temp = line.rstrip("\n").split("\t")
        plist.append(list(map(float,temp)))
#OUTPUT FORMATTER TO FIX FLOATS
float_formatter = "{:.10f}".format
np.set_printoptions(formatter={'float_kind':float_formatter})
#FUNCTIONS
bruteforcepmatrix = brute_force()
improvedpmatrix = my_imp()

#PLOTTING
if plotting:
    plotpoints(qlist,plist,bruteforcepmatrix,improvedpmatrix)