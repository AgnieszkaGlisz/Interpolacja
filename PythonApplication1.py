import numpy as np
import matplotlib.pyplot as plt
import copy
import pandas as pd

def load(name):
    x = pd.read_csv(name, header=None, usecols=[0]).values
    y = pd.read_csv(name, header=None, usecols=[1]).values
    l = x.size - 1

    return x, y, l

def lagrange(x, y, x0, n):
    Y0 = 0
    for i in range (0, n):
        prod = y[i]
        for j in range (0, n):   
            if i!= j:
                prod *=(x0-x[j])/(x[i] - x[j])
        Y0 += prod
    return Y0

def createLagrangeFigure(xp, yp, n, distance):
    nodeNr = int((xp[n])/distance)
    x = [] 
    y = []  
    x[0] = xp[0]
    y[0] = yp[0]
    
    for i in range (1, nodeNr):
           x[i] = x[i-1]+distance
           y[i] = lagrange(xp, yp, x[i], n)

    #plt.plot(xp, yp, 'navy', label = 'given data', linewidth=1.0)
    #plt.plot(x, y, 'magenta', label = 'interpolation')
    #plt.yscale('symlog')
    #plt.grid(True)
    #plt.legend() 
    #plt.show()

def pivoting(N, b, A):
    x_ref = np.linalg.solve(A, b)
    L = np.eye(N)
    P = np.eye(N)
    U = copy.deepcopy(A)
    

    for k in range(0,  N-1):
        pivot=max(abs(U[k:N,k]))
        for j in range (k, N):
            if(abs(U[j,k])==pivot):
                ind=j
                break

        U[[k,ind],k:N]=U[[ind,k],k:N]
        L[[k,ind],0:k]=L[[ind,k],0:k]
        P[[k,ind],:]=P[[ind,k],:]


        for j in range(k+1, N):
            L[j][k] = U[j][k]/U[k][k]
            U[j][k:N] = U[j][k:N] - L[j][k]*U[k][k:N]

    x = np.linalg.solve(U, np.linalg.solve(L, P.dot(b)))

    return x

def spline(xp, yp, n, distance):
    results = []
    wspA = np.zeros((4*n, 4*n))
    b = np.zeros(4*n)
    r = 0
    h = xp[1]
    for i in range(0,n):
        c = 4*i
        wspA[r,c] = 1
        b[r] = yp[i]
        r+=1
        
        wspA[r, c] = 1
        wspA[r, c + 1] = h
        wspA[r, c + 2] = h **2
        wspA[r, c + 3] = h **3
        b[r] = yp[i+1]
        r+=1

        if i > 0:
            wspA[r, c + 1] = -1
            wspA[r, c - 1] = 3 * h **2
            wspA[r, c - 2] = 2 * h
            wspA[r, c - 3] = 1
            r+=1

            wspA[r, c + 2] = -2
            wspA[r, c - 1] = 6 * h
            wspA[r, c - 2] = 2
            r+=1

    wspA[r, 2] = 2
    r+=1

    wspA[r, 4 * (n-1) + 2] = 2
    wspA[r, 4 * (n-1) + 3] = 6 * (xp[n] - xp[n-1])

    results = pivoting(4*n, b, wspA)
    xn = []
    yn = []

    for i in range (0, n):
        r = 4*i
        for j in range (int(xp[i]), int(xp[i+1]), distance):
            h = j - xp[i]
            f = results[r] + results[r+1]*h + results[r+2]*(h**2) + results[r+3]*(h**3)
            xn.append(j)
            yn.append(f)

    return xn, yn

def createSplineFigure(xp,yp,size,s):
    x, y = spline(xp, yp, size, s)

    plt.plot(xp, yp, 'navy', label = 'given data', linewidth=1.8)
    plt.plot(x, y, 'plum', label = 'interpolation')
    plt.grid(True)
    plt.legend() 
    plt.show()


xp, yp, size = load('Hel_yeah.csv')#load("GlebiaChallengera.csv")
#"Hel_yeah.csv")
createLagrangeFigure(xp, yp, size, 1000)
createSplineFigure(xp,yp, size, 100)
2+2