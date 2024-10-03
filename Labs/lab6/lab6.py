import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

def driver():
    h = 0.01*2.**(-np.arange(0,10))
    f = lambda x:np.cos(x)
    fprime1 = lambda s:(f(s+h)-f(s))/h
    fprime2 = lambda s: (f(s+h)-f(s-h))/(2*h)

    error1 = np.abs(fprime1(np.pi/2)+np.linspace(1,1,10))
    error2 = np.abs(fprime2(np.pi/2)+np.linspace(1,1,10))

    tol = 10e-10
    Nmax = 100

    # print(error1)
    # print(error2)

    f1 = lambda x: 4*x[0]**2+x[1]**2-4
    f2 = lambda x: x[0] + x[1] - np.sin(x[0]-x[1])

    x0 = np.array([1,0])

    x01 = np.array([0.1, 0.1, -0.1])

    [xstar, ier,its] = SlackerNewton(x0,tol,Nmax)
    print('root:',xstar)
    print(f1(xstar))
    print(f2(xstar))
    print('Number if iterations:',its)
    
    [xstar1, ier1,its1] = SlackerNewton1(x01,tol,Nmax)
    print('root:',xstar)
    print('Evaluation at root:',evalF(xstar1))
    print('number of iterations:',its1)

    X = np.array([0,1])
    print(evalJ2(X,0.0005))
    
    





def evalF(x): 

    F = np.zeros(2)
    
    F[0] = 4*x[0]**2+x[1]**2-4
    F[1] = x[0] + x[1] - np.sin(x[0]-x[1])
    return F
    
def evalJ(x): 

    
    J = np.array([[8*x[0], 2*x[1]], 
        [1-np.cos(x[0]-x[1]), 1+np.cos(x[0]-x[1])]]) 
        
    return J

def evalJ2(X,h):
    
    f1 = lambda x: 4*x[0]**2+x[1]**2-4
    

    J = np.array([[(f1(X[0]+h,X[1])-f1(X[0],X[1]))/h]])

    return J

def evalF1(x): 

    F = np.zeros(3)
    
    F[0] = 3*x[0]-math.cos(x[1]*x[2])-1/2
    F[1] = x[0]-81*(x[1]+0.1)**2+math.sin(x[2])+1.06
    F[2] = np.exp(-x[0]*x[1])+20*x[2]+(10*math.pi-3)/3
    return F
    
def evalJ1(x): 

    
    J = np.array([[3.0, x[2]*math.sin(x[1]*x[2]), x[1]*math.sin(x[1]*x[2])], 
        [2.*x[0], -162.*(x[1]+0.1), math.cos(x[2])], 
        [-x[1]*np.exp(-x[0]*x[1]), -x[0]*np.exp(-x[0]*x[1]), 20]])
    return J

def SlackerNewton(x0,tol,Nmax):
    count = 0
    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    J = evalJ(x0)
    Jinv = inv(J)
    denom = 2
    for its in range(Nmax):
        
        
        F = evalF(x0)
        
        if count % denom == 0:
            J = evalJ(x0)
            Jinv = inv(J)
            denom = denom*2
            
        count = count + 1
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier,its]
        
        x0 = x1

    xstar = x1
    ier = 1
    return[xstar,ier,its]   

def SlackerNewton1(x0,tol,Nmax):
    count = 0
    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    J = evalJ1(x0)
    Jinv = inv(J)
    denom = 2
    for its in range(Nmax):
        
        
        F = evalF1(x0)
        
        if count % denom == 0:
            J = evalJ1(x0)
            Jinv = inv(J)
            denom = denom*2
            
        count = count + 1
        x1 = x0 - Jinv.dot(F)
        if (norm(x1-x0) < tol):
            xstar = x1
            ier =0
            return[xstar, ier,its]
        
        x0 = x1

    xstar = x1
    ier = 1
    return[xstar,ier,its]


driver()