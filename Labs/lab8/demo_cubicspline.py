import matplotlib.pyplot as plt
import numpy as np
import math
from numpy.linalg import inv 
from numpy.linalg import norm


def driver():
    
    f = lambda x: np.exp(x)
    a = 0
    b = 1
    
    ''' number of intervals'''
    Nint = 3
    xint = np.linspace(a,b,Nint+1)
    yint = f(xint)

    ''' create points you want to evaluate at'''
    Neval = 100
    xeval =  np.linspace(xint[0],xint[Nint],Neval+1)

#   Create the coefficients for the natural spline    
    (M,C,D) = create_natural_spline(yint,xint,Nint)

#  evaluate the cubic spline     
    yeval = eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D)
    
    
    ''' evaluate f at the evaluation points'''
    fex = f(xeval)
        
    nerr = norm(fex-yeval)
    print('nerr = ', nerr)
    
    plt.figure()    
    plt.plot(xeval,fex,'ro-',label='exact function')
    plt.plot(xeval,yeval,'bs--',label='natural spline') 
    plt.legend
    plt.show()
     
    err = abs(yeval-fex)
    plt.figure() 
    plt.semilogy(xeval,err,'ro--',label='absolute error')
    plt.legend()
    plt.show()
    
def create_natural_spline(yint,xint,N):
    f = lambda x: np.exp(x)
#    create the right  hand side for the linear system
    b = np.zeros(N+1)
#  vector values
    h = np.zeros(N+1)
    i = 1
    h[0] = xint[i]-xint[i-1]  
    for i in range(1,N):
       h[i] = xint[i+1] - xint[i]
       b[i] = (yint[i+1]-yint[i])/h[i] - (yint[i]-yint[i-1])/h[i-1]

#  create the matrix A so you can solve for the M values
    A = np.zeros(N+1,N+1)
    v = np.zeros(N+1)
    v1 = np.zeros(N)
    
    for i in range(N+1):
        if i == 0 or i == N:
            v[i] = 1
        else: 
            v[i] = (h[i]+h[i-1])/3
    A = A + np.diag(v,0)

    for i in range(N):
        if i == N-1:
            v1[i] = 0
        else:
            v1[i] = (h[i])/6
    A = A + np.diag(v1,-1)

    for i in range(N):
        if i == 0:
            v1[i] = 0
        else:
            v1[i] = (h[i])/6
    A = A + np.diag(v1,1)


    

#  Invert A    
    Ainv = inv(A)

# solver for M    
    M  = Ainv*b
    
#  Create the linear coefficients
    C = np.zeros(N)
    D = np.zeros(N)
    for j in range(N):
       C[j] = f(xint[j])/h[j]-h[j]/6*M[j] # find the C coefficients
       D[j] = f(xint[j+1])/h[j]-h[j]/6*M[j+1] # find the D coefficients
    return(M,C,D)
       
def eval_local_spline(xeval,xi,xip,Mi,Mip,C,D):
# Evaluates the local spline as defined in class
# xip = x_{i+1}; xi = x_i
# Mip = M_{i+1}; Mi = M_i

    hi = xip-xi
    

    yeval = (xip-xeval)**3/(6*hi)*Mi+(xeval-xi)**3*Mip/(6*hi)+C*(xip-xeval)+D(xeval-xi)
    return yeval 
    
    
def  eval_cubic_spline(xeval,Neval,xint,Nint,M,C,D):
    
    yeval = np.zeros(Neval+1)
    
    for j in range(Nint):
        '''find indices of xeval in interval (xint(jint),xint(jint+1))'''
        '''let ind denote the indices in the intervals'''
        atmp = xint[j]
        btmp= xint[j+1]
        
#   find indices of values of xeval in the interval
        ind= np.where((xeval >= atmp) & (xeval <= btmp))
        xloc = xeval[ind]

# evaluate the spline
        yloc = eval_local_spline(xloc,atmp,btmp,M[j],M[j+1],C[j],D[j])
#   copy into yeval
        yeval[ind] = yloc

    return(yeval)
           
driver()               


