import numpy as np
import numpy.linalg as la
from numpy.linalg import inv
from numpy.linalg import norm
import matplotlib.pyplot as plt

def driver(): 

    f = lambda x: 1/(1+(10*x)**2)
    
    N = 10
    a = -1
    b = 1
    xint = np.zeros(N)
    for j in range(N):
        h=2/(N-1)
        xj = -1+(j)*h
        xint[j] = xj

    

    yint = f(xint)
    

    
#     ''' Create interpolation nodes'''
#     xint = np.linspace(a,b,N+1)
# #    print('xint =',xint)
#     '''Create interpolation data'''
#     yint = f(xint)
# #    print('yint =',yint)
    
    ''' Create the Vandermonde matrix'''
    V = Vandermonde(xint,N-1)
#    print('V = ',V)

    ''' Invert the Vandermonde matrix'''    
    Vinv = inv(V)
#    print('Vinv = ' , Vinv)
    
    ''' Apply inverse to rhs'''
    ''' to create the coefficients'''
    coef = Vinv @ yint
    
#    print('coef = ', coef)

# No validate the code
    Neval = 1000    
    xeval = np.linspace(a,b,Neval+1)
    yeval = eval_monomial(xeval,coef,N-1,Neval)

# exact function
    yex = f(xeval)
    
    err =  norm(yex-yeval) 
    print('err = ', err)
    

    plt.plot(xeval,yeval)
    plt.plot(xeval,yex)
    plt.show()
    print(xint)
    return

def  eval_monomial(xeval,coef,N,Neval):

    yeval = coef[0]*np.ones(Neval+1)
    
#    print('yeval = ', yeval)
    
    for j in range(1,N+1):
      for i in range(Neval+1):
#        print('yeval[i] = ', yeval[i])
#        print('a[j] = ', a[j])
#        print('i = ', i)
#        print('xeval[i] = ', xeval[i])
        yeval[i] = yeval[i] + coef[j]*xeval[i]**j

    return yeval

   
def Vandermonde(xint,N):

    V = np.zeros((N+1,N+1))
    
    ''' fill the first column'''
    for j in range(N+1):
       V[j][0] = 1.0

    for i in range(1,N+1):
        for j in range(N+1):
           V[j][i] = xint[j]**i

    return V     

driver()    
