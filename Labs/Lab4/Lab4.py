import numpy as np

def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    N = np.zeros((Nmax,1))

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          N[count-1] = xstar
          return [N,ier]
        
       x0 = x1
       N[count-1] = x1

    xstar = x1
    ier = 1
    return [N, ier]

f = lambda x: (10/(x+4))**0.5
x0 = 1.5
tol = 10**(-10)

[N,ier] = fixedpt(f,x0,tol,50)

print(N)