import numpy as np
import matplotlib.pyplot as plt
    
def driver():


     f1 = lambda x: x-4*np.sin(2*x)-3
     X = np.linspace(-2,5,1000)
     y1 = f1(X)
     plt.plot(X,y1)
     plt.show()    

     Nmax = 1000
     tol = 1e-10
     #to find the zero of f1, we need to find the fixed point of x-f1
     f2 = lambda x: -np.sin(2*x)+5*x/4-3/4
     

# test f1 '''
     x0 = -0.5
     [xstar,ier,count] = fixedpt(f2,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('value at f(xstar):',f1(xstar))
     print('Error message reads:',ier)
     print('Number of iterations required:',count)
     x0 = 3
     [xstar,ier,count] = fixedpt(f2,x0,tol,Nmax)
     print('the approximate fixed point is:',xstar)
     print('value at f(xstar):',f1(xstar))
     print('Error message reads:',ier)
     print('Number of iterations required:',count)


# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier,count]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier,count]
    

driver()