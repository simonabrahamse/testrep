import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

def driver():
    Nmax =100
    tol = 1e-10
    x0 = np.array([1,1,1])
    f = lambda x: x[0]**2+4*x[1]**2+4*x[2]**2 - 16
    fx = lambda x: 2*x[0]
    fy = lambda x: 8*x[1]
    fz = lambda x: 8*x[2]

    j = 0
    X = np.zeros(Nmax)

    while j<Nmax:
        
        F =f(x0)
        Fx = fx(x0)
        Fy = fy(x0)
        Fz = fz(x0)
        x = x0[0]
        y= x0[1]
        z= x0[2]
        d = F/(fx(x0)**2+fy(x0)**2+fz(x0)**2)
        x1 = np.array([x-d*Fx,y-d*Fy,z-d*Fz])
        X[j] = norm(x0-x1)
        
        if(norm(x0-x1)<tol):
            xstar = x1
            print('Number of iterations taken:',j)
            print(xstar)
            print(f(x1))
            its = j
            j=Nmax
        j = j+1
        x0 = x1

    X = X[0:its+1]
    [_lambda3,alpha3] = compute_order(X,xstar)
    print('Order of convergence:',alpha3)


    

def compute_order(x,xstar):
   """
   approximates the order of convergence given:
   x: array of iterate values
   xstar: fixed point/solution
   """
   diff1 = x[1::]

   diff2 = x[0:-1]

   fit  = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)
   _lambda = np.exp(fit[1])
   alpha = fit[0]
   print(f"lambda is {_lambda}")
   return [_lambda,alpha]
    

        
           


if __name__ == '__main__':
    # run the drivers only if this is called from the command line
    driver()