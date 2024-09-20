import numpy as np

def compute_order(x,xstar):
   """
   approximates the order of convergence given:
   x: array of iterate values
   xstar: fixed point/solution
   """
   diff1 = np.abs(x[1::]-xstar)

   diff2 = np.abs(x[0:-1]-xstar)

   fit  = np.polyfit(np.log(diff2.flatten()),np.log(diff1.flatten()),1)
   _lambda = np.exp(fit[1])
   alpha = fit[0]
   print(f"lambda is {_lambda}")
   return [_lambda,alpha]

def aitkens(x):
   xn = x[:-2]
   xn1 = x[1:-1]
   xn2 = x[2::]
   xhat = xn - (xn1-xn)**2/(xn2-2*xn1+xn)

   return xhat

def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    N = []

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          N.append(x1)
          return [np.array(N),xstar,ier]
        
       x0 = x1
       N.append(x1)

    xstar = x1
    ier = 1
    return [np.array(N), xstar, ier]

f = lambda x: (10/(x+4))**0.5
x0 = 1.5
tol = 10**(-10)

[N, xstar, ier] = fixedpt(f,x0,tol,50)

xhat = aitkens(N)

[_lambda,alpha] = compute_order(xhat,xstar)
print(_lambda)
print(alpha)