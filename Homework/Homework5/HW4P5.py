import numpy as np
import matplotlib.pyplot as plt
import pylab

def driver():
   
   f = lambda x: x**6-x-1
   fp = lambda x: 6*x**5-1
   x0 = 2
   x1 = 1
   tol = 10e-10
   Nmax = 1000

   [p1,pstar,info,it] = newton(f,fp,x0,tol,Nmax)
   print('pstar:',pstar)
   print('f(pstar) =',f(pstar))

#    print('iterations used:',it)
#    print('error message:',info)
   
   [pstar2,ier,p2,count] = secant(f,x0,x1,tol,Nmax)
   print('pstar:',pstar2)
   print('f(pstar) =',f(pstar2))
   print('error message:',ier)
   print('Number of iterations:',count)

   X1 = np.abs(p1[1:it +2] - pstar)
   Y1 = np.abs(p1[0:it+1] - pstar)

   fig = plt.figure()
   ax = fig.add_subplot(1,1,1)
   line, = ax.plot(X1,Y1,label="Original")
   plt.title("Rate of convergence graph")
   plt.xlabel("Xn error")
   plt.ylabel("Xn+1 error")
   
   ax.set_xscale('log')
   ax.set_yscale('log')
   pylab.show()

   X2 = np.abs(p2[1:count] - pstar2)
   Y2 = np.abs(p1[0:count-1] - pstar2)
   

   fig = plt.figure()
   ax = fig.add_subplot(1,1,1)
   line, = ax.plot(X2,Y2,label="Original")
   plt.title("Rate of convergence graph")
   plt.xlabel("Xn error")
   plt.ylabel("Xn+1 error")
   
   ax.set_xscale('log')
   ax.set_yscale('log')
   pylab.show()








def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1)
  p[0] = p0
  for it in range(Nmax):
      
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          
          pstar = p1
          info = 0
          
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

def secant(f,x0,x1,tol,Nmax):
    count = 0
    p = np.zeros(Nmax+1)
    p[0] = x0
    if f(x1) - f(x0) == 0:
        ier = 1
        pstar = x1
        return [pstar,ier,p,count]
    count =1
    while count < Nmax:
        
        x2 = x1 -f(x1)*(x1-x0)/(f(x1)-f(x0))
        if np.abs(x2-x1) < tol:
            pstar = x2
            ier=0
            
            return [pstar,ier,p,count]
        
        # update
        x0 = x1
        x1 = x2
        p[count] = x2
        count = count + 1
        
        if np.abs(f(x1)-f(x0))==0:
            ier = 1
            pstar = x2
            return [pstar,ier,p,count]
    pstar = x2
    ier = 1
    return [pstar,ier,p,count]


    
   


driver()               
