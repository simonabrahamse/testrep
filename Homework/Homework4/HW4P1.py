import numpy as np
from scipy import special
import matplotlib.pyplot as plt

def driver():

    Ti = 20
    Ts = -15
    alpha = .138*10**(-6)
    tol = 10**(-13)

# a = 0
# b = 0.15

    txt = np.zeros(100)

    t = 60*24*60

    x = np.linspace(0,0.15,100)

#Determine approximately where the root lies
    # Tactual = (Ti-Ts)*special.erf(x/(2*np.sqrt(alpha*t)))+Ts

    # plt.plot(x,Tactual)
    # plt.plot(x,txt)
    # plt.xlabel('depth (m)')
    # plt.ylabel('Temperature after 60 days (C)')
    # plt.show()

# use routines    
    f = lambda x: (Ti-Ts)*special.erf(x/(2*np.sqrt(alpha*t)))+Ts
    fp = lambda x: (Ti-Ts)*(2/np.sqrt(np.pi))*np.exp(x/(2*np.sqrt(alpha*t)))
    a = 0
    b = 0.15

#    f = lambda x: np.sin(x)
#    a = 0.1
#    b = np.pi+0.1

    tol = 1e-13

    [astar,ier,count] = bisection(f,a,b,tol)
    print('the approximate root is',astar)
    print('the error message reads:',ier)
    print('f(astar) =', f(astar))
    print('number of iterations:',count)

    

    Nmax = 1000
    p0 = 0.08

    [p,pstar,info,it] = newton(f,fp,p0,tol, Nmax)
    print('the approximate root is', '%16.16e' % pstar)
    print('the error message reads:', '%d' % info)
    print('Number of iterations:', '%d' % it)




# define routines
def bisection(f,a,b,tol):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 
    count = 1
    fa = f(a)
    fb = f(b)
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier,count]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier,count]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier,count]

    count = 0
    d = 0.5*(a+b)
    while (abs(d-a)> tol):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier,count]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    return [astar, ier,count]
      
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

driver()               


 