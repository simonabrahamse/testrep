# import libraries
import numpy as np

def driver():

# use routines    
    f = lambda x: np.exp(x**2+7*x-30)-1
    fp = lambda x: (2*x+7) * np.exp(x**2+7*x-30)
    fpp = lambda x: (2*x+7)**2 * np.exp(x**2+7*x-30) + 2 * np.exp(x**2+7*x-30)
    a = 2
    b = 4.5
    Nmax = 1000

#    f = lambda x: np.sin(x)
#    a = 0.1
#    b = np.pi+0.1

    tol = 1e-7

    [p,pstar,info,it] = hybrid(f,fp,fpp,a,b,tol,Nmax)
    print('the approximate root is',pstar)
    print('the error message reads:',info)
    print('f(astar) =', f(pstar))
    print('Number of iterations required is:',it)




# define routines
def hybrid(f,fp,fpp,a,b,tol,Nmax):
    
    p = np.zeros(Nmax+1)
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      info   - error message
#            - info = 1 => Failed
#            - info = 0 == success

#     first verify there is a root we can find in the interval 

    fa = f(a)
    fb = f(b)
    if (fa*fb>0):
       info = 1
       astar = a
       p = f(astar)
       
       return [p,pstar,info,it]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      info =0
      p = f(astar)
      return [p,pstar,info,it]

    if (fb ==0):
      astar = b
      info = 0
      p = f(astar)
      return [p,pstar,info,it]

    count = 1
    d = 0.5*(a+b)
    test = 1 - (fp(d)**2-f(d)*fpp(d)/(fp(d))**2)
    while (np.abs(test) > 1):
      fd = f(d)
      
      if (fd ==0):
        astar = d
        info = 0
        return [p,pstar,info,it]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
      test = f(d)*fpp(d)/(fp(d)**2)
    
    p0 = d
    
    

    p = np.zeros(Nmax+1)
    p[0] = p0
    for it in range(Nmax):
        
        p1 = p0-f(p0)/fp(p0)
        p[it+1] = p1
        if (abs(p1-p0) < tol):
            pstar = p1
            info = 0
            it = count + it
            
            return [p,pstar,info,it]
        p0 = p1
    pstar = p1
    info = 1

    it = count + it

    return [p,pstar,info,it]


    
      
driver()               

