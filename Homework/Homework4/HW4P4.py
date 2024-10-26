# import libraries
import numpy as np
        
def driver():
#f = lambda x: (x-2)**3
#fp = lambda x: 3*(x-2)**2
#p0 = 1.2

  f = lambda x: np.exp(3*x) -27*x**6+27*x**4*np.exp(x)-9*x**2*np.exp(2*x)
  fp = lambda x: 3 * np.exp(3 * x) + (-18 * x**2 - 18 * x) * np.exp(2 * x) + (27 * x**4 + 108 * x**3) * np.exp(x) - 162 * x**5
  fpp = lambda x: 9 * np.exp(3 * x) + (-36 * x**2 - 72 * x - 18) * np.exp(2 * x) + (27 * x**4 + 216 * x**3 + 324 * x**2) * np.exp(x) - 810 * x**4

  g = lambda x: f(x)/fp(x)
  gp = lambda x: (fp(x)**2-f(x)*fpp(x))/(fp(x)**2)
  p0 = 3

  Nmax = 1000
  tol = 10e-5

  (p1,pstar1,info1,it1) = newton(f,fp,p0,tol, Nmax)
  print('the approximate root is', '%16.16e' % pstar1)
  print('the error message reads:', '%d' % info1)
  print('Number of iterations:', '%d' % it1)

  p1 = p1[0:it1]

  [_lambda1,alpha1] = compute_order(p1,pstar1)
  print('Order of convergence:',alpha1)

  (p2,pstar2,info2,it2) = newton(g,gp,p0,tol, Nmax)
  print('the approximate root is', '%16.16e' % pstar2)
  print('the error message reads:', '%d' % info2)
  print('Number of iterations:', '%d' % it2)
  
  p2 = p2[0:it2]

  [_lambda2,alpha2] = compute_order(p2,pstar2)
  print('Order of convergence:',alpha2)

  (p3,pstar3,info3,it3) = newton_mod(f,fp,p0,tol, Nmax,2)
  print('the approximate root is', '%16.16e' % pstar3)
  print('the error message reads:', '%d' % info3)
  print('Number of iterations:', '%d' % it3)

  p3 = p3[0:it3]

  [_lambda3,alpha3] = compute_order(p3,pstar3)
  print('Order of convergence:',alpha3)








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


def newton_mod(f,fp,p0,tol,Nmax,m):
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
      p1 = p0-m*f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

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



driver()