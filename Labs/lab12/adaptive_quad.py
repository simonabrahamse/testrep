# get lgwts routine and numpy
import numpy as np
# adaptive quad subroutines
# the following three can be passed
# as the method parameter to the main adaptive_quad() function

def trap(a,b,f,N):
    h = (b-a)/N
    I1 = 0
    for i in range(1,N):
        xj = a+i*h
        I1 = I1 + 2*f(xj)
    return (h/2)*(f(a)+f(b)+I1)
    


def simpson(a,b,f,N):
    h = (b-a)/N
    I1 = 0
    I2 = 0
    for i in range (1,int(N/2)):
        xj = a + 2*i*h
        I1 = I1 + 2*f(xj)
    for i in range(0,int(N/2+1)):
        xj = a+(2*i+1)*h
        I2 = I2 + 4*f(xj)
    

    return (f(a) + I1 + I2 + f(b))*(h/3)

def eval_composite_trap(M,a,b,f):
  """
  put code from prelab with same returns as gauss_quad
  you can return None for the weights
  """

def eval_composite_simpsons(M,a,b,f):
  """
  put code from prelab with same returns as gauss_quad
  you can return None for the weights
  """

def eval_gauss_quad(M,a,b,f):
  """
  Non-adaptive numerical integrator for \int_a^b f(x)w(x)dx
  Input:
    M - number of quadrature nodes
    a,b - interval [a,b]
    f - function to integrate
  
  Output:
    I_hat - approx integral
    x - quadrature nodes
    w - quadrature weights

  Currently uses Gauss-Legendre rule
  """
  x,w = lgwt(M,a,b)
  I_hat = np.sum(f(x)*w)
  return I_hat,x,w

def adaptive_quad(a,b,f,tol,M,method):
  """
  Adaptive numerical integrator for \int_a^b f(x)dx
  
  Input:
  a,b - interval [a,b]
  f - function to integrate
  tol - absolute accuracy goal
  M - number of quadrature nodes per bisected interval
  method - function handle for integrating on subinterval
         - eg) eval_gauss_quad, eval_composite_simpsons etc.
  
  Output: I - the approximate integral
          X - final adapted grid nodes
          nsplit - number of interval splits
  """
  # 1/2^50 ~ 1e-15
  maxit = 50
  left_p = np.zeros((maxit,))
  right_p = np.zeros((maxit,))
  s = np.zeros((maxit,1))
  left_p[0] = a; right_p[0] = b;
  # initial approx and grid
  s[0],x,_ = method(M,a,b,f)
  # save grid
  X = []
  X.append(x)
  j = 1
  I = 0
  nsplit = 1
  while j < maxit:
    # get midpoint to split interval into left and right
    c = 0.5*(left_p[j-1]+right_p[j-1]);
    # compute integral on left and right spilt intervals
    s1,x,_ = method(M,left_p[j-1],c,f); X.append(x)
    s2,x,_ = method(M,c,right_p[j-1],f); X.append(x)
    if np.max(np.abs(s1+s2-s[j-1])) > tol:
      left_p[j] = left_p[j-1]
      right_p[j] = 0.5*(left_p[j-1]+right_p[j-1])
      s[j] = s1
      left_p[j-1] = 0.5*(left_p[j-1]+right_p[j-1])
      s[j-1] = s2
      j = j+1
      nsplit = nsplit+1
    else:
      I = I+s1+s2
      j = j-1
      if j == 0:
        j = maxit
  return I,np.unique(X),nsplit

