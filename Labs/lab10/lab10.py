import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import math
from scipy.integrate import quad

def driver():
    w=lambda x: 1.
    f = lambda x: 1/(1+x**2)
    a=-1
    b=1
    n=2
    
    top = lambda x: eval_top(j,f,w,eval_legendre(j,x))
    bottom = lambda x: eval_bottom(j,w,eval_legendre(j,x))

    aj = quad(top,a,b)/quad(bottom,a,b)

    quad()




def eval_legendre(n,x):
   phi = np.zeros(n+1)
   phi[0] = 1
   phi[1] = x
   for i in range(2,n+1):
      phi[i] = 1/(n+1)*((2*n+1)*x*phi[n-1]-n*phi[n-2])
   
   return phi

def eval_top(j,f,w,phi,x):
    top = phi[j-1] * f(x)*w(x)
    return top


def eval_bottom(j,w,phi,x):
    bottom = phi[j-1]**2*w(x)
    return bottom



driver()

