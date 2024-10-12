import matplotlib.pyplot as plt
import numpy as np
import math
import time
from numpy.linalg import inv
from numpy.linalg import norm

def driver():
    Nmax = 100
    x0= np.array([0,-1,-1])
    tol = 1e-6
    t = time.time()
    for j in range(50):
        [xstar,gval,ier] = SteepestDescent(x0,tol,Nmax)
    elapsed = time.time() - t

    print("Steepest Descent: the solution is:",xstar)
    print("Steepest Descent: g evaluated at this point is ", gval)
    # print("Steepest Descent: error is", ier)
    print('Steepest Descent: Time per run:', elapsed / 50)
    
    t = time.time()
    for j in range(50):
        [xstar,ier,its] = Newton(x0,tol,Nmax)
    elapsed = time.time() - t
    print('Newton: the error message reads:',ier)
    print('Newton: number of iterations is:',its)
    # print('Function evaluated at xstar =',evalF(xstar))
    print('Newton: Time per run:', elapsed / 50)

    tol_comb = 5e-2
    t = time.time()
    for j in range(50):
        [xstar, ier, its] = Combined(x0, tol_comb, tol, Nmax)
    elapsed = time.time() - t
    tol_comb = 5e-2

    print('Combined: the solution is:', xstar)
    print('Combined: Error =', ier)
    # print('Function evaluated at xstar =',evalF(xstar))
    print('Combined: Time per run:', elapsed / 50)


###########################################################

#functions:

def evalF(x):

    F = np.zeros(3)
    F[0] = x[0] +np.cos(x[0]*x[1]*x[2])-1.
    F[1] = (1 -x[0])**(0.25) + x[1] +0.05*x[2]**2 -0.15*x[2]-1
    F[2] = -x[0]**2-0.1*x[1]**2 +0.01*x[1]+x[2] -1

    return F


def evalJ(x):
    J = np.array([[1 - x[1] * x[2] * np.sin(x[0] * x[1] * x[2]), -x[0] * x[2] * np.sin(x[0] * x[1] * x[2]), -x[0] * x[1] * np.sin(x[0] * x[1] * x[2])],
                [-0.25 * (1 - x[0])**(-0.75) if (1 - x[0]) > 0 else float('nan'), 1, 0.1 * x[2] - 0.15],
                [-2 * x[0], -0.2 * x[1] + 0.01, 1]
               ])
    return J


def evalg(x):

    F = evalF(x)
    g = np.dot(F,F)
    return g


def eval_gradg(x):
    F = evalF(x)
    J = evalJ(x)  

    gradg = 2*J.T.dot(F)
    return gradg



###############################

### steepest descent code 


def SteepestDescent(x,tol,Nmax):
    for its in range(Nmax):
        g1 = evalg(x)
        z = eval_gradg(x)
        z0 = norm(z)

        if z0 == 0:
            print("zero gradient")

        z = z/z0
        alpha1 = 0
        alpha3 = 1
        dif_vec = x - alpha3*z
        g3 = evalg(dif_vec)


        while g3>=g1:
            alpha3 = alpha3/2
            dif_vec = x - alpha3*z
            g3 = evalg(dif_vec)

           

        if alpha3<tol:
            print("no likely improvement")
            ier = 0
            return [x,g1,ier]
  
        alpha2 = alpha3/2
        dif_vec = x - alpha2*z
        g2 = evalg(dif_vec)

        h1 = (g2 - g1)/alpha2
        h2 = (g3-g2)/(alpha3-alpha2)
        h3 = (h2-h1)/alpha3

        alpha0 = 0.5*(alpha2 - h1/h3)
        dif_vec = x - alpha0*z
        g0 = evalg(dif_vec)


        if g0<=g3:
            alpha = alpha0
            gval = g0

        else:
            alpha = alpha3
            gval =g3

        x = x - alpha*z

        if abs(gval - g1)<tol:
            ier = 0
            return [x,gval,ier]

    print('max iterations exceeded')    
    ier = 1        
    return [x,g1,ier]


def Combined(x0,tol_comb,tol_final,Nmax):
    [x_steep, g_steep, ier_steep] = SteepestDescent(x0,tol_comb,Nmax)
    [x_newton, ier_newton, its_newton] = Newton(x_steep, tol_final, Nmax)

    its = its_newton
    return [x_newton, 0, its]

def Newton(x0, tol, Nmax):
    for its in range(Nmax):
        J = evalJ(x0)
        Jinv = inv(J)
        F = evalF(x0)
        x1 = x0 - Jinv.dot(F)
        if norm(x1 - x0) < tol:
            xstar = x1
            ier = 0
            return [xstar, ier, its]
        x0 = x1
    xstar = x1
    ier = 1
    return [xstar, ier, its]


if __name__ == '__main__':
  # run the drivers only if this is called from the command line
  driver()        