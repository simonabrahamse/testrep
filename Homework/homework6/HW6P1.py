import numpy as np
import math
import time
from numpy.linalg import inv
from numpy.linalg import norm


def driver():
    x0 = np.array([1, 1]) 
  
    Nmax = 100
    tol = 1e-10
    # Newton's Method
    t = time.time()
    for j in range(50):
        [xstar, ier, its] = Newton(x0, tol, Nmax)
    elapsed = time.time() - t
    print('Newton Result:', xstar)

    print('Newton: Error message:', ier)
    print('Newton: Time per run:', elapsed / 50)
    print('Newton: Iterations:', its)

   # Lazy Newton Method

    x0 = [1,-1]
    t = time.time()
    for j in range(20):
        [xstar, ier, its] = LazyNewton(x0, tol, Nmax)
    elapsed = time.time() - t
    print('Lazy Newton Result:', xstar)
    print('Lazy Newton: Error message:', ier)
    print('Lazy Newton: Time per run:', elapsed / 20)
    print('Lazy Newton: Iterations:', its)

   

    # Broyden's Method
    x0 = [0,0]
    t = time.time()
    for j in range(20):
        [xstar, ier, its] = Broyden(x0, tol, Nmax)
    elapsed = time.time() - t
    print('Broyden Result:', xstar)
    print('Broyden: Error message:', ier)
    print('Broyden: Time per run:', elapsed / 20)
    print('Broyden: Iterations:', its)


    


###########################################################

#functions:


def evalF(x):
    F = np.zeros(2)
    F[0] = x[0]**2 + x[1]**2 - 4  # f(x, y)\
    F[1] = math.exp(x[0]) + x[1] - 1  # g(x, y)
    return F


def evalJ(x):
    J = np.array([[2 * x[0], 2 * x[1]],
                  [math.exp(x[0]), 1]])
    return J

###########################################################

#Newton:

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


###########################################################

#LazyNewton:


def LazyNewton(x0, tol, Nmax):
    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):
        
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


###########################################################

#Broyden:


# def Broyden(x0, tol, Nmax):
#     A0 = evalJ(x0)
#     F = evalF(x0)
#     A = np.linalg.inv(A0)
#     s = -A.dot(F)
#     xk = x0 + s
#     for its in range(Nmax):
#         w = F
#         F = evalF(xk)
#         y = F - w
#         z = -A.dot(y)
#         p = np.dot(s, z)
#         u = np.dot(s, A)
#         A += np.outer(s + z, u) / p
#         s = -A.dot(F)
#         xk = xk + s
#         if norm(s) < tol:
#             xstar = xk
#             ier = 0
#             return [xstar, ier, its]
#     xstar = xk
#     ier = 1
#     return [xstar, ier, its]

def Broyden(x0, tol, Nmax):
    
    n = len(x0)  # Dimension of the system
    x = x0       # Current guess
    B = np.eye(n)  # Initial Jacobian approximation (identity matrix)
    
    for k in range(Nmax):
        # Evaluate the function at the current guess
        Fx = evalF(x)
        
        # Check if the function is close enough to zero
        if np.linalg.norm(Fx) < tol:
            ier = 0
            return [x,ier,k]
        
        # Solve B_k Δx_k = -F(x_k)
        delta_x = np.linalg.solve(B, -Fx)
        
        # Update x
        x_new = x + delta_x
        
        # Evaluate the function at the new guess
        Fx_new = evalF(x_new)
        
        # Compute y_k = F(x_new) - F(x)
        y = Fx_new - Fx
        
        # Broyden's update for the Jacobian approximation
        delta_x_T = delta_x.reshape(-1, 1)  # Column vector
        B = B + (y - B @ delta_x).reshape(-1, 1) @ delta_x_T.T / np.dot(delta_x, delta_x)
        
        # Update x
        x = x_new
    
    print('Maximum number of iterations reached.')
    ier = 1
    return [x,ier,k]


if __name__ == '__main__':
    driver()

