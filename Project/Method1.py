import numpy as np
from scipy.integrate import quad
import math
import matplotlib.pyplot as plt
import time

def driver():
    # part a)

    gx = lambda x: np.e**x/np.sqrt(x) #actual function
    f = lambda x: np.e**x #top part of fucntion
    a2 = 0 # location of the discontinuity
    alpha = 0.5 # order of the polynomial on the bottom

    a=0
    b=1

    N = 5000 #Number of simpson nodes
    n = 10 #number of taylor nodes

    time1 = time.time()
    for i in range(100):
        I = simpsonWeakSing(a,b,f,N,a2,n,alpha)
    time2 = time.time()
    tim = (time2-time1)/100
    print(I)
    print('Time take to perform first method:',tim)



    time1 = time.time()
    for i in range(100):
        I2,err = quad(gx,a,b,limit=2)
    time2 = time.time()
    tim = (time2-time1)/100
    print(I2)
    print('Time take to perform Gaussian Quadrature:',tim)

    # Question 2
    actual =2.92531
    n = 1
    N= 100
    I = simpsonWeakSing(a,b,f,N,a2,n,alpha)
    print('Order 1 Taylor Expansion error:',np.abs(actual-I))

    n = 2
    I = simpsonWeakSing(a,b,f,N,a2,n,alpha)
    print('Order 2 Taylor Expansion error:',np.abs(actual-I))

    n=3
    I = simpsonWeakSing(a,b,f,N,a2,n,alpha)
    print('Order 3 Taylor Expansion error:',np.abs(actual-I))

    n=4
    I = simpsonWeakSing(a,b,f,N,a2,n,alpha)
    print('Order 4 Taylor Expansion error:',np.abs(actual-I))

    n=10
    I = simpsonWeakSing(a,b,f,N,a2,n,alpha)
    print('Order 10 Taylor Expansion error:',np.abs(actual-I))

    # Question 3



    
    
    

    
def Gx(fx,a,n,alpha,eval,lowbound):
    # fx is the top function
    # a is the position of the discontinuity
    # n is the order of the taylor expansion
    # eval is the point of evaluation
    # lowbound is the lower bound of the interval
    if eval == a:
        return 0
    Tn = taylor_approximation(fx,lowbound,n)
    Gx = (fx(eval)- Tn(eval))/((eval-a)**alpha)
    return Gx

def Sx(fx,a,n,alpha,eval,lowbound):
    # fx is the top function
    # a is the position of the discontinuity
    # n is the order of the taylor expansion
    # eval is the point of evaluation
    # lowbound is the lower bound of the interval
    
    Tn = taylor_approximation(fx,lowbound,n)
    Sx = Tn(eval)/((eval-a)**alpha)
    return Sx

def simpsonWeakSing(a,b,f,N,a2,n,alpha):
    # f is the top function
    # a is the lower bound of the interval
    # b is the upper bound
    # N is the order of the simpson approx
    # a2 is the position of the discontinuity
    # n is the order of the taylor expansion
    # alpha is the order of the bottom polynomial``
    h = (b-a)/N
    I1 = 0
    I2 = 0
    for i in range (1,int(N/2)):
        xj = a + 2*i*h
        I1 = I1 + 2*Gx(f,a2,n,alpha,xj,a)
    for i in range(0,int(N/2)):
        xj = a+(2*i+1)*h
        I2 = I2 + 4*Gx(f,a2,n,alpha,xj,a)
    I = (f(a) + I1 + I2 + f(b))*(h/3)

    I1 = 0
    I2 = 0
    for i in range (1,int(N/2)):
        xj = a + 2*i*h
        I1 = I1 + 2*Sx(f,a2,n,alpha,xj,a)
    for i in range(0,int(N/2)):
        xj = a+(2*i+1)*h
        I2 = I2 + 4*Sx(f,a2,n,alpha,xj,a)
    I = I + (f(a) + I1 + I2 + f(b))*(h/3)
    return I

def simpson(a,b,f,N):
    h = (b-a)/N
    I1 = 0
    I2 = 0
    for i in range (1,int(N/2)):
        xj = a + 2*i*h
        I1 = I1 + 2*f(xj)
    for i in range(0,int(N/2)):
        xj = a+(2*i+1)*h
        I2 = I2 + 4*f(xj)
    

    return (I1 + I2 + f(b))*(h/3)


def taylor_approximation(f, a, n):

    def factorial(x):
        if x == 0:
            return 1
        else:
            return x * factorial(x - 1)

    def taylor_term(x, k):
        return (f(a)**(k) / factorial(k)) * (x - a)**k

    def taylor_approx(x):
        result = 0
        for k in range(n + 1):
            result += taylor_term(x, k)
        return result

    return taylor_approx


driver()