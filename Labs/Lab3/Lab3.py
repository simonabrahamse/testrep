import numpy as np
import matplotlib.pyplot as plt

def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]

def bisection(a,b,fx,err):
    er = err +1
    if fx(a)*fx(b) > 0:
        print('inputs are of the same sign')
        return
    maxiter = np.log(err/(b-a))
    i = 0
    while abs(er) >= err and i < 15:
        d = 0.5*(a+b)
        i = i+1
        fa = fx(a)
        fb = fx(b)
        fd = fx(d)
        if fa*fd > 0:
            a = d
            #fa = fd
            er = (d-b)/d
        else:
            b = d
            #fb = fd
            er = (d-a)/d
    return d

# Part a
f = lambda x : x**2*(x-1)
a = 0.5
b = 2
x1 = bisection(a,b,f,0.01)
print(x1)
# Part b
a = -1
b = 0.5
x2 = bisection(a,b,f,0.01)
print(x2)
# part c
a= -1
b = 2
x3 = bisection(a,b,f,0.01)
print(x3)

print('Problem 2')
#part a 
g = lambda x : (x-1)*(x-3)*(x-5)
a = 0
b= 2.4
e = 10**-5
print('part a: x =',bisection(a,b,g,e))
# part b
g = lambda x : (x-1)**2*(x-3)
a = 0
b = 2
print('part b: x =',bisection(a,b,g,e))
# part c
g = lambda x : np.sin(x)
a = 0
b = .1
print('part c1: x =',bisection(a,b,g,e))
a = 0.5
b = 3*np.pi/4
print('part c2: x =',bisection(a,b,g,e))

print('Problem 3')
x_0 = 1
tol = 10**-10
# part a 
# f = lambda x : x*(1+(7-x**5)/(x**2))**3
# print(f(7**(1/5)))

# print(fixedpt(f,x_0,tol,100))

# part b
# f = lambda x : x-(x**5-7)/(x**2)
# print(f(7**(1/5)))

# print(fixedpt(f,x_0,tol,100))

# part c
f = lambda x : x-(x**5-7)/(5*x**4)
print(f(7**(1/5)))
print('part c: [x,ier] =',fixedpt(f,x_0,tol,10000))

# part d
f = lambda x : x-(x**5-7)/(12)
print(f(7**(1/5)))
print('part d: [x,ier] =',fixedpt(f,x_0,tol,10000))
