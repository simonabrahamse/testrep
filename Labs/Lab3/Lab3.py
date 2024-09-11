import numpy as np
import matplotlib.pyplot as plt

def fx(x):
    fx = (x**2)*(x-1)
    return fx

def bisection(a,b,err):
    er = err +1
    # if fx(a)*fx(b) > 0:
    #     print('inputs are of the same sign')
    #     return 0
    maxiter = np.log(err/(b-a))
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

ans = bisection(0.5,2,0.001)
print(ans)
ans = bisection(-1,0.5,.01)
print(ans)
ans = bisection(-1,2,.01)
print(ans)