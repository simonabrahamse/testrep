import numpy as np

def fx(x):
    y=np.e**x
    return y-1

x = 9.999999995000000*10**(-10)

y= fx(x)
print(y)

f = lambda x: x + 1/2*x**2 

print(f(x))