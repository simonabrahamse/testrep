import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

def driver():
    f = lambda x: 1/(1+100*x**2)

    for i in range(3,21):
        x,y,xfine,yfine,p = monomial(f,i)
        plt.plot(x,y,'o')
        plt.plot(xfine,yfine)
        plt.plot(xfine,p)
        plt.show()

def monomial(f,N):
    x = np.linspace(-1,1,N)
    y = f(x)


    V = np.vander(x,increasing=True)

    c=np.linalg.solve(V,y)

    xfine = np.linspace(-1,1,1001)
    yfine = f(xfine)

    p = np.polyval(c[::-1],xfine)

    return x,y,xfine,yfine,p

driver()