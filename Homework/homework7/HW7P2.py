import numpy as np
import matplotlib.pyplot as plt

def driver():
    f = lambda x: 1/(1+100*x**2)
       
    for i in range(3,21):
        xfine, yfine, poly, x, y = lagrange(f,i)
        plt.plot(x, y, 'o')
        plt.plot(xfine, poly)
        plt.plot(xfine, yfine)
        plt.show()



def p(xeval,w,x,y):
    n = np.sum(w/(xeval-x)*y)
    d = np.sum(w/(xeval-x))
    c = n/d
    return c


def lagrange(f,N):
    x= np.linspace(-1,1,N)
    y = f(x)

    w = np.ones(N)
    for j in range(N):
        for i in range(N):
            if i != j:
                w[j] /= (x[j]-x[i])

    xfine = np.linspace(-1, 1, 1001)
    yfine = f(xfine)
    # poly = np.zeros(1001)
    # for i in range(1001):
    #       poly[i] = p(xfine[i])
    poly = np.array([p(xi,w,x,y) for xi in xfine])

    return xfine, yfine, poly, x, y


driver()