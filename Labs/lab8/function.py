import numpy as np






def subroutine(x1,fx1,x2,fx2,alpha):

    f = lambda x: fx1 + (fx2-fx1)/(x2-x2)*(x-x1)

    return f(alpha)