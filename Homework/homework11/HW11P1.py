import numpy as np
import scipy
import scipy.integrate

def driver():

    f = lambda x: 1/(1+x**2)
    a = -5
    b= 5
    Ns = 108
    NT = 1291

    Iact,err = scipy.integrate.quad(f,a,b,epsabs=10**-6)

    I = simpson(a,b,f,Ns)
    I1 = trap(a,b,f,NT)
    print('value from scipy:',Iact)
    print('Simpsons error:',np.abs(Iact-I))
    print('Trapezoidal error:', np.abs(I1-Iact))



def trap(a,b,f,N):
    h = (b-a)/N
    I1 = 0
    for i in range(1,N):
        xj = a+i*h
        I1 = I1 + 2*f(xj)
    return (h/2)*(f(a)+f(b)+I1)
    


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
    

    return (f(a) + I1 + I2 + f(b))*(h/3)


driver()