import numpy as np
import scipy
import scipy.integrate

def driver():

    f = lambda t: t**-2*np.cos(t**-1)/(t**-3)
    a = 0
    b= 1
    Ns = 4


    Iact,err = scipy.integrate.quad(f,a,b,epsabs=10**-6)

    I = simpson(a,b,f,Ns)
    print('value from scipy:',Iact)
    print('Simpsons estimate:',I)
    print('Simpsons error:',np.abs(Iact-I))



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


driver()