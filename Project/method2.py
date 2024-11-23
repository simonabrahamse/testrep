import numpy as np
from scipy.integrate import quad

def driver():
    gx1 = lambda x: x*np.log(x)
    gx2 = lambda x: -x*np.log(x)
    gx3 = lambda x: np.cos(x)*np.log(np.sin(x))

    nodes2 = np.array([1.825748064736159,-1.325748064736159])
    nodes6 = np.array([4.967362978287758,-1.620501504859126*10,2.585153761832639*10,
                       -2.222599466791883*10,9.930104998037539,-1.817995878141594])
    nodes10  = np.array([7.832432020568779,-4.565161670374749*10,1.452168846354677*100,
                         -2.901348302886379*100,3.870862162579900*100,-3.523821383570681*100,
                         2.172421547519342*100,-8.707796087382991*10,2.053584266072635*10, -2.166984103403823])

    a = 0
    b = np.pi
    N = 10
    I2,err = quad(gx3,a,b)

    I = krlog(gx1,a,b,N,nodes10)

    print(I)
    print(I2)


def krlog(gx,a,b,N,nodes):
    h = (b-a)/N
    I1 = 0
    for i in range(2,N):
        xj = a+i*h
        I1 = I1 + 2*gx(xj)

    I1 += (h/2)*(I1)
    I2 = 0
    I3 = 0
    for i in range(1,N+1):
        if i > 0:
            # print(i)
            I2 += h*nodes[i-1]*gx((i)*h)
            I3 += h*nodes[i-1]*gx(b+(i)*h)
    # for i in range(-N,0):
    #     if i < 0:
    #         if b+i*h != 0:
    #             I3 += gx(b+(i)*h)
    return I1+I2+I3

driver()