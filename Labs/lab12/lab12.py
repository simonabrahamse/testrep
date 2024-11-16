import numpy as np

def driver():

    f = lambda x: x**2
    a = -1
    b= 1
    N = 20

    I = simpson(a,b,f,N)
    I1 = trap(a,b,f,N)
    print(I)
    print(I1)



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
    for i in range(0,int(N/2+1)):
        xj = a+(2*i+1)*h
        I2 = I2 + 4*f(xj)
    

    return (f(a) + I1 + I2 + f(b))*(h/3)

    
# def CompTrap(a,b,f,n):
#     h = (b-a)/n
#     xnode = a+np.arange(0,n+1)*h
#     I_trap = h*f(xnode[0])*1/2
#     for j in range(1,n):
#         I_trap = I_trap+h*f(xnode[j])
#         I_trap= I_trap + 1/2*h*f(xnode[n])
#     return I_trap

# def CompSimp(a,b,f,n):
#     h = (b-a)/n
#     xnode = a+np.arange(0,n+1)*h
#     I_simp = f(xnode[0])
#     nhalf = n/2
#     for j in range(1,int(nhalf)+1):
#         # even part
#         I_simp = I_simp+2*f(xnode[2*j])
#         # odd part
#         I_simp = I_simp +4*f(xnode[2*j-1])
#         I_simp= I_simp + f(xnode[n])
#         I_simp = h/3*I_simp
#     return I_simp



driver()