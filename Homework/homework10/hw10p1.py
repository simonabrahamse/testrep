import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt

def driver():
    A = np.array([
        [1, 0, 0, 0, 0, 0, 0],        
        [0, 1, 0, 0, -1, 0, 0],        
        [0, 0, 1, 0, 0, -1, 0],         
        [0, 0, 0, 1, 0, 0, -1/6],    
        [0, 0, 0, 0, 1, 0, 0],       
        [0, 0, 0, 0, 0, 1, 0],       
        [0, 0, 0, 0, 0, 0, 1/120]   
    ])

    A_b = np.array([
        [1, 0, 0, 0, 0, 0, 0],      
        [0, 1, 0, -1, 0, 0, 0],            
        [0, 0, 1, 0, -1, 0, 0],             
        [0, 0, 0, 1, 0, -1/6, 0],   
        [0, 0, 0, 0, 1, 0, -1/24],  
        [0, 0, 0, 0, 0, 1, 0],      
        [0, 0, 0, 0, 0, 0, 1/120]   
    ])


    A_c = np.array([
        [1, 0, 0, 0, 0, 0, 0],           
        [0, 1, 0, -1, 0, 0, 0],           
        [0, 0, 1, 0, -1, 0, 0],            
        [0, 0, 0, 1, 0, 0, -1/6],      
        [0, 0, 0, 0, 1, 0, 0],        
        [0, 0, 0, 0, 0, 1, 0],        
        [0, 0, 0, 0, 0, 0, 1/120]     
    ])

    b = np.array([0, 1, 0, -1/6, 0, 1/120, 0])

    
    coeffsa = solve(A, b)
    print(coeffsa)
    coeffsb = solve(A_b,b)
    coeffsc = solve(A_c,b)

    x = np.linspace(0,5,100)

    y = np.sin(x)
    P = pade_approx(x,coeffsa)
    P2 = pade_approxb(x,coeffsb)
    P3 = pade_approxc(x,coeffsc)
    T = x-x**6/6+x**5/120

    
    plt.plot(x, np.abs(y-P),'o', label='Pade approximation error')
    plt.plot(x, np.abs(y-P2),'o', label='Pade approximation error part b')
    plt.plot(x, np.abs(y-P3),'o', label='Pade approximation error part c')
    plt.plot(x, np.abs(y-T), 'o', label='Maclaurin Polynomial error')

    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.yscale('log')
    plt.title('Pade approximation accuracy vs. Maclaurin')
    plt.legend()
    plt.grid(True)
    plt.show()




def pade_approx(x, coeffs):
    a0, a1, a2, a3, b1, b2, b3 = coeffs
    numerator = a0 + a1 * x + a2 * x**2 + a3 * x**3
    denominator = 1 + b1 * x + b2 * x**2 + b3 * x**3
    return numerator / denominator

def pade_approxb(x,coeffsb):
    a0, a1, a2, b1, b2, b3, b4 = coeffsb
    numerator = a0 + a1 * x + a2 * x**2
    denominator = 1 + b1 * x + b2 * x**2 + b3*x**3 + b4*x**4
    return numerator / denominator

def pade_approxc(x,coeffsc):
    a0, a1, a2, a3, a4, b1, b2 = coeffsc
    numerator = a0 + a1 * x + a2 * x**2 + a3 * x**3 + a4*x**4
    denominator = 1 + b1 * x + b2 * x**2 
    return numerator / denominator


driver()