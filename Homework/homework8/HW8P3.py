import numpy as np
import matplotlib.pyplot as plt


def driver():
    N=10
    f = lambda x: np.sin(10*x)
    nodes = np.linspace(0,2*np.pi,20)
    y_nodes = f(nodes)
    M = tridiagonal(N,f,nodes)

    h = nodes[1]-nodes[0]

    a = f(nodes)
    b = np.zeros(N)
    for i in range(N-1):
        b[i] = (f(nodes[i+1])-f(nodes[i]))/h-h/6*(M[i+1]+2*M[i])
    c = np.zeros(N)
    for i in range(N):
        c[i] = M[i]/2
    d = np.zeros(N)
    for i in range(N-1):
        d[i] = (M[i+1]-M[i])/(6*h)

    x = np.linspace(0,2*np.pi,1000)
    xm = np.delete(x,-1)
    S = np.zeros(999)
    for i in range(999):
        S[i] = evaluate_spline(x[i],nodes,a,b,c,d)

    plt.plot(x, f(x), label='Original Function', color='black', linewidth=2)
    plt.plot(xm, S, label='Periodic Cubic Spline', linestyle='--')

    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Periodic Cubic Spline Interpolation with N = 10')
    plt.legend()
    plt.grid(True)
    plt.show()

def evaluate_spline(x, nodes, a, b, c, d):
    # Find the interval index i where x is located
    i = int((10*x)/(2*np.pi))
    print(i)
    if i < 0:
        i = 0
    elif i >= len(a):
        i = len(a) - 1

    # Calculate the difference (x - x_i)
    dx = x - nodes[i]

    # Evaluate the polynomial for the interval [x_i, x_{i+1}]
    S_x = a[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3
    return S_x



def tridiagonal(N,f,nodes):
    M = np.zeros((N,N))
    diag1 = np.zeros(N) +4
    diag1[0] = 2
    diag1[-1]=2
    diag2 = np.ones(N-1)

    M = M+np.diag(diag1)
    M = M+np.diag(diag2,-1)
    M = M+np.diag(diag2,1)
    
    M[0,N-1] = 1
    M[N-1,0] = 1
    
    h = nodes[1]-nodes[0]

    R = np.zeros(N)
    for i in range(1,N):
        R[i] = (f(nodes[i+1])-f(nodes[i]))/h-(f(nodes[i])-f(nodes[i-1]))/h

    Ms = np.linalg.solve(M,R)



    return Ms


driver()