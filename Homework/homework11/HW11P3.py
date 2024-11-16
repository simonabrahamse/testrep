import numpy as np
from numpy.linalg import solve

A = np.array([
    [1, 1, 1],
    [1, 2**(3/2),4**(3/2)],
    [1, 2**2,4**2]
])
b = np.array([1, 0, 0])

coefficients = solve(A,b)
print(coefficients)