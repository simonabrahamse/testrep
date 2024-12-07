import numpy as np

# Define the matrix
A = np.array([
    [12, 10, 4],
    [10, 8, -5],
    [4, -5, 3]
], dtype=float)

def householder_tridiagonalize(A):
    n = A.shape[0]
    for k in range(n - 2):
        # Extract the vector to zero out below diagonal
        x = A[k+1:, k]
        norm_x = np.linalg.norm(x)
        v = x.copy()
        v[0] += np.sign(x[0]) * norm_x
        v /= np.linalg.norm(v)

        # Form the Householder matrix
        H_k = np.eye(n)
        H_k[k+1:, k+1:] -= 2.0 * np.outer(v, v)
        print(H_k)
        # Apply similarity transformation
        A = H_k @ A @ H_k
    return A

# Transform the matrix
tridiagonal_matrix = householder_tridiagonalize(A)
print("Tridiagonal form:\n", tridiagonal_matrix)
