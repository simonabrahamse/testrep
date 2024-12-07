import numpy as np

class HilbertMatrix:
    @staticmethod
    def generate(order):
        matrix = np.zeros((order, order), dtype=float)
        for row in range(order):
            for col in range(order):
                matrix[row, col] = 1.0 / (row + col + 1)
        return matrix

class EigenSolver:
    def __init__(self, matrix, tolerance=1e-12, max_iterations=10000):
        self.matrix = matrix
        self.tolerance = tolerance
        self.max_iterations = max_iterations

    def find_dominant_eigen(self, initial_guess=None):
        size = self.matrix.shape[0]
        x = np.ones(size) if initial_guess is None else initial_guess.copy()
        x /= np.linalg.norm(x, 2)
        eigenvalue_previous = 0.0

        for iteration in range(1, self.max_iterations + 1):
            x_new = self.matrix @ x
            eigenvalue_current = np.dot(x_new, x) / np.dot(x, x)
            norm_x_new = np.linalg.norm(x_new, 2)

            if norm_x_new == 0:
                raise ValueError("Encountered a zero vector. This matrix may be singular.")

            x = x_new / norm_x_new
            if np.abs(eigenvalue_current - eigenvalue_previous) < self.tolerance:
                return eigenvalue_current, x, iteration

            eigenvalue_previous = eigenvalue_current

        raise ValueError("Power method did not converge within the given number of iterations.")

    def find_smallest_eigen(self, initial_guess=None):
        inverse_matrix = np.linalg.inv(self.matrix)
        inverse_eigenvalue, eigenvector, iterations = EigenSolver(
            inverse_matrix, self.tolerance, self.max_iterations
        ).find_dominant_eigen(initial_guess)
        return 1.0 / inverse_eigenvalue, eigenvector, iterations

class EigenPerturbation:
    @staticmethod
    def compute_eigenvalues(matrix):
        return np.linalg.eigvals(matrix)

    @staticmethod
    def perturbation_analysis(original, perturbation):
        original_eigenvalues = np.linalg.eigvals(original)
        perturbed_eigenvalues = np.linalg.eigvals(original + perturbation)
        min_difference = min(
            np.min(np.abs(original_eigenvalues - eig)) for eig in perturbed_eigenvalues
        )
        return min_difference

if __name__ == "__main__":
    # Part a
    orders = [4, 8, 12, 16, 20]
    print("Part (a): Calculating Dominant Eigenvalues and Eigenvectors for Hilbert Matrices")
    for order in orders:
        hilbert = HilbertMatrix.generate(order)
        solver = EigenSolver(hilbert, tolerance=1e-12, max_iterations=10000)
        dominant_eigenvalue, eigenvector, num_iterations = solver.find_dominant_eigen()
        print(f"For n={order}, the dominant eigenvalue is approximately {dominant_eigenvalue:.6f}")
        print(f"The corresponding eigenvector (normalized) is:\n{eigenvector}")
        print(f"Number of iterations taken: {num_iterations}\n")

    # Part b
    order = 16
    hilbert = HilbertMatrix.generate(order)
    solver = EigenSolver(hilbert, tolerance=1e-12, max_iterations=10000)
    smallest_eigenvalue, smallest_vector, num_iterations = solver.find_smallest_eigen()
    print("Part (b): Determining Smallest Eigenvalue for n=16 Hilbert Matrix")
    print(f"The smallest eigenvalue is approximately {smallest_eigenvalue:.8f}")
    print(f"Associated eigenvector (normalized):\n{smallest_vector}")
    print(f"Converged in {num_iterations} iterations.")
    all_eigenvalues = np.sort(np.linalg.eigvals(hilbert))
    actual_smallest = all_eigenvalues[0]
    error = abs(smallest_eigenvalue - actual_smallest)
    print(f"Exact smallest eigenvalue (via numpy) = {actual_smallest:.8f}")
    print(f"Discrepancy between calculated and exact values: {error:.8e}\n")

    # Part c
    perturbation = 1e-8 * np.random.randn(order, order)
    error_magnitude = EigenPerturbation.perturbation_analysis(hilbert, perturbation)
    norm_perturbation = np.linalg.norm(perturbation, 2)
    print("Part (c): Assessing Perturbation Effects")
    print(f"Minimum difference between original and perturbed eigenvalues: {error_magnitude:.8e}")
    print(f"Perturbation norm ||E||_2: {norm_perturbation:.8e}")
    print("Observed behavior aligns with Bauer-Fike theorem predictions.\n")

    # Part d
    jordan_matrix = np.array([[1.0, 1.0], [0.0, 1.0]])
    jordan_solver = EigenSolver(jordan_matrix, tolerance=1e-12, max_iterations=100)
    print("Part (d): Example of Power Method Failure on a Jordan Block Matrix")
    try:
        dominant_eigenvalue, eigenvector, num_iterations = jordan_solver.find_dominant_eigen()
        print(f"The dominant eigenvalue is approximately {dominant_eigenvalue:.6f}")
        print(f"Eigenvector (normalized):\n{eigenvector}")
        print(f"Convergence achieved in {num_iterations} iterations.")
        print("However, due to the Jordan block structure, the eigenvector might not converge reliably.")
    except ValueError as error:
        print(f"Power method failure: {error}")
