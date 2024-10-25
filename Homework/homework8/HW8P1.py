import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange, CubicSpline, PchipInterpolator

# Define the function
def f(x):
    return 1 / (1 + x**2)

# Define equispaced nodes
def equispaced_nodes(a, b, n):
    return np.linspace(a, b, n)

# Range for interpolation
x_range = np.linspace(-5, 5, 1000)

# Part 1(a): Lagrange Interpolation
def lagrange_interpolation(n):
    nodes = equispaced_nodes(-5, 5, n)
    values = f(nodes)
    poly = lagrange(nodes, values)
    return poly(x_range)

# Part 1(b): Hermite Interpolation (using PchipInterpolator for Hermite-like behavior)
def hermite_interpolation(n):
    nodes = equispaced_nodes(-5, 5, n)
    values = f(nodes)
    hermite_poly = PchipInterpolator(nodes, values)
    return hermite_poly(x_range)

# Part 1(c): Natural Cubic Spline
def natural_cubic_spline(n):
    nodes = equispaced_nodes(-5, 5, n)
    values = f(nodes)
    spline = CubicSpline(nodes, values, bc_type='natural')
    return spline(x_range)

# Part 1(d): Clamped Cubic Spline
def clamped_cubic_spline(n):
    nodes = equispaced_nodes(-5, 5, n)
    values = f(nodes)
    # Derivative estimates at endpoints (approximate using centered difference)
    f_prime_start = (f(nodes[1]) - f(nodes[0])) / (nodes[1] - nodes[0])
    f_prime_end = (f(nodes[-1]) - f(nodes[-2])) / (nodes[-1] - nodes[-2])
    spline = CubicSpline(nodes, values, bc_type=((1, f_prime_start), (1, f_prime_end)))
    return spline(x_range)

# Plot results for n = 5, 10, 15, 20
n_values = [5, 10, 15, 20]
for n in n_values:
    plt.figure(figsize=(12, 8))
    plt.plot(x_range, f(x_range), label='Original Function', color='black', linewidth=2)

    # Perform each interpolation and plot
    plt.plot(x_range, lagrange_interpolation(n), label=f'Lagrange (n={n})', linestyle='--')
    plt.plot(x_range, hermite_interpolation(n), label=f'Hermite (n={n})', linestyle=':')
    plt.plot(x_range, natural_cubic_spline(n), label=f'Natural Cubic Spline (n={n})', linestyle='-.')
    plt.plot(x_range, clamped_cubic_spline(n), label=f'Clamped Cubic Spline (n={n})', linestyle='-')

    plt.title(f'Interpolation Comparison for n = {n}')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.show()
