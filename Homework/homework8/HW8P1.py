import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange, CubicSpline, PchipInterpolator

def driver():
    n_values = [5, 10, 15, 20]
    for n in n_values:
        plt.figure(figsize=(12, 8))
        plt.plot(x_range, f(x_range), label='Original Function', color='black', linewidth=2)

        lg = lagrange_interpolation(n)
        hermite = hermite_interpolation(n)
        cubic = natural_cubic_spline(n)
        clamped = clamped_cubic_spline(n)

        # Perform each interpolation and plot
        plt.plot(x_range, lg, label=f'Lagrange (n={n})', linestyle='--')
        plt.plot(x_range, hermite, label=f'Hermite (n={n})', linestyle=':')
        plt.plot(x_range, cubic, label=f'Natural Cubic Spline (n={n})', linestyle='-.')
        plt.plot(x_range, clamped, label=f'Clamped Cubic Spline (n={n})', linestyle='-')

        plt.title(f'Interpolation Comparison for n = {n}')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.legend()
        plt.grid(True)
        plt.show()

        plt.figure(figsize=(12, 8))

        plt.plot(x_range, np.abs(lg-f(x_range)),'o', label=f'Lagrange (n={n})')
        plt.plot(x_range, np.abs(hermite-f(x_range)),'o',label=f'Hermite (n={n})')
        plt.plot(x_range, np.abs(cubic-f(x_range)),'o', label=f'Natural Cubic Spline (n={n})')
        plt.plot(x_range, np.abs(clamped-f(x_range)),'o', label=f'Clamped Cubic Spline (n={n})')

        plt.title(f'Interpolation error Comparison for n = {n}')
        plt.xlabel('x')
        plt.ylabel('error at each point')
        # plt.ylim(-0.5/n,0.5/n)
        plt.yscale('log')
        plt.legend()
        plt.grid(True)
        plt.show()





# Range for interpolation
x_range = np.linspace(-5, 5, 1000)
f = lambda x: 1/(1+x**2)

# Part 1(a): Lagrange Interpolation
def lagrange_interpolation(n):
    nodes = np.linspace(-5, 5, n)
    values = f(nodes)
    poly = lagrange(nodes, values)
    return poly(x_range)

# Part 1(b): Hermite Interpolation (using PchipInterpolator for Hermite-like behavior)
def hermite_interpolation(n):
    nodes = np.linspace(-5, 5, n)
    values = f(nodes)
    hermite_poly = PchipInterpolator(nodes, values)
    return hermite_poly(x_range)

# Part 1(c): Natural Cubic Spline
def natural_cubic_spline(n):
    nodes = np.linspace(-5, 5, n)
    values = f(nodes)
    spline = CubicSpline(nodes, values, bc_type='natural')
    return spline(x_range)

# Part 1(d): Clamped Cubic Spline
def clamped_cubic_spline(n):
    nodes = np.linspace(-5, 5, n)
    values = f(nodes)
    
    f_prime = lambda x: -1/((1+x**2)**2)

    f_prime_start = f_prime(nodes[0])
    f_prime_end = f_prime(nodes[-1])
    spline = CubicSpline(nodes, values, bc_type=((1, f_prime_start), (1, f_prime_end)))
    return spline(x_range)

driver()