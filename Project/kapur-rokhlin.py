import numpy as np
from scipy.integrate import quad

def driver():
    # Regularized test functions to handle scalars and arrays
    gx1 = lambda x: np.where(x > 0, x * np.log(x), 0)
    gx2 = lambda x: np.where(x > 0, -x * np.log(x), 0)
    gx3 = lambda x: np.where(np.sin(x) > 0, np.cos(x) * np.log(np.sin(x)), 0)
    gx3transform = lambda x: np.where(x > 0, np.pi*np.cos(np.pi*x) * np.log(x), 0)

    # Kapur-Rokhlin weights for different orders
    nodes2 = np.array([1.825748064736159, -1.325748064736159])
    nodes6 = np.array([4.967362978287758, -16.20501504859126, 25.85153761832639,
                       -22.22599466791883, 9.930104998037539, -1.817995878141594])
    nodes10 = np.array([7.832432020568779, -45.65161670374749, 145.2168846354677,
                        -290.1348302886379, 387.08621625799, -352.3821383570681,
                        217.2421547519342, -87.07796087382991, 20.53584266072635, -2.166984103403823])

    a = 0
    b = 1
    N2 = 2
    N6 = 6
    N10 = 10

    # Reference integral using scipy's quad
    Iactual, err = quad(lambda x: x * np.log(x) if x > 0 else 0, a, b)

    # Kapur-Rokhlin quadrature
    I2 = krlog(gx1,a,b,N2,nodes2)
    I6 = krlog(gx1,a,b,N6,nodes6)
    I10 = krlog(gx1, a, b, N10, nodes10)

    print(f"Computed integral using Kapur-Rokhlin quadrature with 2 nodes: {I2}")
    print(f"Computed integral using Kapur-Rokhlin quadrature with 6 nodes: {I6}")
    print(f"Computed integral using Kapur-Rokhlin quadrature with 10 nodes: {I10}")
    print(f"Reference integral using quad: {Iactual}")


    # Reference integral using scipy's quad
    Iactual, err = quad(lambda x: -x * np.log(x) if x > 0 else 0, a, b)

    # Kapur-Rokhlin quadrature
    I2 = krlog(gx2,a,b,N2,nodes2)
    I6 = krlog(gx2,a,b,N6,nodes6)
    I10 = krlog(gx2, a, b, N10, nodes10)

    print(f"Computed integral using Kapur-Rokhlin quadrature with 2 nodes: {I2}")
    print(f"Computed integral using Kapur-Rokhlin quadrature with 6 nodes: {I6}")
    print(f"Computed integral using Kapur-Rokhlin quadrature with 10 nodes: {I10}")
    print(f"Reference integral using quad: {Iactual}")

    # Reference integral using scipy's quad
    a = 0
    b = np.pi

    Iactual, err = quad(lambda x: np.cos(x) * np.log(np.sin(x)) if x > 0 else 0, a, b)

    a = 0
    b= 1

    # Kapur-Rokhlin quadrature
    I2 = krlog(gx3transform,a,b,N2,nodes2)
    I6 = krlog(gx3transform,a,b,N6,nodes6)
    I10 = krlog(gx3transform, a, b, N10, nodes10)

    print(f"Computed integral using Kapur-Rokhlin quadrature with 2 nodes: {I2}")
    print(f"Computed integral using Kapur-Rokhlin quadrature with 6 nodes: {I6}")
    print(f"Computed integral using Kapur-Rokhlin quadrature with 10 nodes: {I10}")
    print(f"Reference integral using quad: {Iactual}")







def krlog(gx, a, b, N, nodes):
    """
    Kapur-Rokhlin quadrature for integrals of the form
    ∫_a^b f(x) + g(x) log(x) dx.
    """
    h = (b - a) / N
    x_points = np.linspace(a, b, N + 1)

    # Trapezoidal rule contribution
    I1 = 0.5 * (gx(x_points[0]) + gx(x_points[-1])) + np.sum(gx(x_points[1:-1]))
    I1 *= h

    # Kapur-Rokhlin corrections
    I2 = 0  # Lower-end correction
    I3 = 0  # Upper-end correction
    num_correction_nodes = len(nodes)

    for l in range(1, num_correction_nodes + 1):
        if l <= N:  # Only apply corrections if they fall within the range
            I2 += nodes[l - 1] * gx(a + l * h)
            I3 += nodes[l - 1] * gx(b - l * h)

    I2 *= h
    I3 *= h

    return I1 + I2 + I3


driver()
