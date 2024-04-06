import sympy as sp
import math
from sympy.utilities.lambdify import lambdify
from colors import bcolors
x = sp.symbols('x')
err = sp.symbols('err')


def trapezoidal_rule(f, a, b, n):
    """
    This function calculates the approximation of an integral according to the trapezoidal method
    :param f: The function on which the integral will be performed
    :param a: Lower bound
    :param b: Upper bound
    :param n: The number of sections to be divided
    :return: Integral of the function in the given field
    """
    second_diff = sp.diff(f, x, 2)  # Second derivative to the function that defined by 'x' for error estimate

    f = lambdify(x, f)
    h = (b - a) / n
    T = f(a) + f(b)
    integral = 0.5 * T  # Initialize with endpoints

    for i in range(1, n):
        x_i = a + i * h
        integral += f(x_i)

    integral *= h

    # Calculating the error estimate
    err = b # estimated error for upper bound
    second_diff = lambdify(x, second_diff)
    Error = 1/12 * h**2 * (b - a) * second_diff(err)

    print("The estimated error for err = ", err, " is: ", Error)

    return integral


if __name__ == '__main__':
    f = sp.exp(x ** 2)
    result = trapezoidal_rule(f, 0, 1, 2)
    print(bcolors.OKBLUE,"Approximate integral:", result, bcolors.ENDC)
