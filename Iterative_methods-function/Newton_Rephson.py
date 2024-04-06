import numpy as np
import sympy as sp
import math

from sympy import sympify
from sympy.utilities.lambdify import lambdify

x = sp.symbols('x')


def max_steps(a, b, err):
    s = int(np.floor(- np.log2(err / (b - a)) / np.log2(2) - 1))
    return s

def newton_rephson(P, test_range, epsilon):
    """
    This function calculates polynom root (P) by using Newton's formula: xr+1= xr-(f(x) / f'(x))
    :param P: Polynom
    :type P: function
    :param test_range: A range of values in which we will check if a root exists
    :type test_range: list with 2 values
    :param epsilon: Another stopping condition - when the difference of x between the 2 iterations is less than or equal
     to epsilon
    :type epsilon:float
    :return: None
    """

    p_derivative = sp.diff(P, x)  # p_derivative = p'(x)
    p_derivative = lambdify(x, p_derivative)  # convert to an equivalent NumPy function to perform calculations
    P = lambdify(x, P) # convert to an equivalent NumPy function to perform calculations
    max_step = 0
    if len(test_range) == 0:
        x0 = 0.5
        max_step = 50

    elif len(test_range) == 1:
        x0 = test_range[0]


    else:
        x0 = (test_range[-1] + test_range[0]) / 2  # x0 is the first value that are guessed and placed in a Newton's formula
        max_step = max_steps(test_range[0], test_range[-1], epsilon)
    count_iteration = 0

    if P(x0) == 0:
        print('The root is ', x0, '. The number of iterations is ', count_iteration)
        return

    if p_derivative(x0) == 0:
        test_range[0] += epsilon
        x0 = (test_range[-1] + test_range[0]) / 2
    x_root = x0 - ((P(x0)) / (p_derivative(x0)))  # x value of the next iteration, Newton's formula


    """if abs(p_derivative(x0)) < 1:
        print('The convergence condition is not met, therefore in the field:', test_range, 'there is no roots')
        return"""

    for i in range(max_step):
        old_x0 = x0
        x0 = x_root
        x_root = x0 - ((P(x0)) / (p_derivative(x0)))  # x value of the next iteration, Newton's formula
        if (round(P(x_root), 10000) == 0) or (abs(x0 - old_x0) <= epsilon):
            print('The root is ', round(x_root, 6), '\n f( ', x_root, ')= ', round(P(x_root), 10),
                  '\nThe number of iterations is ', i+1)
            return


    print("For:", P, " in the field", test_range, "there is no root")



if __name__ == '__main__':
    epsilon = 0.0001
    #try number 1
    f = x**3 - x - 1
    domain = [-2, 2]
    """successes
    
    try number 2
    f = -2 * sp.log(x) + 2
    domain = []
    newton_rephson(f, domain, epsilon)
    successes"""
    """f = input("Enter polynom: ")
    f = sympify(f)
    f = -x**3 + sp.sin(x**2) -1
    domain = [-2, 2]"""
    newton_rephson(f, domain, 0.0001)


