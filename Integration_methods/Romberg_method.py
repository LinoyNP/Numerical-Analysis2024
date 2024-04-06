import numpy as np
from colors import bcolors

def romberg_integration1(func, a, b, n, epsilon):
    """
    This function implements Romberg method by predefined level of accuracy (and not by a number of steps only)
    :param func: The function to be integrated.
    :param a: The lower limit of integration.
    :type a: float
    :param b: The upper limit of integration.
    :type b: float
    :param n:  The number of iterations (higher value leads to better accuracy).
    :type n: int
    :param epsilon: The level of accuracy of the error for which we will stop performing an integral approximation calculation
    :type epsilon: float
    Return:The approximate definite integral of the function over [a, b].
    :rtype: list of lists - 2D matrix
    """
    h = b - a
    R = []
    temp = []
    temp.append(0.5 * h * (func(a) + func(b)))  # T(1))
    R.append(temp)
    h/=2
    temp = []
    temp.append(0.5 * R[0][0] + h * func(a + h))  # T(0.5)
    temp.append(temp[0] + (1 / 3) * (temp[0] - R[0][0]))  # T(1, 0.5)
    R.append(temp)
    i = 2
    real_steps = 0
    while abs(R[i-1][-1] - R[i-2][-1]) > epsilon and i < n:
        temp = []
        h /= 2

        sum_term = 0
        for k in range(1, 2 ** i, 2):
            sum_term += func(a + k * h)

        temp.append(0.5 * R[i - 1][0] + h * sum_term)

        for j in range(1, i+1):
            temp.append(temp[j - 1] + (temp[j - 1] - R[i - 1][j - 1]) / ((4 ** j) - 1))
        R.append(temp)
        real_steps = i
        i += 1
    print("The number of steps to perform in order to approximate an integral with an error of:%f is: %d" % (epsilon, real_steps))
    return R[-1][-1]
def romberg_integration(func, a, b, n):
    """
    Romberg Integration

    Parameters:
    func (function): The function to be integrated.
    a (float): The lower limit of integration.
    b (float): The upper limit of integration.
    n (int): The number of iterations (higher value leads to better accuracy).

    Returns:
    float: The approximate definite integral of the function over [a, b].
    """
    h = b - a
    R = np.zeros((n, n), dtype=float)

    R[0, 0] = 0.5 * h * (func(a) + func(b))

    for i in range(1, n):
        h /= 2
        sum_term = 0

        for k in range(1, 2 ** i, 2):
            sum_term += func(a + k * h)

        R[i, 0] = 0.5 * R[i - 1, 0] + h * sum_term

        for j in range(1, i + 1):
            R[i, j] = R[i, j - 1] + (R[i, j - 1] - R[i - 1, j - 1]) / ((4 ** j) - 1)

    return R[n - 1, n - 1]


def f(x):
    return 1/(2+x ** 4)


if __name__ == '__main__':

    a = 0
    b = 1
    n = 5
    integral = romberg_integration(f, a, b, n)

    epsilon = 0.0001  #Accuracy of 5 digits after the point
    integral_epsilon = romberg_integration1(f, a, b, n, epsilon)
    """print( f" Division into n={n} sections ")
    print(bcolors.OKBLUE, f"Approximate integral in range [{a},{b}] is {integral}", bcolors.ENDC)
    print( f" Division into n={n} sections ")"""
    print(bcolors.OKBLUE, f"Approximate integral in range [{a},{b}] is {integral_epsilon}", bcolors.ENDC)
