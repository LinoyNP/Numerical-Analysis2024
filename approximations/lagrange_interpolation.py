from colors import bcolors
import sympy as sp
from sympy.utilities.lambdify import lambdify

x_variable = sp.symbols('x_variable')

def lagrange_interpolation(x_data, y_data, x):
    """
    Lagrange Interpolation
    This function create polynom and interpolate according to the Lagrange
    Parameters:
    x_data (list): List of x-values for data points.
    y_data (list): List of y-values for data points.
    x (float): The x-value where you want to evaluate the interpolated polynomial.

    Returns:
    float: The interpolated y-value at the given x.
    """
    n = len(x_data)
    result = 0.0
    L = []  # L will contain lists that will contain the part of polynom- pn(x)=L[0]*y0
    #  + L[1]*y1 +...+ L[n]*yn

    for i in range(n): # build polynom
        term = temp = y_data[i]
        for j in range(n):
            if i != j:
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
                temp *= (x_variable - x_data[j]) / (x_data[i] - x_data[j])
        L.append(temp)
        result += term

    P = 0  # Complete polynom
    for i in range(len(L)):
        P += L[i]
    print("The polynomial is: ", P)
    P = lambdify(x_variable, P)
    print(P(x))
    return result

if __name__ == '__main__':

    x_data = [1, 2, 5]
    y_data = [1, 0, 2]
    x_interpolate = 3  # The x-value where you want to interpolate
    y_interpolate = lagrange_interpolation(x_data, y_data, x_interpolate)
    print(bcolors.OKBLUE, "\nInterpolated value at x =", x_interpolate, "is y =", y_interpolate, bcolors.ENDC)


