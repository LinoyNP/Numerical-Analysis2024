"""
Date: 8.4.24
Group: Chaya Mizrachi ID: 214102584,
        Yael Siman-Tov ID:325181295,
        Linoy Nisim Pur ID: 324029685
Source Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
GitHub of this project: https://github.com/LinoyNP/Numerical-Analysis2024.git
Name:Linoy Nisim Pur ID:324029685
"""
from Integration_methods.Romberg_method import *
from Integration_methods.Trapezoidal_method import *
from approximations.polynomial_interpolation import *
from approximations.cubicSpline import *
import sympy as sp
from sympy.utilities.lambdify import lambdify
x = sp.symbols('x')
def main():
    print('''
        \tDate: 8.4.24
        \tGroup: Chaya Mizrachi ID: 214102584, Yael Siman-Tov ID:325181295, Linoy Nisim Pur ID: 324029685
        \tGit: https://github.com/LinoyNP/Numerical-Analysis2024.git
        \tName: Linoy Nisim Pur''')
    #Part1 - Section 8
    table_points = [(1.2, -3.5), (1.3, -3.69), (1.4, 0.9043), (1.5, 1.1293), (1.6, 2.3756)]
    xa = 1.35
    xb = 1.55

    print(bcolors.OKBLUE, "----------------- polynomial interpolation -----------------\n", bcolors.ENDC)
    print(bcolors.OKBLUE, "Finding an approximation with polynomial interpolation, to the point: ", bcolors.ENDC, xa, '\n')
    fa = polynomialInterpolation(table_points, xa)
    print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n",
          bcolors.ENDC)
    print(bcolors.OKBLUE, "Finding an approximation with polynomial interpolation, to the point: ", bcolors.ENDC, xb, '\n')
    fb = polynomialInterpolation(table_points, xb)
    print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n",
          bcolors.ENDC)


    # Part2- Section 6
    f = (2*(x**2) + sp.cos(2 * (sp.exp(-2 * x)))) / (2 * x**3 + x**2 - 6)
    f = lambdify(x, f)
    epsilon = 0.0001  # Accuracy of 5 digits after the point
    approximation = romberg_integration1(f, round(fa, 1), round(fb, 1), 20, epsilon)
    print("The integration approximation is: ", round(approximation, 5))

if __name__ == '__main__':
    main()