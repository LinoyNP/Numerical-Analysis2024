import numpy as np
import math as m
from Iterative_methods_matrix.Jacobi import jacobi_iterative
from sympy import *

x = Symbol('x')

def cubicSpline(table_points, point, y1_derivative, y2_derivative):
    table_size = len(table_points)
    #matrix = [[0 for x in range(table_size)] for y in range(table_size)]  # Initialize the matrix to zeros
    matrix = list()
    matrix.append(list())
    h = list()
    for i in range(table_size - 1):
        h.append(table_points[i + 1][0] - table_points[i][0])  # h = xi - xi-1
    #matrix[0][0] = matrix[]
    # first line
    for j in range(table_size):
        if j == 0:
            matrix[0].append(h[0]/3)
        elif j == 1:
            matrix[0].append(h[0]/6)
        else:
            matrix[0].append(0)

    for i in range(1, table_size-1):
        matrix.append(list())
        for j in range(table_size):
            if j == i-1: matrix[i].append(h[i-1]/6)
            elif j == i: matrix[i].append((h[i-1] + h[i])/3)
            elif j == i+1: matrix[i].append(h[i]/6)
            else:
                matrix[i].append(0)

    #last line:
    matrix.append(list())
    for row in range(table_size):
        if row == table_size-2:
            matrix[table_size-1].append(h[-1] / 6)
        elif row == table_size-1:
            matrix[table_size-1].append(h[-1] / 3)
        else:
            matrix[table_size-1].append(0)

    d = list()
    for i in range(table_size):
        if i == 0:  # d0=(y1-y0)/h0 - y'(0)
            d.append(((table_points[1][1]-table_points[0][1])/h[0]) - y1_derivative)
        elif i == table_size-1:  # dn= y'(n)- (yn-yn-1)/h0
            d.append(y2_derivative - ((table_points[table_size-1][1] - table_points[table_size-2][1])/h[0]))
        else:
            d.append(((table_points[i + 1][1] - table_points[i][1]) / h[i]) - ((table_points[i][1] - table_points[i - 1][1]) / h[i - 1]))

    print("matrix: " + str(matrix))
    print("vector b: " + str(d))

    # get m vector
    result = np.zeros_like(d, dtype=np.double)
    print("\nJacobi middle results: ")
    M = (jacobi_iterative(matrix, d, result))
    print("\nvector M: " + str(list(map(float, M))))

    # find si(x)
    S = 0
    for i in range(table_size-1):
        S = table_points[i+1][1]*(x - table_points[i][0])/h[i] - table_points[i][1] * (x - table_points[i+1][0])/h[i]
        S += M[i+1]/6 * (((x - table_points[i][0])**3)/h[i] - (h[i] * (x -table_points[i][0])))
        S -= M[i]/6 * (((x - table_points[i+1][0])**3)/h[i] - (h[i] * (x -table_points[i+1][0])))

    # find the location of x0:
    loc = 0
    for i in range(1, len(f)):
        if x0 < f[i][0] and x0 > f[i - 1][0]:
            loc = i
            break

    if loc == 0:
        print("no range found for x0")
        return

    S = 0
    S = f[loc][1] * (x - f[loc-1][0]) / h[loc-1] - f[loc-1][1] * (x - f[loc][0]) / h[loc-1]
    S += M[loc] / 6 * (((x - f[loc-1][0]) ** 3) / h[loc-1] - (h[loc-1] * (x - f[loc-1][0])))
    S -= M[loc] / 6 * (((x - f[loc][0]) ** 3) / h[loc-1] - (h[loc-1] * (x - f[loc][0])))

    print("\nx0 between f(x" + str(loc - 1) + ") = " + str(f[loc - 1][0]) + " and f(x" + str(loc) + ") = " + str(
        f[loc][0]) + " so:")
    print("s" + str(loc - 1) + "(" + str(x0) + ") = " + str(float(S.subs(x, x0))))




def natural_cubic_spline(f, x0):
    h = list()
    for i in range(len(f) - 1):
        h.append(f[i + 1][0] - f[i][0])

    g = list()
    g.append(0)  # g0
    for i in range(1, len(f) - 1):
        g.append(h[i] / (h[i] + h[i - 1]))
    g.append(0)  # gn

    m = list()
    m.append(0)
    for i in range(1, len(f)):
        m.append(1 - g[i])

    d = list()
    d.append(0)  # d0=0
    for i in range(1, len(f) - 1):
        d.append((6 / (h[i - 1] + h[i])) * (((f[i + 1][1] - f[i][1]) / h[i]) - ((f[i][1] - f[i - 1][1]) / h[i - 1])))
    d.append(0)  # dn

    # building the matrix
    mat = list()

    # first row
    mat.append(list())
    mat[0].append(2)
    for j in range(len(f) - 1):
        mat[0].append(0)

    for i in range(1, len(f) - 1):
        mat.append(list())
        for j in range(len(f)):
            if j == i - 1:  # put miu
                mat[i].append(m[i])
            elif j == i:
                mat[i].append(2)
            elif j == i + 1:  # put lambda
                mat[i].append(g[i])
            else:
                mat[i].append(0)

    # last row
    mat.append(list())
    for j in range(len(f) - 1):
        mat[len(f) - 1].append(0)
    mat[len(f) - 1].append(2)

    print("matrix: " + str(mat))
    print("vector b: " + str(d))

    # get m vector
    result = np.zeros_like(d, dtype=np.double)
    print("\nJacobi middle results: ")
    M = (jacobi_iterative(mat, d, result))
    print("\nvector M: " + str(list(map(float, M))))

    """# find S:
    for loc in range(1, len(f)):
        s = (((f[loc][0] - x) ** 3) * M[loc - 1] + ((x - f[loc - 1][0]) ** 3) * M[loc]) / (6 * h[loc - 1])
        s += (((f[loc][0] - x) * f[loc - 1][1]) + ((x - f[loc - 1][0]) * f[loc][1])) / h[loc - 1]
        s -= (((f[loc][0] - x) * M[loc - 1] + (x - f[loc - 1][0]) * M[loc]) * h[loc - 1]) / 6
        print("s" + str(loc - 1) + "(x) = " + str(s))"""
    S = 0
    for i in range(len(f) - 1):
        S = f[i + 1][1] * (x - f[i][0]) / h[i] - f[i][1] * (x - f[i + 1][0]) / h[i]
        S += M[i + 1] / 6 * (((x - f[i][0]) ** 3) / h[i] - (h[i] * (x - f[i][0])))
        S -= M[i] / 6 * (((x - f[i + 1][0]) ** 3) / h[i] - (h[i] * (x - f[i + 1][0])))
        print("s" + str(i) + "(x) = " + str(S))

    # find the location of x0:
    loc = 0
    for i in range(1, len(f)):
        if x0 < f[i][0] and x0 > f[i - 1][0]:
            loc = i
            break

    if loc == 0:
        print("no range found for x0")
        return
    """s = (((f[loc][0] - x) ** 3) * M[loc - 1] + ((x - f[loc - 1][0]) ** 3) * M[loc]) / (6 * h[loc - 1])
    s += (((f[loc][0] - x) * f[loc - 1][1]) + ((x - f[loc - 1][0]) * f[loc][1])) / h[loc - 1]
    s -= (((f[loc][0] - x) * M[loc - 1] + (x - f[loc - 1][0]) * M[loc]) * h[loc - 1]) / 6"""

    S = 0
    S = f[loc][1] * (x - f[loc-1][0]) / h[loc-1] - f[loc-1][1] * (x - f[loc][0]) / h[loc-1]
    S += M[loc] / 6 * (((x - f[loc-1][0]) ** 3) / h[loc-1] - (h[loc-1] * (x - f[loc-1][0])))
    S -= M[loc] / 6 * (((x - f[loc][0]) ** 3) / h[loc-1] - (h[loc-1] * (x - f[loc][0])))

    print("\nx0 between f(x" + str(loc - 1) + ") = " + str(f[loc - 1][0]) + " and f(x" + str(loc) + ") = " + str(
        f[loc][0]) + " so:")
    print("s" + str(loc - 1) + "(" + str(x0) + ") = " + str(float(S.subs(x, x0))))




if __name__ == '__main__':
    f = [(1, 1), (2, 2), (3, 1), (4, 1.5), (5, 1)]
    x0 = 4.5
    """f = [(0,0), (m.pi/6, 0.5), (m.pi/4, 0.7072), (m.pi/2, 1)]
    x0 = m.pi/3"""
    y_0_deri = 1
    y_n_deri = 0

    print("func: " + str(f))
    print("x0 = " + str(x0) + "\n")
    #natural_cubic_spline(f, x0)
    cubicSpline(f, x0, y_0_deri, y_n_deri)