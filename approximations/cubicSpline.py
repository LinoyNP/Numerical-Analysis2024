import numpy as np

def cubicSpline(table_points, point, type = False):
    table_size = len(table_points)
    matrix = [[0 for x in range(table_size)] for y in range(table_size)]
    if type: # this is natural cubicSpline
        #matrix[0][0] = matrix[]
        for i in range(1, table_size):
            h0 = table_points[i][0] - table_points[i-1][0]
            h1 = table_points[i+1][1] - table_points[i][1]
            for j in range(i-1, i+2):
                if i != table_size-1:  # this is not last line
                    matrix[i][j].append(1/6)

    else:
        h = table_points[1][0] - table_points[0][0]
        matrix[0][0] = 1/3 * h
        matrix[0][1] = 1/6 * h
        for i in range(1, table_size):

            h = table_points[i+1][0] - table_points[i][0]

from Iterative_methods_matrix.Jacobi import jacobi_iterative
from sympy import *

x = Symbol('x')


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

    # find S:
    for loc in range(1, len(f)):
        s = (((f[loc][0] - x) ** 3) * M[loc - 1] + ((x - f[loc - 1][0]) ** 3) * M[loc]) / (6 * h[loc - 1])
        s += (((f[loc][0] - x) * f[loc - 1][1]) + ((x - f[loc - 1][0]) * f[loc][1])) / h[loc - 1]
        s -= (((f[loc][0] - x) * M[loc - 1] + (x - f[loc - 1][0]) * M[loc]) * h[loc - 1]) / 6
        print("s" + str(loc - 1) + "(x) = " + str(s))

    # find the location of x0:
    loc = 0
    for i in range(1, len(f)):
        if x0 < f[i][0] and x0 > f[i - 1][0]:
            loc = i
            break

    if loc == 0:
        print("no range found for x0")
        return

    s = (((f[loc][0] - x) ** 3) * M[loc - 1] + ((x - f[loc - 1][0]) ** 3) * M[loc]) / (6 * h[loc - 1])
    s += (((f[loc][0] - x) * f[loc - 1][1]) + ((x - f[loc - 1][0]) * f[loc][1])) / h[loc - 1]
    s -= (((f[loc][0] - x) * M[loc - 1] + (x - f[loc - 1][0]) * M[loc]) * h[loc - 1]) / 6

    print("\nx0 between f(x" + str(loc - 1) + ") = " + str(f[loc - 1][0]) + " and f(x" + str(loc) + ") = " + str(
        f[loc][0]) + " so:")
    print("s" + str(loc - 1) + "(" + str(x0) + ") = " + str(float(s.subs(x, x0))))


if __name__ == '__main__':
    f = [(1, 1), (2, 2), (3, 1), (4, 1.5), (5, 1)]
    x0 = 4.5

    print("func: " + str(f))
    print("x0 = " + str(x0) + "\n")
    natural_cubic_spline(f, x0)