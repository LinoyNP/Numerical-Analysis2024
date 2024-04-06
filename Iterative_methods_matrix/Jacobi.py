import numpy as np
from numpy.linalg import norm

from colors import bcolors
from matrix_utility import is_diagonally_dominant, DominantDiagonalFix, is_square_matrix

"""
Performs Jacobi iterations to solve the line system of equations, Ax=b, 
starting from an initial guess, ``x0``.

Terminates when the change in x is less than ``tol``, or
if ``N`` [default=200] iterations have been exceeded.

Receives 5 parameters:
    1.  a, the NxN matrix that method is being performed on.
    2.  b, vector of solution. 
    3.  X0,  the desired initial guess.
        if x is None, the initial guess will bw determined as a vector of 0's.
    4.  TOL, tolerance- the desired limitation of tolerance of solution's anomaly.
        if tolerance is None, the default value will set as 1e-16.
    5.  N, the maxim number of possible iterations to receive the most exact solution.
        if N is None, the default value will set as 200.

Returns variables:
    1.  x, the estimated solution

"""

def jacobi_iterative(A, b, X0, TOL=1e-16, N=200):
    n = len(A)
    k = 1

    if is_diagonally_dominant(A):

        print('Matrix is diagonally dominant - preforming jacobi algorithm\n')
        print( "Iteration" + "\t\t\t".join([" {:>12}".format(var) for var in ["x{}".format(i) for i in range(1, len(A) + 1)]]))
        print("-----------------------------------------------------------------------------------------------")

        while k <= N:
            x = np.zeros(n, dtype=np.double)
            for i in range(n):
                sigma = 0
                for j in range(n):
                    if j != i:
                        sigma += A[i][j] * X0[j]
                x[i] = (b[i] - sigma) / A[i][i]

            print("{:<15} ".format(k) + "\t\t".join(["{:<15} ".format(val) for val in x]))

            if norm(x - X0, np.inf) < TOL:
                return tuple(x)

            k += 1
            X0 = x.copy()

        print("Maximum number of iterations exceeded. try solve it with another way")
        return tuple(x)
    else:
        L = list()
        D = list()
        U = list()
        for i in range(n):
            L.append(list())
            D.append(list())
            U.append(list())
            for j in range(A.shape[1]):  # Building L, D, U
                if j < i:
                    L[i].append(A[i][j])
                    D[i].append(0)
                    U[i].append(0)
                elif j == i:
                    L[i].append(0)
                    D[i].append(A[i][j])
                    U[i].append(0)
                else:
                    L[i].append(0)
                    D[i].append(0)
                    U[i].append(A[i][j])
        #U_inverse = np.linalg.inv(np.array(U))
        D_inverse = np.linalg.inv(np.array(D))
        G = np.multiply(-1, D_inverse)
        G = np.dot(G, np.add(L, U))
        if norm(G, np.inf) > 1:  # np.inf = max norm
            print("The condition for convergence is not met and therefore cannot be solved by Jacobi's method")
            return None

        while(k < N):
            x = np.add(np.dot(G, X0), np.dot(D_inverse, b))  # x= G*xr + D^-1*b
            if norm(x - X0, np.inf) < TOL:
                return tuple(x)
            X0 = x.copy()
            k+=1


if __name__ == "__main__":

    A = np.array([[3, -1, 1], [0, 1, -1], [1, 1, -2]])
    #A = np.arange(9).reshape(3,3) #[1,2,3][4,5,6,][7,8,9]
    print(A)
    b = np.array([4, -1, -3])

    x = np.zeros_like(b, dtype=np.double)
    solution = jacobi_iterative(A, b, x)
    if solution is not None:
        print(bcolors.OKBLUE,"\nApproximate solution:", solution)
