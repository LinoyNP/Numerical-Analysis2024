import numpy as np

from colors import bcolors
from inverse_matrix.matrix_utility import swap_rows_elementary_matrix, row_addition_elementary_matrix, matrix_multiply
""" instead of calculating the Ax=b with inverse matrix-> x=A^-1b(because not always it is easy to calculate A^-1) 
so we can calculate this equation with LU factorization- 
A=LU
U is 
L is the multy of elementary matrix inverse that need to to arrive at an upper triangular matrix"
"""

def lu(A):

        N = len(A)
        L = np.eye(N) # Create an identity matrix of size N x N

        for i in range(N):
    
            # Partial Pivoting: Find the pivot row with the largest absolute value in the current column
            pivot_row = i
            v_max = A[pivot_row][i]
            for j in range(i + 1, N):
                if abs(A[j][i]) > abs(v_max):
                    v_max = A[j][i]
                    pivot_row = j

            # if a principal diagonal element is zero,it denotes that matrix is singular,
           ## # and will lead to a division-by-zero later.
            if 0 <= A[pivot_row][i] <= 1e-10:
                raise ValueError("can't perform LU Decomposition")

            # Swap the current row with the pivot row
            if pivot_row != i:
                e_matrix = swap_rows_elementary_matrix(N, i, pivot_row)
                #print(f"elementary matrix for swap between row {i} to row {pivot_row} :\n {e_matrix} \n")
                A = np.dot(e_matrix, A)
                #A= matrix_multiply(e_matrix, A)
                """print(f"The matrix after elementary operation :\n {A}")
                print(bcolors.OKGREEN,"---------------------------------------------------------------------------", bcolors.ENDC)"""

            for j in range(i + 1, N):
    
                #  Compute the multiplier
                m = -A[j][i] / A[i][i]
                e_matrix = row_addition_elementary_matrix(N, j, i, m)
                e_inverse = np.linalg.inv(e_matrix)
                L = np.dot(L, e_inverse) #L=l^-1 * l^-2 *... l^-n
                A = np.dot(e_matrix, A)
                #A = matrix_multiply(e_matrix, A)
                #print(f"elementary matrix to zero the element in row {j} below the pivot in column {i} :\n {e_matrix} \n")
                #print(f"The matrix after elementary operation :\n {A}")
                #print(bcolors.OKGREEN,"---------------------------------------------------------------------------", bcolors.ENDC)

        U = A[:, :N] # matrix A without vector b
        return L, U


def lu_solve(A_b):

    L, U = lu(A_b)
    print("LU is: ", matrix_multiply(L, U))
    inverse_L = np.linalg.inv(L)
    inverse_U = np.linalg.inv(U)
    b = np.array(A_b)[:, -1]
    result = np.dot(inverse_U, inverse_L)
    result = np.dot(result, b)
    print(bcolors.HEADER, "The solution by LU: x= u^-1 * L*-1 * b is ")
    for x in result:
        print("{:.6f}".format(x), )
    print(bcolors.ENDC)
    return result
