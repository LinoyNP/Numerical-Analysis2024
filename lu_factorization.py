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
            v_max = abs(A[pivot_row][i])
            for j in range(i + 1, N):
                if abs(A[j][i]) > v_max:
                    v_max = A[j][i]
                    pivot_row = j

            # if a principal diagonal element is zero,it denotes that matrix is singular,
           ## # and will lead to a division-by-zero later.
            if A[pivot_row][i] == 0:
                raise ValueError("can't perform LU Decomposition")

            # Swap the current row with the pivot row
            if pivot_row != i:
                e_matrix = swap_rows_elementary_matrix(N, i, pivot_row)
                print(f"elementary matrix for swap between row {i} to row {pivot_row} :\n {e_matrix} \n")
                A = np.dot(e_matrix, A)
                #A= matrix_multiply(e_matrix, A)
                print(f"The matrix after elementary operation :\n {A}")
                print(bcolors.OKGREEN,"---------------------------------------------------------------------------", bcolors.ENDC)

            for j in range(i + 1, N):
    
                #  Compute the multiplier
                m = -A[j][i] / A[i][i]
                e_matrix = row_addition_elementary_matrix(N, j, i, m)
                e_inverse = np.linalg.inv(e_matrix)
                L = np.dot(L, e_inverse) #L=l^-1 * l^-2 *... l^-n
                #A = np.dot(e_matrix, A)
                A = matrix_multiply(e_matrix, A)
                #print(f"elementary matrix to zero the element in row {j} below the pivot in column {i} :\n {e_matrix} \n")
                #print(f"The matrix after elementary operation :\n {A}")
                #print(bcolors.OKGREEN,"---------------------------------------------------------------------------", bcolors.ENDC)

        U = A
        return L, U


# function to calculate the values of the unknowns
def backward_substitution(mat):
    N = len(mat)
    x = np.zeros(N)  # An array to store solution

    # Start calculating from last equation up to the first
    for i in range(N - 1, -1, -1):

        x[i] = mat[i][N-1]

        # Initialize j to i+1 since matrix is upper triangular
        for j in range(i + 1, N):
            x[i] -= mat[i][j] * x[j]

        #x[i] = (x[i] / mat[i][i])

    return x

def lu_solve(A_b):
    #check if A_b is invertible matrix
    if np.linalg.det(A_b)==0:
        print("The matrix is singular")
        return
    L, U = lu(A_b)
    #print(matrix_multiply(L, U))
    """if (L == 0) and (U == 1):
        print("Not possible to solve LU")
        return"""
    print("Lower triangular matrix L:\n", L)
    print("Upper triangular matrix U:\n", U)

    #find matrix solution- vector b
    result = backward_substitution(U)
    print(bcolors.OKBLUE,"\nSolution for the system:")
    for x in result:
        print("{:.6f}".format(x))


"""if __name__ == '__main__':


    A_b = [[1, 2, 3],
           [4, 5, 6],
           [7, 8, 9]]
    [[1, -1, 2, -1],
        [2, -2, 3, -3],
        [1, 1, 1, 0],
        [1, -1, 4, 3]]
    A_b=[[1, 2, 3, 4, 5],
    [2, 3, 4, 5, 1],
    [8, 8, 8, 8, 1],
    [24, 15, 22, 1, 8]]


    lu_solve(A_b)"""