""" Functions to quiz 1
Date:19/02/2024
Group: Chaya Mizrachi ID:, Yael Siman-Tov ID:, Linoy Nisim Pur ID: 
Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
Name:
"""

import numpy as np
from numpy.linalg import norm, inv
from matrix_utility import row_addition_elementary_matrix, backward_substitution,partial_pivoting,swap_rows_elementary_matrix,matrix_multiply,print_matrix,swap_row
from condition_of_linear_equations import  norm
from colors import bcolors


def gaussianElimination(mat):
    """if np.linalg.det(mat)==0:
        print("The matrix is singular")
        return"""
    N = len(mat)

    singular_flag = forward_substitution(mat)


    if singular_flag != -1:

        if mat[singular_flag][N]:
            return "Singular Matrix (Inconsistent System)"
        else:
            return "Singular Matrix (May have infinitely many solutions)"

    # if matrix is non-singular: get solution to system using backward substitution
    # אם המטריצה אינה סינגולרית: קבל פתרון למערכת באמצעות התאמה לאחור

    return backward_substitution(mat)


def forward_substitution(mat):
    N=len(mat)
    A=mat.copy()
    for k in range (N):
        result1 = partial_pivoting(A, k, N)
        if isinstance(result1, str):  # Check if partial pivoting returns an error message
            print(result1)
            return k  # Matrix is singular
        print("Matrix after pivoting: ", result1)
        A = result1  # Update A with the result of partial pivoting

        for i in range(k + 1, N):
          m = -A[i][k] / A[k][k] #the multiple
          B=row_addition_elementary_matrix(N, i, k,  m)
          C=np.dot(B, A)
          #C=matrix_multiply(B, A)

          """print('Elementary matrix:')
          print_matrix(B)
          print('*')
          print('Original matrix:')
          print_matrix(A)
          print('=')
          print('Result matrix:')
          print_matrix(C)
          print("------------------------------------------------------------------")"""

          A = C
    mat[:] = A.tolist()  # Update the original matrix with the modified one
    return -1



"""if __name__ == '__main__':

    print('''   Date:
    Group: Chaya Mizrachi ID: , Yael Siman-Tov ID:, Linoy Nisim Pur ID: 
    Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
    Name:''')


    A_b = [[1, 0.5, (1/3), 1],
        [0.5, (1/3), (1/4), 0],
        [(1/3),(1/4), (1/5), 0],]

    A= [[1, 0.5, (1 / 3)],
           [0.5, (1 / 3), (1 / 4)],
           [(1 / 3), (1 / 4), (1 / 5)], ]

    result = gaussianElimination(A_b)
    if isinstance(result, str):
        print(result)
    else:
        print(bcolors.OKBLUE,"\nSolution for the system:")
        for x in result:
            print("{:.6f}".format(x)) # דיוק של 6 ספרות אחרי הנקודה (נקודה צפה)
 """
