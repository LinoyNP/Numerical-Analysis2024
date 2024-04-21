"""
Date: 21.4.24
Group: Chaya Mizrachi ID: 214102584,
        Yael Siman-Tov ID:325181295,
        Linoy Nisim Pur ID: 324029685
Source Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
GitHub of this project: https://github.com/LinoyNP/Numerical-Analysis2024.git
Name:Linoy Nisim Pur ID:324029685
"""
import numpy as np

from condition_of_linear_equations import norm
from lu_factorization import lu_solve
def main():
    print('''
            \tDate: 21.4.24
            \tGroup: Chaya Mizrachi ID: 214102584, Yael Siman-Tov ID:325181295, Linoy Nisim Pur ID: 324029685
            \tGit: https://github.com/LinoyNP/Numerical-Analysis2024.git
            \tName: Linoy Nisim Pur''')
    # Section 1
    numQ = 1 # the number of question
    matrix_b = [[1, 1/2, 1/3, 1], [1/2, 1/3, 1/4, 0], [1/3, 1/4, 1/5, 0]]
    lu_solve(matrix_b)
    matrix = np.array(matrix_b)[:, :-1]
    matrix_norm = norm(matrix)
    print("Matrix's norm + question numer is: ", matrix_norm + numQ)

    # Section 2
    print('\n\n')
    numQ = 2
    matrix_b = [[5, 1, 10, 3/2], [10, 8, 1, -7], [4, 10, -5, 2]]
    matrix = np.array(matrix_b)[:, :-1] # matrix without vector b
    lu_solve(matrix_b)
    matrix_norm = norm(matrix)
    print("Matrix's norm + question numer is: ", matrix_norm + numQ)
if __name__ == '__main__':
    main()