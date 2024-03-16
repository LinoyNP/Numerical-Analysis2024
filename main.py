"""
Date:
Group: Chaya Mizrachi ID: 214102584,
        Yael Siman-Tov ID:325181295,
        Linoy Nisim Pur ID: 324029685
Source Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
GitHub of this project:
Name:Linoy Nisim Pur ID:324029685
"""
from gussianElimination import gaussianElimination
from lu_factorization import lu_solve
from colors import bcolors
def main():
    print('''
    \tDate:
    \tGroup: Chaya Mizrachi ID: 214102584, Yael Siman-Tov ID:325181295, Linoy Nisim Pur ID: 324029685
    \tGit: 
    \tName: Linoy Nisim Pur''')
    Matrix=[[1, -1, 2, -1],
        [2, -2, 3, -3],
        [1, 1, 1, 0],
        [1, -1, 4, 3]]
    """[[1, 2, 3, 4, 5],
            [2, 3, 4, 5, 1],
            [8, 8, 8, 8, 1],
            [24, 15, 22, 1, 8]]
    result = gaussianElimination(Matrix)
    if isinstance(result, str):
        print(result)
    else:
        print(bcolors.OKBLUE, "\nSolution for the system:")
        for x in result:
            print("{:.6f}".format(x))  # דיוק של 6 ספרות אחרי הנקודה (נקודה צפה)"""

    lu_solve(Matrix)
    

if __name__ == '__main__':
    main()