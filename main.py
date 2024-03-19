"""
Date: 18.3.24
Group: Chaya Mizrachi ID: 214102584,
        Yael Siman-Tov ID:325181295,
        Linoy Nisim Pur ID: 324029685
Source Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
GitHub of this project: https://github.com/LinoyNP/Quiz2.git
Name:Linoy Nisim Pur ID:324029685
"""
from gussianElimination import gaussianElimination
from lu_factorization import lu_solve ,lu
from Iterative_methods.bisection_method import bisection_method
from colors import bcolors

def main():
    print('''
    \tDate: 18.3.24
    \tGroup: Chaya Mizrachi ID: 214102584, Yael Siman-Tov ID:325181295, Linoy Nisim Pur ID: 324029685
    \tGit: https://github.com/LinoyNP/Quiz2.git
    \tName: Linoy Nisim Pur''')
    Matrix=[[2, 3, 4, 5, 6, 70],
            [-5, 3, 4, -2, 3, 20],
            [4, -5, -2, 2, 6, 26],
            [4, 5, -1, -2, -3, -12],
            [5, 5, 3, -3, 5, 37]]
    
    result = gaussianElimination(Matrix)
    if isinstance(result, str):
        print(result)
    else:
        print(bcolors.OKBLUE, "\nSolution for the system by Gussian Elimination :")
        for x in result:
            print("{:.6f}".format(x))  # דיוק של 6 ספרות אחרי הנקודה (נקודה צפה)"""""

    SOL = lu_solve(Matrix)

    #Section B
    #Note: The solution did not come out right - but for the sake of the examiner's completeness I will solve section b with the solution I received
    f = lambda x: ( SOL[4] * x ** 3 + SOL[0] * x ** 2 + SOL[2]) / (SOL[2] * x -SOL[4])
    roots = bisection_method(f, 1, 3)
    print(bcolors.OKBLUE, f"\nThe equation f(x) has an approximate root at x = {roots}", bcolors.ENDC, )

    

if __name__ == '__main__':
    main()