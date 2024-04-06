import numpy as np

from colors import bcolors
from matrix_utility import *
from gussianElimination import *
from Iterative_methods_matrix.gauss_seidel import *





def solveMatrix(matrixA,vectorb):
    detA = np.linalg.det(matrixA)
    #print(bcolors.YELLOWBG, "\nDET(A) = ", detA)

    if detA != 0:
        #print("CondA = ", Cond(matrixA, InverseMatrix(matrixA, vectorb)), bcolors.ENDC)
        print(bcolors.OKBLUE, "\nnon-Singular Matrix - Perform GaussJordanElimination",bcolors.ENDC)
        for line in range(len(matrixA)):  # Adding vectorb to matrix
            matrixA[line].append(vectorb[line][0])
        #print(matrixA)
        result = gaussianElimination(matrixA)
        #print("The result is: ", np.array(result))
        return result
    else:
        X0 = np.zeros_like(vectorb)
        result = gauss_seidel(matrixA, vectorb, X0, 1e-16, 200)
        vector_result =[]
        for element in result:
            vector_result.append(element)
        print("The result is: ", np.array(vector_result))
        return vector_result
        """print("Singular Matrix - Perform LU Decomposition\n")
        print("Matrix U: \n")
        print(np.array(UMatrix(matrixA, vectorb)))
        print("\nMatrix L: \n")
        print(np.array(LMatrix(matrixA, vectorb)))
        print("\nMatrix A=LU: \n")
        result = matrix_multiply(LMatrix(matrixA, vectorb), UMatrix(matrixA, vectorb)) #change
        print(np.array(result))
        return result"""


def polynomialInterpolation(table_points, x):
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points] # Makes the initial matrix

    b = [[point[1]] for point in table_points]

    print(bcolors.OKBLUE, "The matrix obtained from the points: ", bcolors.ENDC,'\n', np.array(matrix))
    print(bcolors.OKBLUE, "\nb vector: ", bcolors.ENDC,'\n',np.array(b))
    matrixSol = solveMatrix(matrix, b)

    result = sum([matrixSol[i] * (x ** i) for i in range(len(matrixSol))])
    print(bcolors.OKBLUE, "\nThe polynom:", bcolors.ENDC)
    print('P(X) = '+'+'.join([ '('+str(matrixSol[i])+') * x^' + str(i) + ' ' for i in range(len(matrixSol))])  )
    print(bcolors.OKGREEN, f"\nThe Result of P(X={x}) is:", bcolors.ENDC)
    print(result)
    #return result


if __name__ == '__main__':

    """table_points = [(0, 0), (1, 0.8415), (2, 0.9093), (3, 0.1411), (4, -0.7568), (5, -0.9589), (6, -0.2794)]
    x = 1.28"""
    table_points = [(1, 3), (2, 4), (3, -1)]
    x = 1.5
    print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------\n", bcolors.ENDC)
    print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
    print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x,'\n')
    polynomialInterpolation(table_points, x)
    print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n", bcolors.ENDC)
