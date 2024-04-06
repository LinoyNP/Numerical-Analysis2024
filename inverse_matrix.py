from colors import bcolors
from matrix_utility import row_addition_elementary_matrix, scalar_multiplication_elementary_matrix
import numpy as np

"""
Function that find the inverse of non-singular matrix
The function performs elementary row operations to transform it into the identity matrix. 
The resulting identity matrix will be the inverse of the input matrix if it is non-singular.
 If the input matrix is singular (i.e., its diagonal elements become zero during row operations), it raises an error.
"""

"מוצאת את המטריצה ההופכית למטריצה שנתונה באמצעות פעולות שורה אלמנטריות."
def inverse(matrix):
    print( f"=================== Finding the inverse of a non-singular matrix using elementary row operations ===================\n {matrix}\n")
    if np.linalg.det(matrix)==0:
        raise ValueError("The matrix is singular")
    if matrix.shape[0] != matrix.shape[1]:  # בדיקה אם המטריצה היא ריבועית
        raise ValueError("Input matrix must be square.")
    n = matrix.shape[0]    # קביעת גודל המטריצה (n) כגודל השורות של המטריצה

    identity = np.identity(n)#  יצירת מטריצת יחידה בגודל המתאים למטריצה
    for i in range(n):#ביצוע פעולות שורה על המטריצה כך שתהפוך למטריצת יחידה
        if matrix[i, i] == 0: # בדיקה אם איבר האלכסון הוא אפס - אם כן, המטריצה היא סינגולרית
            raise ValueError("Matrix is singular, cannot find its inverse.")

        if matrix[i, i] != 1: # אם איבר האלכסון הוא שונה מ־1, נצטרך לעשות פעולות כדי להפוך אותו ל־1
            scalar = 1.0 / matrix[i, i] # חישוב הסקלר כך שיהיה אפשר לכפול אותו ולקבל 1
            elementary_matrix = scalar_multiplication_elementary_matrix(n, i, scalar)# יצירת מטריצה אלמנטרית שתכפיל את השורה הנוכחית כך שאיבר האלכסון יהיה 1

            print(f"elementary matrix to make the diagonal element 1 :\n {elementary_matrix} \n")
            matrix = np.dot(elementary_matrix, matrix) # כפל המטריצה האלמנטרית במטריצה המקורית
            print(f"The matrix after elementary operation :\n {matrix}")
            print("------------------------------------------------------------------------------------------------------------------")
            identity = np.dot(elementary_matrix, identity) # כפל המטריצה האלמנטרית במטריצת היחידה

            # אפס את האיברים מעל ומתחת לאיבר האלכסון
            for j in range(n):
                if i != j:
                    scalar = -matrix[j, i]  # חישוב ערך הסקלר על ידי הכפלת האיבר במיקום (j, i) במטריצה ב־(-1)
                    elementary_matrix = row_addition_elementary_matrix(n, j, i, scalar)
                    print(f"elementary matrix for R{j+1} = R{j+1} + ({scalar}R{i+1}):\n {elementary_matrix} \n")  # יצירת מטריצה אלמנטרית שתוסיף לשורה הנוכחית את השורה האחרת, כך שנוכל לאפס את האיבר הנוכחי

                    # כפל המטריצה האלמנטרית במטריצה המקורית
                matrix = np.dot(elementary_matrix, matrix)
                print(f"The matrix after elementary operation :\n {matrix}")
                print("------------------------------------------------------------------------------------------------------------------")
                # כפל המטריצה האלמנטרית במטריצת היחידה
                identity = np.dot(elementary_matrix, identity)

    return identity

if __name__ == '__main__':
    Matrix=np.array([[2, 3, 4, 5, 6],
            [-5, 3, 4, -2, 3],
            [4, -5, -2, 2, 6],
            [4, 5, -1, -2, -3],
            [5, 5, 3, -3, 5]])
    print(inverse(Matrix))
