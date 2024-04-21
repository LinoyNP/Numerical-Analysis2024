"""
Date:
Group: Chaya Mizrachi ID: 214102584, Yael Siman-Tov ID:325181295, Linoy Nisim Pur ID: 324029685
Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
Name:
"""
import numpy as np
from inverse_matrix import inverse
from colors import bcolors
from matrix_utility import print_matrix
"""
על מנת לחשב את טיב ההצגה יש לחשב את:
 * הנורמה של מטריצה A
 * מטריצה הפיכה 
 * נורמה של מטריצה הפיכה
 ונחזיר את הכפל בין המטריצות
"""
# find norm of matrix A
def norm(mat):
    size = len(mat)  # size of matrix
    max_row = 0  # In max_row will be the max sum between rows- יכיל את סכום השורה הגבוה ביותר
    for row in range(size):
        sum_row = 0  # In sum_row will be the sum of each row -סכום כל שורה
        for col in range(size):  # On this for make sum to one row- בלולאה זו סוכמים את השורה
            sum_row += abs(mat[row][col])  # Adding element to sum - מוסיפים את הערך המוחלט של האיבר לסכום השורה בה הוא נמצא
        if sum_row > max_row:  # check if the sum of current row is bigger than sum of previous row- בודקים האם הסכום שחישבנו עכשיו גדול יותר מהסכום השמור
            max_row = sum_row  # enter to max_row the biggest sum- עדכון סכום השורה הגבוה
    return max_row  #return max sum of row-מחזירים את סכום השורה הגבוה ביותר
def condition_number(A):
    # Step 1: Calculate the max norm (infinity norm) of A -חישוב הנורמה של המטריצה ע"י פונקציה
    norm_A = norm(A)

    # Step 2: Calculate the inverse of A - חישוב הפונקציה ההפיכה ע"י פונקציה
    A_inv = inverse(A)

    # Step 3: Calculate the max norm of the inverse of A- חישוב הנורמה של המטריצה ההפיכה ע"י שלחית המטריצה ההפיכה לפונקציה
    norm_A_inv = norm(A_inv)

    # Step 4: Compute the condition number -חישוב הקונד של המטריצה
    cond = norm_A * norm_A_inv  # Multiply the norm of A with the norm of reverse A- מכפלת הנורמה של המטריצה עם הנורמה של המטריצה ההפיכה

    """print(bcolors.OKBLUE, "A:", bcolors.ENDC)
    print_matrix(A)

    print(bcolors.OKBLUE, "inverse of A:", bcolors.ENDC)
    print_matrix(A_inv)

    print(bcolors.OKBLUE, "Max Norm of A:", bcolors.ENDC, norm_A, "\n")

    print(bcolors.OKBLUE, "max norm of the inverse of A:", bcolors.ENDC, norm_A_inv)"""

    return cond

def main():
    A = np.array([[2, 1.7, -2.5],
                  [1.24, -2, -0.5],
                  [3, 0.2, 1]])
    """A = np.array([[1, 2, 3],
                [2, 3, 4],
                [3, 4, 6]])"""
    cond = condition_number(A)

    print(bcolors.OKGREEN, "\n condition number: ", cond, bcolors.ENDC)
if __name__ == '__main__':
    main()

