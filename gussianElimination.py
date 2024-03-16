""" Functions to quiz 1
Date:19/02/2024
Group: Chaya Mizrachi ID: 214102584, Yael Siman-Tov ID:325181295, Linoy Nisim Pur ID: 324029685
Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
Name:
"""

import numpy as np
from numpy.linalg import norm, inv
from inverse_matrix.matrix_utility import row_addition_elementary_matrix, scalar_multiplication_elementary_matrix,partial_pivoting,swap_rows_elementary_matrix,matrix_multiply,print_matrix
from condition_of_linear_equations import  norm
from colors import bcolors


"""
*) זוהי הפונקציה העיקרית שמתזמרת את האלמנציה של גאוס.
זה מתחיל בקריאה לפעולה forward_substitution(mat), המבצעת את שלבי האלמנציה .
*) אם forward_substitution מחזירה אינדקס k, זה אומר שהמטריצה היא סינגולרית , 
והפונקציה מחזירה הודעה המציינת אם המערכת לא עקבית (כלומר שאין לה פיתרון ) או עשויה להכיל אינסוף פתרונות.
*) אם המטריצה אינה סינגולרית (זה אומר שייש לה פיתרון יחיד ), 
*) היא קוראת לפונקציה backward_substitution(mat) כדי להשיג את הפתרון של המערכת.( את ווקטור הנעלם X אותו אנו מחפשים שמכיל את ערכי המשתנים).
**) (תזכורת: מטריצה סינגולרית הינה מטריצה שאין לה פיתרון או שייש לה אינסוף פיתרונות)

"""
def gaussianElimination(mat):
    if np.linalg.det(mat)==0:
        print("The matrix is singular")
        return
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

""" 
*) פונקציית עזר להחלפת שתי שורות במטריצה. 
*) הפונקציה מקבלת את המטריצה (mat) ואת שני האינדקסים של השורות שיש להחליף (i וְ-j).

"""
# function for elementary operation of swapping two rows
def swap_row(mat, i, j):
    N = len(mat)
    for k in range(N + 1):
        temp = mat[i][k]
        mat[i][k] = mat[j][k]
        mat[j][k] = temp


"""
*)פונקציה זו מטרתה לבצע תת-הקמה קדימה (Forward Substitution) במטריצה מורכבת ממערך משוואות לינאריות (מערכת לינארית).
*) הוא מתחיל בבדיקת "Partial Pivoting"( כלומר רק מחליף בין שורות) - מוצאת שורת הפיבוט (pivot_row) בעמודה הנוכחית שמכילה את הערך המקסימלי בעמודה.
*) אם אחד מהאיברים הראשיים באלכסון שווים לאפס, זה מציין שהמטריצה היא מטריצה סינגולרית (Singular) והפונקציה מחזירה את האינדקס k של השורה הסינגולרית.
 *) אחרת, הוא מחליף את השורה הנוכחית עם  pivot_row  ומבטל את החלק המשולש התחתון של המטריצה 
על ידי הפחתת כפולות של שורת הציר( כלומר הופך את כול העמודה של מתחת לpivot לאפסים).
**) פונקציה זו מבצעת את שלבי האלמנציה של גאוס שהופכים את המטריצה למטריצה משולשת עליונה כך שאיברי הציר- (pivot )
 הם יהיו הכי מקסימליים מבין אותם ערכיים שעומדים באותו עמודה של אותו pivot ( אבל רק בשורות התחתונות לשורה של ה pivot הנ"ל)
"""

"""
יוצרת מטריצה משולשת עליונה
"""
def forward_substitution(mat):
    N=len(mat)
    A=mat.copy()
    for k in range (N):
        result1 = partial_pivoting(A, k, N)
        if isinstance(result1, str):  # Check if partial pivoting returns an error message
            print(result1)
            return k  # Matrix is singular

        A = result1  # Update A with the result of partial pivoting

        for i in range(k + 1, N):
          m = -A[i][k] / A[k][k] #the multiple
          B=row_addition_elementary_matrix(N, i, k,  m)
          C=matrix_multiply(B, A)

          print('Elementary matrix:')
          print_matrix(B)
          print('*')
          print('Original matrix:')
          print_matrix(A)
          print('=')
          print('Result matrix:')
          print_matrix(C)
          print("------------------------------------------------------------------")

          A = C
    mat[:] = A.tolist()  # Update the original matrix with the modified one
    return -1




"""
def forward_substitution(mat):
    
    N = len(mat)
    for k in range(N):

        # Partial Pivoting: Find the pivot row with the largest absolute value in the current column
        # Pivoting חלקי: מצא את השורה הפיבוטית עם הערך המרבי בעמודה הנוכחית

        pivot_row = k
        v_max = mat[pivot_row][k]
        for i in range(k + 1, N):
            if abs(mat[i][k]) > v_max:
                v_max = mat[i][k]
                pivot_row = i

        # if a principal diagonal element is zero,it denotes that matrix is singular,
        # and will lead to a division-by-zero later.
        # אם אחד מהאיברים הראשיים באלכסון הוא אפס, זה מציין שהמטריצה היא סינגולרית,
        # ויביא לחלוקה באפס מאוחר יותר.

        if not mat[k][pivot_row]:# mat[k][pivot_row]==0
            return k  # Matrix is singular

        # Swap the current row with the pivot row
        if pivot_row != k:
            swap_row(mat, k, pivot_row)
        # End Partial Pivoting

        for i in range(k + 1, N):


            #  Compute the multiplier
            #  חישוב "הכופל" 'm'
            m = mat[i][k] / mat[k][k]

            # subtract fth multiple of corresponding kth row element
            # החסר מכפלה שלישית מהאיבר התואם בשורה ה-k

            for j in range(k + 1, N + 1):
            #for j in range(k , N + 1):
                mat[i][j] -= mat[k][j] * m

            # filling lower triangular matrix with zeros
            # מילוי של מטריצה תלת-משולבת עם אפסים

            mat[i][k] = 0

    return -1
"""
"""
*) פונקציה זו מבצעת תת-הקמה לאחור (Backward Substitution) במערך משוואות לינאריות, לאחר שהמטריצה כבר עברה תת-הקמה קדימה.
*) כלומר הפונקציה מקבלת את המטריצה (mat) ומחשבת את ערכי המשתנים (unknowns) באמצעות חישובים אחוריים.
*) מתחילה משורה האחרונה ועוברת על המשוואות לאחור, מחשבת את ערכי המשתנים, ומכניסה אותם למערך הפתרונות X.

"""
# function to calculate the values of the unknowns
def backward_substitution(mat):
    N = len(mat)
    x = np.zeros(N)  # An array to store solution [0,0,....0]

    # Start calculating from last equation up to the first
    # התחל חישובים מהמשוואה האחרונה ועד הראשונה

    for i in range(N - 1, -1, -1):

        x[i] = mat[i][N-1]

        # Initialize j to i+1 since matrix is upper triangular
        #מכיוון שהמטריצה היא עליונה משולבת אתחול של j ל-i+1

        for j in range(i + 1, N):
            x[i] -= mat[i][j] * x[j]

        #x[i] = (x[i] / mat[i][i])

    return x

"""
*) הבלוק if __name__ == '__main__': הוא ריצת התוכנית הראשית.
*) הוא מגדיר מטריצה של מערכת משוואות  (A_b)- זוהי מטריצת המקדמים(A) עם ווקטור התוצאה (b)  .
*) זה קורא gaussianElimination(A_b) כדי לפתור את מערכת המשוואות הלינאריות.
*) אם התוצאה ( result) היא מחרוזת, היא מדפיסה את ההודעה המציינת שהיא סינגולרית מחוסר עקביות( שאין לה פיתרון) או שייש אינסוף פיתרונות.
*) אחרת result אמור להיות רשימה X שמכילה את ווקטור הנעלם ( את ערכי המשתנים הלא ידועים של המטריצה )
 , הוא מדפיס את הפתרון (הרשימה X) בפורמט המוגדר.
"""


if __name__ == '__main__':

    print('''   Date:19/02/2024
    Group: Chaya Mizrachi ID: 214102584, Yael Siman-Tov ID:325181295, Linoy Nisim Pur ID: 324029685
    Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
    Name:Yael Siman Tov''')


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
        num=norm(A)+1


        print("the norm of the matrix A plus 1 is :", num)
