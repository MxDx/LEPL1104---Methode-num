import numpy as np

def exercice_1(X0, U0, X1, U1):
    b = (U1 - U0) / (X1 - X0)     
    a = U0 - b*X0
    return a, b

a = [[1, 3, 9], [1, 1, 1], [1, 2, 4]]
b = [2, 1, 1]

np.linalg.solve(a, b)