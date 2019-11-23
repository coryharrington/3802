import numpy as np
import copy as cp
from math import *

def construct_a(a, g_s, g_g):
    for i in range(n):
        for j in range(n):
            if i == 0:
                if j == 0:
                    a[i][j] = g_s[0]
                elif j == 1:
                    a[i][j] = -g_s[0]
                else:
                    a[i][j] = 0
            elif i < (n - 1):
                if j < (i - 1):
                    a[i][j] = 0
                elif j == (i - 1):
                    a[i][j] = g_s[i - 1]
                elif j == i:
                    a[i][j] = -1 * (g_s[i - 1] + g_s[i] + g_g[i - 1])
                elif j == (i + 1):
                    a[i][j] = g_s[i]
            else:
                if j < (n - 2):
                    a[i][j] = 0
                elif j == (n - 2):
                    a[i][j] = g_s[i - 1]
                else:
                    a[i][j] = -1 * (g_s[i - 1] + g_g[i - 1])

def gauss_sidel(n, a, x, b, tol):
    L = cp.deepcopy(a)
    U = cp.deepcopy(a)
    for i in range(n):
        for j in range(n):
            if j > i:
                L[i][j] = 0
                U[i][j] = a[i][j]
            else:
                L[i][j] = a[i][j]
                U[i][j] = 0
    x_new = cp.deepcopy(x)
    c = 0
    max_iter = 10000
    while((c == 0 or np.linalg.norm(x_new - x) > tol) and c < max_iter):
        x = cp.deepcopy(x_new)
        x_new = np.linalg.inv(L).dot(b - U.dot(x))
        c += 1
    return x_new

if __name__ == "__main__":
    n = 100
    # nxn matrix
    a = np.array([np.zeros(n) for i in range(n)])
    #r_s = [400]
    #r_g = [50]
    #g_s = [1/400]
    #g_g = [1/50]
    #for i in range(1, n):
    #    r_s.append(r_s[i - 1] / 4)
    #    g_s.append(1 / r_s[i])
    #    r_g.append(r_g[i - 1] / 2)
    #    g_g.append(1 / r_g[i])
    g_s = [1 for i in range(n)]
    g_g = [1 for i in range(n)]
    construct_a(a, g_s, g_g)
    print("a, is : ", a, "\n\n")
    b = np.zeros(n)
    b[0] = 1
    x = np.linalg.solve(a, b)
    print("x, analytical is : ", x)
    y = gauss_sidel(n, a, np.zeros(n), b, 0.1)
    print("x, numerical is : ", y)
    print("rms is : ", np.linalg.norm(y - x))
    g = (1 + sqrt(5)) / 2
    print("golden ratio approximation percent error is : ", abs(y[0] - g) / g)

