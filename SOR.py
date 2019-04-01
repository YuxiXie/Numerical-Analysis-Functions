import math
import numpy as np
from numpy import matlib
from numpy import mat
from numpy import linalg


def initialize(n):
    A = matlib.identity(n, dtype=float)
    b = np.zeros(n)
    for i in range(n):
        b[i] = math.pi
        A[i, i] = 2 * (i + 1)
        if (i >= 2):
            A[i, i-2] = 0.5 * (i + 1)
            if (i >= 4):
                A[i, i-4] = 0.25 * (i + 1)
        if (i < n - 2):
            A[i, i+2] = 0.5 * (i + 1)
            if (i < n - 4):
                A[i, i+4] = 0.25 * (i + 1)
    return (A, b)


def SOR(n, TOL, N, x0, w):
    (A, b) = initialize(n)
    x = np.zeros(n)
    k = 0
    while (k <= N):
        err = 0
        for i in range(n):
            y = 0.0
            for j in range(n):
                if (j < i):
                    y += A[i, j] * x[j]
                elif (j > i):
                    y += A[i, j] * x0[j]
            x[i] = (1 - w) * x0[i] + (w / A[i, i]) * (-y + b[i])
            e = math.fabs(x[i] - x0[i])
            if (e > err):
                err = e
        if (err < TOL):
            b = np.transpose(mat(b))
            y = linalg.inv(A).dot(b)
            f = 0
            for i in range(n):
                e = math.fabs(x[i] - y[i, 0])
                if (e > f):
                    f = e
            k += 1
            # print('k = ' + str(k))
            print('x = ')
            print(x)
            # print(('error1 = ' + ("%.15f" %err)).ljust(30) + ('error2 = ' + ("%.15f" %f)).ljust(30))
            return
        k += 1
        x0 = x
    print("Maximum number of iterations exceeded")

if __name__ == "__main__":
    x0 = np.zeros(80)
    SOR(80, 0.00001, 10, x0, 1.25)
