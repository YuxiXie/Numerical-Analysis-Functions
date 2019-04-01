import math
import numpy as np
from numpy import matlib


def obtain_A():
    N = int(input())
    A = matlib.zeros((N, N))
    for i in range(N):
        row = input()
        row = row.strip()
        row = row.split(' ')
        for j in range(N):
            A[i, j] = float(row[j])
    return (A, N)


def sum_off_diagonal(A, N):
    sum = 0.0
    for i in range(N):
        for j in range(N):
            if (i != j):
                sum += math.fabs(A[i, j])
    return sum


def sgn(x):
    if (x > 0):
        return 1
    elif (x < 0):
        return -1
    else:
        return 0


def Jacobi(A, N, TOL = 0.000000000005, max_it = 100):
    it = 1
    while (it <= max_it):
        for i in range(N - 1):
            j = i + 1
            for k in range(j):
                P = matlib.eye(N, N)
                if (A[j, k] != A[k, j]):
                    c = 2 * A[j, k] * sgn(A[j, j] - A[k, k])
                    b = math.fabs(A[j, j] - A[k, k])
                    P[j, j] = P[k, k] = math.sqrt((1 + b / math.sqrt(c*c + b*b)) / 2)
                    P[k, j] = c / (2 * P[j, j] + math.sqrt(c*c + b*b))
                    P[j, k] = - P[k, j]
                else:
                    P[j, j] = P[k, k] = math.sqrt(1/2)
                    P[k, j] = P[j, j]
                    P[j, k] = - P[k, j]
                A = np.transpose(P).dot(A).dot(P)
        if (sum_off_diagonal(A, N) < TOL):
            print('it = ' + str(it) + ',    A = ')
            for r in range(N):
                row = ' '
                for c in range(N):
                    row += (' ' + str(A[r, c])).rjust(28)
                print(row)
            return A
        it += 1
    print("Maximum number of iterations exceeded")
    return

if __name__ == "__main__":
    (A, N) = obtain_A()
    A = Jacobi(A, N)
    for i in range(N):
        print('miu%d = %.11f'%(i+1, A[i, i]))
