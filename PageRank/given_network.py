import math
import numpy as np


################################   Initialize   ################################
def obtain_A(N):
    A = np.zeros((N, N))
    for i in range(N):
        row = input()
        row = row.strip()
        row = row.split(' ')
        for j in range(N):
            A[i, j] = float(row[j])
    return A

def obtain_n(A, N):
    n = []
    for i in range(N):
        sum = 0.0
        for j in range(N):
            sum += A[i, j]
        n.append(sum)
    return n

def obtain_G(A, n, q, N):
    G = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            G[i, j] = q / N + (1 - q) * A[j, i] / n[j]
    return G

def infinite_norm(x, n):
    norm = 0.0
    p = 0
    for i in range(n):
        xp = math.fabs(x[i, 0])
        if (xp > norm):
            norm = xp
            p = i
    return norm

def sum_vec(x, n):
    sum = 0.0
    for i in range(n):
        sum += x[i, 0]
    return sum

def power_method(G, x, N, TOL = 0.00005, max_it = 100):
    x0 = np.zeros((N, 1))
    k = 1
    x = x / infinite_norm(x, N)
    while (k <= max_it):
        y = G.dot(x)
        miu = infinite_norm(y, N)
        y = y / miu
        if (int(miu * 100000) < 5):
            print("G has the eigenvalue 0, select a new x and restart")
            return (0, x0)
        err = infinite_norm(x - y, N)
        x = y
        if (err < TOL):
            x = x / sum_vec(x, N)
            print(('k = ' + str(k)).ljust(10) + ('miu = ' + str(miu)).ljust(30) + ('err = ' + str(err)).ljust(30))
            print('p = ')
            print(x)
            return (miu, x)
        k += 1
    print("Maximum number of iterations exceeded")
    return (0, x0)


###################   Prove that G is a stochastic matrix   ####################
N = 15
q = 0.15
A = obtain_A(N)
n = obtain_n(A, N)
G = obtain_G(A, n, q, N)
for j in range(N):
    sum = 0.0
    for i in range(N):
        sum += G[i, j]
    print('sum of the %dth column:  %f' %(j, sum))

#################   Verify the given dominant eigenvector p   ##################
p = np.zeros((N, 1))
for i in range(N):
    p[i, 0] = float(input())
Gp = G.dot(p)
err = p - Gp
print(err)

######################   Change the jump probability q   #######################
x0 = np.ones((N, 1))
q = 0
G = obtain_G(A, n, q, N)
(lumda, p1) = power_method(G, x0, N)
q = 0.5
G = obtain_G(A, n, q, N)
(lumda, p2) = power_method(G, x0, N)

#####################   Improve the page rank of Page 7   ######################
q = 0.15
A[2 - 1, 7 - 1] = 2
A[12 - 1, 7 - 1] = 2
n = obtain_n(A, N)
G = obtain_G(A, n, q, N)
(lumda, p) = power_method(G, x0, N)

##############################   Remove Page 10   ##############################
N = 14
x0 = np.ones((N, 1))
B = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        if (i < 9 and j < 9):
            B[i, j] = A[i, j]
        elif (i < 9):
            B[i, j] = A[i, j + 1]
        elif (j < 9):
            B[i, j] = A[i + 1, j]
        else:
            B[i, j] = A[i + 1, j + 1]
n = obtain_n(B, N)
G = obtain_G(B, n, q, N)
(lumda, p) =power_method(G, x0, N)
