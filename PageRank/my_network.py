import math
import numpy as np
import matplotlib.pyplot as plt


def obtain_A():
    A = np.load('adjacency_matrix0.npy')
    return A

def obtain_n(A, N):
    n = []
    for i in range(N):
        sum = 0.0
        for j in range(N):
            sum += A[i, j]
        n.append(sum)
    return n

def obtain_m(A, N):
    m = []
    for i in range(N):
        sum = 0.0
        for j in range(N):
            sum += A[j, i]
        m.append(sum)
    return m

def obtain_G(A, n, q, N):
    G = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            G[i, j] = q / N
            if (A[j, i] > 0):
                G[i, j] += (1 - q) * A[j, i] / n[j]
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

def power_method(G, x, N, TOL = 0.000005, max_it = 100):
    x0 = np.zeros((N, 1))
    k = 1
    while (k <= max_it):
        y = G.dot(x)
        miu = infinite_norm(y, N)
        y /= miu
        if (int(miu * 1000000) < 5):
            print("G has the eigenvalue 0, select a new x and restart")
            return (0, x0)
        err = infinite_norm(x - y, N)
        x = y
        if (err < TOL):
            x = x / sum_vec(x, N)
            print(('k = ' + str(k)).ljust(10) + ('miu = ' + str(miu)).ljust(30) + ('err = ' + str(err)).ljust(30))
            print(x)
            return (miu, x)
        k += 1
    print("Maximum number of iterations exceeded")
    return (0, x0)

def sort_PR(PR, N):
    pages = []
    for i in range(N):
        p = (PR[i, 0], i)
        pages.append(p)
    pages = sorted(pages, key=lambda x:x[0], reverse=True)
    return pages

def display_result(pages, n, N, filename):
    id_list = []
    rank_list = []
    pr_list = []
    innum_list = []
    fpr = open(filename, 'w')
    for i in range(N):
        id_list.append(pages[i][1])
        rank_list.append(i + 1)
        pr_list.append(pages[i][0])
        innum_list.append(n[pages[i][1]])
        fpr.write(('id = ' + str(pages[i][1])).ljust(20) + ('PR = ' + str(pages[i][0])).ljust(40) + '\n')
    fpr.close()

    plt.plot(rank_list, pr_list)
    plt.ylabel('PageRank')
    plt.xlabel('Ranking of pages')
    plt.show()

    cValue = 'r'
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.ylabel('Number of inlinks')
    plt.xlabel('Ranking of pages')
    ax1.scatter(rank_list, innum_list, c = cValue, marker = 'o')
    plt.show()
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    plt.ylabel('Number of inlinks')
    plt.xlabel('PageRank')
    ax2.scatter(pr_list, innum_list, c = cValue, marker = 'o')
    plt.show()
    fig = plt.figure()
    ax3 = fig.add_subplot(111)
    plt.ylabel('Ranking of pages')
    plt.xlabel('Page ID')
    ax3.scatter(id_list, rank_list, c = cValue, marker = 'o')
    plt.show()


################################   Initialize   ################################
N = 30000
q = 0.15
A = obtain_A()
n = obtain_n(A, N)
G = obtain_G(A, n, q, N)

##############################   Get page ranks   ##############################
pr = np.ones((N, 1))
(lumda, PR) = power_method(G, pr, N)
pages = sort_PR(PR, N)
m = obtain_m(A, N)
display_result(pages, m, N, 'PR_result.txt')

######################   Change the jump probability q   #######################
q = 0
G = obtain_G(A, n, q, N)
(lumda, PR1) = power_method(G, pr, N)
pages1 = sort_PR(PR1, N)
m = obtain_m(A, N)
display_result(pages1, m, N, 'PR1_result.txt')

q = 0.5
G = obtain_G(A, n, q, N)
(lumda, PR2) = power_method(G, pr, N)
pages2 = sort_PR(PR2, N)
m = obtain_m(A, N)
display_result(pages2, m, N, 'PR2_result.txt')

# #####################   Improve the page rank of Page x   ######################
x = 1060
q = 0.15
B = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        B[i, j] = A[i, j]
        if (j == x - 1):
            B[i, j] *= 100
n1 = obtain_n(B, N)
G = obtain_G(B, n1, q, N)
(lumda, PR3) = power_method(G, pr, N)
pages3 = sort_PR(PR3, N)
m = obtain_m(B, N)
display_result(pages3, m, N, 'PR3_result.txt')

# ##############################   Remove Page y   ###############################
y = 26875
N -= 1
pr = np.ones((N, 1))
C = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        if (i < y - 1 and j < y - 1):
            C[i, j] = A[i, j]
        elif (i < y - 1):
            C[i, j] = A[i, j + 1]
        elif (j < y - 1):
            C[i, j] = A[i + 1, j]
        else:
            C[i, j] = A[i + 1, j + 1]
n2 = obtain_n(C, N)
G = obtain_G(C, n2, q, N)
(lumda, PR4) = power_method(G, pr, N)
pages4 = sort_PR(PR4, N)
m = obtain_m(C, N)
display_result(pages4, m, N, 'PR4_result.txt')
