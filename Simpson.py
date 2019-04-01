import math
import matplotlib.pyplot as plt
import numpy as np


def Simpson(a, b, n):
    h = (b - a) / n
    f0 = f1 = 0.0
    S = f0 + f1
    for i in range(int(n/2 - 1)):
        j = (i + 1) * 2
        x = a + h * j
        S += 2 * math.sqrt(x) * math.log(x)
    for i in range(int(n/2)):
        j = (i + 1) * 2 - 1
        x = a + h * j
        S += 4 * math.sqrt(x) * math.log(x)
    S *= h / 3
    return S


def Err_function(a, b):
    N = np.arange(4, 1000, 2, dtype=int)
    H = (b - a) / N
    S = []
    for n in N:
        h = (b - a) / n
        s = Simpson(a, b, n)
        S += [s]
    S = np.asarray(S)
    E = S + 4/9
    plt.plot(H, E)
    plt.xlabel('h')
    plt.ylabel('Error: Simpson_result - (-4/9)')
    plt.show()

    plt.plot(N, E)
    plt.xlabel('n')
    plt.ylabel('Error: Simpson_result - (-4/9)')
    plt.show()
    return

if __name__ == "__main__":
    Err_function(0, 1)
