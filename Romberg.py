import math
import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

a = 0.0
b = 1.0


def Romberg(k, j, R, F):
    if (F[k-1, j-1]):
        return R[k-1, j-1]
    if (k == 1):
        h = b - a
        f0 = f1 = 0.0
        r = (h / 2) * (f0 + f1)
    elif (j == 1):
        r = 0.0
        h = (b - a) / math.pow(2, k-1)
        n = math.pow(2, k-2)
        for i in range(int(n)):
            x = a + (2 * i + 1) * h
            f = math.sqrt(x) * math.log(x)
            r += f
        r *= 2 * h
        r += Romberg(k-1, 1, R, F)
        r /= 2
    else:
        r = Romberg(k, j-1, R, F) + (Romberg(k, j-1, R, F) - Romberg(k-1, j-1, R, F)) / (math.pow(4, j-1) - 1)
    F[k-1, j-1] = 1
    R[k-1, j-1] = r
    return r


def Err_function(a, b, n):
    N = np.arange(1, n + 1, 1, dtype=int)
    R = np.matlib.zeros((n, n))
    F = R
    r = Romberg(n, n, R, F)
    E = R + 4/9
    fig = figure()
    ax = Axes3D(fig)
    K = np.arange(1, n + 1, 1, dtype=int)
    print(K)
    J = np.arange(1, n + 1, 1, dtype=int)
    K, J = np.meshgrid(K, J)
    ax.set_xlabel('k')
    ax.set_ylabel('j')
    ax.set_zlabel('Error: R(k, j)(x) - f(x)')
    ax.plot_surface(K, J, E, cmap='hot')
    show()
    return

if __name__ == "__main__":
    Err_function(0, 1, 20)
