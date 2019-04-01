import math
import numpy as np
from numpy import matlib
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D


f = - 1.5 / 1.04

def g(x, y):
    g = x * (6 - x) + y * (5 - y)
    return g

def PEFD(a, b, c, d, n, m, TOL = 0.0000000001, max_it = 1000):
    h = (b - a) / n
    k = (d - c) / m
    x = np.arange(a, b + h, h, dtype = float)
    y = np.arange(c, d + k, k, dtype = float)
    w = matlib.zeros((n + 1, m + 1))
    lumda = math.pow(h/k, 2)
    miu = 2 * (1 + lumda)
    it = 1
    while (it <= max_it):
        z = (-h*h*f + g(a, y[m-1]) + lumda*g(x[1], d) + lumda*w[1, m-2] + w[2, m-1]) / miu
        norm = math.fabs(z - w[1, m-1])
        w[1, m-1] = z
        for l in range(n - 3):
            i = l + 2
            z = (-h*h*f + lumda*g(x[i], d) + w[i-1, m-1] + w[i+1, m-1] + lumda*w[i, m-2]) / miu
            if (math.fabs(w[i, m-1] - z) > norm):
                norm = math.fabs(w[i, m-1] - z)
            w[i, m-1] = z
        z = (-h*h*f + g(b, y[m-1]) + lumda*g(x[n-1], d) + w[n-2, m-1] + lumda*w[n-1, m-2]) / miu
        if (math.fabs(w[n-1, m-1] - z) > norm):
            norm = math.fabs(w[n-1, m-1] - z)
        w[n-1, m-1] = z
        for l in range(m - 3):
            j = m - 2 - l
            z = (-h*h*f + g(a, y[j]) + lumda*w[1, j+1] + lumda*w[1, j-1] + w[2, j]) / miu
            if (math.fabs(w[1, j] - z) > norm):
                norm = math.fabs(w[1, j] - z)
            w[1, j] = z
            for p in range(n - 3):
                i = p + 2
                z = (-h*h*f + w[i-1, j] + lumda*w[i, j+1] + w[i+1, j] + lumda*w[i, j-1]) / miu
                if (math.fabs(w[i, j] - z) > norm):
                    norm = math.fabs(w[i, j] - z)
                w[i, j] = z
            z = (-h*h*f + g(b, y[j]) + w[n-2, j] + lumda*w[n-1, j+1] + lumda*w[n-1, j-1]) / miu
            if (math.fabs(w[n-1, j] - z) > norm):
                norm = math.fabs(w[n-1, j] - z)
            w[n-1, j] = z
        z = (-h*h*f + g(a, y[1]) + lumda*g(x[1], c) + lumda*w[1, 2] + w[2, 1]) / miu
        if (math.fabs(w[1, 1] - z) > norm):
            norm = math.fabs(w[1, 1] - z)
        w[1, 1] = z
        for l in range(n - 3):
            i = l + 2
            z = (-h*h*f + lumda*g(x[i], c) + w[i-1, 1] + lumda*w[i, 2] + w[i+1, 1]) / miu
            if (math.fabs(w[i, 1] - z) > norm):
                norm = math.fabs(w[i, 1] - z)
            w[i, 1] = z
        z = (-h*h*f + g(b, y[1]) + lumda*g(x[n-1], c) + w[n-2, 1] + lumda*w[n-1, 2]) /miu
        if (math.fabs(w[n-1, 1] - z) > norm):
            norm = math.fabs(w[n-1, 1] - z)
        w[n-1, 1] = z
        if (norm <= TOL):
            for i in range(n - 1):
                for j in range(m - 1):
                    print('%.1f    %.10f    %.10f' %(x[i + 1], y[j + 1], w[i + 1, j + 1]))
            return (x, y, w)
        it += 1
    print("Maximum number of iterations exceeded")
    return


def main():
    a = 0.0
    b = 6.0
    c = 0.0
    d = 5.0
    n = int((b - a) / 0.5)
    m = int((d - c) * 3)

    (x, y, w) = PEFD(a, b, c, d, n, m)

    h = (b - a) / n
    k = (d - c) / m
    X = np.arange(a + h, b, h, dtype = float)
    Y = np.arange(c + h, d, k, dtype = float)
    W = matlib.zeros((m - 1, n - 1))
    for i in range(m - 1):
        for j in range(n - 1):
            W[i, j] = w[j + 1, i + 1]

    fig = figure()
    ax = Axes3D(fig)
    X, Y = np.meshgrid(X, Y)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('w')
    ax.plot_surface(X, Y, W, cmap = plt.cm.winter)
    show()


if __name__ == "__main__":
    main()
