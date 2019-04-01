import math
import matplotlib.pyplot as plt
import numpy as np


# 求解T1(0,1,h)和T2(0,1,h)
def Trapezoidal(n, a, b):
    h = (b - a)/n
    a0 = 1.25 * math.exp(2) - 0.75
    a1 = 2/h - 1
    a2 = 1.0
    b0 = (-n/2 + 1.25) * math.exp(2) - (n/2 + 0.75)
    b1 = -(2*n + 1)
    b2 = 2/h + 1
    for i in range(n - 1):
        x = a + h * (i + 1)
        a0 += 2 * x*math.exp(2*x) + (x*x - x/2)*math.exp(2) - (x*x + x/2)
        a1 -= 2 * x
        a2 += 2 * x*x
        b0 += 2 * math.exp(2*x) + x*math.exp(2) - x
        b2 += 2 * x
    T1 = (a0*b2 - a2*b0) / (a1*b2 - a2*b1)
    T2 = (a0*b1 - a1*b0) / (a2*b1 - a1*b2)
    return (T1, T2)


# 输入n和x，求解利用复合梯形公式得到的y(x)以及误差
def SolveFunction(n, a, b, X):
    T = Trapezoidal(n, a, b)
    T1 = T[0]
    T2 = T[1]
    k = (math.exp(2)-1)/2 - T2
    b = T1 - (math.exp(2)+1)/4
    Y0 = np.exp(2 * X)
    Y1 = np.exp(2 * X) + k * X + b
    Err1 = np.subtract(Y1, Y0)
    Err2 = np.divide(Err1, Y0)
    return (Y0, Y1, Err1, Err2)


# 命令行打印结果
def PrintResult(a, b):
    X = np.array([0, 0.25, 0.75, 1])
    for i in range(4):
        n = 8 * (i + 1)
        print('n = ' + str(n) + ':')
        R = SolveFunction(n, a, b, X)
        Y0 = R[0]
        Y1 = R[1]
        Err1 = R[2]
        Err2 = R[3]
        j = 0
        for y in Y0:
            x = X[j]
            f_real = y
            f_estimated = Y1[j]
            abs_err = Err1[j]
            rel_err = Err2[j]
            j += 1
            s = 'x = ' + repr(x)
            print(s.rjust(15), repr(f_real).rjust(25), repr(f_estimated).rjust(25), end = ' ')
            print(repr(abs_err).rjust(30), repr(rel_err).rjust(25), end = ' ')
            print('\n')
        print('\n')
    return


# 作图：函数图像以及误差曲线
def Function_view(a, b):
    x = np.linspace(0, 1, 1000)
    f = np.exp(2 * x)
    F = []
    Err1 = []
    Err2 = []
    for i in range(4):
        n = (i + 1) * 8
        r = SolveFunction(n, a, b, x)
        F += [r[1]]
        Err1 += [r[2]]
        Err2 += [r[3]]
    plt.plot(x, f, label='real')
    plt.plot(x, F[0], label='estimated: n = 8')
    plt.plot(x, F[1], label='estimated: n = 16')
    plt.plot(x, F[2], label='estimated: n = 24')
    plt.plot(x, F[3], label='estimated: n = 32')
    plt.xlabel('t')
    plt.ylabel('Function y')
    plt.legend()
    plt.show()

    plt.plot(x, Err1[0], label='estimated: n = 8')
    plt.plot(x, Err1[1], label='estimated: n = 16')
    plt.plot(x, Err1[2], label='estimated: n = 24')
    plt.plot(x, Err1[3], label='estimated: n = 32')
    plt.xlabel('t')
    plt.ylabel('Error: estimated_y(t) - y(t)')
    plt.legend()
    plt.show()

    plt.plot(x, Err2[0], label='estimated: n = 8')
    plt.plot(x, Err2[1], label='estimated: n = 16')
    plt.plot(x, Err2[2], label='estimated: n = 24')
    plt.plot(x, Err2[3], label='estimated: n = 32')
    plt.xlabel('t')
    plt.ylabel('Error: [estimated_y(t) - y(t)]/y(t)')
    plt.legend()
    plt.show()
    return

if __name__ == "__main__":
    a = 0.0
    b = 1.0
    PrintResult(a, b)
    Function_view(a, b)
