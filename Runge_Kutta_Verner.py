import math


def Runge_Kutta_Verner(a, b, alpha, TOL, hmax, hmin):
    t = a
    w = alpha
    h = hmax
    flag = 1
    while (flag):
        x = w / t
        k1 = h * x * (1 - x)
        x = (w + k1/6) / (t + h/6)
        k2 = h * x * (1 - x)
        x = (w + 4*k1/75 + 16*k2/75) / (t + 4*h/15)
        k3 =h * x * (1 - x)
        x = (w + 5*k1/6 - 8*k2/3 + 5*k3/2) / (t + 2*h/3)
        k4 = h * x * (1 - x)
        x = (w - 165*k1/64 + 55*k2/6 - 425*k3/64 + 85*k4/96) / (t + 5*h/6)
        k5 =h * x * (1 - x)
        x = (w + 12*k1/5 - 8*k2 + 4015*k3/612 - 11*k4/36 + 88*k5/255) / (t + h)
        k6 = h * x * (1 - x)
        x = (w - 8263*k1/15000 + 124*k2/75 - 643*k3/680 - 81*k4/250 + 2484*k5/10625) / (t + h/15)
        k7 = h * x * (1 - x)
        x = (w + 3501*k1/1720 - 300*k2/43 + 297275*k3/52632 - 319*k4/2322 + 24068*k5/84065 + 3850*k7/26703) / (t + h)
        k8 = h * x * (1 - x)
        wi = w + 13*k1/160 + 2375*k3/5984 + 5*k4/16 + 12*k5/85 + 3*k6/44
        ei = math.fabs(k1/160 + 125*k3/17952 - k4/144 + 12*k5/1955 + 3*k6/44 - 125*k7/11592 - 43*k8/616)
        R = ei / math.fabs(wi)
        if (R <= TOL):
            t = t + h
            zi = w + 3*k1/40 + 875*k3/2244 + 23*k4/72 + 264*k5/1955 + 125*k7/11592 + 43*k8/616
            w = zi
            y = t/(1 + math.log(t))
            e = math.fabs((zi - y) / y)
            print(('t = ' + str(t)).ljust(30) + ' ' + ('w = ' + str(w)).ljust(30) + ' ' + ('h = ' + str(h)).ljust(30))
            # print(('error1 = ' + ("%.15f" %ei)).ljust(30) + ('error2 = ' + ("%.15f" %e)).ljust(30) + '\n')
        deta = 0.8 * math.pow(TOL/R, 1/6)
        if (deta < 0.2):
            h *= 0.2
        elif (deta > 2):
            h *= 2
        else:
            h *= deta
        if (h > hmax):
            h = hmax
        if (t >= b):
            flag = 0
        elif (t + h > b):
            h = b - t
        elif (h < hmin):
            flag = 0
            print("minimum h exceeded")
    return

if __name__ == "__main__":
    Runge_Kutta_Verner(1, 4, 1, 0.000001, 0.5, 0.05)
