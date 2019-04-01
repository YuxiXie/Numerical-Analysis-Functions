#include<iostream>

using namespace std;

double Newton_backward(double *y, double * F) {
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= i; j++) {
            F[i * 4 + j] = (F[i * 4 + j - 1] - F[i * 4 + j - 5]) / (y[i] - y[i - j]);
        }
    }
    double y0 = 0;
    double a = 1;
    double p = F[3 * 4 + 0];
    for (int i = 1; i <= 3; i++) {
        a *= y0 - y[4 - i];
        p += a * F[3 * 4 + i];
    }
    return p;
}

int main() {
    double ex[4] = {0};
    double y[4] = {0};
    double F[4 * 4] = {0};

    for (int i = 0; i < 4; i++) {
        cin >> F[i * 4] >> ex[i];
        y[i] = F[i * 4] - ex[i];
    }

    double p = Newton_backward(y, F);
    cout << "p = " <<  p << endl;
}
