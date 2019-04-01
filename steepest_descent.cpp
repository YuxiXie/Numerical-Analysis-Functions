#include <iostream>
#include <cmath>

using namespace std;

void steepest_descend(int n, double *x, double TOL, int N) {
    int k = 1;
    double g1;
    double g2;
    double g3;
    double *z = new double[n + 1];
    double z_norm;
    while (k <= N) {
        g1 = 100 * pow(pow(x[1], 2) - x[2], 2) + pow(1 - x[1], 2);
        z[1] = 400 * pow(x[1], 3) - 400 * x[1] * x[2] + 2 * x[1] - 2;
        z[2] = 200 * x[2] - 200 * pow(x[1], 2);
        z_norm = sqrt(pow(z[1], 2) + pow(z[2], 2));
        if (z_norm == 0) {
            cout << "Zero gradient" << endl;
            cout  << "( " << x[1] << ", " << x[2] << " )" << endl;
            cout << g1 << endl << k << endl;
            return ;
        }
        z[1] /= z_norm;
        z[2] /= z_norm;
        double alpha1 = 0;
        double alpha3 = 1;
        double x1 = x[1] - alpha3 * z[1];
        double x2 = x[2] - alpha3 * z[2];
        g3 = 100 * pow(pow(x1, 2) - x2, 2) + pow(1 - x1, 2);
        while (g3 >= g1) {
            alpha3 /= 2;
            x1 = x[1] - alpha3 * z[1];
            x2 = x[2] - alpha3 * z[2];
            g3 = 100 * pow(pow(x1, 2) - x2, 2) + pow(1 - x1, 2);
            if (alpha3 < TOL/2) {
                cout << "No likely improvement" << endl << "( ";
                cout  << "( " << x[1] << ", " << x[2] << " )" << endl;
                cout << g1 << endl << k << endl;
                return ;
            }
        }
        double alpha2 = alpha3 / 2;
        x1 = x[1] - alpha2 * z[1];
        x2 = x[2] - alpha2 * z[2];
        g2 = 100 * pow(pow(x1, 2) - x2, 2) + pow(1 - x1, 2);
        double h1 = (g2 - g1) / alpha2;
        double h2 = (g3 - g2) / (alpha3 - alpha2);
        double h3 = (h2 - h1 ) / alpha3;
        double alpha0 = (alpha2 - h1 / h3) / 2;
        x1 = x[1] - alpha0 * z[1];
        x2 = x[2] - alpha0 * z[2];
        double g0 = 100 * pow(pow(x1, 2) - x2, 2) + pow(1 - x1, 2);
        double alpha, g;
        if (g0 < g3) {
            alpha = alpha0;
            g = g0;
        } else {
            alpha = alpha3;
            g = g3;
        }
        x[1] = x[1] - alpha * z[1];
        x[2] = x[2] - alpha * z[2];
        if (fabs(g - g1) < TOL) {
            cout << "( " << x[1] << ", " << x[2] << " )" << endl;
            cout << g << endl << k << endl;
            return ;
        }
        k++;
    }
    delete [] z;
    cout << "Maximum iterations exceeded" << endl;
    return ;
}



int main() {
    double x[3] = {0, 0, 0};
    steepest_descend(2, x, 0.005, 100);
    return 0;
}
