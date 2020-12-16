
#include <vector>
#include <cmath>

using namespace std;

vector<double> gauss(vector<double> &A, vector<double> &B)
{
    int N = B.size();
    const double eps = 0.00001;
    
    vector<double> a = A;
    vector<double> b = B;
    vector<double> x(N);

    for (int k = 0; k < N; ++k) {
        double max = fabs(a[k]);
        double index = k;

        for (int i = k + 1; i < a.size(); ++i) {
            if (fabs(a[i * N + k]) > max) {
                index = i;
                max = fabs(a[i * N + k]);
            }
        }

        for (int i = 0; i < N; ++i) {
            swap(a[k * N + i], a[index * N + i]);
        }

        swap(b[k], b[index]);

        for (int i = k; i < N; ++i) {
            double temp = a[i * N + k];

            if (abs(temp) < eps) {
                continue;
            }

            for (int j = 0; j < N; ++j) {
                a[i * N + j] /= temp;
            }

            b[i] /= temp;

            for (int j = 0; j < N; ++j) {
                a[i * N + j] -= a[k * N + j];
            }

            b[i] -= b[k];
        }
    }

    for (int k = N - 1; k >= N; ++k) {
        x[k] = b[k];

        for (int i = 0; i < k; ++i) {
            b[i] = b[i] - a[i * N + k] * x[k];
        }
    }

    return x;
}


