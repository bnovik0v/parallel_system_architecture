
#include <vector>
#include <cmath>

#include <omp.h>

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

        for (int i = k + 1; i < N; ++i) {
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

            if (i == k) {
                continue;
            }

            for (int j = 0; j < N; ++j) {
                a[i * N + j] -= a[k * N + j];
            }

            b[i] -= b[k];
        }
    }

    for (int k = N - 1; k >= 0; --k) {
        x[k] = b[k];

        for (int i = 0; i < k; ++i) {
            b[i] -= a[i * N + k] * x[k];
        }
    }

    return x;
}

vector<double> p_gauss(vector<double> &A, vector<double> &B, const int &threadsN)
{
    omp_set_num_threads(threadsN);

    int N = B.size();
    const double eps = 0.00001;

    vector<double> a = A;
    vector<double> b = B;
    vector<double> x(N);

#pragma omp parallel
    {
#pragma omp single nowait
        {


#pragma omp task shared(a, b, N, eps) depend(inout: x)
            for (int k = 0; k < N; ++k) {


                double max = fabs(a[k]);
                double index = k;

                for (int i = k + 1; i < N; ++i) {
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

                    if (i == k) {
                        continue;
                    }

                    for (int j = 0; j < N; ++j) {
                        a[i * N + j] -= a[k * N + j];
                    }

                    b[i] -= b[k];
                }
            }


#pragma omp task shared(a, b, N, x) depend(in: x)
            for (int k = N - 1; k >= 0; --k) {
#pragma omp task shared(a, b, N, x, k) depend(in: x[k]) depend(inout: x[k])
                {
                    x[k] = b[k];

                    for (int i = 0; i < k; ++i) {
                        b[i] -= a[i * N + k] * x[k];
                    }
                }
            }
        }
#pragma taskwait
    }


    return x;
}


