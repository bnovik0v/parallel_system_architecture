#include <cmath>
#include <omp.h>

using namespace std;


void gauss(double *&a, double *&b, double *&x, const int &N)
{
    const double eps = 0.00001;

    for (int k = 0; k < N; ++k) {
        double max = fabs(a[k]);
        int index = k;

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
}


void p_gauss(double *&a, double *&b, double *&x, const int &N, const int &threadsN)
{
    omp_set_num_threads(threadsN);

    const double eps = 0.00001;

#pragma omp parallel default(none) shared(a, b, x, N, eps)
    {
#pragma omp single nowait
        {
            for (int k = 0; k < N; ++k) {
                double max = fabs(a[k]);
                int index = k;

#pragma omp taskloop shared(a, N, index, max)
                for (int i = k + 1; i < N; ++i) {
                    if (fabs(a[i * N + k]) > max) {
#pragma omp critical
                        {
                            index = i;
                            max = fabs(a[i * N + k]);
                        }
                    }
                }

#pragma omp taskloop shared(a, N, index)
                for (int i = 0; i < N; ++i) {
                    swap(a[k * N + i], a[index * N + i]);
                }

#pragma omp task shared(b, index)
                swap(b[k], b[index]);


                for (int i = k; i < N; ++i) {
                    double temp = a[i * N + k];

                    if (abs(temp) < eps) {
                        continue;
                    }

#pragma omp taskloop shared(a, i, temp)
                    for (int j = 0; j < N; ++j) {
                        a[i * N + j] /= temp;
                    }

#pragma omp task shared(b, i, temp)
                    b[i] /= temp;

                    if (i == k) {
                        continue;
                    }

#pragma omp taskloop shared(a, i, k, temp)
                    for (int j = 0; j < N; ++j) {
                        a[i * N + j] -= a[k * N + j];
                    }

#pragma omp task shared(b, i, k, temp)
                    b[i] -= b[k];
                }
            }

            for (int k = N - 1; k >= 0; --k) {
#pragma omp task shared(x, k) depend(out:x[k])
                x[k] = b[k];


                for (int i = 0; i < k; ++i) {
#pragma omp task shared(a, b, x, k) depend(in:x[k])
                    b[i] -= a[i * N + k] * x[k];
                }
            }
        }
#pragma omp taskwait
    }
}