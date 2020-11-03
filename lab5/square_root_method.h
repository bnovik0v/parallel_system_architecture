#pragma once

#include <cmath>

#include <omp.h>


void squareRootMethod(double *&A, double *&b, double *&x, const int &N) {
    auto *G = new double[N * N]; //int sizeG = N * (N + 1) / 2;
    auto *y = new double[N];

    G[0] = sqrt(A[0]);
    y[0] = b[0] / G[0];

    for (int j = 1; j < N; ++j) {
        G[j * N] = A[j * N] / G[0];
    }

    for (int i = 1; i < N; ++i) {
        for (int j = 1; j <= i; ++j) {
            double sumG = 0;

            for (int k = 0; k < i; ++k) {
                sumG += G[i * N + k] * G[j * N + k];
            }

            if (i == j) {
                G[i * N + i] = sqrt(A[i * N + i] - sumG);
            } else {
                G[j * N + i] = (A[j * N + i] - sumG) / G[i * N + i];
            }
        }

        double sumGY = 0;

        for (int k = 0; k < i; ++k) {
            sumGY += G[i * N + k] * y[k];
        }

        y[i] = (b[i] - sumGY) / G[i * N + i];
    }

    x[N - 1] = y[N - 1] / G[(N - 1) * N + (N - 1)];

    for (int i = N - 2; i >= 0; --i) {
        double sumGY = 0;

        for (int k = i + 1; k < N; ++k) {
            sumGY += G[k * N + i] * x[k];
        }

        x[i] = (y[i] - sumGY) / G[i * N + i];
    }

    delete[] G;
    delete[] y;
}


void squareRootMethodParallel(double *&A, double *&b, double *&x, const int &N, const int &threadsN) {
    omp_set_num_threads(threadsN);

    auto *G = new double[N * N]; //int sizeG = N * (N + 1) / 2;
    auto *y = new double[N];

#pragma omp parallel
    {
#pragma omp single nowait
        {
#pragma omp task shared(A, G) depend(inout: G[0])
            {
                G[0] = sqrt(A[0]);
            }

            for (int i = 1; i < N; ++i) {
#pragma omp task shared(A, G) depend(in: G[0]) depend(inout: G[i * N])
                {
                    G[i * N] = A[i * N] / G[0];
                }
            }

            for (int i = 1; i < N; ++i) {
#pragma omp task shared(A, G) depend(in: G[i * N + i - 1]) depend(inout: G[i * N + i])
                {
                    double sumG = 0;

                    for (int k = 0; k < i; ++k) {
                        sumG += G[i * N + k] * G[i * N + k];
                    }

                    G[i * N + i] = sqrt(A[i * N + i] - sumG);
                }

                for (int j = 1; j < i; ++j) {
#pragma omp task shared(A, G) depend(in: G[j * N + i - 1], G[i * N + i - 1], G[i * N + i]) depend(inout: G[j * N + i])
                    {
                        double sumG = 0;

                        for (int k = 0; k < i; ++k) {
                            sumG += G[i * N + k] * G[j * N + k];
                        }

                        G[j * N + i] = (A[j * N + i] - sumG) / G[i * N + i];
                    }
                }
            }

#pragma omp task shared(b, G, y) depend(in: G[0]) depend(inout: y[0])
            {
                y[0] = b[0] / G[0];
            }

            for (int i = 1; i < N; ++i) {
#pragma omp task shared(b, G, y) depend(in: G[i * N + i - 1], y[i - 1]) depend(inout: y[i])
                {
                    double sumGY = 0;

                    for (int k = 0; k < i; ++k) {
                        sumGY += G[i * N + k] * y[k];
                    }

                    y[i] = (b[i] - sumGY) / G[i * N + i];
                }
            }

#pragma omp task shared(x, G, y) depend(in: y[N - 1]) depend(inout: x[N - 1])
            {
                x[N - 1] = y[N - 1] / G[(N - 1) * N + (N - 1)];
            }


            for (int i = N - 2; i >= 0; --i) {
#pragma omp task shared(x, G, y) depend(in: x[i+1]) depend(inout: x[i])
                {
                    double sumGY = 0;

                    for (int k = i + 1; k < N; ++k) {
                        sumGY += G[k * N + i] * x[k];
                    }

                    x[i] = (y[i] - sumGY) / G[i * N + i];
                }
            }
        }
#pragma taskwait
    }

    delete[] G;
    delete[] y;
}