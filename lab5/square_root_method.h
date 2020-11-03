#pragma once

#include <cmath>

#include <omp.h>


void squareRootMethod(double *&A, double *&b, double *&x, const int &N) {
    double *G = new double[N * N]; //int sizeG = N * (N + 1) / 2;
    double *y = new double[N];

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

    double *G = new double[N * N]; //int sizeG = N * (N + 1) / 2;
    double *y = new double[N];

#pragma omp parallel default(none) shared(A, b, x, N, G, y, threadsN)
    {
#pragma omp single
        {
            G[0] = sqrt(A[0]);
            y[0] = b[0] / G[0];
        }

#pragma omp for schedule(static, N/threadsN)
        for (int j = 1; j < N; ++j) {
            G[j * N] = A[j * N] / G[0];
        }

        for (int i = 1; i < N; ++i) {
            for (int j = 1; j <= i; ++j) {
                double sumG = 0;

#pragma for reduction(+ : sumG) schedule(static, N/threadsN)
                for (int k = 0; k < i; ++k) {
                    sumG += G[i * N + k] * G[j * N + k];
                }

#pragma omp single
                {
                    if (i == j) {
                        G[i * N + i] = sqrt(A[i * N + i] - sumG);
                    } else {
                        G[j * N + i] = (A[j * N + i] - sumG) / G[i * N + i];
                    }
                }
            }

            double sumGY = 0;

#pragma for reduction(+ : sumGY) schedule(static, N/threadsN)
            for (int k = 0; k < i; ++k) {
                sumGY += G[i * N + k] * y[k];
            }

#pragma omp single nowait
            y[i] = (b[i] - sumGY) / G[i * N + i];
        }

#pragma omp single
        {
            x[N - 1] = y[N - 1] / G[(N - 1) * N + (N - 1)];
        }

        for (int i = N - 2; i >= 0; --i) {
            double sumGY = 0;

#pragma for reduction(+ : sumGY) schedule(static, N/threadsN)
            for (int k = i + 1; k < N; ++k) {
                sumGY += G[k * N + i] * x[k];
            }

#pragma omp single
            x[i] = (y[i] - sumGY) / G[i * N + i];
        }
    }

    delete[] G;
    delete[] y;
}