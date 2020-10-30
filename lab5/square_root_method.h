#pragma once

#include <cmath>

#include <omp.h>


void squareRootMethod(double *&A, double *&b, double *&x, const int &N)
{
    double *G = nullptr;
    G = new double [N * N]; //int sizeG = N * (N + 1) / 2;

    G[0] = sqrt(A[0]);

    for (int j = 1; j < N; ++j) {
        G[j * N] = A[j * N] / G[0];
    }

    for (int i = 1; i < N; ++i) {
        double sumG = 0;

        for (int k = 0; k < i; ++k) {
            sumG += G[i * N + k] * G[i * N + k];
        }

        G[i * N + i] = sqrt(A[i * N + i] - sumG);
    }

    for (int j = 2; j < N; ++j) {
        for (int i = 1; i < j; ++i) {
            double sumG = 0;

            for (int k = 0; k < i; ++k) {
                sumG += G[i * N + k] * G[j * N + k];
            }

            G[j * N + i] = (A[j * N + i] - sumG) / G[i * N + i];
        }
    }

    double *y = nullptr;
    y = new double [N];

    y[0] = b[0] / G[0];

    for (int i = 1; i < N; ++i) {
        double sumGY = 0;

        for (int k = 0; k < i; ++k) {
            sumGY += G[i * N + k] * y[k];
        }

        y[i] = (b[i] - sumGY) / G[i * N + i];
    }

    x[N-1] = y[N-1] / G[(N - 1) * N + (N - 1)];

    for (int i = N - 2; i >= 0; --i) {
        double sumGY = 0;

        for (int k = i + 1; k < N; ++k) {
            sumGY += G[k * N + i] * x[k];
        }

        x[i] = (y[i] - sumGY) / G[i * N + i];
    }

    delete [] G;
    delete [] y;
}

#include <iostream>

void squareRootMethodParallel(double *&A, double *&b, double *&x, const int &N, const int &threadsN)
{
    double *G = nullptr;
    G = new double[N * N]; //int sizeG = N * (N + 1) / 2;

    double *y = nullptr;
    y = new double[N];

    #pragma omp parallel default(shared) num_threads(threadsN)
    {
        #pragma omp single
        G[0] = sqrt(A[0]);

        #pragma omp for schedule(static)
        for (int j = 1; j < N; ++j) {
            G[j * N] = A[j * N] / G[0];
        }

        // можно снаружи
        for (int i = 1; i < N; ++i) {
            double sumG = 0;

            #pragma omp parallel for reduction(+:sumG) schedule(static)
            for (int k = 0; k < i; ++k) {
                sumG += G[i * N + k] * G[i * N + k];
            }

            #pragma omp single
            G[i * N + i] = sqrt(A[i * N + i] - sumG);
        }

        // только внутри
        for (int j = 2; j < N; ++j) {
            for (int i = 1; i < j; ++i) {
                double sumG = 0;

                #pragma omp parallel for reduction(+:sumG) schedule(static)
                for (int k = 0; k < i; ++k) {
                    sumG += G[i * N + k] * G[j * N + k];
                }

                #pragma omp single
                G[j * N + i] = (A[j * N + i] - sumG) / G[i * N + i];
            }
        }

        //!!!!!!!!!!!!!!!!!!!!!!!!!
        #pragma omp single
        y[0] = b[0] / G[0];

        // только внутри
        for (int i = 1; i < N; ++i) {
            double sumGY = 0;

            #pragma omp parallel for reduction(+:sumGY) schedule(static)
            for (int k = 0; k < i; ++k) {
                sumGY += G[i * N + k] * y[k];
            }

            #pragma omp single
            y[i] = (b[i] - sumGY) / G[i * N + i];
        }

        //!!!!!!!!!!!!!!!!!!!!!!!!

        #pragma omp single
        x[N - 1] = y[N - 1] / G[(N - 1) * N + (N - 1)];

        for (int i = N - 2; i >= 0; --i) {
            double sumGY = 0;

            #pragma omp parallel for reduction(+:sumGY) schedule(static)
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