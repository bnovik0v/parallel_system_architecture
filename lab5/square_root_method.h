#pragma once

#include <iostream>
using namespace std;

#include <cmath>


void squareRootMethod(double *&A, double *&b, double *&x, const int &N)
{
    double *G = nullptr;
    G = new double [N * N]; //int sizeG = N * (N + 1) / 2;

    G[0] = sqrt(A[0]);

    for (int j = 0; j < N; ++j) {
        for (int i = j + 1; i < N; ++i) {
            G[j * N + i] = 0;
        }
    }

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