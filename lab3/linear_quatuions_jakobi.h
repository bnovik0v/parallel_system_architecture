#pragma once

#include <cstring>
#include <cmath>

#include <iostream>

const double EPS = 0.001;

void prepareX(double *&x, double *&A, double *&b, const int N) {
    for (int i = 0; i < N; ++i) {
        x[i] = b[i] / A[i * N + i];
    }
}

double calcRelativeX(double *&a, const int &rowN, double *&x, double &b, const int &N) {
    double rel_x = 0;

    for (int i = 0; i < N; ++i)
        if (i != rowN)
            rel_x += a[i] * x[i];

    return (b - rel_x) / a[rowN];
}

double calcAccuracy(double *&x1, double *&x0, const int &N) {
    double norm = 0;

    for (int i = 0; i < N; ++i)
        norm += (x1[i] - x0[i]) * (x1[i] - x0[i]);

    return sqrt(norm);
}

void yakobi(double *&A, double *&x, double *&b, const int &N) {
    auto *curX = new double[N];
    auto *prevX = new double[N];

    prepareX(curX, A, b, N);

    do {
        memcpy(prevX, curX, sizeof(double) * N);

        for (int i = 0; i < N; ++i) {
            double *a = A + i * N;
            curX[i] = calcRelativeX(a, i, prevX, b[i], N);
        }
    } while (calcAccuracy(curX, prevX, N) > EPS);

    memcpy(x, curX, sizeof(double) * N);

    delete[] curX;
    delete[] prevX;
}

void yakobi_parallel(double *&A, double *&x, double *&b, const int &N, const int &rank, const int &size) {
    auto *curX = new double[N];
    auto *prevX = new double[N];

    if (rank == 0)                                          // process 0 prepares curX values
        prepareX(curX, A, b, N);

    MPI_Bcast(curX, N, MPI_DOUBLE, 0, MPI_COMM_WORLD); // process 0 sends the prepared value to all

    int taskAmountForProc = N / size;
    int partOfA = N * N / size;
    int partOfB = N / size;
    int partOfX = N / size;

    auto *bufA = new double[partOfA];
    auto *bufB = new double[partOfB];
    auto *bufX = new double[partOfB];

    MPI_Scatter(A, partOfA, MPI_DOUBLE, bufA, partOfA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, partOfB, MPI_DOUBLE, bufB, partOfB, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    auto *acc = new double;                                // accuracy of iteration

    do {
        memcpy(prevX, curX, sizeof(double) * N);

        for (int i = 0; i < taskAmountForProc; ++i) {
            double *a = bufA + i * N;
            int rowN = i + rank * taskAmountForProc;
            bufX[i] = calcRelativeX(a, rowN, prevX, bufB[i], N);
        }

        MPI_Allgather(bufX, partOfX, MPI_DOUBLE, curX, partOfX, MPI_DOUBLE, MPI_COMM_WORLD);
        // gather curX from all in bufX
        // processes and send to them

        if (rank == 0) // only rank 0 calculates accuracy
            *acc = calcAccuracy(curX, prevX, N);

        MPI_Bcast(acc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // process 0 sends accuracy to others

        MPI_Barrier(MPI_COMM_WORLD);                        // wait until all processes to be here
    } while (*acc > EPS);

    if (rank == 0) // process 0 forms the result
        memcpy(x, curX, sizeof(double) * N);

    delete[] bufA;
    delete[] bufB;
    delete[] bufX;
    delete[] curX;
    delete[] prevX;
    delete acc;
}