#pragma once

#include <cstring>
#include <cmath>

const double EPS = 0.001;

double calcRelativeX(double *& a, const int & rowN, double *& x, double & b, const int & N)
{
    double rel_x = 0;

    for (int i = 0; i < N; ++i)
        if (i != rowN)
            rel_x += a[i] * x[i];

    return (b - rel_x) / a[rowN];
}

double calcAccuracy(double *& x1, double *& x0, const int & N)
{
    double norm = 0;

    for (int i = 0; i < N; ++i)
        norm += (x1[i] - x0[i]) * (x1[i] - x0[i]);

    return sqrt(norm);
}

void yakobi(double *& A, double *& x, double *& b, const int & N)
{
    auto * curX = new double [N];
    auto * prevX = new double [N];

    memcpy(curX, x, sizeof(double) * N);

    do {
        memcpy(prevX, curX, sizeof(double) * N);

        for (int i = 0; i < N; ++i) {
            double *a = A + i * N;
            curX[i] = calcRelativeX(a, i, prevX, b[i], N);
        }
    } while(calcAccuracy(curX, prevX, N) > EPS);

    memcpy(x, curX, sizeof(double) * N);

    delete [] curX;
    delete [] prevX;
}

void yakobi_parallel(double *& A, double *& x, double *& b, const int & N, const int & size, const int & rank)
{
    int taskAmountForProc = N / size;

    int partOfA = N * N / size;
    auto * bufA = new double [partOfA];
    MPI_Scatter(A, partOfA, MPI_DOUBLE, bufA, partOfA, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int partOfB = N / size;
    auto * bufB = new double [partOfB];
    MPI_Scatter(b, partOfB, MPI_DOUBLE, bufB, partOfB, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int partOfX = N / size;
    auto * bufX = new double [partOfB];

    auto * curX = new double [N];
    memcpy(curX, x, sizeof(double) * N);

    auto *prevX = new double[N];

    do {
        memcpy(prevX, curX, sizeof(double) * N);

        for (int i = 0; i < taskAmountForProc; ++i) {
            double *a = A + i * N;
            int rowN = i + rank * taskAmountForProc;
            bufX[i] = calcRelativeX(a, rowN, prevX, b[i], N);
        }

        MPI_Allgather(bufX, partOfX, MPI_DOUBLE, curX, partOfX, MPI_DOUBLE, MPI_COMM_WORLD);
    } while(calcAccuracy(curX, prevX, N) > EPS);

    memcpy(x, curX, sizeof(double) * N);

    delete [] bufA;
    delete [] bufB;
    delete [] bufX;
    delete [] curX;
    delete [] prevX;
}