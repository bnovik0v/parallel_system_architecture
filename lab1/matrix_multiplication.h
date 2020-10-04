#pragma once

#include <mpi/mpi.h>

void linearMatrixMultiplication(double *&A,  double *&B, double *&C, const int &n)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void parallelMatrixMultiplication(double **A,  double **B, double **C, const int &n)
{
    double *bufA = new double[4];
    double *bufB = new double[4];
    MPI_Scatter(A, n, MPI_DOUBLE, bufA, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    delete [] bufA;
    delete [] bufB;
}