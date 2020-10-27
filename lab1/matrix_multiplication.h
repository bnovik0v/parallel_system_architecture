#pragma once

#include <mpi/mpi.h>

void matrixMultiplication(double *&A, double *&B, double *&C, const int &n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                C[i * n + j] += A[i * n + k] * B[j * n + k];
}

void parallelMatrixMultiplication(double *&A, double *&B, double *&C, const int &n, const int &procRank,
                                  const int &procAmount) {
    int part = n * n / procAmount;
    int taskAmount = n / procAmount;

    auto *bufA = new double[part];
    auto *bufB = new double[part];
    auto *bufC = new double[part];

    MPI_Scatter(A, part, MPI_DOUBLE, bufA, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, part, MPI_DOUBLE, bufB, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < procAmount; ++i) {
        for (int j = 0; j < taskAmount; ++j)
            for (int k = 0; k < taskAmount; ++k) {
                double temp = 0;
                for (int l = 0; l < n; ++l)
                    temp += bufA[j * n + l] * bufB[k * n + l];
                bufC[((i + procRank) % procAmount) * taskAmount + j * n + k] = temp;
            }

        MPI_Sendrecv_replace(bufB, part, MPI_DOUBLE, ((procRank - 1) < 0 ? procAmount - 1 : procRank - 1),
                             0, (procRank + 1) % procAmount, 0, MPI_COMM_WORLD, nullptr);
    }

    MPI_Gather(bufC, part, MPI_DOUBLE, C, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    delete[] bufA;
    delete[] bufB;
    delete[] bufC;
}