#pragma once

void linearMatrixMultiplication(double **A,  double **B, double **C, const int &n)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                C[i][j] += A[i][k] * B[k][j];
}

void parallelMatrixMultiplication(double **A,  double **B, double **C, const int &n)
{

}