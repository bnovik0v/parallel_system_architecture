#include <iostream>
#include <cstdlib>
#include <fstream>

#include "mpi_init.h"
#include "matrix_multiplication.h"
#include "timer.h"

using namespace std;

// Matrix
void generate2DMatrix(double *&M, const int &n);

void randFill2DMatrix(double *&M, const int &n);

void free2DMatrix(double *&M);

void show2DMatrix(double *&M, const int &n);

int main(int argc, char **argv) {
    int size = 0, rank = -1;
    mpi_setup(argc, argv, size, rank);

    cout << "rank=" << rank << " size=" << size << endl;

    std::ofstream out;
    if (rank == 0) {
        out.open("matrix_mult_p" + to_string(size) + ".csv", ios::trunc);
        out << "d,t" << endl;
    }

    double *A = nullptr,
            *B = nullptr,
            *C = nullptr;

    Timer timer;

    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 5; ++j) {
            int N = 100 * i + 100;

            if (N % size != 0) continue;

            if (rank == 0) {
                generate2DMatrix(A, N);
                generate2DMatrix(B, N);
                generate2DMatrix(C, N);

                randFill2DMatrix(A, N);
                randFill2DMatrix(B, N);
            }

            MPI_Barrier(MPI_COMM_WORLD);

            if (rank == 0)
                timer.reset();

            if (size == 1)
                matrixMultiplication(A, B, C, N);
            else
                parallelMatrixMultiplication(A, B, C, N, rank, size);

            if (rank == 0) {
                out << N << "," << timer.elapsed() << endl;

                free2DMatrix(A);
                free2DMatrix(B);
                free2DMatrix(C);
            }
        }

    mpi_endup();

    return 0;
}

void generate2DMatrix(double *&M, const int &n) {
    M = new double[n * n];
    for (int i = 0; i < n * n; ++i) {
        M[i] = 0;
    }
}

void randFill2DMatrix(double *&M, const int &n) {
    for (int i = 0; i < n * n; ++i)
        M[i] = rand() % 10;
}

void free2DMatrix(double *&M) {
    delete[] M;
    M = nullptr;
}

void show2DMatrix(double *&M, const int &n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            cout << M[i * n + j] << " ";
        cout << endl;
    }
    cout << endl;
}