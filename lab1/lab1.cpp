#include <iostream>
#include <cstdlib>

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
    double *A  = nullptr, *B  = nullptr, *C = nullptr;

    //MPI_INIT
    mpi_setup(argc, argv, size, rank);

    Timer timer;
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 5; ++j) {
            int N = 100*i+100;

            if (N % size != 0) continue;

            generate2DMatrix(A, N);
            generate2DMatrix(B, N);
            generate2DMatrix(C, N);

            randFill2DMatrix(A, N);
            randFill2DMatrix(B, N);

            timer.reset();

            //matrixMultiplication(A, B, C, N);
            parallelMatrixMultiplication(A, B, C, N, rank, size);

            if (rank == 0) cout << "d:" << N << ", t:" << timer.elapsed() << endl;

            free2DMatrix(A);
            free2DMatrix(B);
            free2DMatrix(C);
        }

    /*if (rank==0) {
        //show2DMatrix(A, N);
        //show2DMatrix(B, N);
        //show2DMatrix(C, N);
    }*/

    //MPI_FINALIZE
    mpi_endup();

    return 0;
}

void generate2DMatrix(double *&M, const int &n)
{
    M = new double[n*n];
    for (int i = 0; i < n * n; ++i) {
        M[i] = 0;
    }
}

void randFill2DMatrix(double *&M, const int &n)
{
    for (int i = 0; i < n * n; ++i)
            M[i] = rand() % 10;
}

void free2DMatrix(double *&M)
{
    delete [] M;
    M = nullptr;
}

void show2DMatrix(double *&M, const int &n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            cout << M[i * n + j] << " ";
        cout << endl;
    }
    cout << endl;
}