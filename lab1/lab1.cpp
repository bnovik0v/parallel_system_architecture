#include <iostream>
#include <cstdlib>
#include "timer.h"

#include "matrix_multiplication.h"

using namespace std;

// MPI
void setup(int &argc, char **argv, int &size, int &rank);
void endup();

// Matrix
void generate2DMatrix(double *&M, const int &n);
void randFill2DMatrix(double *&M, const int &n);
void free2DMatrix(double *&M);
void show2DMatrix(double *&M, const int &n);

int main(int argc, char **argv) {
    int size = 0, rank = -1;
    double *A  = nullptr, *B  = nullptr, *C = nullptr;

    //MPI_INIT
    setup(argc, argv, size, rank);

    Timer timer;
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 5; ++j) {
            int N = 10*i+10;

            if (N % size != 0) continue;

            generate2DMatrix(A, N);
            generate2DMatrix(B, N);
            generate2DMatrix(C, N);

            randFill2DMatrix(A, N);
            randFill2DMatrix(B, N);

            timer.reset();

            //linearMatrixMultiplication(A, B, C, N);
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
    endup();

    return 0;
}

void setup(int &argc, char **argv, int &size, int &rank)
{
    int res;

    res = MPI_Init(&argc, &argv);
    if (res != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Init failed\n");
        exit(0);
    }

    res = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (res != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Comm_size failed\n");
        exit(0);
    }

    res = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (res != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Comm_rank failed\n");
        exit(0);
    }
}

void endup() {
    int res = MPI_Finalize();
    if (res != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Finalize failed\n");
        exit(0);
    }
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