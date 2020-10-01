#include <iostream>
#include <cstdlib>
#include <mpi/mpi.h>

#include "matrix_multiplication.h"

using namespace std;

// MPI
void setup(int &argc, char **argv, int &size, int &rank);

void endup();

//
double ** generate2DMatrix(const int &n);
void randFill2DMatrix(double **M, const int &n);
void free2DMatrix(double **M, const int &n);
void show2DMatrix(double **M, const int &n);

int main(int argc, char **argv) {
    //MPI_INIT
    int size, rank;
    setup(argc, argv, size, rank);


    const int N = 4;
    double **A = nullptr;
    double **B = nullptr;
    double **C = nullptr;

    A = generate2DMatrix(N);
    B = generate2DMatrix(N);
    C = generate2DMatrix(N);

    randFill2DMatrix(A, N);
    randFill2DMatrix(B, N);

    linearMatrixMultiplication(A, B, C, N);
    parallelMatrixMultiplication(A, B, C, N);

    show2DMatrix(A, N);
    show2DMatrix(B, N);
    show2DMatrix(C, N);

    free2DMatrix(A, N);
    free2DMatrix(B, N);
    free2DMatrix(C, N);

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

double ** generate2DMatrix(const int &n)
{
    auto **M = new double*[n];
    for (int i = 0; i < n; ++i) {
        M[i] = new double[n];
        for (int j = 0; j < n; ++j)
            M[i][j] = 0;
    }

    return M;
}

void randFill2DMatrix(double **M, const int &n)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M[i][j] = rand();
}

void free2DMatrix(double **M, const int &n)
{
    for (int i = 0; i < n; ++i)
       delete []  M[i];

    delete [] M;
}

void show2DMatrix(double **M, const int &n)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            cout << M[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}