#include <iostream>
#include <cstdlib>

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
    //MPI_INIT
    int size, rank;
    setup(argc, argv, size, rank);

    const int N = 4;
    double *A = nullptr;
    double *B = nullptr;
    double *C = nullptr;

    generate2DMatrix(A, N);
    generate2DMatrix(B, N);
    generate2DMatrix(C, N);

    randFill2DMatrix(A, N);
    randFill2DMatrix(B, N);

    double *bufA = new double[4];
    MPI_Scatter(A, 4, MPI_DOUBLE, bufA, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    cout << "PR:" << rank << endl;
    for (int i = 0; i < 4; i++)
        cout << bufA[i] << " ";
    cout << endl << endl;


    //linearMatrixMultiplication(A, B, C, N);
    parallelMatrixMultiplication(A, B, C, N);


    /*show2DMatrix(A, N);
    show2DMatrix(B, N);
    show2DMatrix(C, N);
*/
    free2DMatrix(A);
    free2DMatrix(B);
    free2DMatrix(C);

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