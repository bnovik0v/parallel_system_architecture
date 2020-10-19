#include <iostream>
#include <cstdlib>

#include "mpi_init.h"
#include "linear_quatuions_jakobi.h"
#include "timer.h"

using namespace std;

void fillMatrix(double *& matrix, const int & N);

int main(int argc, char **argv)
{
    int size, rank;
    mpi_setup(argc, argv, size, rank);

    const int N = 12;

    cout << "rank: " << rank << " size: " << size << endl;


    auto *A = new double [N * N],
         *x = new double [N],
         *b = new double [N];

    fillMatrix(A, N * N);
    fillMatrix(b, N);

    for (int i = 0; i < N; ++i) {
        A[i * N + i] = A[i * N + i] * 1000; // матрица должна преобладать на диагонали
        x[i] = b[i] / A[i * N + i]; // подготовка
    }


    yakobi(A, x, b, N);
    //yakobi_parallel(A, x, b, N, rank, size);

    if (rank == 0)
        for (int i = 0; i < N; ++i) {
            cout << x[i] << " ";
        }


    delete [] A;
    delete [] x;
    delete [] b;

    mpi_endup();

    return 0;
}

void fillMatrix(double *& matrix, const int &N)
{
    for (int i = 0; i < N; ++i)
        matrix[i] = (rand() % 10) + 1;
}



