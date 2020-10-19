
#include <cstdlib>
#include <iostream>

#include "mpi_init.h"
#include "linear_quatuions_jakobi.h"
#include "timer.h"

void fillMatrix(double *& matrix, const int & N);


int main(int argc, char **argv)
{
    int size, rank;
    mpi_setup(argc, argv, size, rank);




    const int N = 4;

    double *A = new double [N * N],
           *x = new double [N],
           *b = new double [N];

    fillMatrix(A, N * N);
    fillMatrix(b, N);

    for (int i = 0; i < N; ++i) {
        A[i * N + i] = A[i * N + i] * 10;
        x[i] = b[i] / A[i * N + i];
    }



    yakobi(A, x, b, N);




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



