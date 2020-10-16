
#include <cstdlib>
#include <iostream>

#include "linear_quatuions_jakobi.h"

void fillMatrix(double *& matrix, const int & N);


int main()
{
    const int N = 4;

    double *A = new double [N * N],
           *x = new double [N],
           *b = new double [N];

    fillMatrix(A, N * N);
    fillMatrix(b, N);

    for (int i = 0; i < N; ++i) {
        A[i * N + i] = 10 + (rand() % 90);
        x[i] = b[i] / A[i * N + i];
    }

    for (int i = 0; i < N*N; ++i) {
        std::cout << A[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << b[i] << " ";
    }

    yakobi(A, x, b, N);


    for (int i = 0; i < N; ++i) {
        std::cout << x[i] << " ";
    }


    delete [] A;
    delete [] x;
    delete [] b;

    return 0;
}

void fillMatrix(double *& matrix, const int &N)
{
    for (int i = 0; i < N; ++i)
            matrix[i] = rand() % 10;
}



