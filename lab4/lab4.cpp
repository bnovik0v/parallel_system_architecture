//
// Created by master on 16.12.2020.
//
#include <iostream>

#include "gauss_method.h"
#include "timer.h"


void fill_array(double *&arr, const int &N)
{
    srand(0);
    for (int i = 0; i < N; ++i) {
        arr[i] = (rand() % 10) + 1;
    }
}

void show_matrix(double *&arr, const int &N, const int &M)
{
    for (int i = 0; i < N * M; ++i) {
        cout << arr[i] << " ";

        if ((i + 1) % N == 0)
        {
            cout << endl;
        }
    }
    cout << endl;
}

int main()
{
    int N = 50;


    double *A = new double[N*N];
    double *B = new double[N];
    double *X = new double[N];

    fill_array(A, N * N);
    fill_array(B, N);

    show_matrix(A, N, N);
    show_matrix(B, N, 1);

    p_gauss(A, B, X, N, 4);

    show_matrix(X, N, 1);

    delete [] A;
    delete [] B;
    delete [] X;

    return 0;
}