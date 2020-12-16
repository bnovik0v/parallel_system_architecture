//
// Created by master on 16.12.2020.
//
#include <iostream>

#include "gauss_method.h"


void fill_array(vector<double> &arr)
{
    for (double & el : arr) {
        el = (rand() % 10) + 1;
    }
}

void show_matrix(vector<double> &arr, const int &N)
{
    for (int i = 0; i < arr.size(); ++i) {
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
    int N = 6;
    vector<double> A(N*N);
    vector<double> B(N);

    fill_array(A);
    fill_array(B);

    show_matrix(A, N);
    show_matrix(B, N);

    auto X = gauss(A, B);

    show_matrix(X, N);

    return 0;
}