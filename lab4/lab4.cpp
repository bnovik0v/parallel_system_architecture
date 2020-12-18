
#include <iostream>
#include <fstream>

#include "gauss_method.h"
#include "timer.h"


void fill_array(double *&arr, const int &N)
{
    srand(42);
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

void demo(const int &N)
{
    auto *A = new double[N * N];
    auto *B = new double[N];
    auto *X = new double[N];

    fill_array(A, N * N);
    fill_array(B, N);

    show_matrix(A, N, N);
    show_matrix(B, N, 1);

    gauss(A, B, X, N);

    show_matrix(X, N, 1);

    ///////////////////////////////////////

    fill_array(A, N * N);
    fill_array(B, N);

    p_gauss(A, B, X, N, 4);

    show_matrix(X, N, 1);

    delete[] A;
    delete[] B;
    delete[] X;
}

int main()
{
    //demo(60);


    Timer timer;

    std::ofstream out;
    out.open("gauss100.csv", std::ios::trunc);
    out << "d,p1,p2,p3,p4" << std::endl;

    int A = 10, B = 100, STEP = 10, MAX_NUM_THREADS = 4, NUM_ITER = 30;


    for (int N = A; N <= B; N+=STEP) {
        out << N;
        cout << "estimation of " <<  N << "x" << N << " matrix" << endl;

        for (int num_threads = 1; num_threads <= MAX_NUM_THREADS; ++num_threads) {
            double total_time = 0;



            for (int iter = 0; iter < NUM_ITER; ++iter) {
                auto *A = new double[N * N];
                auto *B = new double[N];
                auto *X = new double[N];


                fill_array(A, N * N);
                fill_array(B, N);

                timer.reset();

                if (num_threads > 1)
                    p_gauss(A, B, X, N, num_threads);
                else
                    gauss(A, B, X, N);

                total_time += timer.elapsed();

                delete[] A;
                delete[] B;
                delete[] X;
            }

            out << "," << total_time / NUM_ITER;
        }
        out << endl;
    }

    out.close();

    return 0;

}