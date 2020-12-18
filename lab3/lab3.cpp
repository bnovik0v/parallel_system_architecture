#include <iostream>
#include <fstream>
#include <cstdlib>

#include "mpi_init.h"
#include "linear_quatuions_jakobi.h"
#include "timer.h"

using namespace std;

void fillMatrix(double *&matrix, const int &N, const int &M);

void fillWithDiagonalDominance(double *&matrix, const int &N, const int &M);

int main(int argc, char **argv) {
    int size, rank;
    mpi_setup(argc, argv, size, rank);
    cout << "rank: " << rank << " size: " << size << endl;

    std::ofstream out;
    if (rank == 0) {
        out.open("solv_linear_eq_p" + to_string(size) + ".csv", ios::trunc);
        out << "d,t" << endl;
    }

    Timer timer;

    double *A = nullptr, *x = nullptr, *b = nullptr;

    for (int N = size; N <= 1000; ++N) {                // ~1000 dimensions to check

        if (N % size != 0)                              // dimension needs to be % to size
            continue;

        for (int i = 0; i < 5; ++i) {                   // repeat it 5 times

            if (rank == 0) {                            // only process 0 has values
                A = new double[N * N];
                x = new double[N];
                b = new double[N];

                fillMatrix(A, N, N);
                fillMatrix(b, N, 1);

                fillWithDiagonalDominance(A, N, N);  // matrix needs to be diagonal dominant
            }

            MPI_Barrier(MPI_COMM_WORLD);                // barrier for processes to start in one time

            if (rank == 0)
                timer.reset();

            if (size == 1)
                yakobi(A, x, b, N);
            else
                yakobi_parallel(A, x, b, N, rank, size);

            if (rank == 0) {                            // process 0 outputs the results
                out << N << "," << timer.elapsed() << endl;

                delete[] A;
                delete[] x;
                delete[] b;
            }
        }
    }

    if (rank == 0)
        out.close();

    mpi_endup();

    return 0;
}

void fillWithDiagonalDominance(double *&matrix, const int &N, const int &M) {
    double sum = 0;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            if (i != j)
                sum += matrix[i * N + j];       // sum of all elements except diagonal ones

    for (int i = 0; i < N; ++i)
        matrix[i * N + i] += sum;               // to make diagonal dominant
}

void fillMatrix(double *&matrix, const int &N, const int &M) {
    for (int i = 0; i < N * M; ++i)
        matrix[i] = (rand() % 10) + 1;          // fill matrix with random values
}





