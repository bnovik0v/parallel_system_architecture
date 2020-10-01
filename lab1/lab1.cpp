#include <iostream>
#include <mpi/mpi.h>

using namespace std;

int main(int argc, char **argv) {
    int size, rank;
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

    fprintf(stdout, "Hello World from rank %d of %d~\n", rank, size);

    res = MPI_Finalize();
    if (res != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Finalize failed\n");
        exit(0);
    }
}