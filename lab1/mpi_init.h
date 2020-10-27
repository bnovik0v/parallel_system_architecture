#pragma once

#include <mpi.h>

void mpi_setup(int &argc, char **argv, int &size, int &rank) {
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

void mpi_endup() {
    int res = MPI_Finalize();
    if (res != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Finalize failed\n");
        exit(0);
    }
}