#include <stdio.h>
#include <omp.h>

int main() {
    int id;
#pragma omp parallel private(id)
    {
        id = omp_get_thread_num();
        printf("Hello World! by thread %d\n", id);
    }
    return 0;
}