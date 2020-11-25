#pragma once

namespace qs {

    int _partition(double *&array, int start, int end) {
        int temp;
        int marker = start;

        for (int i = start; i <= end; ++i) {
            if (array[i] < array[end]) {
                temp = array[marker];
                array[marker] = array[i];
                array[i] = temp;
                ++marker;
            }
        }

        temp = array[marker];
        array[marker] = array[end];
        array[end] = temp;

        return marker;
    }

    void _qsort(double *&array, int start, int end) {
        if (start >= end) {
            return;
        }
        int pivot = _partition(array, start, end);
        _qsort(array, start, pivot - 1);
        _qsort(array, pivot + 1, end);
    }

    void qsort(double *&array, const int &N) {
        _qsort(array, 0, N - 1);
    }

    void _pqsort(double *&array, int start, int end) {
        if (start >= end) {
            return;
        }
        int pivot = _partition(array, start, end);

#pragma omp task
        _pqsort(array, start, pivot - 1);
#pragma omp task
        _pqsort(array, pivot + 1, end);
    }

    void pqsort(double *&array, const int &N, const int &numThreads) {
#pragma omp parallel default(none) shared(array, N) num_threads(numThreads)
        {
#pragma omp single nowait
            {
                _pqsort(array, 0, N - 1);
            }
#pragma taskwait
        }
    }

}