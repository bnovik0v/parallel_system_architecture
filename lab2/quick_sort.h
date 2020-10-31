#pragma once

#include <omp.h>

namespace myalg {
    int partition(std::vector<int> &vec, int low, int high) {
        int i = low;
        int j = high;
        int pivot = low; //(i + j) / 2;

        while (i < j) {
            while (vec[i] <= vec[pivot]) i++;
            while (vec[j] > vec[pivot]) j--;
            if (i < j) {
                std::swap(vec[i], vec[j]);
            }
        }
        std::swap(vec[pivot], vec[j]);

        return j;
    }

    void qsort(std::vector<int> &vec, int low, int high) {
        if (low < high) {
            int p = partition(vec, low, high);
            qsort(vec, low, p - 1);
            qsort(vec, p + 1, high);
        }
    }

    void qsort(std::vector<int> &vec) {
        qsort(vec, 0, vec.size());
    }

    void pqsort(std::vector<int> &vec, int low, int high) {
        if (low < high) {
            int p = partition(vec, low, high);
#pragma omp task shared(vec) firstprivate(low, high)
            pqsort(vec, low, p - 1);
#pragma omp task shared(vec) firstprivate(low, high)
            pqsort(vec, p + 1, high);
        }
    }


    void pqsort(std::vector<int> &vec) {
#pragma omp parallel shared(vec)
        {
#pragma omp single nowait
            pqsort(vec, 0, vec.size());
        }
    }
}

