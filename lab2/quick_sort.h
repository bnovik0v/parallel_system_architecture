#pragma once

#include <omp.h>

namespace myalg {
    int __partition__(std::vector<int> &vec, int low, int high) {
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

    void __qsort__(std::vector<int> &vec, int low, int high) {
        if (low < high) {
            int p = __partition__(vec, low, high);
            __qsort__(vec, low, p - 1);
            __qsort__(vec, p + 1, high);
        }
    }

    void qsort(std::vector<int> &vec) {
        __qsort__(vec, 0, vec.size());
    }

    void __pqsort__(std::vector<int> &vec, int low, int high) {
        if (low < high) {
            int p = __partition__(vec, low, high);
            #pragma omp task shared(vec) firstprivate(low, high)
                __pqsort__(vec, low, p - 1);
            #pragma omp task shared(vec) firstprivate(low, high)
                __pqsort__(vec, p + 1, high);
        }
    }


    void pqsort(std::vector<int> &vec) {
        #pragma omp parallel shared(vec) default(none)
        {
            #pragma omp single nowait
                __pqsort__(vec, 0, vec.size());
        }
    }
}

