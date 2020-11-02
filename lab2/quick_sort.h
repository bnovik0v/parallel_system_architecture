#pragma once

#include <algorithm>

#include <omp.h>


namespace myalg {
    int partition(std::vector<int> &vec, int low, int high) {
        int i = low;
        int j = high;
        int pivot = vec[(i + j) / 2]; //(i + j) / 2;

        while (i <= j) {
            while (vec[i] < pivot) i++;
            while (vec[j] > pivot) j--;

            if (i >= j) {
                break;
            }

            std::swap(vec[i++], vec[j--]);
        }

        return j;
    }

    void qsort(std::vector<int> &vec, int low, int high) {
        if (low < high) {
            int p = partition(vec, low, high);
            qsort(vec, low, p);
            qsort(vec, p + 1, high);
        }
    }

    void qsort(std::vector<int> &vec) {
        qsort(vec, 0, vec.size() - 1);
    }

    void pqsort(std::vector<int> &vec, int low, int high) {
        if (low < high) {
            int p = partition(vec, low, high);
#pragma omp task default(none) shared(vec) firstprivate(low, high, p)
            pqsort(vec, low, p);
#pragma omp task default(none)  shared(vec) firstprivate(low, high, p)
            pqsort(vec, p + 1, high);
        }
    }


    void pqsort(std::vector<int> &vec) {
#pragma omp parallel default(none) shared(vec)
        {
#pragma omp single nowait
            pqsort(vec, 0, vec.size() - 1);
        }
    }
}

