#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "quick_sort.h"


std::vector<int> init_vector(const int &size) {
    std::vector<int> vec(size);
    for (int &i : vec) {
        i = (int) random() % size;
    }
    return vec;
}

void print_vector(const std::vector<int> &vec) {
    for (int i : vec) {
        std::cout << i << " ";
    }
}

void sort_demo(void (*sort_func)(std::vector<int> &), const int &size) {
    std::vector<int> vec = init_vector(size);

    std::cout << "vector : ";
    print_vector(vec);

    sort_func(vec);

    std::cout << "\tsorted : ";
    print_vector(vec);

    std::cout << std::endl;
}

clock_t estimate_time(void (*sort_func)(std::vector<int> &), const int &size, const int &count = 100) {
    double time = 0;
    for (int i = 0; i < count; i++) {
        std::vector<int> vec = init_vector(size);
        clock_t start = clock();
        sort_func(vec);
        time += (double) (clock() - start);
    }
    return time / count;
}

int main() {
    int num_threads = 4;

    omp_set_num_threads(num_threads);

    sort_demo(&myalg::qsort, 10);
    sort_demo(&myalg::pqsort, 10);

    std::cout << std::endl;

    std::cout << "seq : " << estimate_time(myalg::qsort, 25000, 50) << " ticks" << std::endl;
    std::cout << "par : " << estimate_time(myalg::pqsort, 25000, 50) << " ticks" << std::endl;
}