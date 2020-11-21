#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <omp.h>

#include "shell_sort.h"


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
    srand(0);
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
        srand(0);
        std::vector<int> vec = init_vector(size);

        clock_t start = clock();

        sort_func(vec);

        time += (double) (clock() - start);
    }
    return time / count;
}

int main() {
    sort_demo(&myalg::shell_sort, 30);
    sort_demo(&myalg::p_shell_sort, 30);

    std::cout << std::endl;

    std::cout << "seq : " << estimate_time(&myalg::shell_sort, 500, 100) << " ticks" << std::endl;
    omp_set_num_threads(2);
    std::cout << "par 2: " << estimate_time(&myalg::p_shell_sort, 500, 100) << " ticks" << std::endl;
    omp_set_num_threads(4);
    std::cout << "par 4: " << estimate_time(&myalg::p_shell_sort, 500, 100) << " ticks" << std::endl;
}