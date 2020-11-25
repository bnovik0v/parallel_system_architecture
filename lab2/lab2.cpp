#include <iostream>
#include <fstream>

#include "quick_sort.h"
#include "timer.h"

void fillArray(double *&array, const int &N)
{
    srand(0);

    for (int i = 0; i < N; ++i) {
        array[i] = rand() % N;
    }
}

void printArray(double *&array, const int &N)
{
    for (int i = 0; i < N; ++i) {
        std::cout << array[i] << " ";
    }
    std::cout << std::endl;
}


int main() {

    const int N = 10;




    Timer timer;

    std::ofstream out;
    out.open("qsort.csv", std::ios::trunc);
    out << "d,p1,p2,p4,p8" << std::endl;



    for (int N = 100; N <= 50000; N+=100)
    {
        auto *array = new double[N];

        out << N;

        for (int numThreads = 1; numThreads <= 8; numThreads*=2) {
            double fullTime = 0;
            for (int i = 0; i < 5; ++i) {
                fillArray(array, N);

                timer.reset();
                if (numThreads == 1) {
                    qsort(array, N);
                }
                else {
                    pqsort(array, N, numThreads);
                }
                fullTime += timer.elapsed();
            }
            out << ","  << fullTime / 5;
        }

        delete [] array;

        out << std::endl;
    }

    return 0;
}