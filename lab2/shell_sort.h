#pragma once

namespace myalg {

    void shell_sort(std::vector<int> &vec) {
        int length = vec.size();

        for (int h = length / 2; h > 0; h = h / 2) {
            for (int i = 0; i < h; i++) {
                for (int f = h + i; f < length; f = f + h) {
                    int j = f;
                    while (j > i && vec[j - h] > vec[j]) {
                        std::swap(vec[j], vec[j-h]);
                        j = j - h;
                    }
                }
            }
        }
    }

    void p_shell_sort(std::vector<int> &vec) {
        int length = vec.size();
        for (int h = length / 2; h > 0; h = h / 2) {
            #pragma omp parallel for shared(vec, length, h) default(none)
            for (int i = 0; i < h; i++) {
                for (int f = h + i; f < length; f = f + h) {
                    int j = f;
                    while (j > i && vec[j - h] > vec[j]) {
                        std::swap(vec[j], vec[j-h]);
                        j = j - h;
                    }
                }
            }
        }
    }
}