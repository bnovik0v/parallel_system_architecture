#pragma once

#include <cstring>
#include <cmath>


double calcRelativeX(double *& a, const int & rown, double *& x, double & b, const int & N)
{
    double rel_x = 0;

    for (int i = 0; i < N; ++i)
        if (i != rown)
            rel_x += a[i] * x[i];

    return (b - rel_x) / a[rown];
}

double calcAccuracy(double *& x1, double *& x0, const int & N)
{
    double norm = 0;

    for (int i = 0; i < N; ++i)
        norm += (x1[i] - x0[i]) * (x1[i] - x0[i]);

    return sqrt(norm);
}

void yakobi(double *& A, double *& x, double *& b, const int & N)
{
    const double EPS = 0.001;

    auto * x_cur = new double [N];
    memcpy(x_cur, x, sizeof(double) * N);

    auto * x_prev = new double [N];

    double *a = nullptr;

    do {
        memcpy(x_prev, x_cur, sizeof(double) * N);

        for (int i = 0; i < N; ++i) {
            a = A + i * N;
            x_cur[i] = calcRelativeX(a, i, x_prev, b[i], N);
        }
    } while(calcAccuracy(x_cur, x_prev, N) >= EPS);

    memcpy(x, x_cur, sizeof(double) * N);

    delete [] x_cur;
    delete [] x_prev;
}

void yakobi_parallel(double *& A, double *& x, double *& b, const int & N)
{

}