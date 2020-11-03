#include <cstdlib>
#include <iostream>

#include "square_root_method.h"

#include "timer.h"

using namespace std;


void fillSymmetricMatrixWithRandValues(double *&M, const int &N);
void fillVectWithRandValues(double *&V, const int &N);
void addDiagonalDominanceToMatrix(double *&M, const int &N);

void showMatrix(double *&M, const int &N);
void showVect(double *&V, const int &N);


int main() {

    double *A = nullptr, *b = nullptr, *x = nullptr;

    const int N = 2000;

    A = new double [N * N];
    b = new double [N];
    x = new double [N];

    fillSymmetricMatrixWithRandValues(A, N); // maxtrix should be symmetric
    addDiagonalDominanceToMatrix(A, N); // matrix should be positive determined
    fillVectWithRandValues(b, N);

    //showMatrix(A, N);
    //showVect(b, N);

    Timer time;

    time.reset();
    squareRootMethod(A, b, x, N);
    cout << time.elapsed() << endl;
    //showVect(x, N);


    time.reset();
    squareRootMethodParallel(A, b, x, N, 2);
    cout << time.elapsed() << endl;
    //showVect(x, N);


    delete [] A;
    delete [] b;
    delete [] x;

    return 0;
}

void showVect(double *&V, const int &N)
{
    for (int i = 0; i < N; ++i) {
        cout << V[i] << " ";
    }
    cout << endl;
    cout << endl;
}

void showMatrix(double *&M, const int &N)
{
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << M[i * N + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void fillSymmetricMatrixWithRandValues(double *&M, const int &N) {
    const int VALUE_LIMIT = 10;
    const int VALUE_SHIFT = 1;

    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            double val = (rand() % VALUE_LIMIT) + VALUE_SHIFT;
            M[i * N + j] = M[j * N + i] = val;
        }
    }
}

void fillVectWithRandValues(double *&V, const int &N) {
    const int VALUE_LIMIT = 10;
    const int VALUE_SHIFT = 1;

    for (int i = 0; i < N; ++i) {
        V[i] = (rand() % VALUE_LIMIT) + VALUE_SHIFT;
    }
}

void addDiagonalDominanceToMatrix(double *&M, const int &N)
{
    double sum = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                sum += M[i * N + j];
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        M[i * N + i] += sum;
    }
}