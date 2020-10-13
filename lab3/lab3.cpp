
#include <cstdlib>

void fillMatrix(double *& matrix, const int &N, const int &M);

void yakobi(double *& A, double *& x, double *& b, const int &N);

int main()
{

    double *A = new double [16],
           *x = new double [4],
           *b = new double [4];

    fillMatrix(A, 4, 4);
    fillMatrix(b, 4, 1);





    delete [] A;
    delete [] x;
    delete [] b;

    return 0;
}

void fillMatrix(double *& matrix, const int &N, const int &M)
{
    for (int i = 0; i < N * M; ++i)
            matrix[i] = rand();
}



void yakobi(double *& A, double *& x, double *& b, const int &N)
{
    // e = 0.001
    //
    // x_next = nullptr
    // x_cur = x


    //do
    // x_next = A
    //
    //while norm(x_n - x_n-1) >= e







}