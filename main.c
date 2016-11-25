#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mkl.h>
#include <string.h>

#include "util.h"
//#include "cblas.h"

void dgetf2_nopiv2(int m, int n, double *A, int lda)
{
    int N = min(m, n);
    for (int k = 0; k < N; ++k) {
        cblas_dscal(N-k-1, 1.0/A[k*lda+k], A+k*lda+k+1, 1);
        cblas_dger(CblasColMajor, m-k-1, n-k-1, -1.0,
                   A+k*lda+k+1, 1, A+(k+1)*lda+k, lda,
                   A+(k+1)*lda+k+1, lda);
    }
}

void dgetrf_nopiv(int N, double *A, int lda, int b/*lock_size*/)
{
    assert( (N % b) == 0 );
    int NT = N / b;
    
    for (int k = 0; k < NT; ++k) {

        double *LU = A+b*(k*lda+k);
        dgetf2_nopiv2(N-k*b, b, LU, lda);
        cblas_dtrsm(CblasColMajor, CblasRight,
                    CblasLower, CblasNoTrans, CblasNonUnit,
                    b, N-k*b, -1.0, LU, lda,
                    A+b*k*lda, lda);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    N-k*b, N-k*b, b, -1.0,
                    A+b*(k*lda+k+1), lda,
                    A+b*((k+1)*lda+k), lda,
                    1.0, A+b*((k+1)*lda+k+1), lda);
    }
}

int main(int argc, char *argv[])
{
    srand(time(NULL)+(long)&argc);
    int N = 10;
    double *A = tdp_matrix_new(N, N);
    tdp_matrix_rand(N, N, A, -10.0, 10.0);

    tdp_matrix_print(N, N, A, N, stdout);
    dgetf2_nopiv2(N, N, A, N);
    tdp_matrix_print(N, N, A, N, stdout);
    
    return EXIT_SUCCESS;
}
