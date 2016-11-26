#include <mkl.h>
#include <assert.h>

#include "getrf.h"
#include "util.h"

void dgetf2_nopiv(int m, int n, double *A, int lda)
{
    int N = min(m, n);
    for (int k = 0; k < N; ++k) {
        cblas_dscal(N-k-1, 1.0/A[k*lda+k], A+k*lda+k+1, 1);
        cblas_dger(CblasColMajor, m-(k+1), n-(k+1), -1.0,
                   A+k*lda+k+1, 1, A+(k+1)*lda+k, lda,
                   A+(k+1)*lda+k+1, lda);
    }
}

void dgetrf_nopiv(int N, double *A, int lda, int b/*block_size*/)
{
    assert( (N % b) == 0 );
    int NT = N / b;


    for (int k = 0; k < NT; ++k) {
        dgetf2_nopiv(N-k*b, b, A+b*(k*lda+k), lda);
        cblas_dtrsm(CblasColMajor, CblasLeft,
                    CblasLower, CblasNoTrans, CblasUnit,
                    b, N-(k+1)*b,
                    1.0, A+b*(k*lda+k), lda,
                    A+b*((k+1)*lda+k), lda);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    N-(k+1)*b, N-(k+1)*b, b, -1.0,
                    A+b*(k*lda+k+1), lda,
                    A+b*((k+1)*lda+k), lda,
                    1.0, A+b*((k+1)*lda+k+1), lda);
    }
}
