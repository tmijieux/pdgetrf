#include <assert.h>

#include "util.h"
#include "getrf.h"


// factorisation LU "scalaire"
void
dgetf2_nopiv(int m, int n, double *A, int lda)
{
    int N = min(m, n);
    for (int k = 0; k < N; ++k) {
        cblas_dscal(m-k-1, 1.0/A[k*lda+k], A+k*lda+k+1, 1);
        cblas_dger(CblasColMajor, m-k-1, n-k-1, -1.0,
                   A+k*lda+k+1, 1, A+(k+1)*lda+k, lda,
                   A+(k+1)*lda+k+1, lda);
    }
}

// factorisation LU "bloc"
void
dgetrf_nopiv(int N, double *A, int lda, int b/*lock_size*/)
{
    assert( (N % b) == 0 );
    int NT = N / b;
    printf("block_size=%d\n", b);
    printf("NT=%d\n", NT);

    for (int k = 0; k < NT; ++k) {
        dgetf2_nopiv(N-k*b, b, A+b*(k*lda+k), lda);
        cblas_dtrsm(CblasColMajor, CblasLeft,
                    CblasLower, CblasNoTrans, CblasUnit,
                    b, N-b*(k+1), 1.0, A+b*(k*lda+k), lda,
                    A+b*((k+1)*lda+k), lda);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    N-b*(k+1), N-b*(k+1), b, -1.0,
                    A+b*(k*lda+k+1), lda,
                    A+b*((k+1)*lda+k), lda,
                    1.0, A+b*((k+1)*lda+k+1), lda);
    }
}

// General Matrix Solve Vector "scalaire" (solve Ax=b avec b vecteur)
void
dgesv2(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
       const int N, double *A, const int lda, double *X, const int incX)
{
    assert( order == CblasColMajor );
    assert( TransA == CblasNoTrans );

    dgetf2_nopiv(N, N, A, lda);

    // L y = B ("descente")
    cblas_dtrsv(CblasColMajor, CblasLower,
                CblasNoTrans, CblasUnit,
                N, A, lda, X, incX);

    // U x = y ("remontée")
    cblas_dtrsv(CblasColMajor, CblasUpper,
                CblasNoTrans, CblasNonUnit,
                N, A, lda, X, incX);
}

// general matrix solve vector "bloc" (solve Ax=b avec b vecteur)
void
dgesv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
      const int N, double *A, const int lda, double *X, const int incX)
{
    assert( order == CblasColMajor );
    assert( TransA == CblasNoTrans );

    dgetrf_nopiv(N, A, lda, 5);

    // L y = B ("descente")
    cblas_dtrsv(CblasColMajor, CblasLower,
                CblasNoTrans, CblasUnit,
                N, A, lda, X, incX);

    // U x = y ("remontée")
    cblas_dtrsv(CblasColMajor, CblasUpper,
                CblasNoTrans, CblasNonUnit,
                N, A, lda, X, incX);
}
