#include <assert.h>

#include "getrf.h"
#include "incblas.h"
#include "gesv.h"

// General Matrix Solve Vector "scalaire (blas2)" (solve Ax=b avec b vecteur)
void tdp_dgesv2(const CBLAS_ORDER order,
                const CBLAS_TRANSPOSE TransA,
                const int64_t N, double *A, const int64_t lda,
                double *X, const int64_t incX)
{
    assert( order == CblasColMajor );
    assert( TransA == CblasNoTrans );

    tdp_dgetf2_nopiv(N, N, A, lda);

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
void tdp_dgesv(const CBLAS_ORDER order,
               const CBLAS_TRANSPOSE TransA,
               const int64_t N, double *A, const int64_t lda,
               double *X, const int64_t incX, const int64_t block_size)
{
    assert( order == CblasColMajor );
    assert( TransA == CblasNoTrans );
    assert( (N % block_size) == 0 );

    tdp_dgetrf_nopiv(N, A, lda, block_size);

    // L y = B ("descente")
    cblas_dtrsv(CblasColMajor, CblasLower,
                CblasNoTrans, CblasUnit,
                N, A, lda, X, incX);

    // U x = y ("remontée")
    cblas_dtrsv(CblasColMajor, CblasUpper,
                CblasNoTrans, CblasNonUnit,
                N, A, lda, X, incX);
}

void tdp_pdgesv(const int64_t N, double *A, const int64_t lda,
                double *X, tdp_trf_dist *dist, tdp_proc *proc)
{

}
