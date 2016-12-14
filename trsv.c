#include <assert.h>
#include <mpi.h>

#include "trsv.h"
#include "incblas.h"
#include "proc.h"

/*
 * aucun parallelisme entre les noeuds
 * due à la répartition en colonne des matrices:
 * l'élimination ne peut être fait que
 * par le processus qui a trouvé l'inconnue.
 */

static void pdtrsv_upper_nonunit(
    tdp_trf_dist *dist, tdp_proc *proc,
    const int64_t N, int64_t b, const double *A,
    const int64_t lda, double *X, const int64_t incX)
{
    int64_t NB = N / b;
    for (int64_t k = NB-1; k >= 0; --k) {
        if (dist->block_owner[k] == proc->rank) {
            int64_t K = dist->block_idx[k];
            cblas_dtrsv(CblasColMajor, CblasUpper,
                        CblasNoTrans, CblasNonUnit,
                        b, A+b*(K*lda+k), lda, X+k*b, incX);

            //eliminer les inconnus:
            cblas_dgemv(CblasColMajor, CblasNoTrans,
                        k*b, b,
                        -1.0, A+b*(K*lda), lda,
                        X+k*b, incX,
                        1.0, X, incX  );
        }
        MPI_Bcast(X, N, MPI_DOUBLE, dist->block_owner[k], MPI_COMM_WORLD);
    }
}

void pdtrsv_lower_unit(
    tdp_trf_dist *dist, tdp_proc *proc,
    const int64_t N, int64_t b, const double *A,
    const int64_t lda, double *X, const int64_t incX)
{
    int64_t NB = N / b;
    for (int64_t k = 0; k < NB; ++k) {
        if (dist->block_owner[k] == proc->rank) {
            int64_t K = dist->block_idx[k];
            cblas_dtrsv(CblasColMajor, CblasLower,
                        CblasNoTrans, CblasUnit,
                        b, A+b*(K*lda+k), lda, X+k*b, incX);

            //eliminer les inconnus:
            cblas_dgemv(CblasColMajor, CblasNoTrans,
                        N-(k+1)*b, b,
                        -1.0, A+b*(K*lda+k+1), lda,
                        X+k*b, incX,
                        1.0, X+(k+1)*b, incX  );
        }
        MPI_Bcast(X, N, MPI_DOUBLE, dist->block_owner[k], MPI_COMM_WORLD);
    }
}

void tdp_pdtrsv(
    const CBLAS_ORDER order, const CBLAS_UPLO Uplo,
    const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
    const int64_t N, int64_t b, const double *A,
    const int64_t lda, double *X, const int64_t incX,
    tdp_trf_dist *dist, tdp_proc *proc)
{
    assert( (N % b) == 0 );
    assert( order == CblasColMajor );
    assert( TransA == CblasNoTrans );

    if (Uplo == CblasUpper) {
        assert( Diag == CblasNonUnit );

        pdtrsv_upper_nonunit(dist, proc, N, b, A, lda, X, incX);
    } else {
        assert( Uplo == CblasLower );
        assert( Diag == CblasUnit );

        pdtrsv_lower_unit(dist, proc, N, b, A, lda, X, incX);
    }
}
