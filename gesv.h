#ifndef TDP_GESV_H
#define TDP_GESV_H

#include "incblas.h"

// general matrix solve vector "bloc" (solve Ax=b avec b vecteur)
void tdp_dgesv2_nopiv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                      const int64_t N, double *A, const int64_t lda,
                      double *X, const int64_t incX);


// General Matrix Solve Vector "scalaire" (solve Ax=b avec b vecteur)
void tdp_dgesv_nopiv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                     const int64_t N, double *A, const int64_t lda,
                     double *X, const int64_t incX, const int64_t block_size);

void tdp_pdgesv_nopiv(const int64_t N, double *A, const int64_t lda,
                      double *X, int64_t incX, const int64_t b,
                      tdp_trf_dist *dist, tdp_proc *proc);

#endif // TDP_GESV_H
