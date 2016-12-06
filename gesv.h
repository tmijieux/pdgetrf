#ifndef TDP_GESV_H
#define TDP_GESV_H

#include "incblas.h"

// General Matrix Solve Vector "scalaire" (solve Ax=b avec b vecteur)
void tdp_dgesv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
               const int N, double *A, const int lda,
               double *X, const int incX, const int block_size);

// general matrix solve vector "bloc" (solve Ax=b avec b vecteur)
void tdp_dgesv2(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA,
                const int N, double *A, const int lda,
                double *X, const int incX);

#endif // TDP_GESV_H
