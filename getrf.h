#ifndef GETRF_H
#define GETRF_H

#ifdef MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif


void dgetf2_nopiv(int m, int n, double *A, int lda);
void dgetrf_nopiv(int N, double *A, int lda, int b/*block_size*/);


// General Matrix Solve Vector "scalaire" (solve Ax=b avec b vecteur)
void dgesv2(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
            const int N, double *A, const int lda, double *X, const int incX);

// general matrix solve vector "bloc" (solve Ax=b avec b vecteur)
void dgesv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
           const int N, double *A, const int lda, double *X, const int incX);

#endif // GETRF_H
