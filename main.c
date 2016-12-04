#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#ifdef MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif

#include "util.h"

// factorisation LU "scalaire"
static void
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
static void
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
static void
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
static void
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


void test_dgesv_1(void) // test solve Ax=b "bloc"
{
    // initialization:
    int const N = 10;
    double *A = tdp_matrix_new(N, N);
    double *X = tdp_vector_new(N);
    for (int i = 0; i < N-1; ++i) {
        A[i] = 1.0;
        A[(i+1)*N+i] = 1.0 + (double) i;
        X[i] = 2.0 + (double) i;
    }
    A[N-1] = 1.0;
    X[N-1] = 1.0;

    // initial "check"
    printf("A:\n");
    tdp_matrix_print(N, N, A, N, stdout);
    printf("b:\n");
    tdp_vector_print(N, X, stdout);

    // solve:
    dgesv(CblasColMajor, CblasNoTrans, N, A, N, X, 1);

    printf("\"LU\":\n");
    tdp_matrix_print(N, N, A, N, stdout);
    printf("solution:\n");
    tdp_vector_print(N, X, stdout);
}


void test_dgesv2_1(void) // test solve Ax=b "scalaire"
{
    // initialization:
    int const N = 10;
    double *A = tdp_matrix_new(N, N);
    double *X = tdp_vector_new(N);
    for (int i = 0; i < N-1; ++i) {
        A[i] = 1.0;
        A[(i+1)*N+i] = 1.0 + (double) i;
        X[i] = 2.0 + (double) i;
    }
    A[N-1] = 1.0;
    X[N-1] = 1.0;

    // initial "check"
    printf("A:\n");
    tdp_matrix_print(N, N, A, N, stdout);
    printf("b:\n");
    tdp_vector_print(N, X, stdout);

    // solve:
    dgesv2(CblasColMajor, CblasNoTrans, N, A, N, X, 1);

    printf("\"LU\":\n");
    tdp_matrix_print(N, N, A, N, stdout);
    printf("solution:\n");
    tdp_vector_print(N, X, stdout);
}

/* #define TEST(type, symbol)                              \ */
/*     do {                                                \ */
/*         printf("Testing \"%s\" %s.\n", #type, #symbol);	\ */
/*         test_##type(symbol);				\ */
/*     }while(0) */

/* #define BENCH(type, symbol)                             \ */
/*     do {						\ */
/*         /\* printf("\nBenching %s:\n", #symbol);  *\/	\ */
/*         printf("%s\n", #symbol);			\ */
/*         bench_##type(symbol);                           \ */
/*     }while(0) */

#define TEST(type)                              \
    do {                                        \
        printf("Testing \"%s\".\n", #type);	\
        test_##type();				\
    }while(0)


int main(int argc, char *argv[])
{
    srand(time(NULL)+(long)&argc);

    TEST(dgesv_1);
    TEST(dgesv2_1);

    return EXIT_SUCCESS;
}
