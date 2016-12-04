#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include "util.h"
#include "getrf.h"

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
