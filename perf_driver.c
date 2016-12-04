#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include "util.h"
#include "perf/perf.h"
#include "getrf.h"

void test_dgetf2_nopiv(void)
{
    for (int N = 100; N <= 2500; N += 100) {
        double *A = tdp_matrix_new(N, N);
        tdp_matrix_rand(N, N, A, -1.0, 1.0);

        perf_t p1, p2;
        perf(&p1);

        dgetf2_nopiv(N, N, A, N);

        perf(&p2);
        perf_diff(&p1, &p2);

        printf("%d %lu %g\n", N, PERF_MICRO(p2), PERF_MFLOPS(p2, (2.0/3.0)*CUBE(N)));
        free(A);
    }
}


void test_1_dgetrf_nopiv(int N, int b)
{
    double *A = tdp_matrix_new(N, N);
    tdp_matrix_rand(N, N, A, -1.0, 1.0);

    perf_t p1, p2;
    perf(&p1);

    dgetrf_nopiv(N, A, N, b);

    perf(&p2);
    perf_diff(&p1, &p2);

    printf("%d %d %lu %g\n", N, b, PERF_MICRO(p2), PERF_MFLOPS(p2, (2.0/3.0)*CUBE(N)));
    free(A);
}

void test_dgetrf_nopiv(void)
{
    for (int N = 100; N <= 2000; N += 100)
        test_1_dgetrf_nopiv(N, N /10);

    test_1_dgetrf_nopiv(10000, 1000);
}

void test_dgetrf_nopiv_block_size(int N)
{
    for (int b = 25; b < 250; b+=25) {
        if ((N % b) == 0)
            test_1_dgetrf_nopiv(N, b);
    }
}

#define TEST(type)                              \
    do {                                        \
        printf("# test \"%s\"\n", #type);	\
        test_##type();				\
    }while(0)


int main(int argc, char *argv[])
{
    (void) argv;

    srand(time(NULL)+(long)&argc);

    TEST(dgetf2_nopiv);
    TEST(dgetrf_nopiv);
    test_dgetrf_nopiv_block_size(5000);

    return EXIT_SUCCESS;
}
