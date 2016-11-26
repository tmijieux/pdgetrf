#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include "util.h"
#include "perf/perf.h"
#include "getrf.h"

int main(int argc, char *argv[])
{
    (void) argv;
    srand(time(NULL)+(long)&argc);
    int N = 10;
    double *A = tdp_matrix_new(N, N);
    tdp_matrix_rand(N, N, A, -10.0, 10.0);

    tdp_matrix_print(N, N, A, N, stdout);
    dgetf2_nopiv(N, N, A, N);
    tdp_matrix_print(N, N, A, N, stdout);
    
    return EXIT_SUCCESS;
}
