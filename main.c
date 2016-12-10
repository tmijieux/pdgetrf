#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <unistd.h>

#include "util.h"
#include "getrf.h"
#include "gesv.h"
#include "proc.h"
#include "perf/perf.h"

int main(int argc, char *argv[])
{
    srand(time(NULL)+(long)&argc);
    MPI_Init(NULL, NULL);

    tdp_proc proc;
    tdp_trf_dist dist;

    MPI_Comm_size(MPI_COMM_WORLD, &proc.group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc.rank);

    int N = 20000, b = 200;
    char hostname[1024];
    gethostname(hostname, 1024);
    hostname[1023] = 0;
    printf("#rank=%d; group_size=%d; '%s'\n", proc.rank, proc.group_size, hostname);

    printf("dist start\n");
    tdp_trf_dist_snake(&dist, N, b, &proc);
    printf("dist end\nalloc start\n");
    
    double *A = tdp_matrix_new(N, b*dist.local_block_count);
    printf("r=%d BC=%ld, alloc end\nA=%p size=%luMB\n", proc.rank, dist.local_block_count, A,
           (N*b*dist.local_block_count*8UL) / (1024*1024));
    printf("rand START\n");
    tdp_matrix_rand(N, b*dist.local_block_count, A, -1.0, 1.0);
    printf("rand END\n");


    //double t = MPI_Wtime();
    perf_t p1, p2;
    perf(&p1);
    tdp_pdgetrf_nopiv(N, A, N, b, &dist, &proc);
    //t = MPI_Wtime() - t;
    perf(&p2);
    perf_diff(&p1, &p2);
    
    uint64_t micro = PERF_MICRO(p2);
    uint64_t max = 0;
    MPI_Reduce(&micro, &max, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!proc.rank) {
        printf("time=%lu µs || %lu s\n", max, max / 1000000UL);
        printf("MFLOPS=%g\n",PERF_MFLOPS2( max, (2.0/3.0)*CUBE(N)));
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
