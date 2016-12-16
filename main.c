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

typedef struct tdp_trf_time_ {
    perf2_t dist_time;
    perf2_t alloc_time;
    perf2_t rand_time;
    perf2_t compute_time;
    perf2_t total_time;
} tdp_trf_time;

static void print_proc_info(tdp_proc *proc)
{
    char hostname[1024];
    gethostname(hostname, 1024);
    hostname[1023] = 0;
    printf("#rank=%d; group_size=%d; node hostname='%s'\n",
           proc->rank, proc->group_size, hostname);
}

static void trf_rand_matrix(tdp_proc *proc, int N, int b, tdp_trf_time *time)
{
    perf2_t p1, p2, p3, p4;
    tdp_trf_dist dist;

    perf2(&p1);
    tdp_trf_dist_snake(&dist, N, b, proc);
    perf2(&time->dist_time);

    perf2(&p2);
    double *A = tdp_matrix_new(N, b*dist.local_block_count);
    perf2(&time->alloc_time);

    printf("r=%d BC=%ld A=%p size=%luMB\n",
           proc->rank, dist.local_block_count, A,
           (N*b*dist.local_block_count*8UL) / (1024*1024));

    perf2(&p3);
    tdp_matrix_rand(N, b*dist.local_block_count, A, -1.0, 1.0);
    perf2(&time->rand_time);

    perf2(&p4);
    //int64_t info;
    //int64_t *ipiv = malloc(sizeof*ipiv * N);

    //tdp_pdgetrf(N, A, N, b,ipiv, &dist, proc, &info);

    tdp_pdgetrf_nopiv(N, A, N, b, &dist, proc);

    perf2(&time->compute_time);

    // ----------------
    time->total_time = time->compute_time;

    perf2_diff(&p1, &time->dist_time);
    perf2_diff(&p2, &time->alloc_time);
    perf2_diff(&p3, &time->rand_time);
    perf2_diff(&p4, &time->compute_time);
    perf2_diff(&p1, &time->total_time);
}

#define REDUCE_PRINT_TIME(v_, t_) do {                          \
        uint64_t tmp_ = PERF2_MICRO((t_)-> PASTE(v_, _time));    \
        MPI_Reduce(&tmp_, &(v_), 1, MPI_UNSIGNED_LONG,          \
                   MPI_MAX, 0, MPI_COMM_WORLD);                 \
        if (!rank) {                                            \
            printf(#v_"_time=%lu µs || %g s\n",                 \
                   v_, (double) v_ / 1000000UL);                \
        }                                                       \
    } while(0)

static void print_time(int rank, tdp_trf_time *time, int64_t N)
{
    uint64_t compute, rand, dist, alloc, total;
    if (!rank)
        puts("");
    REDUCE_PRINT_TIME(dist, time);
    REDUCE_PRINT_TIME(alloc, time);
    REDUCE_PRINT_TIME(rand, time);
    REDUCE_PRINT_TIME(compute, time);
    REDUCE_PRINT_TIME(total, time);

    if (!rank) {
        puts("");
        printf("N=%ld\n", N);
        printf("MFLOPS=%g\n", PERF_MFLOPS2( compute, (2.0/3.0)*CUBE(N)));
        printf("GFLOPS=%g\n", PERF_MFLOPS2( compute,
                                            (2.0/3.0)*CUBE(N)) / 1000.0);
        printf("TFLOPS=%g\n", PERF_MFLOPS2( compute,
                                            (2.0/3.0)*CUBE(N)) / 1000000.0);
        puts("");
    }
}

static void
tdp_proc_init(tdp_proc *proc)
{
    MPI_Comm_size(MPI_COMM_WORLD, &proc->group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc->rank);
}

int main(int argc, char *argv[])
{
    perf2_t i1, i2;
    srand(time(NULL)+(long)&argc);

    tdp_proc proc;
    tdp_trf_time time;

    perf2(&i1);
    MPI_Init(NULL, NULL);
    perf2(&i2);
    perf2_diff(&i1, &i2);

    tdp_proc_init(&proc);
    fprintf(stderr, "init_time[%d]=%lu µs\n\n", proc.rank, PERF2_MICRO(i2));
    print_proc_info(&proc);

    if (argc < 3) {
        fprintf(stderr, "usage: %s MatrixSize BlockSize\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* int b = 200;// 30, 40, 50, 80, 100, 125, 160, 200, 250, 400, 500, 625 */
    /* int N = 14000; */
    int64_t b = strtoul(argv[2], NULL, 10);
    int64_t N = strtoul(argv[1], NULL, 10);
    
    trf_rand_matrix(&proc, N, b, &time);

    print_time(proc.rank, &time, N);

    perf_t f1, f2;
    perf(&f1);
    MPI_Finalize();
    perf(&f2);
    perf_diff(&f1, &f2);
    fprintf(stderr, "fini_time[%d]=%lu µs\n\n", proc.rank, PERF_MICRO(f2));

    return EXIT_SUCCESS;
}
