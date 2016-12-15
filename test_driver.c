#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "util.h"
#include "getrf.h"
#include "gesv.h"

static void init_test_matrices(int N, double *A, double *X)
{
   for (int i = 0; i < N-1; ++i) {
        A[i] = 1.0;
        A[(i+1)*N+i] = 1.0 + (double) i;
        X[i] = 2.0 + (double) i;
    }
    A[N-1] = 1.0;
    X[N-1] = 1.0;
}

void test_dgesv_nopiv_1(void) // test solve Ax=b "bloc"
{
    // initialization:
    int const N = 10;
    double *A = tdp_matrix_new(N, N);
    double *X = tdp_vector_new(N);
    init_test_matrices(N, A, X);

    // initial "check"
    printf("A:\n");
    tdp_matrix_print(N, N, A, N, stdout);

    printf("b:\n");
    tdp_vector_print(N, X, stdout);

    // solve:
    tdp_dgesv_nopiv(CblasColMajor, CblasNoTrans, N, A, N, X, 1, N/2);

    printf("\"LU\":\n");
    tdp_matrix_print(N, N, A, N, stdout);
    printf("solution:\n");
    tdp_vector_print(N, X, stdout);
}

void test_dgesv2_nopiv_1(void) // test solve Ax=b "scalaire"
{
    // initialization:
    int const N = 10;
    double *A = tdp_matrix_new(N, N);
    double *X = tdp_vector_new(N);
    init_test_matrices(N, A, X);

    // initial "check"
    printf("A:\n");
    tdp_matrix_print(N, N, A, N, stdout);
    printf("b:\n");
    tdp_vector_print(N, X, stdout);

    // solve:
    tdp_dgesv2_nopiv(CblasColMajor, CblasNoTrans, N, A, N, X, 1);

    printf("\"LU\":\n");
    tdp_matrix_print(N, N, A, N, stdout);
    printf("solution:\n");
    tdp_vector_print(N, X, stdout);
}

static void dist_snake_init_test(
    tdp_proc *proc, tdp_trf_dist *dist, int64_t N, int64_t b,
    double *A,  double *X, int lda)
{
    assert( (N % b) == 0 );
    int64_t NB = N / b;

    if (dist->block_owner[0] == proc->rank) {
        for (int i = 0; i < N; ++i)
            A[i] = 1.0;
        for (int i = 1; i < b; ++i)
            A[i*lda+i-1] = i;
    }
    for (int k = 1; k < NB; ++k) {
        if (dist->block_owner[k] == proc->rank) {
            int64_t K = dist->block_idx[k];
            double *B = A+(K*b)*lda+k*b-1;

            for (int i = 0; i < b; i++)
                B[i*lda+i] = (double) k*b+i;
        }
    }

    for (int i = 0; i < N-1; ++i)
        X[i] = 2.0 + i;
    X[N-1] = 1.0;
}

static void print_distributed_matrix(
    tdp_trf_dist *dist, tdp_proc *proc, int64_t N, int64_t b, double *A)
{
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < proc->group_size; ++i) {
        if (proc->rank == i) {
            printf("rank=%d:\n", proc->rank);
            tdp_matrix_print(N, b*dist->local_block_count, A, N, stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void test_pdgesv_nopiv(tdp_proc *proc)
{
    int N = 10, b = 2;

    tdp_trf_dist dist;
    tdp_trf_dist_snake(&dist, N, b, proc);

    double *A = tdp_matrix_new(N, b*dist.local_block_count);
    double *X = tdp_vector_new(N);

    dist_snake_init_test(proc, &dist, N, b, A, X, N);
    if (!proc->rank) {
        printf("B:\n");
        tdp_vector_print(N, X, stdout);
    }

    print_distributed_matrix(&dist, proc, N, b, A);
    tdp_pdgesv_nopiv(N, A, N, X, 1, b, &dist, proc);

    if (!proc->rank)
        puts("GETRF:");
    print_distributed_matrix(&dist, proc, N, b, A);

    if (!proc->rank) {
        printf("solution:\n");
        tdp_vector_print(N, X, stdout);
    }
}

#define TEST(type)                              \
    do{                                         \
        printf("Testing \"%s\".\n", #type);	\
        test_##type();				\
    }while(0)

static void tdp_proc_init(tdp_proc *proc)
{
    MPI_Comm_size(MPI_COMM_WORLD, &proc->group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc->rank);
}

void test_dgetf2(void) // test solve Ax=b "scalaire"
{
    // initialization:
    int const N = 10;
    double *A = tdp_matrix_new(N, N);
    double *X = tdp_vector_new(N);
    init_test_matrices(N, A, X);
    tdp_matrix_print(N, N, A, N, stdout);

    int64_t info;
    int64_t *ipiv = malloc(sizeof*ipiv * N);
    tdp_dgetf2(N, N, A, N, ipiv, &info);
    tdp_matrix_print(N, N, A, N, stdout);

    puts("ipiv: ");
    for (int64_t i = 0; i < N; ++i) {
        printf("%ld ", ipiv[i]);
    }
    puts("");
}

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    srand(time(NULL)+(long)&argc);

    /* TEST(dgesv_1); */
    /* TEST(dgesv2_1); */

    tdp_proc proc;
    tdp_proc_init(&proc);

    //test_pdgesv_nopiv(&proc);
    test_dgetf2();

    MPI_Finalize();
    return EXIT_SUCCESS;
}


