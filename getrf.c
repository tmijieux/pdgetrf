#include <assert.h>
#include <mpi.h>
#include <string.h>

#include "util.h"
#include "getrf.h"
#include "incblas.h"
#include "proc.h"

// factorisation LU "scalaire (BLAS2)"
void tdp_dgetf2_nopiv(int64_t m, int64_t n, double *A, int64_t lda)
{
    int64_t N = min(m, n);
    for (int64_t k = 0; k < N; ++k) {
        cblas_dscal(m-k-1, 1.0/A[k*lda+k], A+k*lda+k+1, 1);
        cblas_dger(CblasColMajor, m-k-1, n-k-1, -1.0,
                   A+k*lda+k+1, 1, A+(k+1)*lda+k, lda,
                   A+(k+1)*lda+k+1, lda);
    }
}

// factorisation LU "bloc (BLAS3)"
void tdp_dgetrf_nopiv(int64_t N, double *A, int64_t lda, int64_t b/*lock_size*/)
{
    assert( (N % b) == 0 );
    int64_t NT = N / b;
    /* printf("block_size=%d\n", b); */
    /* printf("NT=%d\n", NT); */

    for (int64_t k = 0; k < NT; ++k) {
        tdp_dgetf2_nopiv(N-k*b, b, A+b*(k*lda+k), lda);
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

static void
compute_prev_block(int64_t NB, int64_t rank, tdp_trf_dist *dist)
{
    for (int64_t i = 0; i < NB; ++i) {
        int64_t j = 0;
        if (dist->block_owner[i] != rank) {
            j = i;
            while (j >= 0 && dist->block_owner[j] != rank)
                -- j;
            if (j >= 0)
                dist->block_idx[i] = dist->block_idx[j];
            else
                dist->block_idx[i] = -1;
        }
    }
}

void tdp_trf_dist_snake2(tdp_trf_dist *dist, int64_t N, int64_t b, tdp_proc *proc)
{
    assert( (N % b) == 0 );

    int64_t NB = N / b;
    int64_t *owner = malloc(sizeof*owner * NB);
    int64_t *block_idx = malloc(sizeof*owner * NB);
    int64_t k = 0;
    int64_t inc = 1;
    int64_t owned_count = 0;

    int64_t *idx = calloc(proc->group_size, sizeof idx[0]);
    for (int64_t i = 0; i < NB; ++i) {
        owner[i] = k;
        block_idx[i] = idx[k] ++;

        if (k == proc->rank)
            ++ owned_count;

        if (proc->group_size == 1) {
            k = proc->rank;
            inc = 0;
        } else if (k == 0)
            inc = 1;
        else if (k == proc->group_size-1)
            inc = -1;
        k = k+inc;
    }

    #ifdef DEBUG
    if (!proc->rank) {
        for (int64_t i = 0; i < NB; ++i)
            fprintf(stderr, "#owner[%ld] = %ld\n", i, owner[i]);
        printf("#owned=%ld\n", owned_count);
    }
    #endif
    free(idx);

    dist->block_owner = owner;
    dist->block_idx = block_idx;
    dist->local_block_count = owned_count;

    compute_prev_block(NB, proc->rank, dist);
}

void tdp_trf_dist_snake(
    tdp_trf_dist *dist, int64_t N, int64_t b, tdp_proc *proc)
{
    assert( (N % b) == 0 );

    int64_t NB = N / b;
    int64_t *owner = malloc(sizeof*owner * NB);
    int64_t *block_idx = malloc(sizeof*owner * NB);
    int64_t k = 0;
    int64_t owned_count = 0;

    int64_t *idx = calloc(proc->group_size, sizeof idx[0]);
    bool direct = true;

    for (int64_t i = 0; i < NB; ++i) {
        if (direct) {
            owner[i] = k;
            block_idx[i] = idx[k] ++;
            if (k == proc->rank)
                ++ owned_count;
            ++ k;
            if (k == proc->group_size) {
                direct = false;
                --k;
            }
        } else {
            owner[i] = k;
            block_idx[i] = idx[k] ++;
            if (k == proc->rank)
                ++ owned_count;
            -- k;
            if (k == -1) {
                direct = true;
                ++ k;
            }
        }
    }

    #ifdef DEBUG
    if (!proc->rank) {
        for (int64_t i = 0; i < NB; ++i)
            fprintf(stderr, "#owner[%ld] = %ld\n", i, owner[i]);
        printf("#owned=%ld\n", owned_count);
    }
    #endif
    free(idx);

    dist->block_owner = owner;
    dist->block_idx = block_idx;
    dist->local_block_count = owned_count;

    compute_prev_block(NB, proc->rank, dist);
}


static MPI_Datatype block_type;
static MPI_Datatype tmp_block_type;

static void setup_mpi_types(int64_t b, int64_t ld_block, int64_t ld_tmp)
{
    MPI_Datatype tmp;
    MPI_Aint lb, ext;

    // block
    MPI_Type_vector(b, b, ld_block, MPI_DOUBLE, &tmp);
    MPI_Type_get_extent(tmp, &lb, &ext);
    ext = b*sizeof(double);
    MPI_Type_create_resized(tmp, lb, ext, &block_type);
    MPI_Type_commit(&block_type);

    // tmp
    MPI_Type_vector(b, b, ld_tmp, MPI_DOUBLE, &tmp);
    MPI_Type_get_extent(tmp, &lb, &ext);
    ext = b*sizeof(double);
    MPI_Type_create_resized(tmp, lb, ext, &tmp_block_type);
    MPI_Type_commit(&tmp_block_type);
}

static void tdp_pdgetrf_nopiv_impl(
    int64_t N, double *A, int64_t lda, int64_t b/*block size*/,
    tdp_trf_dist *dist, tdp_proc *proc   )
{
    setup_mpi_types(b, lda, N);

    assert( (N % b) == 0 );
    const int64_t NB = N / b; /* nombre de bloc */

    double *tmp = tdp_matrix_new(N, b);

    for (int64_t k = 0; k < NB; ++k) {
        bool me = dist->block_owner[k] == proc->rank;
        int64_t K = dist->block_idx[k];
        int64_t L = dist->local_block_count - (K+1);
        #ifdef DEBUG
        printf("r=%d; k=%ld; K=%ld; L=%ld, NB=%ld;%s\n",
               proc->rank, k, K, L, NB, me ? " me!" : "");
        #endif
        double *trf_ptr = me ? A+b*(K*lda+k) : tmp+b*k;
        int64_t trf_ld = me ? lda : N;
        double *col_ptr = me ? A+b*(K*lda+k+1) : tmp+b*(k+1);
        int64_t col_ld = me ? lda : N;
        double *line_ptr = A+b*((K+1)*lda+k);
        double *block_ptr = A+b*((K+1)*lda+k+1);

        if (me)
            // DGETF(k,k)
            tdp_dgetf2_nopiv(b, b, trf_ptr, trf_ld);

        //printf("r=%d checkpoint A\n", proc->rank);
        // BCAST( A(k,k)<k> )

        MPI_Datatype btype = me ? block_type : tmp_block_type;
        MPI_Bcast(trf_ptr, 1, btype,
                  dist->block_owner[k], MPI_COMM_WORLD);

         if (me) {
            cblas_dtrsm(CblasColMajor, CblasRight,
                        CblasUpper, CblasNoTrans, CblasNonUnit,
                        N-b*(k+1), b, 1.0, trf_ptr, trf_ld, col_ptr, col_ld);
            // DTRSM SUB_COL(k)
            //printf("r=%d checkpoint B.1\n", proc->rank);
        } else {
            cblas_dtrsm(CblasColMajor, CblasLeft,
                        CblasLower, CblasNoTrans, CblasUnit,
                        b, b*L, 1.0, trf_ptr, trf_ld, line_ptr, lda);
            // DTRSM RIGHT_LINE(proc->rank)
            //printf("r=%d checkpoint B.2\n", proc->rank);
        }

        // BCAST(SUB_COL<k>);
        MPI_Bcast(col_ptr, NB-(k+1), btype,
                  dist->block_owner[k], MPI_COMM_WORLD);
        //printf("r=%d checkpoint C\n", proc->rank);

        if (me) {
            cblas_dtrsm(CblasColMajor, CblasLeft,
                        CblasLower, CblasNoTrans, CblasUnit,
                        b, b*L, 1.0, trf_ptr, trf_ld, line_ptr, lda);
            // DTRSM RIGHT_LINE(k)
            //printf("r=%d checkpoint D\n", proc->rank);
        }

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    N-b*(k+1), b*L, b, -1.0,
                    col_ptr, col_ld, line_ptr, lda, 1.0, block_ptr, lda);
        // DGEMM( SUB_COL<k> x RIGHT_LINE )
        //printf("r=%d checkpoint E\n", proc->rank);

        #ifdef DEBUG
        MPI_Barrier(MPI_COMM_WORLD);
        if (proc->rank == 0) {
            printf("\n________________\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
    }
    free(tmp);
}

static void tdp_pdgetrf_nopiv_impl2(
    int64_t N, double *A, int64_t lda, int64_t b/*block size*/,
    tdp_trf_dist *dist, tdp_proc *proc   )
{
    setup_mpi_types(b, lda, N);

    assert( (N % b) == 0 );
    const int64_t NB = N / b; /* nombre de bloc */

    double *tmp = tdp_matrix_new(N, b);

    for (int64_t k = 0; k < NB; ++k) {
        bool me = dist->block_owner[k] == proc->rank;
        int64_t K = dist->block_idx[k];
        int64_t L = dist->local_block_count - (K+1);
        #ifdef DEBUG
        printf("r=%d; k=%ld; K=%ld; L=%ld, NB=%ld;%s\n",
               proc->rank, k, K, L, NB, me ? " me!" : "");
        #endif
        double *trf_ptr = me ? A+b*(K*lda+k) : tmp+b*k;
        int64_t trf_ld = me ? lda : N;
        double *col_ptr = me ? A+b*(K*lda+k+1) : tmp+b*(k+1);
        int64_t col_ld = me ? lda : N;
        double *line_ptr = A+b*((K+1)*lda+k);
        double *block_ptr = A+b*((K+1)*lda+k+1);

        if (me)
            // DGETF(k,k)
            tdp_dgetf2_nopiv(N-k*b, b, trf_ptr, trf_ld);

        //printf("r=%d checkpoint A\n", proc->rank);
        // BCAST( A(k,k)<k> )

        MPI_Datatype btype = me ? block_type : tmp_block_type;
        MPI_Bcast(trf_ptr, NB-k, btype,
                  dist->block_owner[k], MPI_COMM_WORLD);

        if (L > 0) {
            cblas_dtrsm(CblasColMajor, CblasLeft,
                        CblasLower, CblasNoTrans, CblasUnit,
                        b, b*L, 1.0, trf_ptr, trf_ld, line_ptr, lda);
            // DTRSM RIGHT_LINE(k)
            //printf("r=%d checkpoint D\n", proc->rank);

            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        N-b*(k+1), b*L, b, -1.0,
                        col_ptr, col_ld, line_ptr, lda, 1.0, block_ptr, lda);
            // DGEMM( SUB_COL<k> x RIGHT_LINE )
            //printf("r=%d checkpoint E\n", proc->rank);
        }

        #ifdef DEBUG
        MPI_Barrier(MPI_COMM_WORLD);
        if (proc->rank == 0) {
            printf("\n________________\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
    }
    free(tmp);
}

void tdp_pdgetrf_nopiv(int64_t N, double *A,
                       int64_t lda, int64_t b/*block size*/,
                       tdp_trf_dist *dist, tdp_proc *proc  )
{
    assert( A != NULL );

    if (proc->group_size == 1)
        tdp_dgetrf_nopiv(N, A, lda, b); // not distributed
    else
        tdp_pdgetrf_nopiv_impl(N, A, lda, b, dist, proc);
}



/* --------------------------------------------- */
// WITH PIVOT

// factorisation LU "scalaire (BLAS2)"
void tdp_dgetf2(int64_t m, int64_t n, double *A, int64_t lda,
                int64_t ipiv[m], int64_t *info)
{
    int64_t N = min(m, n);
    *info = 0L;

    for (int64_t k = 0; k < N; ++k) {
        int64_t p = k + cblas_idamax(m-k, A+k*lda+k, 1);
        ipiv[k] = p;
        if (A[p*lda+p] != 0) {
            if (p != k)
                cblas_dswap(n, A+k, lda, A+k+p, lda);

            cblas_dscal(m-k-1, 1.0/A[k*lda+k], A+k*lda+k+1, 1);
        } else if (*info == 0) {
            *info = k;
        }
        cblas_dger(CblasColMajor, m-k-1, n-k-1, -1.0,
                   A+k*lda+k+1, 1, A+(k+1)*lda+k, lda,
                   A+(k+1)*lda+k+1, lda);
    }

}


void tdp_pdgetrf( int64_t N, double *A, int64_t lda, int64_t b,
		  int64_t ipiv[N], tdp_trf_dist *dist,
		  tdp_proc *proc,int64_t *info)
{
    setup_mpi_types(b, lda, N);

    assert( (N % b) == 0 );
    const int64_t NB = N / b; /* nombre de bloc */

    double *tmp = tdp_matrix_new(N, b);

    for (int64_t k = 0; k < NB; ++k) {
      bool me = dist->block_owner[k] == proc->rank;
      int64_t K = dist->block_idx[k];
      int64_t L = dist->local_block_count - (K+1);
#ifdef DEBUG
      printf("r=%d; k=%ld; K=%ld; L=%ld, NB=%ld;%s\n",
	proc->rank, k, K, L, NB, me ? " me!" : "");
#endif
      double *trf_ptr = me ? A+b*(K*lda+k) : tmp+b*k;
      int64_t trf_ld = me ? lda : N;
      double *col_ptr = me ? A+b*(K*lda+k+1) : tmp+b*(k+1);
      int64_t col_ld = me ? lda : N;
      double *line_ptr = A+b*((K+1)*lda+k);
      double *block_ptr = A+b*((K+1)*lda+k+1);

      if (me)
	// DGETF(k,k)
	tdp_dgetf2(N-k*b, b, trf_ptr, trf_ld, ipiv+k*b, info);

      //printf("r=%d checkpoint A\n", proc->rank);
      // BCAST( A(k,k)<k> )

      MPI_Datatype btype = me ? block_type : tmp_block_type;
      MPI_Bcast(trf_ptr, NB-k, btype,
		dist->block_owner[k], MPI_COMM_WORLD);
      MPI_Bcast(ipiv+k*b, b, MPI_LONG,
		dist->block_owner[k], MPI_COMM_WORLD);

      if (L > 0) {
	for(int i=0; i<b; ++i){
	  ipiv[k*b+i]+=k*b;
	  if (ipiv[k*b+i] != k*b+i)
	    cblas_dswap(N, A+ipiv[k*b+i], lda, A+k*b+i, lda);
	}
	cblas_dtrsm(CblasColMajor, CblasLeft,
		    CblasLower, CblasNoTrans, CblasUnit,
		    b, b*L, 1.0, trf_ptr, trf_ld, line_ptr, lda);
	// DTRSM RIGHT_LINE(k)
	//printf("r=%d checkpoint D\n", proc->rank);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		    N-b*(k+1), b*L, b, -1.0,
		    col_ptr, col_ld, line_ptr, lda, 1.0, block_ptr, lda);
	// DGEMM( SUB_COL<k> x RIGHT_LINE )
	//printf("r=%d checkpoint E\n", proc->rank);
      }

#ifdef DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
      if (proc->rank == 0) {
      printf("\n________________\n");
    }
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    free(tmp);

}
