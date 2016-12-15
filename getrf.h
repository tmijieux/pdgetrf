#ifndef GETRF_H
#define GETRF_H

#include <stdint.h>
#include "proc.h"

void tdp_dgetf2_nopiv(int64_t m, int64_t n, double *A, int64_t lda);
void tdp_dgetrf_nopiv(int64_t N, double *A, int64_t lda, int64_t b/*block_size*/);
void tdp_pdgetrf_nopiv(int64_t N, double *A,
                       int64_t lda, int64_t b/*block size*/,
                       tdp_trf_dist *dist, tdp_proc *proc  );
void tdp_trf_dist_snake(
    tdp_trf_dist *dist, int64_t N, int64_t b, tdp_proc *proc);

void tdp_dgetf2(int64_t m, int64_t n, double *A, int64_t lda,
                int64_t ipiv[m], int64_t *info);

#endif // GETRF_H
