#ifndef GETRF_H
#define GETRF_H


void dgetf2_nopiv(int m, int n, double *A, int lda);
void dgetrf_nopiv(int N, double *A, int lda, int b/*block_size*/);

#endif // GETRF_H
