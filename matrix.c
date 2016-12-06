#include "matrix.h"

void LA_matrix_init(LA_matrix *m, int M, int N)
{
    m->buf = tdp_matrix_new(M, N);
    m->M = M;
    m->N = N;
    m->ld = M;
}

void LA_vector_init(LA_vector *v, int N)
{
    v->N = N;
    v->inc = 1;
    v->buf = tdp_vector_new(N);
}

void LA_vector_init_full(LA_vector *v, int N, int inc, double *buf)
{
    v->N = N;
    v->inc = inc;
    v->buf = buf;
}


void LA_matrix_mult(
    const LA_matrix *m1, const LA_matrix *m2, LA_matrix *out)
{
    assert( m1->N == m2->M );

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m1->M, m2->N, m1->N,
                1.0, m1->buf, m1->ld,
                m2->buf, m2->ld, 0.0,
                out->buf, out->ld);
}

void LA_matrix_submatrix(const LA_matrix *m, LA_matrix *sub,
                         int i, int j, int Msub, int Nsub)
{
    assert( i >= 0 );
    assert( j >= 0 );
    assert( Msub >= 0 );
    assert( Nsub >= 0 );
    assert( i + Msub < m->M );
    assert( j + Nsub < m->N );

    sub->buf = m->buf + j * m->ld + i;
    sub->ld = m->ld;
    sub->M = Msub;
    sub->N = Nsub;
}

void LA_matrix_solve(
    LA_matrix *A, LA_vector *BX, const inst block_size)
{
    assert( A->N == A->M );
    tdp_dgetrf_nopiv(N, A->buf, A->ld, block_size);

    // L y = B ("descente")
    cblas_dtrsv(CblasColMajor, CblasLower,
                CblasNoTrans, CblasUnit,
                A->N, A->buf, A->ld, BX->buf, BX->inc);

    // U x = y ("remontée")
    cblas_dtrsv(CblasColMajor, CblasUpper,
                CblasNoTrans, CblasNonUnit,
                A->N, A->buf, A->ld, BX->buf, BX->inc);
}
