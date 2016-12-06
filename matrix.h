#ifndef LA_MATRIX_H
#define LA_MATRIX_H

typedef struct LA_matrix {
    int M, N;
    int ld;
    double *buf;
} LA_matrix;

typedef struct LA_vector {
    int N;
    int inc;
    double *buf;
} LA_vector;

void LA_matrix_init(LA_matrix *m, int M, int N);
void LA_matrix_mult(
    const LA_matrix *m1, const LA_matrix *m2, LA_matrix *out);
void LA_matrix_submatrix(const LA_matrix *m, LA_matrix *sub,
                         int i, int j, int Msub, int Nsub);
void LA_matrix_solve(
    LA_matrix *A, LA_vector *BX, const inst block_size);

void LA_vector_init(LA_vector *v, int N);
void LA_vector_sub(LA_vector *v, int N, int inc, double *buf);

#endif // LA_MATRIX_H
