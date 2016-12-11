#ifndef TDP_UTIL_H
#define TDP_UTIL_H

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <immintrin.h> //AVX
#include <assert.h>
#include <omp.h>

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef max
#define max(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifdef PASTE
#undef PASTE
#endif
#define PASTE(X, Y) PASTE2(X, Y)
#define PASTE2(X, Y) X ## Y

#ifdef DEQUAL
#undef DEQUAL
#endif
#define DEQUAL(X_, Y_, EPS_) (fabs((X_) - (Y_)) < (EPS_))

#ifdef SQUARE
#undef SQUARE
#endif
#define SQUARE(X_) ((uint64_t)(X_)*(uint64_t)(X_))

#ifdef CUBE
#undef CUBE
#endif
#define CUBE(X_) ((uint64_t)(X_)*(uint64_t)(X_)*(uint64_t)(X_))

#ifdef __FMA__
#define MM256_FMADD_PD(a, b, c) _mm256_fmadd_pd(a, b, c)
#else
#define MM256_FMADD_PD(a, b, c) _mm256_add_pd(_mm256_mul_pd(a, b), c)
#endif

double *tdp_matrix_new(int64_t m/*rows*/, int64_t n/*columns*/);
double *tdp_avx256_aligned_matrix_new(int64_t m/*rows*/, int64_t n/*columns*/);
void tdp_matrix_zero(int64_t m/*rows*/, int64_t n/*columns*/, double *mat);
void tdp_matrix_one(int64_t m/*row*/, int64_t n/*column*/,
                    double value, double *mat, int64_t lda/*leading dimension*/);
void tdp_matrix_fill(int64_t m/*row*/, int64_t n/*column*/,
                     double value, double *mat, int64_t lda/*leading dimension*/);

void tdp_matrix_print(int64_t m/*row*/, int64_t n/*column*/,
                      double *mat, int64_t lda/*leading dimension*/,
                      FILE *outstream);
void tdp_matrix_rand(int64_t m/*rows*/, int64_t n/*columns*/,
                     double *mat, double min, double max);

double *tdp_vector_new(int64_t m);
void tdp_vector_rand(int64_t m, double min, double max, double *v);
void tdp_vector_one(int64_t m, double value, double *v);
void tdp_vector_zero(int64_t m, double *v);
void tdp_vector_print(int64_t m, double *v, FILE *out);

void tdp_print_cache_size(void);
uint64_t tdp_get_cache_size(int64_t id);
double *tdp_cache_garbage(void);

#endif // TDP_UTIL_H
