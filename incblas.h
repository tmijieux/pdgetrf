#ifndef TDP_INCBLAS_H
#define TDP_INCBLAS_H

#ifdef MKL
#include <mkl.h>
typedef CBLAS_LAYOUT CBLAS_ORDER;
#else
#include <cblas.h>
#endif

#endif // TDP_INCBLAS_H
