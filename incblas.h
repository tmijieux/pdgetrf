#ifndef TDP_INCBLAS_H
#define TDP_INCBLAS_H

#ifdef TDP_USE_MKL
#include <mkl.h>
typedef CBLAS_LAYOUT CBLAS_ORDER;
#else
#include <cblas.h>
#endif

#endif // TDP_INCBLAS_H
