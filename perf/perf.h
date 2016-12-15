#ifndef PERF_H
#define PERF_H

#include <sys/time.h>
#include <stdint.h>

typedef double perf_t;

void perf(perf_t * p);
void perf_diff(const perf_t * begin, perf_t * end);
void perf_printh(const perf_t * p);
void perf_printmicro(const perf_t * p);
double perf_mflops(const perf_t * p, const uint64_t nb_op);

/* #define PERF_MICRO(p) ((uint64_t) (&p)->tv_usec + ((&p)->tv_sec * 1000000UL)) */
/* #define PERF_MFLOPS(p, nb_op) ( (double)nb_op / PERF_MICRO(p) ) */

#define PERF_MICRO(p) ((uint64_t)(p * 1000000UL))
#define PERF_MFLOPS(p, nb_op) ( (double)nb_op / PERF_MICRO(p) )
#define PERF_MFLOPS2(p, nb_op) ( (double)nb_op / (uint64_t)p )


typedef struct timeval perf2_t;

void perf2(perf2_t * p);
void perf2_diff(const perf2_t * begin, perf2_t * end);
void perf2_printh(const perf2_t * p);
void perf2_printmicro(const perf2_t * p);
double perf2_mflops(const perf2_t * p, const uint64_t nb_op);

#define PERF2_MICRO(p) ((uint64_t) (&p)->tv_usec + ((&p)->tv_sec * 1000000UL))
#define PERF2_MFLOPS(p, nb_op) ( (double)nb_op / PERF2_MICRO(p) )
#define PERF2_MFLOPS2(p, nb_op) ( (double)nb_op / (uint64_t)p )


#endif // PERF_H
