#ifndef PERF_H
#define PERF_H

#include <sys/time.h>
#include <stdint.h>

typedef struct timeval perf_t;

void perf(perf_t * p);
void perf_diff(const perf_t * begin, perf_t * end);
void perf_printh(const perf_t * p);
void perf_printmicro(const perf_t * p);
double perf_mflops(const perf_t * p, const uint64_t nb_op);

#define PERF_MICRO(p) ((uint64_t) (&p)->tv_usec + ((&p)->tv_sec * 1000000UL))
#define PERF_MFLOPS(p, nb_op) ( (double)nb_op / PERF_MICRO(p) )
#define PERF_MFLOPS2(p, nb_op) ( (double)nb_op / p )

#endif // PERF_H
