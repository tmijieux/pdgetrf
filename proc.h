#ifndef TDP_PROC_H
#define TDP_PROC_H

#include <stdint.h>

typedef struct tdp_proc_ {
    int rank;
    int group_size;

} tdp_proc;

typedef struct tdp_trf_dist_ {
    int64_t *block_owner; // block_owner[x] = rank of process owning block x

    int64_t *block_idx;
    // block_idx[x] = index of block x, when block_owner[k] == rank
    // and block_idx[x] = index of first owned block before block x, otherwise

    int64_t local_block_count;
} tdp_trf_dist;

#endif // TDP_PROC_H
