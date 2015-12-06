#ifndef _BUILD_DE_BWT_
#define _BUILD_DE_BWT_
#include <stdint.h>
#include "rest_index.h"

typedef struct {
    uint64_t C[5];    // count of A/C/G/T/$
    uint64_t bwt_len;

    uint32_t *bwt;    // bwt char

    uint32_t unipath_n; // unipath's count
    uint32_t unipath_pos_n; // unipath pos' count (one unipath may have multi-pos)
    uint64_t *u_pos;  // offset of unipath
} bwt_t;


int build_de_bwt(hash_idx h_idx, char *prefix);

#endif
