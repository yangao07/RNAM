#ifndef _BUILD_DE_BWT_
#define _BUILD_DE_BWT_
#include <stdint.h>
#include "rest_index.h"

int build_de_bwt(char *prefix, hash_idx *h_idx, de_bwt_t *de_idx);
void de_bwt_free(de_bwt_t *de_idx);

#endif
