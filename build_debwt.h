#ifndef _BUILD_DE_BWT_
#define _BUILD_DE_BWT_
#include <stdint.h>
#include "rest_index.h"

int build_debwt(char *prefix, hash_idx *h_idx, debwt_t *de_idx);
void debwt_free(debwt_t *de_idx);
void push_debwt_bwt(uint8_t bwt_nt, debwt_t *db_idx);

#endif
