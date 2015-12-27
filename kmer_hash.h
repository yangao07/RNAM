#ifndef _KMER_HASH_H_
#define _KMER_HASH_H_
#include "rest_index.h"

#define nt_A 0
#define nt_C 1
#define nt_G 2
#define nt_T 3
#define nt_N 4

void hash_init_num(hash_idx *h_idx);
void hash_init_node(hash_idx *h_idx);
void hash_free(hash_idx *h_idx);

int kmer_tol_count(const char *prefix, hash_idx *h_idx);
int kmer_gen(const char *prefix, hash_idx *h_idx);
int spe_kmer_gen(const char *prefix, hash_idx *h_idx, debwt_t *de_idx);
int kmer_merge(hash_idx h_idx, debwt_t *de_idx);

#endif
