#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "build_de_bwt.h"
#include "rest_index.h"
#include "kmer_hash.h"
#include "kstring.h"

void de_bwt_init(de_bwt_t *de_idx)
{
    de_idx->bwt_l = 0;
    de_idx->str_m = 1000000;
    de_idx->bwt_str = (kstring_t*)calloc(1, sizeof(kstring_t));

    de_idx->n_unipath = 0;
}

void de_bwt_free(de_bwt_t *de_idx)
{
    free(de_idx->bwt_str->s);
    free(de_idx->bwt_str);
}

int build_de_bwt(char *prefix, hash_idx *h_idx, de_bwt_t *de_idx)
{
    fprintf(stderr, "[build_de_bwt] Building de Bruijn-BWT index for genome ...\n");

    hash_init(h_idx); 
    de_bwt_init(de_idx);
    // FIRST run for generating k-mer
    //   update in/out flag
    //   update bwt_char
    //NEW: one-run for the whole generating
    kmer_gen(prefix, h_idx, de_idx);
    // SECOND run for generating unipath
    //   extra hash_num&node for $-contained kmer
    //   pos of unipath
    //   update bwt_char of head and tail kmer 
    //   uid and offset of kmer
    //kmer_re_gen(prefix, h_idx, de_idx);
    // merge normal kmer and $-contained kmer
    kmer_merge(*h_idx, de_idx);
    hash_free(h_idx);
    int i; for (i=0; i<de_idx->bwt_str->l; ++i) fprintf(stdout, "%c\n", de_idx->bwt_str->s[i]);
    
    fprintf(stderr, "[build_de_bwt] Building de Bruijn-BWT index done!\n");
    return 0;
}
