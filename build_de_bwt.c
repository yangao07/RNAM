#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "build_de_bwt.h"
#include "rest_index.h"
#include "kmer_hash.h"

int build_de_bwt(hash_idx h_idx, char *prefix)
{
    fprintf(stderr, "[build_de_bwt] Building de Bruijn-BWT index for genome ...\n");

    hash_init(&h_idx);
    // FIRST run for generating k-mer
    //   update in/out flag
    //   update bwt_char
    kmer_gen(prefix, &h_idx);
    // SCOUND for generating unipath
    //   extra hash_num&node for $-contained kmer
    //   pos of unipath
    //   update bwt_char of head and tail kmer 
    //   uid and offset of kmer
    //unipath_gen(h_idx);
    // merge normal kmer and $-contained kmer
    //kmer_merge(h_idx, s_num, s_node, bwt_s);
    
    hash_free(&h_idx);
    fprintf(stderr, "[build_de_bwt] Building de Bruijn-BWT index done!\n");
    return 0;
}
