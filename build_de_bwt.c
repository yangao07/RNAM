#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "build_de_bwt.h"
#include "rest_index.h"
#include "kmer_hash.h"
#include "kstring.h"
#include "utils.h"

void de_bwt_init(de_bwt_t *de_idx, hash_idx h_idx)
{
    de_idx->bwt_l = h_idx.kmer_real_count + h_idx.hp.k * h_idx.uni_tol_count;
    de_idx->bwt_str = (kstring_t*)calloc(1, sizeof(kstring_t));

    de_idx->n_unipath = h_idx.uni_tol_count;
    de_idx->n_offset = h_idx.off_tol_count;

    de_idx->uni_id = (debwt_count_t*)malloc(de_idx->bwt_l * sizeof(debwt_count_t));
    de_idx->uni_offset_c = (debwt_count_t*)malloc((de_idx->n_unipath+1) * sizeof(debwt_count_t));
    de_idx->uni_offset = (debwt_count_t*)malloc(de_idx->n_offset * sizeof(debwt_count_t));
}

void de_bwt_free(de_bwt_t *de_idx)
{
    free(de_idx->bwt_str->s); free(de_idx->bwt_str);
    free(de_idx->uni_id); free(de_idx->uni_offset); free(de_idx->uni_offset_c);
}

void print_de_bwt(de_bwt_t de_idx)
{
    stdout_printf("BWT_STR:\n%s\n", de_idx.bwt_str->s);
    stdout_printf("UID:\n");
    uint32_t i, j;
    /*
    for (i = 0; i < de_idx.bwt_l-de_idx.n_unipath; ++i) {
        fprintf(stdout, "# %d:\t%d\n", i+1, de_idx.uni_id[i]); 
    }*/
    stdout_printf("Uni_offset_c:\n");
    for (i = 0; i < de_idx.n_unipath; ++i) stdout_printf("U_c %d:\t%d\n", i+1, de_idx.uni_offset_c[i+1]);
    stdout_printf("Uni_offset:\n");
    for (i = 0; i < de_idx.n_unipath; ++i) {
        stdout_printf("U %d:\t", i+1);
        for (j = de_idx.uni_offset_c[i]; j < de_idx.uni_offset_c[i+1]; ++j)
            stdout_printf("%d ", de_idx.uni_offset[j]);
        stdout_printf("\n");
    }
}

/*
#define bwt_pre_pos(db, i) ()

void de_bwt_cal_sa(de_bwt_t *db_idx, int intv)
{
    int uni_id_i, uni_id;
    int i;

    // calculate uni_id
    for (i = 0; i < uni_n; ++i) {
        uni_id_i = i, uni_id = i;
        while () {
            if (uni_id_i % intv == 0) db_idx->uni_id[uni_id_i/intv] = uni_id;
            uni_id_i = bwt_pre_pos(db_idx, uni_id_i);
        }
        if (uni_id_i % intv == 0) db_idx->uni_id[uni_id_i/intv] = uni_id;
    }
}

// calculate 
void de_bwt_update_bwt(de_bwt_t *db_idx)
{
    int n_occ, i, k;
    uint64_t c[4];

    c[0] = c[1] = c[2] = c[3] = 0;

    uint64_t *buf;

    for (i = k = 0; i < db_idx->bwt_l; ++i) {
        if (i % occ_interval = 0) {
            memcpy(buf+k, c, sizeof(uint64_t) * 4);
            k += sizeof(uint64_t);
        }
        if (i % 21 == 0) buf[k++] = db_idx->bwt[i/21]; // 64/3 = 21...1
        ++c[bwt_nt(db_idx, i)];
    }
    // last
    memcpy(buf+k, c, sizeof(uint64_t) * 4);
    // update bwt
    free(db_idx->bwt); db_idx->bwt = buf;
}
*/

int build_de_bwt(char *prefix, hash_idx *h_idx, de_bwt_t *de_idx)
{
    err_printf("[build_de_bwt] Building de Bruijn-BWT index for genome ...\n");

    hash_init_num(h_idx); 
    // FIRST run for counting the total number of kmer
    kmer_tol_count(prefix, h_idx);
    // SECOND run for generating k-mer and countint spe-kmer
    //   update in/out flag
    //   update bwt_char
    //   count the total number of spe-kmer
    kmer_gen(prefix, h_idx);
    // THIRD run for generating unipath
    //   extra hash_num&node for $-contained kmer
    //   pos of unipath
    //   update bwt_char of head and tail kmer 
    //   uid and offset of kmer
    de_bwt_init(de_idx, *h_idx);
    spe_kmer_gen(prefix, h_idx, de_idx);
    // merge normal kmer and $-contained kmer
    kmer_merge(*h_idx, de_idx);
    hash_free(h_idx);
    print_de_bwt(*de_idx);
    
    err_printf("[build_de_bwt] Building de Bruijn-BWT index done!\n");
    return 0;
}
