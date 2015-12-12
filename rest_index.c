#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <zlib.h>
#include "rest_index.h"
#include "build_de_bwt.h"
#include "bntseq.h"
#include "utils.h"
//#include "error.h"

int rest_index_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   rest index [option] <ref.fa>\n");
    fprintf(stderr, "                    bulid index for <ref.fa>\n\n");
    fprintf(stderr, "Option:  \n");
    fprintf(stderr, "         -k [INT]     Length of kmer to construct de Bruijn graph. [Def=22]\n");
    fprintf(stderr, "         -s [INT]     Length of first level's hashed sequence. [Def=12]\n");
    fprintf(stderr, "\n");
    return 0;
}

// hash-node: (uint64_t)
// [1-18][19--20][21-23][24--36][37-64]
//  9mer  in/out  bwt_c  offset   uid
//              for tmp: (next_c)
//                       [34--36]
void hash_init_idx_para(hash_idx *h)
{
    h->hp.k = 22;
    h->hp.k_m = 0xfffffffffff; // 44
    h->hp.hash_k = 13;
    h->hp.hash_n = 26;
    h->hp.hash_m = 0x3ffffff; // 26
    h->hp.hash_size = pow(4, 13);
    
    h->hp.remn_k = 9;
    h->hp.remn_n = 18;
    h->hp.remn_m = 0x3ffff; // 18
    h->hp.remn_ni = 46;

    h->hp.in_ni = 45;   // in
    h->hp.out_ni = 44;  // out
    h->hp.inout_m = 0x1;    // 1

    h->hp.bwt_char_ni = 41; 
    h->hp.char_m = 0x7;     // 3

    h->hp.next_char_ni = 38;

    h->hp.uni_ni = 0;
    h->hp.uni_m = 0x1;

    /*
    h->hp.offset_ni = 28;
    h->hp.offset_m = 0x1fff;// 13

    h->hp.uid_ni = 0;
    h->hp.uid_m= 0xfffffff; // 28
    */
}

#ifdef __LIT__
void hash_init_idx32_para(hash_idx *h)
{
    h->hp.k = 4;
    h->hp.k_m = 0xff; // 8
    h->hp.hash_k = 2;
    h->hp.hash_n = 4;
    h->hp.hash_m = 0xf; // 4
    h->hp.hash_size = pow(4, 2);
    
    h->hp.remn_k = 2;
    h->hp.remn_n = 4;
    h->hp.remn_m = 0xf; // 4
    h->hp.remn_ni = 28;

    h->hp.in_ni = 27;   // in
    h->hp.out_ni = 26;  // out
    h->hp.inout_m = 0x1;    // 1

    h->hp.bwt_char_ni = 23; 
    h->hp.char_m = 0x7;     // 3
    h->hp.next_char_ni = 20;

    h->hp.sk_ni = 1;
    h->hp.sk_n = 5;
    h->hp.sk_m = 0x1f; // 5

    h->hp.uni_ni = 0;  
    h->hp.uni_m = 0x1;
}

void hash_reset_idx_para(hash_idx *h)
{
    if (h->hp.k != 4) {
        h->hp.k_m = 0;
        uint8_t i=0;
        while (i < h->hp.k << 1) {
            h->hp.k_m <<= 1;
            h->hp.k_m |= 1;
            i++;
        }
        h->hp.hash_k = h->hp.k - 2;
        h->hp.hash_n = h->hp.hash_k << 1;
        h->hp.hash_m = 0;
        i=0;
        while (i < h->hp.hash_n) {
            h->hp.hash_m <<= 1;
            h->hp.hash_m |= 1;
            i++;
        }
        h->hp.hash_size = pow(4, h->hp.hash_k);
    }
}
#else
// hash-node: (uint32_t)
// [1-18][19--20][21-23][24--26][27-31][32]
//  9-mer in/out  bwt_c  next_c  k-len  uni
//                               k-len: for special-kmer, len<k (2^5=32)
void hash_init_idx32_para(hash_idx *h)
{
    h->hp.k = 22;
    h->hp.k_m = 0xfffffffffff; // 44
    h->hp.hash_k = 13;
    h->hp.hash_n = 26;
    h->hp.hash_m = 0x3ffffff; // 26
    h->hp.hash_size = pow(4, 13);
    
    h->hp.remn_k = 9;
    h->hp.remn_n = 18;
    h->hp.remn_m = 0x3ffff; // 18
    h->hp.remn_ni = 14;

    h->hp.in_ni = 13;   // in
    h->hp.out_ni = 12;  // out
    h->hp.inout_m = 0x1;    // 1

    h->hp.bwt_char_ni = 9; 
    h->hp.char_m = 0x7;     // 3
    h->hp.next_char_ni = 6;

    h->hp.sk_ni = 1;
    h->hp.sk_n = 5;
    h->hp.sk_m = 0x1f; // 5
    
    h->hp.uni_ni = 0;
    h->hp.uni_m = 0x1;
}

void hash_reset_idx_para(hash_idx *h)
{
    if (h->hp.k != 22) {
        h->hp.k_m = 0;
        uint8_t i=0;
        while (i < h->hp.k << 1) {
            h->hp.k_m <<= 1;
            h->hp.k_m |= 1;
            i++;
        }
        h->hp.hash_k = h->hp.k - 9;
        h->hp.hash_n = h->hp.hash_k << 1;
        h->hp.hash_m = 0;
        i=0;
        while (i < h->hp.hash_n) {
            h->hp.hash_m <<= 1;
            h->hp.hash_m |= 1;
            i++;
        }
        h->hp.hash_size = pow(4, h->hp.hash_k);
    }
}
#endif

int rest_index(int argc, char *argv[])
{
    char *prefix=0; int c;

    hash_idx h_idx; de_bwt_t de_idx;

    //hash_init_idx_para(&h_idx);
    hash_init_idx32_para(&h_idx);

    while ((c = getopt(argc, argv, "k:")) >= 0)
    {
        switch (c)
        {
            case 'k': h_idx.hp.k = atoi(optarg); break;
            default: return rest_index_usage();
        }
    }
    if (optind + 1 > argc) return rest_index_usage();
    hash_reset_idx_para(&h_idx);
    prefix = strdup(argv[optind]);

    /*{ // generate for&rev.pac
        gzFile fp = xzopen(prefix, "r");
        fprintf(stderr, "[%s] Pack genome FASTA ... ", __func__);
        uint64_t l_pac = bns_fasta2bntseq(fp, prefix, 0);
        fprintf(stderr, "done!\n");
        err_gzclose(fp);
    }*/
    { // generate de Bruijn graph and BWT
        build_de_bwt(prefix, &h_idx, &de_idx);
        de_bwt_free(&de_idx);
    }
    free(prefix);
    return 0;
}
