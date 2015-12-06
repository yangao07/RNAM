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
//     [1-18][19--20][21-23][24--36][37-64]
//      9mer  in/out  bwt_c  offset   uid
//                  for tmp: (next_c)
//                           [34--36]
void hash_init_idx_para(hash_idx *h)
{
    h->hp.k = 22;
    h->hp.k_m = 0xfffffffffff; // 44
    h->hp.hash_k = 13;
    h->hp.hash_size = pow(4, 13);
    
    h->hp.remn_k = 9;
    h->hp.remn_n = 18;
    h->hp.remn_ni = 46;

    h->hp.in_ni = 45;   // in
    h->hp.out_ni = 44;  // out
    h->hp.inout_m = 0x1;    // 1
    h->hp.bwt_char_ni = 41; 
    h->hp.char_m = 0x7;     // 3
    h->hp.offset_ni = 28;
    h->hp.offset_m = 0x1fff;// 13
    h->hp.uid_ni = 0;
    h->hp.uid_m= 0xfffffff; // 28
}

void hash_init_idx32_para(hash_idx *h)
{
    h->hp.k = 22;
    h->hp.k_m = 0xfffffffffff; // 44
    h->hp.hash_k = 13;
    h->hp.hash_size = pow(4, 13);
    
    h->hp.remn_k = 9;
    h->hp.remn_n = 18;
    h->hp.remn_ni = 12;

    h->hp.in_ni = 11;   // in
    h->hp.out_ni = 10;  // out
    h->hp.inout_m = 0x1;    // 1
    h->hp.bwt_char_ni = 7; 
    h->hp.char_m = 0x7;     // 3
    h->hp.offset_ni = 0;
    //h->hp.offset_m = 0x1fff;// 13
    //h->hp.uid_ni = 0;
    //h->hp.uid_m= 0xfffffff; // 28
}

void hash_reset_idx_para(hash_idx *h)
{
    if (h->hp.k != 22) {
        h->hp.k_m = 1;
        uint8_t i=0;
        while (i < h->hp.k) {
            h->hp.k_m <<= 1;
            h->hp.k_m |= 1;
            i++;
        }
        h->hp.hash_k = h->hp.k - 9;
        h->hp.hash_size = pow(4, h->hp.hash_k);
    }
}

int rest_index(int argc, char *argv[])
{
    char *prefix=0; int c;
    uint64_t l_pac;

    hash_idx h_idx;

    hash_init_idx_para(&h_idx);
    //hash_init_idx32_para(&h_idx);

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
        l_pac = bns_fasta2bntseq(fp, prefix, 0);
        fprintf(stderr, "done!\n");
        err_gzclose(fp);
    }*/
    { // generate de Bruijn graph and BWT
        build_de_bwt(h_idx, prefix);
    }
    free(prefix);
    return 0;
}
