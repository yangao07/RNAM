#ifndef _REST_INDEX_H_
#define _REST_INDEX_H_
#include <stdint.h>
#include "kstring.h"

int rest_index(int argc, char *argv[]);

typedef struct {
    uint8_t k;          // kmer len
    uint8_t k_n;
    uint64_t k_m;       // 22: 0xfffffffffff

    uint8_t hash_k;     // first level's hash-kmer len 
    uint8_t hash_n;
    uint64_t hash_m;
    uint32_t hash_size; // pow(4, hash_klen)

    uint8_t remn_k;     // remain len
    uint8_t remn_n;     // remain len
    uint64_t remn_m;
    uint8_t remn_ni;    // pos of remain kmer

    uint8_t in_ni, out_ni, inout_m; // pos of in/out flag
    uint8_t bwt_char_ni,   char_m;  // pos of bwt_char
    uint8_t next_char_ni;           // pos of next_char

    uint8_t sk_ni, sk_n, sk_m;      // for special-kmer, len<k

    uint8_t spe_ni;// uni_flag: generate #-kmer
    uint8_t spe_m; // 0x1
} hash_para;

#define kmer_int_t uint64_t
#define hash_int_t uint32_t
#define remn_int_t uint32_t
#define _int_size 64
#define kmer_num_t uint32_t 
//#define kmer_node_t uint64_t 
#define kmer_node_t uint32_t 
//#define _node_size 64
#define _node_size 32 

typedef struct {
    hash_para hp; // hash_para
    
    // [key+1]: num/m of key
    uint64_t *kmer_m; // [0]=0, [i] = cumu_count of kmer
    kmer_num_t *kmer_num; // whole count of existing kmer
    uint64_t kmer_tol_count;
    kmer_node_t *kmer_node;

    uint64_t *skmer_m;   // hash index for special kmer
    kmer_num_t *skmer_num;   
    uint64_t skmer_tol_count;
    kmer_node_t *skmer_node;
} hash_idx;

typedef struct {
    uint64_t bwt_l;       // length of bwt
    uint64_t str_m;
    //char *bwt_str;        // bwt string
    kstring_t *bwt_str;

    uint64_t n_unipath;   // number of unipath
    uint64_t *uni_offset; // offset of each unipath
    uint64_t C[5];        // cumulative count of 'A/C/G/T/#'
} de_bwt_t;

#endif
