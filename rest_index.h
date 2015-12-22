#ifndef _REST_INDEX_H_
#define _REST_INDEX_H_
#include <stdint.h>
#include "kstring.h"

int rest_index(int argc, char *argv[]);

typedef struct {
    uint8_t k, k_n;        // kmer len
    uint64_t k_m;          // 22: 0xfffffffffff

    uint8_t hash_k, hash_n;// first level's hash-kmer len 
    uint64_t hash_m;
    uint32_t hash_size;    // pow(4, hash_k)

    uint8_t uid_n, uid_ni; // for spe-kmer
    uint64_t uid_m;

    uint8_t remn_n, remn_k, remn_ni; // remain len
    uint64_t remn_m;

    uint8_t in_ni, out_ni, inout_m;  // pos of in/out flag

    uint8_t bwt_char_ni, next_char_ni, char_m; // pos of bwt_char/next_char

    uint8_t sk_ni, sk_n, sk_m; // for special-kmer, len<k

    uint8_t spe_ni, spe_m;     // spe_flag: generate #-kmer
} hash_para;

#define kmer_int_t uint64_t
#define hash_int_t uint32_t
#define remn_int_t uint32_t
#define _int_size 64
#define debwt_count_t uint32_t
#define kmer_num_t uint32_t  // for each hashKey
#define kmer_node_t uint32_t 
#define _node_size 32 
#define _snode_size 64
#define ref_offset_t uint32_t // for forward-only
//#define ref_offset_t uint64_t // for forward and reverse
#define skmer_node_t uint64_t

typedef struct {
    hash_para hp;               // hash_para
    
    /* hash index for normal kmer */
    // [key+1]: num/m of key
    debwt_count_t *kmer_c;       // [0]=0, [i] = cumu_count of kmer
    kmer_num_t    *kmer_num;     // whole count of existing kmer
    debwt_count_t kmer_tol_count;
    debwt_count_t kmer_real_count;
    kmer_node_t   *kmer_node;
                 // same to kmer_node, size: kmer_tol_count
    debwt_count_t *uni_offset_c; // (cumulative) number of offsets for each unipath
    debwt_count_t next_offset_c;
                 // same to kmer_num/skmer_num, size: uni_tol_count
    debwt_count_t *uni_offset_n;
                 // smae to kmer_node, size: kmer_tol_count
    debwt_count_t *uni_id;       // UID of normal kmer
    debwt_count_t last_uid, cur_uid, next_uid;
                 // smae to skmer_node
    ref_offset_t  *uni_offset;   // 1-base
    debwt_count_t uni_tol_count; // uni_tol = skmer_tol / (k-1)
    debwt_count_t off_tol_count;

    /* hash index for special kmer */
    debwt_count_t *skmer_c;      
    kmer_num_t    *skmer_num;   
    debwt_count_t skmer_tol_count;
    skmer_node_t  *skmer_node;   // [UID][NODE]
                                //  1-32 33-64 
} hash_idx;

typedef struct {
    debwt_count_t bwt_l;         // length of bwt
    //XXX
    kstring_t     *bwt_str;      // bwt_str

    debwt_count_t n_unipath;     // total count of unipath
    debwt_count_t n_offset;      // total count of offsets
    //XXX
    debwt_count_t *uni_id;
    //XXX
    debwt_count_t *uni_offset_c; // cumulative number of offsets for each unipath
    //XXX
    ref_offset_t  *uni_offset;   // offsets of each unipath
    debwt_count_t C[5];          // cumulative count of 'A/C/G/T/#'
} de_bwt_t;

#endif
