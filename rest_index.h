#ifndef _REST_INDEX_H_
#define _REST_INDEX_H_
#include <stdint.h>

int rest_index(int argc, char *argv[]);

typedef struct {
    uint8_t k;          // kmer len
    uint64_t k_m;       // 22: 0xfffffffffff
    uint8_t hash_k;     // first level's hash-kmer len 
    uint32_t hash_size; // pow(4, hash_klen)

    uint8_t remn_k;     // remain len
    uint8_t remn_n;     // remain len
    uint8_t remn_ni;    // pos of remain kmer

    uint8_t in_ni, out_ni, inout_m; // pos of in/out flag
    uint8_t bwt_char_ni,   char_m;   // pos of bwt_char
    uint8_t offset_ni;     // pos of offset on unipath
    uint32_t offset_m;
    uint8_t uid_ni;        // pos of unipath id
    uint32_t uid_m;
} hash_para;

#define kmer_num_t uint32_t 
#define kmer_node_t uint64_t 
//#define kmer_node_t uint32_t 
#define kmer_int_t uint64_t
#define _node_size 64
//#define _node_size 32 

typedef struct {
    hash_para hp; // hash_para
    kmer_num_t *kmer_num;
    kmer_num_t *kmer_m;
    kmer_node_t **kmer_node;
} 
hash_idx;

#endif
