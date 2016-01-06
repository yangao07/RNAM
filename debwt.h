#ifndef _BUILD_DE_BWT_
#define _BUILD_DE_BWT_
#include <stdint.h>
#include "rest_index.h"

#define debwt_int_t uint64_t
#define _DEBWT_INT_SIZE  64
#define debwt_hash_t uint32_t

// default schema of bwt
// [5 OCC][4 BWT-INT]                                [5 OCC][4 BWT-INT] [...][...] [5 OCC][last n BWT-INT(n<=4)]
//        [4-4-.-4] [4--4-.-4] [4--4-.-4] [4--4-.-4]                                                ..[4][4]
// index:  0|1| 15| 16|17| 31| 32|33| 47| 48|49| 63|                                          ..|64i+n-1|64i+n|
//
// MACRO of DEBWT_PARA
#define _OCC_C 5         // 5; A/C/G/T/#
#define _BWT_N 4
#define _BWT_NT_K 4      // bits of per bwt_char 
#define _BWT_NT_B 2
#define _BWT_NT_M 0x7    // 0/1/2/3/4
#define _BWT_INV 16      // number of bwt char in one bwt_int; 16=64/4
#define _BWT_INV_B 4  
#define _BWT_INV_M 0xf   // 15

#define _OCC_INV 64      //store C[5] for every occ_inv bwt_char; 64=16*4
#define _OCC_INV_B 6
#define _OCC_INV_M 0x3f  // 63

#define _BWT_OCC_B 3     // 9=2^3+1
#define _SA_INV 32
#define _SA_INV_B 5
#define _SA_INV_M 0x1f

//#define _BWT_HASH_K 13

typedef struct {
    uint8_t occ_c;              
    uint8_t bwt_nt_k, bwt_nt_b, bwt_nt_m; 
    uint8_t bwt_inv, bwt_inv_b, bwt_inv_m; //  16,4,0xf (16=64/4=2^4)
    uint8_t occ_inv, occ_inv_b, occ_inv_m; //  64,6,0x3f (64=16*4=2^6)
    uint8_t bwt_occ_b;          // 3; 5+4=9=2^3+1

    uint8_t sa_inv;             // 32
    uint8_t sa_inv_b;           // 5
} debwt_para;

typedef struct {
    //debwt_para dp;
    debwt_count_t bwt_l;         // length of bwt_str(seq_len)
    debwt_count_t bwt_size;      // size of *bwt
    debwt_int_t   *bwt;          // bwt_str & OCC
    debwt_int_t   bwt_unit;      // for push_bwt
    debwt_count_t bwt_i, bwt_k;  // index of bwt_str and bwt
    uint8_t bit_table16[65536];
    //uint8_t bit_table32[536870912]; // pow(2,19)

    debwt_count_t n_unipath;     // total count of unipath
    debwt_count_t n_offset;      // total count of offsets

    debwt_count_t n_sa, n_s_sa;
    debwt_count_t *sa_uid, *sa_u_off; 
    debwt_count_t *s_sa_uid;     // normal SA and special SA for all the '#-bwt_char'

    debwt_count_t *uni_offset_c; // cumulative number of offsets for each unipath
    ref_offset_t  *uni_offset;   // offsets of each unipath
    debwt_count_t C[_OCC_C];     // cumulative count of 'A/C/G/T/#'
                                 // after update: [0]: num of #
                                 //               [1]: num of #, A
                                 //               [2]: num of #, A, C
                                 //               [3]: num of #, A, C, G
                                 //               [4]: num of #, A, C, G, T
    uint8_t bwt_hash_k;
    int bwt_hash_size;
    debwt_count_t *bwt_hash;
    debwt_count_t *bwt_hash_il;  // interval length
} debwt_t;

typedef struct {
    debwt_count_t u, l;
} debwt_intv_t;



int build_debwt(char *prefix, hash_idx *h_idx, debwt_t *de_idx, int rev);
void debwt_free(debwt_t *de_idx);
void push_debwt_bwt(uint8_t bwt_nt, debwt_t *db_idx);
int pac_build_debwt(const char *prefix, debwt_pac_t *db_pac, debwt_count_t l_pac, hash_idx *h_idx, debwt_t *de_idx);

#endif
