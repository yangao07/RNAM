#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include "kmer_hash.h"
#include "rest_index.h"
#include "kseq.h"
#include "utils.h"
#include "string.h"

KSEQ_INIT(gzFile, gzread)

char BWT_STR[9] = "ACGTNNNN";

extern unsigned char nst_nt4_table[256];

void hash_init(hash_idx *h_idx)
{
    h_idx->kmer_num = (kmer_num_t*)_err_calloc(h_idx->hp.hash_size, sizeof(kmer_num_t));
    h_idx->kmer_m = (kmer_num_t*)_err_malloc(h_idx->hp.hash_size *  sizeof(kmer_num_t));
    h_idx->kmer_node = (kmer_node_t**)_err_malloc(h_idx->hp.hash_size * sizeof(kmer_node_t*));

    h_idx->skmer_num = (kmer_num_t*)_err_calloc(h_idx->hp.hash_size, sizeof(kmer_num_t));
    h_idx->skmer_m = (kmer_num_t*)_err_malloc(h_idx->hp.hash_size * sizeof(kmer_num_t));
    h_idx->skmer_node = (kmer_node_t**)_err_malloc(h_idx->hp.hash_size * sizeof(kmer_node_t*));

    uint32_t i;
    for (i = 0; i < h_idx->hp.hash_size; ++i) {
        h_idx->kmer_m[i] = 1;
        h_idx->kmer_node[i] = (kmer_node_t*)_err_malloc(sizeof(kmer_node_t));

        h_idx->skmer_m[i] = 1;
        h_idx->skmer_node[i] = (kmer_node_t*)_err_malloc(sizeof(kmer_node_t));
    }
}

void hash_free(hash_idx *h_idx)
{
    uint32_t i;
    for (i = 0; i < h_idx->hp.hash_size; ++i) {
        free(h_idx->kmer_node[i]);
        free(h_idx->skmer_node[i]);
    }
    free(h_idx->kmer_node); free(h_idx->kmer_num); free(h_idx->kmer_m);
    free(h_idx->skmer_node); free(h_idx->skmer_num); free(h_idx->skmer_m);
}

void cons_kmer(kmer_node_t n, hash_int_t hashKey, uint8_t k, uint8_t remn_n, uint8_t remn_ni)
{
    char *seq = (char*)malloc(k+1);
    kmer_int_t all = hashKey;
    all = (all << remn_n) | (n>>(remn_ni));
    int i;
    for (i = 0; i < k; ++i)
    {
        seq[i] = "ACGTN"[(all >> (2 * (k-i-1)))&0x3];
    }
    seq[k] = 0;
    fprintf(stdout, "%s\n", seq);
    free(seq);
}

void con_kmer(kmer_int_t all, uint8_t k)
{
    char *seq = (char*)malloc(k+1);
    int i;
    for (i = 0; i < k; ++i) {
        seq[i] = "ACGTN"[(all >> (2 * (k-i-1)))&0x3];
    }
    seq[k] = 0;
    fprintf(stderr, "%s %lld\n", seq, (long long)all);
    free(seq);
} 

int get_first_kmer(kseq_t *seq, int *seq_i, uint8_t k, kmer_int_t *kmer, uint8_t *next_nt)
{
    int k_i = 0, i = *seq_i;
    uint8_t nt;
    *kmer = 0;
    while (i < (int)seq->seq.l) {
        nt = nst_nt4_table[(int)seq->seq.s[i++]];
        if (nt != nt_N) {
            (*kmer) = (*kmer) << 2 | nt;
            k_i++;
        } else {
            k_i = 0; *kmer = 0;
        }
        if (k_i == k) break;
    }
    *seq_i = i;
    *next_nt = nst_nt4_table[(int)seq->seq.s[i]];
    return (k_i == k);
}

int get_shift_kmer(kseq_t *seq, int *seq_i, uint8_t k, uint64_t k_m, kmer_int_t *kmer, uint8_t *pre_nt, uint8_t *next_nt)
{
    size_t i = *seq_i;
    while (i < seq->seq.l) {
        i++;
        if (*next_nt != nt_N) {
            *pre_nt = (*kmer) >> ((k-1)<<1);
            (*kmer) = ((*kmer) << 2 | (*next_nt)) & k_m; // 44 
            if (i == seq->seq.l) *next_nt = nt_N;
            else *next_nt = nst_nt4_table[(int)seq->seq.s[i]];
            *seq_i = i;
            return 1;
        } else {
            *seq_i = i;
            if (!get_first_kmer(seq, seq_i, k, kmer, next_nt)) return 0;
            *pre_nt = nt_N;
            return 1;
        }
    }
    return 0;
}

#define first_n_bits(type, x, l, n)  (type)((x) >> ((l)-(n)))
#define _first_n_bits(type, x, l, n) (type)(((x) >> ((l)-(n))) << ((l)-(n)))
#define last_n_bits(type, x, l, n)   (type)(((x) << ((l)-(n))) >> ((l)-(n)))
#define _last_n_bits(type, x, l, n)  (type)((x) << ((l)-(n)))

#define set_kmer_key(n, kmerKey, h) { (*(n)) = (kmerKey << ((h)->hp.remn_ni)); }
#define set_MinEdge(n, h) { (*(n)) |= (1 << ((h)->hp.in_ni)); }
#define set_MoutEdge(n, h) { (*(n)) |= (1 << ((h)->hp.out_ni)); }
#define set_bwt_char(n, last_nt, h) { (*(n)) |= (last_nt << ((h)->hp.bwt_char_ni)); }    // [0:12]: [N:ACGT]
#define set_next_char(n, next_nt, h) { (*(n)) |= (next_nt << ((h)->hp.next_char_ni)); }
#define set_sk_n(n, sk_n, h) { (*(n)) |= ((sk_n) << ((h)->hp.sk_ni)); }
#define set_spe(n, h) {(*(n)) |= (1 << ((h)->hp.spe_ni)); }

#define get_kmer_key(n, h)  ((n) >> ((h).hp.remn_ni))
#define has_MinEdge(n, h) (((n) >> ((h).hp.in_ni)) & ((h).hp.inout_m))
#define has_MoutEdge(n, h) (((n) >> ((h).hp.out_ni)) & ((h).hp.inout_m))
#define get_bwt_char(n, h) (((n) >> ((h).hp.bwt_char_ni)) & ((h).hp.char_m))
#define get_next_char(n, h) (((n) >> ((h).hp.next_char_ni)) & ((h).hp.char_m))
#define get_sk_n(n, h) (((n) >> ((h).hp.sk_ni)) & ((h).hp.sk_m))
#define get_spe(n, h) (((n) >> ((h).hp.spe_ni)) & ((h).hp.spe_m))

kmer_num_t bin_search_kmer(uint8_t *ef, hash_idx h, hash_int_t hashKey, remn_int_t kmerKey)
{
    remn_int_t tmpKey;
    int left=0, right = h.kmer_num[hashKey]-1, middle;
    if (right == -1) { *ef = 0; return 0; }

    if (kmerKey < first_n_bits(remn_int_t, h.kmer_node[hashKey][left], _node_size, h.hp.remn_n)) { *ef = 0; return left; }
    else if (kmerKey > first_n_bits(remn_int_t, h.kmer_node[hashKey][right], _node_size, h.hp.remn_n)) { *ef = 0; return right+1; }

    // binary search for hashInt(hashKey, kmerKey)
    while (left <= right) {
        middle = (left+right)>>1;
        tmpKey = first_n_bits(remn_int_t, h.kmer_node[hashKey][middle], _node_size, h.hp.remn_n);
        if (tmpKey == kmerKey) { *ef = 1; return middle; }
        else if (tmpKey > kmerKey) {
            if ((middle == 0) || (first_n_bits(remn_int_t, h.kmer_node[hashKey][middle-1], _node_size, h.hp.remn_n) < kmerKey)) {
                *ef = 0; return middle;
            } else right = middle-1;
        } else left = middle+1;
    }
    *ef = 0;
    return h.kmer_num[hashKey];
}



/*kmer_num_t bin_spe_search_kmer(hash_idx h, hash_int_t hashKey, remn_int_t kmerKey, uint8_t sk_n)
{
    remn_int_t tmpKey; uint8_t tmp_sk_n;
    int left=0, right = h.skmer_num[hashKey]-1, middle;
    if (right == -1) { return 0; }

    if (kmerKey < first_n_bits(remn_int_t, h.skmer_node[hashKey][left], _node_size, h.hp.remn_n)) { return left; }
    else if (kmerKey > first_n_bits(remn_int_t, h.skmer_node[hashKey][right], _node_size, h.hp.remn_n)) { return right+1; }

    // binary search for hashInt(hashKey, kmerKey)
    while (left <= right) {
        middle = (left+right)>>1;
        tmpKey = first_n_bits(remn_int_t, h.skmer_node[hashKey][middle], _node_size, h.hp.remn_n);
        tmp_sk_n = get_sk_n(h.skmer_node[hashKey][middle], h);
        if (tmpKey > kmerKey || (tmpKey == kmerKey && tmp_sk_n > sk_n)) { // bigger
            if ((middle == 0) || (first_n_bits(remn_int_t, h.skmer_node[hashKey][middle-1], _node_size, h.hp.remn_n) < kmerKey) 
                    || (first_n_bits(remn_int_t, h.skmer_node[hashKey][middle-1], _node_size, h.hp.remn_n) == kmerKey && get_sk_n(h.skmer_node[hashKey][middle-1], h) <= sk_n))
                return middle;
            right = middle-1;
        } else if (tmpKey < kmerKey || (tmpKey == kmerKey && tmp_sk_n <= sk_n))
            left = middle+1;
        else err_fatal_simple("ERROR: same skmer and sk_n!\n");
    }
    return h.skmer_num[hashKey];
}*/

int bin_hit_kmer(kmer_num_t *ki, hash_idx *h, hash_int_t hashKey, remn_int_t kmerKey)
{
    remn_int_t tmpKey;
    int left=0, right = h->kmer_num[hashKey]-1, middle;
    if (right == -1) return 0;

    // binary search for hashInt(hashKey, kmerKey)
    while (left <= right) {
        middle = (left+right)>>1;
        tmpKey = first_n_bits(remn_int_t, h->kmer_node[hashKey][middle], _node_size, h->hp.remn_n);
        if (tmpKey == kmerKey) {
            *ki = middle; return 1;
        } else if (tmpKey > kmerKey)  right = middle-1;
        else left = middle+1;
    }
    return 0;
}

int update_spe_hash(hash_idx *h, hash_int_t hashKey, remn_int_t kmerKey, uint8_t pre_nt, uint8_t sk_n)
{
    kmer_num_t k_i = h->skmer_num[hashKey];
    if (k_i >= h->skmer_m[hashKey]) {
        h->skmer_m[hashKey] <<= 1;
        h->skmer_node[hashKey] = (kmer_node_t*)_err_realloc(h->skmer_node[hashKey], h->skmer_m[hashKey] * sizeof(kmer_node_t));
    }
    // copy node
    //memmove(h->skmer_node[hashKey]+k_i+1, h->skmer_node[hashKey]+k_i, (k_n-1-k_i)*sizeof(kmer_node_t));
    kmer_node_t t = kmerKey;
    set_kmer_key(h->skmer_node[hashKey]+k_i, t, h);
    set_bwt_char(h->skmer_node[hashKey]+k_i, pre_nt, h);
    set_sk_n(h->skmer_node[hashKey]+k_i, sk_n, h);
    ++(h->skmer_num[hashKey]);
    return 1;
}

void pre_gen_spe_kmer(uint8_t pre_nt, hash_int_t pre_hashKey, kmer_num_t pre_k_i, kmer_int_t kmerInt, hash_idx *h, de_bwt_t *db) 
{
    if (pre_nt == nt_N || get_spe(h->kmer_node[pre_hashKey][pre_k_i], *h)) return;

    kmer_int_t s_kmerInt = kmerInt - (kmerInt & 0x3);

    hash_int_t s_hashKey; remn_int_t s_kmerKey;

    int i; uint8_t bwt_nt = pre_nt;
    for (i = 0; i < h->hp.k-1; ++i) {
        s_hashKey = first_n_bits(hash_int_t, s_kmerInt, h->hp.k_n, h->hp.hash_n);
        s_kmerKey = last_n_bits(remn_int_t, s_kmerInt, _int_size, h->hp.remn_n);
        //kmer_num_t k_i = bin_spe_search_kmer(*h, s_hashKey, s_kmerKey, h->hp.k-1-i);
        if (!update_spe_hash(h, s_hashKey, s_kmerKey, bwt_nt, h->hp.k-1-i)) err_fatal_simple("ERROR: hash-update fail(1)!\n");
        bwt_nt = first_n_bits(uint8_t, s_kmerInt, h->hp.k_n, 2);
        s_kmerInt = (s_kmerInt<<2) & h->hp.k_m;
    }
    kputc(BWT_STR[bwt_nt], db->bwt_str);
}

void cur_gen_spe_kmer(kmer_int_t kmerInt, hash_int_t hashKey, kmer_num_t k_i, hash_idx *h, de_bwt_t *db) 
{
    if (get_spe(h->kmer_node[hashKey][k_i], *h)) return;

    kmer_int_t s_kmerInt = kmerInt;
    hash_int_t s_hashKey; remn_int_t s_kmerKey;
    int i; uint8_t bwt_nt;
    for (i = 0; i < h->hp.k-1; ++i) {
        bwt_nt = first_n_bits(uint8_t, s_kmerInt, h->hp.k_n, 2);
        s_kmerInt = (s_kmerInt << 2) & h->hp.k_m;

        s_hashKey = first_n_bits(hash_int_t, s_kmerInt, h->hp.k_n, h->hp.hash_n);
        s_kmerKey = last_n_bits(remn_int_t, s_kmerInt, _int_size, h->hp.remn_n);
        //kmer_num_t k_i = bin_spe_search_kmer(*h, s_hashKey, s_kmerKey, h->hp.k-1-i);
        if (!update_spe_hash(h, s_hashKey, s_kmerKey, bwt_nt, h->hp.k-1-i)) err_fatal_simple("ERROR: hash-update fail(2)!\n");
    }
    kputc(BWT_STR[first_n_bits(uint8_t, s_kmerInt, (h->hp.k)<<1, 2)], db->bwt_str);
}

uint8_t update_hash(kmer_num_t k_i, uint8_t ef, hash_idx *h, de_bwt_t *db, kmer_int_t kmerInt, hash_int_t hashKey, remn_int_t kmerKey, uint8_t pre_nt, hash_int_t pre_hashKey, kmer_num_t pre_k_i, uint8_t next_nt, uint8_t cur_spe_flag)
{
    uint8_t next_spe_flag = 0;
    if (ef) {
        { // in
            // 1-in => m-in
            if (!has_MinEdge(h->kmer_node[hashKey][k_i], *h)) {
                if (pre_nt != get_bwt_char(h->kmer_node[hashKey][k_i], *h) || pre_nt == nt_N) {
                    // spe-kmer for old_pre_kmer and new_pre_kmer(NOT N), set spe-flag
                    uint8_t old_pre_nt = get_bwt_char(h->kmer_node[hashKey][k_i], *h);
                    if (old_pre_nt < nt_N) {
                        kmer_int_t old_pre_kmerInt = old_pre_nt; 
                        old_pre_kmerInt = (old_pre_kmerInt << (h->hp.k_n-2)) | (kmerInt >> 2);

                        hash_int_t old_pre_hashKey = first_n_bits(hash_int_t, old_pre_kmerInt, h->hp.k_n, h->hp.hash_n);
                        remn_int_t old_pre_kmerKey = last_n_bits(remn_int_t, old_pre_kmerInt, _int_size, h->hp.remn_n);
                        kmer_num_t old_pre_k_i;
                        if (!bin_hit_kmer(&old_pre_k_i, h, old_pre_hashKey, old_pre_kmerKey)) err_fatal_core(__func__, "ERROR: Cannot hit in hash-table(1)!\npre_nt: %d\n", old_pre_nt);
                        pre_gen_spe_kmer(old_pre_nt, old_pre_hashKey, old_pre_k_i, kmerInt, h, db);
                    }
                    pre_gen_spe_kmer(pre_nt, pre_hashKey, pre_k_i, kmerInt, h, db); 

                    set_MinEdge(h->kmer_node[hashKey]+k_i, h);
                    set_bwt_char(h->kmer_node[hashKey]+k_i, nt_N, h);
                    // uid, uoffset XXX
                }
            }
            // m-in => m-in
            else pre_gen_spe_kmer(pre_nt, pre_hashKey, pre_k_i, kmerInt, h, db); // spe-kmer for new_pre_kmer(NOT N), set spe-flag
        }
        { // out
            // 1-out => m-out
            if (!has_MoutEdge(h->kmer_node[hashKey][k_i], *h)) {
                if (next_nt != get_next_char(h->kmer_node[hashKey][k_i], *h) || next_nt == nt_N) {
                    // spe-kmer for cur_kmer, set spe-flag
                    cur_gen_spe_kmer(kmerInt, hashKey, k_i, h, db);
                    // set_bwt_char(old_next_kmer, #)
                    uint8_t old_next_nt = get_next_char(h->kmer_node[hashKey][k_i], *h);
                    if (old_next_nt < nt_N) {
                        kmer_int_t old_next_kmerInt = ((kmerInt << 2) & h->hp.k_m) | old_next_nt;
                        hash_int_t old_next_hashKey = first_n_bits(hash_int_t, old_next_kmerInt, h->hp.k_n, h->hp.hash_n);
                        remn_int_t old_next_kmerKey = last_n_bits(remn_int_t, old_next_kmerInt, _int_size, h->hp.remn_n);
                        kmer_num_t old_next_k_i;
                        if (!bin_hit_kmer(&old_next_k_i, h, old_next_hashKey, old_next_kmerKey)) err_fatal_core(__func__, "ERROR: Cannot hit in hash-table(2)!\nnext_nt: %d\n", old_next_nt);
                        set_bwt_char(h->kmer_node[old_next_hashKey]+old_next_k_i, nt_N, h);
                    }
                    next_spe_flag = 1; // set_bwt_char(new_next_kmer, #) ()
                    set_MoutEdge(h->kmer_node[hashKey]+k_i, h);
                    set_next_char(h->kmer_node[hashKey]+k_i, nt_N, h);
                    // uid, uoffset XXX
                }
            }
            // m-out => m-out
            else next_spe_flag = 1; // set_bwt_char(new_next_kmer, #)
        }
    } else {
        kmer_num_t k_n = ++(h->kmer_num[hashKey]);
        if (k_n > h->kmer_m[hashKey]) {
            h->kmer_m[hashKey] <<= 1;
            h->kmer_node[hashKey] = (kmer_node_t*)_err_realloc(h->kmer_node[hashKey], h->kmer_m[hashKey] * sizeof(kmer_node_t));
        }
        // copy node
        memmove(h->kmer_node[hashKey]+k_i+1, h->kmer_node[hashKey]+k_i, (k_n-1-k_i)*sizeof(kmer_node_t));
        kmer_node_t n = kmerKey;
        set_kmer_key(h->kmer_node[hashKey]+k_i, n, h);
        uint8_t bwt_nt = cur_spe_flag ? nt_N : pre_nt;
        set_bwt_char(h->kmer_node[hashKey]+k_i, bwt_nt, h);
        set_next_char(h->kmer_node[hashKey]+k_i, next_nt, h);
        // N-in
        if (bwt_nt == nt_N) {
            // uid, uoffset XXX
        }
        // N-out
        if (next_nt == nt_N) {
            cur_gen_spe_kmer(kmerInt, hashKey, k_i, h, db); // spe-kmer for cur_kmer, set spe-flag
        }
    }
    return next_spe_flag;
}

uint8_t hash_check(hash_idx *h, de_bwt_t *db, kmer_int_t kmerInt, uint8_t pre_nt, hash_int_t pre_hashKey, kmer_num_t pre_k_i, hash_int_t *cur_hashKey, kmer_num_t *cur_k_i, uint8_t next_nt, uint8_t cur_spe_flag)
{
    hash_int_t hashKey; remn_int_t kmerKey, k_i;
    uint8_t equal_f=0;

    hashKey = first_n_bits(hash_int_t, kmerInt, (h->hp.k)<<1, h->hp.hash_n);
    kmerKey = last_n_bits(remn_int_t, kmerInt, _int_size, h->hp.remn_n);

    k_i = bin_search_kmer(&equal_f, *h, hashKey, kmerKey);
#ifdef __DEBUG__
    fprintf(stderr, "k_i %d %c\n", k_i, "NY"[equal_f]);
    con_kmer(kmerInt, h->hp.k);
#endif
    
    *cur_hashKey = hashKey; *cur_k_i = k_i;
    return update_hash(k_i, equal_f, h, db, kmerInt, hashKey, kmerKey, pre_nt, pre_hashKey, pre_k_i, next_nt, cur_spe_flag);
}

// for static sort
uint8_t REMN_N;
uint8_t SK_NI;
uint8_t SK_M;
int skmer_comp(const void *sk1, const void *sk2) 
{ 
    remn_int_t kmerKey1 = first_n_bits(remn_int_t, (*(kmer_node_t*)sk1), _node_size, REMN_N);
    remn_int_t kmerKey2 = first_n_bits(remn_int_t, (*(kmer_node_t*)sk2), _node_size, REMN_N);
    uint8_t sk_n1 = ((*(kmer_node_t*)sk1) >> (SK_NI)) & SK_M;
    uint8_t sk_n2 = ((*(kmer_node_t*)sk2) >> (SK_NI)) & SK_M;
    if (kmerKey1 != kmerKey2) return kmerKey1 - kmerKey2;
    else return sk_n1 - sk_n2;
}

void sort_spe_kmer(hash_idx *h)
{
    REMN_N = h->hp.remn_n;
    SK_NI = h->hp.sk_ni;
    SK_M  = h->hp.sk_m;
    uint32_t i;
    for (i = 0; i < h->hp.hash_size; ++i) {
        if (h->skmer_num[i] > 1)
            qsort(h->skmer_node[i], h->skmer_num[i], sizeof(kmer_node_t), skmer_comp);
    }
}

int kmer_gen(const char *prefix, hash_idx *h_idx, de_bwt_t *db_idx)
{
    kmer_int_t kmerInt;
    gzFile fp; kseq_t *fa_seq;
    fp = xzopen(prefix, "r");
    fa_seq = kseq_init(fp);
    int seq_i; uint8_t pre_nt, next_nt, cur_spe_flag, next_spe_flag;
    hash_int_t pre_hashKey, cur_hashKey; kmer_num_t pre_k_i, cur_k_i;

    while (kseq_read(fa_seq) >= 0) {
        seq_i = 0;
        cur_spe_flag = 0;
        fprintf(stderr, "[%s] Generating and sorting kmer for %s ...\n",__func__, fa_seq->name.s);
        if (get_first_kmer(fa_seq, &seq_i, h_idx->hp.k, &kmerInt, &next_nt)) {        // first kmer
#ifdef __DEBUG__
                fprintf(stderr, "seq_i: %d\n", seq_i);
#endif
            next_spe_flag = hash_check(h_idx, db_idx, kmerInt, nt_N, pre_hashKey, pre_k_i, &cur_hashKey, &cur_k_i, next_nt, cur_spe_flag);
            cur_spe_flag = next_spe_flag;
            pre_hashKey = cur_hashKey; pre_k_i = cur_k_i;
            while (get_shift_kmer(fa_seq, &seq_i, h_idx->hp.k, h_idx->hp.k_m, &kmerInt, &pre_nt, &next_nt)) { // generate kmer based on last kmer
#ifdef __DEBUG__
                fprintf(stderr, "seq_i: %d\n", seq_i);
#endif
                next_spe_flag = hash_check(h_idx, db_idx, kmerInt, pre_nt, pre_hashKey, pre_k_i, &cur_hashKey, &cur_k_i, next_nt, cur_spe_flag);
                cur_spe_flag = next_spe_flag;
                pre_hashKey = cur_hashKey; pre_k_i = cur_k_i;
            }
        }
    }
    
    sort_spe_kmer(h_idx);

    fprintf(stderr, "[%s] Generating and sorting done!\n",__func__);
    // output kmer
#ifdef __DEBUG__
    uint64_t count = 0;
    uint32_t i, j;
    for (i = 0; i < h_idx->hp.hash_size; ++i) {
        count += h_idx->kmer_num[i];
        for (j = 0; j < h_idx->kmer_num[i]; ++j) {
            cons_kmer(h_idx->kmer_node[i][j], i, h_idx->hp.k, h_idx->hp.remn_n, h_idx->hp.remn_ni);
        }
    }
    fprintf(stderr, "[%s] Total normal kmer count: %lld\n",__func__, (long long)count);
#endif
    kseq_destroy(fa_seq); err_gzclose(fp);
    return 0;
}

int kmer_merge(hash_idx h_idx, de_bwt_t *de_idx)
{
    uint64_t i, j, k;
    remn_int_t j_key, k_key;
    for (i = 0; i < h_idx.hp.hash_size; ++i) {
        j = 0; k = 0;
        while (j < h_idx.kmer_num[i] && k < h_idx.skmer_num[i]) {
            j_key = get_kmer_key(h_idx.kmer_node[i][j], h_idx);
            k_key = get_kmer_key(h_idx.skmer_node[i][k], h_idx);
            if (j_key >= k_key) { // '#' is always smaller
                kputc(BWT_STR[get_bwt_char(h_idx.skmer_node[i][k], h_idx)], de_idx->bwt_str);
                k++;
            } else {
                kputc(BWT_STR[get_bwt_char(h_idx.kmer_node[i][j], h_idx)], de_idx->bwt_str);
                j++;
            }
        }
        while (j < h_idx.kmer_num[i]) {
            kputc(BWT_STR[get_bwt_char(h_idx.kmer_node[i][j], h_idx)], de_idx->bwt_str);
            j++;
        }
        while (k < h_idx.skmer_num[i]) {
            kputc(BWT_STR[get_bwt_char(h_idx.skmer_node[i][k], h_idx)], de_idx->bwt_str);
            k++;
        }
    }
#ifdef __DEBUG__
    for (i = 0; i < h_idx.hp.hash_size; ++i) {
        for (j = 0; j < h_idx.kmer_num[i]; ++j) {
            cons_kmer(h_idx.kmer_node[i][j], i, h_idx.hp.k, h_idx.hp.remn_n, h_idx.hp.remn_ni);
            fprintf(stdout, "\t\t%c\n", BWT_STR[get_bwt_char(h_idx.kmer_node[i][j], h_idx)]);
        }
    }
    fprintf(stdout, "\n");
    for (i = 0; i < h_idx.hp.hash_size; ++i) {
        for (j = 0; j < h_idx.skmer_num[i]; ++j) {
            cons_kmer(h_idx.skmer_node[i][j], i, h_idx.hp.k, h_idx.hp.remn_n, h_idx.hp.remn_ni);
            fprintf(stdout, "\t\t%c\n", BWT_STR[get_bwt_char(h_idx.skmer_node[i][j], h_idx)]);
        }
    }
#endif
    return 0;
}
