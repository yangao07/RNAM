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

void hash_init_num(hash_idx *h_idx)
{
    h_idx->kmer_num = (kmer_num_t*)_err_calloc(h_idx->hp.hash_size, sizeof(kmer_num_t));
    h_idx->kmer_c = (debwt_count_t*)_err_calloc(h_idx->hp.hash_size+1, sizeof(debwt_count_t));

    h_idx->skmer_num = (kmer_num_t*)_err_calloc(h_idx->hp.hash_size, sizeof(kmer_num_t));
    h_idx->skmer_c = (debwt_count_t*)_err_calloc(h_idx->hp.hash_size+1, sizeof(debwt_count_t));
}

void hash_init_node(hash_idx *h_idx)
{
    h_idx->kmer_node = (kmer_node_t*)_err_malloc(h_idx->kmer_tol_count * sizeof(kmer_node_t));
    h_idx->uni_offset_c = (debwt_count_t*)_err_calloc(1025, sizeof(debwt_count_t));
    h_idx->uni_node = (kmer_int_t*)_err_calloc(1024, sizeof(kmer_int_t));
    h_idx->uni_n = 0; h_idx->uni_m = 1024, h_idx->uni_t = 0;
}

void hash_init_snode(hash_idx *h_idx)
{
    h_idx->skmer_node = (skmer_node_t*)_err_malloc(h_idx->skmer_tol_count * sizeof(skmer_node_t));

    h_idx->cur_uid = 0; h_idx->next_uid = 0; // 0-base
    h_idx->next_offset_c = 0;
    h_idx->uni_offset_n = (debwt_count_t*)_err_calloc((h_idx->uni_tol_count),  sizeof(debwt_count_t));
    h_idx->uni_offset = (debwt_count_t*)_err_malloc(h_idx->off_tol_count * sizeof(debwt_count_t));
}

void hash_free(hash_idx *h_idx)
{
    free(h_idx->kmer_node); free(h_idx->kmer_num); free(h_idx->kmer_c);
    free(h_idx->skmer_node); free(h_idx->skmer_num); free(h_idx->skmer_c);

    free(h_idx->uni_offset_n); free(h_idx->uni_offset);
    free(h_idx->uni_offset_c); free(h_idx->uni_node);
}

void cons_kmer(kmer_node_t n, hash_int_t hashKey, hash_idx h)
{
    uint8_t k = h.hp.k, remn_n = h.hp.remn_n, remn_m = h.hp.remn_m, remn_ni = h.hp.remn_ni;
    char *seq = (char*)malloc(k+1);
    kmer_int_t all = hashKey;
    all = (all << remn_n) | ((n>>(remn_ni))&remn_m);
    int i;
    for (i = 0; i < k; ++i)
    {
        seq[i] = "ACGTN"[(all >> (2 * (k-i-1)))&0x3];
    }
    seq[k] = 0;
    stdout_printf("%s\n", seq);
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
    stdout_printf("%s %lld\n", seq, (long long)all);
    free(seq);
} 

int get_first_kmer(kseq_t *seq, ref_offset_t *seq_i, uint8_t k, kmer_int_t *kmer, uint8_t *next_nt)
{
    ref_offset_t k_i = 0, i = *seq_i;
    uint8_t nt;
    *kmer = 0;
    while (i < seq->seq.l) {
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

int get_first_kmer_hashKey(kseq_t *seq, int *seq_i, uint8_t k, uint8_t hash_k,  hash_int_t *hashKey, uint8_t *next_hashKey_nt, uint8_t *next_nt)
{
    int hash_k_i = 0, k_i = 0, i = *seq_i;
    uint8_t nt;
    *hashKey = 0;
    while (i < (int)seq->seq.l) {
        nt = nst_nt4_table[(int)seq->seq.s[i++]];
        if (nt != nt_N) {
            if (hash_k_i < hash_k) {
                (*hashKey) = (*hashKey) << 2 | nt;
                hash_k_i++;
            }
            k_i++;
        } else {
            k_i = 0; hash_k_i = 0; *hashKey = 0;
        }
        if (k_i == k) break;
    }
    *seq_i = i;
    *next_hashKey_nt = nst_nt4_table[(int)seq->seq.s[i-k+hash_k]];
    *next_nt = nst_nt4_table[(int)seq->seq.s[i]];
    return (k_i == k);
}

int get_shift_kmer(kseq_t *seq, ref_offset_t *seq_i, uint8_t k, kmer_int_t k_m, kmer_int_t *kmer, uint8_t *pre_nt, uint8_t *next_nt)
{
    ref_offset_t i = *seq_i;
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

int get_shift_kmer_hashKey(kseq_t *seq, int *seq_i, uint8_t k, uint8_t hash_k, kmer_int_t hash_k_m, hash_int_t *hashKey, uint8_t *next_hashKey_nt, uint8_t *next_nt)
{
    size_t i = *seq_i;
    while (i < seq->seq.l) {
        i++;
        if (*next_nt != nt_N) {
            (*hashKey) = ((*hashKey) << 2 | (*next_hashKey_nt)) & hash_k_m;
            if (i == seq->seq.l) *next_nt = nt_N;
            else {
                *next_nt = nst_nt4_table[(int)seq->seq.s[i]];
                *next_hashKey_nt = nst_nt4_table[(int)seq->seq.s[i-k+hash_k]];
            }
            *seq_i = i;
            return 1;
        } else {
            *seq_i = i;
            return get_first_kmer_hashKey(seq, seq_i, k, hash_k, hashKey, next_hashKey_nt, next_nt);
        }
    }
    return 0;
}

#define first_n_bits(type, x, l, n)  (type)((x) >> ((l)-(n)))
#define _first_n_bits(type, x, l, n) (type)(((x) >> ((l)-(n))) << ((l)-(n)))
#define last_n_bits(type, x, l, n)   (type)(((x) << ((l)-(n))) >> ((l)-(n)))
#define _last_n_bits(type, x, l, n)  (type)((x) << ((l)-(n)))

// this is ALWAYS the first operation
#define set_kmer_key(n, kmerKey, h) { (*(n)) = (kmerKey << ((h).hp.remn_ni)); }
#define set_bwt_char(n, last_nt, h) { (*(n)) |= (last_nt << ((h).hp.bwt_char_ni)); }    // [0:1,2]: [N:ACGT]
#define set_next_char(n, next_nt, h) { (*(n)) |= (next_nt << ((h).hp.next_char_ni)); }
#define set_sk_n(n, sk_n, h) { (*(n)) |= ((sk_n) << ((h).hp.sk_ni)); }
#define set_spe(n, h) {(*(n)) |= (1 << ((h).hp.spe_ni)); }
#define set_un_spe(n, h) {(*(n)) -= (1 << ((h).hp.spe_ni)); }

#define get_kmer_key(n, h)  (((n) >> ((h).hp.remn_ni)) & ((h).hp.remn_m))
#define get_bwt_char(n, h) (((n) >> ((h).hp.bwt_char_ni)) & ((h).hp.char_m))
#define get_next_char(n, h) (((n) >> ((h).hp.next_char_ni)) & ((h).hp.char_m))
#define get_sk_n(n, h) (((n) >> ((h).hp.sk_ni)) & ((h).hp.sk_m))
#define get_spe(n, h) (((n) >> ((h).hp.spe_ni)) & ((h).hp.spe_m))

#define set_uni_off_flag(n, h) { (*(n)) |= (1 << ((h).hp.uni_off_flag_ni)); }
#define get_uni_off_flag(n, h) (((n) >> ((h).hp.uni_off_flag_ni)) & ((h).hp.uni_off_m))

//#define set_kmer_uid(n, uid, h) { (*(n)) |= ((uid) << ((h)->hp.uid_ni)); }
//#define get_kmer_uid(n, h) (((n) >> ((h).hp.uid_ni)) & ((h).hp.uid_m))

kmer_num_t bin_search_kmer(uint8_t *ef, hash_idx h, hash_int_t hashKey, remn_int_t kmerKey)
{
    kmer_node_t* node = h.kmer_node + h.kmer_c[hashKey];
    remn_int_t tmpKey;
    int left=0, right = h.kmer_num[hashKey]-1, middle;
    if (right == -1) { *ef = 0; return 0; }

    if (kmerKey < first_n_bits(remn_int_t, node[left], _node_size, h.hp.remn_n)) { *ef = 0; return left; }
    else if (kmerKey > first_n_bits(remn_int_t, node[right], _node_size, h.hp.remn_n)) { *ef = 0; return right+1; }

    // binary search for hashInt(hashKey, kmerKey)
    while (left <= right) {
        middle = (left+right)>>1;
        tmpKey = first_n_bits(remn_int_t, node[middle], _node_size, h.hp.remn_n);
        if (tmpKey == kmerKey) { *ef = 1; return middle; }
        else if (tmpKey > kmerKey) {
            if ((middle == 0) || (first_n_bits(remn_int_t, node[middle-1], _node_size, h.hp.remn_n) < kmerKey)) {
                *ef = 0; return middle;
            } else right = middle-1;
        } else left = middle+1;
    }
    *ef = 0;
    return h.kmer_num[hashKey];
}

debwt_count_t bin_search_kInt(kmer_int_t *node, debwt_count_t n, kmer_int_t kInt)
{
    if (n == 0) return 0;
    kmer_int_t tmp_Int;
    debwt_count_t left = 0, right = n-1, middle;
    if (kInt < node[left]) return left;
    if (kInt > node[right]) return right+1;
   
    while (left <= right) {
        middle = (left+right) >> 1; 
        tmp_Int = node[middle];
        if (tmp_Int > kInt) {
            if ((middle == 0) || node[middle-1] < kInt) return middle;
            else right = middle-1;
        } else left = middle+1;
    }
    return n;
}

int bin_hit_kmer(kmer_num_t *ki, hash_idx *h, hash_int_t hashKey, remn_int_t kmerKey)
{
    remn_int_t tmpKey;
    int left=0, right = h->kmer_num[hashKey]-1, middle;
    if (right == -1) return 0;

    // binary search for hashInt(hashKey, kmerKey)
    debwt_count_t k_ns = h->kmer_c[hashKey];
    while (left <= right) {
        middle = (left+right)>>1;
        tmpKey = first_n_bits(remn_int_t, h->kmer_node[k_ns+middle], _node_size, h->hp.remn_n);
        if (tmpKey == kmerKey) {
            *ki = middle; return 1;
        } else if (tmpKey > kmerKey)  right = middle-1;
        else left = middle+1;
    }
    return 0;
}

int bin_hit_kInt(kmer_int_t *node, debwt_count_t n, debwt_count_t *k_i, kmer_int_t kInt)
{
    if (n == 0) return 0;
    kmer_int_t tmp_Int;
    debwt_count_t left = 0, right = n-1, middle;
    while (left <= right) {
        middle = (left+right) >> 1; 
        tmp_Int = node[middle];
        if (tmp_Int == kInt) {
           *k_i = middle; return 1;
        } else if (tmp_Int > kInt) right = middle-1;
        else left = middle+1;
    }
    return 0;
}

int update_spe_hash(hash_idx *h, hash_int_t hashKey, remn_int_t kmerKey, uint8_t pre_nt, uint8_t sk_n)
{
    kmer_num_t k_i = h->skmer_num[hashKey];
    skmer_node_t t = kmerKey;
    set_kmer_key(h->skmer_node+h->skmer_c[hashKey]+k_i, t, *h);
    set_bwt_char(h->skmer_node+h->skmer_c[hashKey]+k_i, pre_nt, *h);
    set_sk_n(h->skmer_node+h->skmer_c[hashKey]+k_i, sk_n, *h);
    // set_kmer_uid for spe-kmer
    //skmer_node_t _uid = uid;
    //set_kmer_uid(h->skmer_node+h->skmer_c[hashKey]+k_i, _uid, h);
    ++(h->skmer_num[hashKey]);
    return 1;
}

void pre_count_spe_kmer(uint8_t pre_nt, hash_int_t pre_hashKey, kmer_num_t pre_k_i, kmer_int_t kmerInt, hash_idx *h) 
{
    if (pre_nt >= nt_N || get_spe(h->kmer_node[h->kmer_c[pre_hashKey]+pre_k_i], *h)) return;

    kmer_int_t s_kmerInt = kmerInt - (kmerInt & 0x3);
    hash_int_t s_hashKey;

    int i;
    for (i = 0; i < h->hp.k-1; ++i) {
        s_hashKey = first_n_bits(hash_int_t, s_kmerInt, h->hp.k_n, h->hp.hash_n);
        h->skmer_c[s_hashKey+1]++;
        s_kmerInt = (s_kmerInt<<2) & h->hp.k_m;
    }

    set_spe(h->kmer_node+h->kmer_c[pre_hashKey]+pre_k_i, *h);
}

void cur_count_spe_kmer(kmer_int_t kmerInt, hash_int_t hashKey, kmer_num_t k_i, hash_idx *h) 
{
    if (get_spe(h->kmer_node[h->kmer_c[hashKey]+k_i], *h)) return;

    kmer_int_t s_kmerInt = kmerInt;
    hash_int_t s_hashKey;
    int i;
    for (i = 0; i < h->hp.k-1; ++i) {
        s_kmerInt = (s_kmerInt << 2) & h->hp.k_m;
        s_hashKey = first_n_bits(hash_int_t, s_kmerInt, h->hp.k_n, h->hp.hash_n);
        h->skmer_c[s_hashKey+1]++;
    }
    set_spe(h->kmer_node+h->kmer_c[hashKey]+k_i, *h);
}

// sort by kInt
void update_uni_offc_by_int(hash_idx *h, kmer_int_t kInt, int off_c, int ef)
{
    debwt_count_t k_i;
    if (ef) { // equal
        // get k_i
        if (!bin_hit_kInt(h->uni_node, h->uni_n, &k_i, kInt)) err_fatal_simple("ERROR: Cannot hit in uni_node by kInt!\n");
        h->uni_offset_c[k_i+1] += off_c;
    } else {  // not equal
        k_i = bin_search_kInt(h->uni_node, h->uni_n, kInt);
        h->uni_n++;
        if (h->uni_n > h->uni_m) {
            // realloc
            h->uni_m <<= 1;
            h->uni_node = (kmer_int_t*)_err_realloc(h->uni_node, h->uni_m * sizeof(kmer_int_t));
            h->uni_offset_c = (debwt_count_t*)_err_realloc(h->uni_offset_c, (h->uni_m+1) * sizeof(debwt_count_t));
        }
        // memmove node and offset_c
        memmove(h->uni_node+k_i+1, h->uni_node+k_i, (h->uni_n-1-k_i)*sizeof(kmer_int_t));
        memmove(h->uni_offset_c+k_i+2, h->uni_offset_c+k_i+1, (h->uni_n-1-k_i)*sizeof(debwt_count_t));
        // add off_c|node 
        h->uni_node[k_i] = kInt;
        h->uni_offset_c[k_i+1] = off_c;
    }
}

// sort by uid 
void update_uni_offc_by_uid(hash_idx *h, kmer_int_t kInt)
{
    debwt_count_t k_i;
    if (!bin_hit_kInt(h->uni_node+h->uni_t, h->uni_n-h->uni_t, &k_i, kInt)) err_fatal_simple("ERROR: Cannot hit in uni_node by uid!\n");
    debwt_count_t tmp_c = h->uni_offset_c[h->uni_t+k_i+1];
    // memmove node and off_c
    memmove(h->uni_node+h->uni_t+1, h->uni_node+h->uni_t, k_i * sizeof(kmer_int_t));
    memmove(h->uni_offset_c+h->uni_t+2, h->uni_offset_c+h->uni_t+1, k_i * sizeof(debwt_count_t));
    // update 
    h->uni_node[h->uni_t] = kInt;
    h->uni_offset_c[h->uni_t+1] = tmp_c + h->uni_offset_c[h->uni_t];
    h->uni_t++;
}

debwt_count_t get_uid_from_off_c(hash_idx h, kmer_int_t kInt)
{
    debwt_count_t uid;
    for (uid = 0; uid < h.uni_t; ++uid)
        if (h.uni_node[uid] == kInt) return uid;
    err_fatal_simple("ERROR: Cannot hit in uni_node by kInt!\n");
    return 0;
}

uint8_t update_hash(kmer_num_t k_i, uint8_t ef, hash_idx *h,
                    kmer_int_t kInt, hash_int_t hKey, remn_int_t kKey, 
                    uint8_t p_nt, hash_int_t p_hKey, kmer_num_t p_ki, 
                    uint8_t n_nt, uint8_t c_sf)
{
    debwt_count_t k_ns = h->kmer_c[hKey];
    uint8_t next_spe_flag = 0;
    if (ef) {
        { // in
            // NEW
            uint8_t old_pre_nt  = get_bwt_char(h->kmer_node[k_ns+k_i], *h); 
            if (p_nt != old_pre_nt || p_nt == nt_N) {
                // new_pre_kmer => spe_kmer, spe_flag
                pre_count_spe_kmer(p_nt, p_hKey, p_ki, kInt, h);
                // old_pre_kmer => spe_kmer, spe_flag
                if (old_pre_nt < nt_N) {
                    kmer_int_t old_pre_kmerInt = old_pre_nt; 
                    old_pre_kmerInt = (old_pre_kmerInt << (h->hp.k_n-2)) | (kInt >> 2);

                    hash_int_t old_pre_hashKey = first_n_bits(hash_int_t, old_pre_kmerInt, h->hp.k_n, h->hp.hash_n);
                    remn_int_t old_pre_kmerKey = last_n_bits(remn_int_t, old_pre_kmerInt, _int_size, h->hp.remn_n);
                    kmer_num_t old_pre_k_i;
                    if (!bin_hit_kmer(&old_pre_k_i, h, old_pre_hashKey, old_pre_kmerKey)) err_fatal_core(__func__, "ERROR: Cannot hit in hash-table(1)!\npre_nt: %d\n", old_pre_nt);
                    pre_count_spe_kmer(old_pre_nt, old_pre_hashKey, old_pre_k_i, kInt, h);
                    // set bwt_char #
                    set_bwt_char(h->kmer_node+k_ns+k_i, nt_N, *h);
                    // uni_offset_c +2
                    update_uni_offc_by_int(h, kInt, 2, 0);
                } else {
                    // uni_offset_c +1
                    update_uni_offc_by_int(h, kInt, 1, 1);
                }
            }
        }
        { // out
            // NEW
            uint8_t old_next_nt = get_next_char(h->kmer_node[k_ns+k_i], *h);
            if (n_nt != old_next_nt || n_nt == nt_N) {
                // cur_kmer => spe_kmer, spe_flag
                cur_count_spe_kmer(kInt, hKey, k_i, h);
                if (old_next_nt < nt_N) {
                    // set old_next_kmer bwt_char #
                    kmer_int_t old_next_kmerInt = ((kInt << 2) & h->hp.k_m) | old_next_nt;
                    hash_int_t old_next_hashKey = first_n_bits(hash_int_t, old_next_kmerInt, h->hp.k_n, h->hp.hash_n);
                    remn_int_t old_next_kmerKey = last_n_bits(remn_int_t, old_next_kmerInt, _int_size, h->hp.remn_n);
                    kmer_num_t old_next_k_i;
                    if (!bin_hit_kmer(&old_next_k_i, h, old_next_hashKey, old_next_kmerKey)) err_fatal_core(__func__, "ERROR: Cannot hit in hash-table(2)!\nnext_nt: %d\n", old_next_nt);
                    if (get_bwt_char(h->kmer_node[h->kmer_c[old_next_hashKey]+old_next_k_i], *h) < nt_N) {
                        // set_bwt_char(old_next_kmer, #)
                        set_bwt_char(h->kmer_node+h->kmer_c[old_next_hashKey]+old_next_k_i, nt_N, *h);
                        // uni_offset_c for old_next_kmer
                        update_uni_offc_by_int(h, old_next_kmerInt, 1, 0); 
                    }
                    // set next_char
                    set_next_char(h->kmer_node+k_ns+k_i, nt_N, *h);
                }
                // new_next_kmer
                next_spe_flag = 1;
            }
        }
    } else {
        kmer_num_t k_n = ++(h->kmer_num[hKey]);
        if (k_n > h->kmer_c[hKey+1]-h->kmer_c[hKey]) err_fatal_simple("Node num > node_m\n");
        // copy node
        memmove(h->kmer_node+k_ns+k_i+1, h->kmer_node+k_ns+k_i, (k_n-1-k_i)*sizeof(kmer_node_t));
        // set node
        kmer_node_t n = kKey;
        set_kmer_key(h->kmer_node+k_ns+k_i, n, *h);
        uint8_t bwt_nt = c_sf ? nt_N : p_nt;
        set_bwt_char(h->kmer_node+k_ns+k_i, bwt_nt, *h);
        set_next_char(h->kmer_node+k_ns+k_i, n_nt, *h);
        // N-in
        if (bwt_nt == nt_N) {
            // uni_offset_c
            update_uni_offc_by_int(h, kInt, 1, 0);
        }
        // N-out
        if (n_nt == nt_N) {
            cur_count_spe_kmer(kInt, hKey, k_i, h); // spe-kmer for cur_kmer, set spe-flag
        }
    }
    return next_spe_flag;
}

int gen_spe_hash(kmer_num_t k_i, hash_idx *h, de_bwt_t *de_idx, ref_offset_t ref_off, 
                 hash_int_t hashKey, remn_int_t kmerKey, 
                 uint8_t p_nt, hash_int_t pre_hashKey, kmer_num_t pre_k_i)
{
    debwt_count_t k_ns = h->kmer_c[hashKey];
    kmer_node_t k_node = h->kmer_node[k_ns+k_i];
    kmer_int_t kmerInt = hashKey; kmerInt <<= h->hp.remn_n; kmerInt |= kmerKey;
    kmer_int_t s_kmerInt; hash_int_t s_hashKey; remn_int_t s_kmerKey;

    uint8_t bwt_nt; int i;

#ifdef __DEBUG__
    cons_kmer(k_node, hashKey, *h);
    con_kmer(kmerInt, h->hp.k);
#endif
    // uid, uni_offset XXX
    if (get_bwt_char(k_node, *h) >= nt_N) { // new unipath , update (cur_uid and last_uid), U_offset, U_num
        h->last_uid = h->cur_uid;
        if (get_uni_off_flag(k_node, *h)) { // offset_c alread updated
            h->cur_uid = get_uid_from_off_c(*h, kmerInt);
        } else { // need to update uid and offset_c
            h->cur_uid = h->next_uid++;
            update_uni_offc_by_uid(h, kmerInt);
            set_uni_off_flag(h->kmer_node+k_ns+k_i, *h);
        }
        // update u_offset
        h->uni_offset[h->uni_offset_c[h->cur_uid]+h->uni_offset_n[h->cur_uid]++] = ref_off;
        err_printf("ref: %d\n", ref_off);
    } 
                                                     // has been set_spe in the pre run
    if (get_bwt_char(k_node, *h) >= nt_N && p_nt != nt_N && get_spe(h->kmer_node[h->kmer_c[pre_hashKey]+pre_k_i], *h)) {
        set_un_spe(h->kmer_node+h->kmer_c[pre_hashKey]+pre_k_i, *h);
        // last_kmer -> spe_kmer 
        bwt_nt = p_nt;
        s_kmerInt = kmerInt - (kmerInt & 0x3);
        for (i = 0; i < h->hp.k-1; ++i) {
            s_hashKey = first_n_bits(hash_int_t, s_kmerInt, (h->hp.k)<<1, h->hp.hash_n);
            s_kmerKey = last_n_bits(remn_int_t, s_kmerInt, _int_size, h->hp.remn_n);
            if (!update_spe_hash(h, s_hashKey, s_kmerKey, bwt_nt, h->hp.k-1-i)) err_fatal_simple("ERROR: hash-update fail(1)!\n");

            bwt_nt = first_n_bits(uint8_t, s_kmerInt, (h->hp.k)<<1, 2);
            s_kmerInt = (s_kmerInt<<2) & h->hp.k_m;
        }
        // generate bwt_char for ^#-kmer
        kputc(BWT_STR[bwt_nt], de_idx->bwt_str);
    } 
    if (get_next_char(k_node, *h) >= nt_N && get_spe(h->kmer_node[k_ns+k_i], *h)) {
        set_un_spe(h->kmer_node+k_ns+k_i, *h);
        // cur_kmer -> unipath
        s_kmerInt = kmerInt;
        for (i = 0; i < h->hp.k-1; ++i) {
            bwt_nt = first_n_bits(uint8_t, s_kmerInt, (h->hp.k)<<1, 2);
            s_kmerInt = (s_kmerInt << 2) & h->hp.k_m;
            s_hashKey = first_n_bits(hash_int_t, s_kmerInt, (h->hp.k)<<1, h->hp.hash_n);
            s_kmerKey = last_n_bits(remn_int_t, s_kmerInt, _int_size, h->hp.remn_n);
            if (!update_spe_hash(h, s_hashKey, s_kmerKey, bwt_nt, h->hp.k-1-i)) err_fatal_simple("ERROR: hash-update fail(2)!\n");
        }
        // generate bwt_char for ^#-kmer
        kputc(BWT_STR[first_n_bits(uint8_t, s_kmerInt, (h->hp.k)<<1, 2)], de_idx->bwt_str);
    }
    return 1;    
}

uint8_t hash_check(hash_idx *h, kmer_int_t kInt, 
                   uint8_t p_nt, hash_int_t p_hK, kmer_num_t p_ki, 
                   hash_int_t *c_hK, kmer_num_t *c_ki, uint8_t n_nt, uint8_t c_sf)
{
    hash_int_t hashKey; remn_int_t kmerKey, k_i;
    uint8_t equal_f=0;

    hashKey = first_n_bits(hash_int_t, kInt, h->hp.k_n, h->hp.hash_n);
    kmerKey = last_n_bits(remn_int_t, kInt, _int_size, h->hp.remn_n);

    k_i = bin_search_kmer(&equal_f, *h, hashKey, kmerKey);
#ifdef __DEBUG__
    stdout_printf("k_i %d %c\n", k_i, "NY"[equal_f]);
    con_kmer(kInt, h->hp.k);
#endif
    
    *c_hK = hashKey; *c_ki = k_i;
    return update_hash(k_i, equal_f, h, kInt, hashKey, kmerKey, p_nt, p_hK, p_ki, n_nt, c_sf);
}

kmer_num_t hash_spe_check(hash_idx *h, de_bwt_t *de_idx, ref_offset_t ref_off, kmer_int_t kmerInt, uint8_t p_nt, hash_int_t pre_hashKey, kmer_num_t pre_k_i, hash_int_t *hashKey)
{
    remn_int_t kmerKey;
    kmer_num_t k_i;

    *hashKey = first_n_bits(hash_int_t, kmerInt, (h->hp.k)<<1, h->hp.hash_n);
    kmerKey = last_n_bits(remn_int_t, kmerInt, _int_size, h->hp.remn_n);

    if (!bin_hit_kmer(&k_i, h, *hashKey, kmerKey)) err_fatal_simple("ERROR: Cannot hit in hash-table!\n");
    if (!gen_spe_hash(k_i, h, de_idx, ref_off, *hashKey, kmerKey, p_nt, pre_hashKey, pre_k_i)) err_fatal_simple("ERROR: hash-update fail!\n");
    return k_i;
}

uint8_t REMN_N; uint64_t REMN_M; uint8_t SK_NI, SK_M;
int skmer_comp(const void *sk1, const void *sk2) 
{ 
    remn_int_t kmerKey1 = (first_n_bits(remn_int_t, (*(skmer_node_t*)sk1), _node_size, REMN_N)) & REMN_M;
    remn_int_t kmerKey2 = (first_n_bits(remn_int_t, (*(skmer_node_t*)sk2), _node_size, REMN_N)) & REMN_M;
    uint8_t sk_n1 = ((*(skmer_node_t*)sk1) >> (SK_NI)) & SK_M;
    uint8_t sk_n2 = ((*(skmer_node_t*)sk2) >> (SK_NI)) & SK_M;
    if (kmerKey1 != kmerKey2) return kmerKey1 - kmerKey2;
    else return sk_n1 - sk_n2;
}

void sort_spe_kmer(hash_idx *h)
{
    REMN_N = h->hp.remn_n;
    REMN_M = h->hp.remn_m;
    SK_NI = h->hp.sk_ni;
    SK_M  = h->hp.sk_m;
    uint32_t i;
    for (i = 0; i < h->hp.hash_size; ++i) {
        if (h->skmer_num[i] > 1) {
            qsort(h->skmer_node+h->skmer_c[i], h->skmer_num[i], sizeof(skmer_node_t), skmer_comp);
        }
    }
}

/* First-run :                       *
 * 1. count the total number of kmer */
int kmer_tol_count(const char *prefix, hash_idx *h_idx)
{
    gzFile fp; kseq_t *fa_seq;
    fp = xzopen(prefix, "r"); fa_seq = kseq_init(fp);
    hash_int_t hashKey; int seq_i; uint8_t next_hashKey_nt, next_nt;

    while (kseq_read(fa_seq) >= 0) {
        seq_i = 0;
        err_printf("[%s] Counting kmer for %s ...\n",__func__, fa_seq->name.s);
        if (get_first_kmer_hashKey(fa_seq, &seq_i, h_idx->hp.k, h_idx->hp.hash_k, &hashKey, &next_hashKey_nt, &next_nt)) {
            h_idx->kmer_c[hashKey+1]++; // kmer_c[0] = 0;
            while (get_shift_kmer_hashKey(fa_seq, &seq_i, h_idx->hp.k, h_idx->hp.hash_k, h_idx->hp.hash_m, &hashKey, &next_hashKey_nt, &next_nt))
                h_idx->kmer_c[hashKey+1]++; // kmer_c[0] = 0;
        }
    }
    debwt_count_t count=0; uint32_t i;
    for (i = 1; i < h_idx->hp.hash_size+1; ++i) {
        count += h_idx->kmer_c[i];
        h_idx->kmer_c[i] = count;
    }
    err_printf("[%s] Total kmer count: %lld\n",__func__, (long long)count);
    h_idx->kmer_tol_count = count;
    hash_init_node(h_idx);
    kseq_destroy(fa_seq); err_gzclose(fp);
    return 0;
}

/* Second-run :                                      *
 * 1. generat and sort normal kmer                   *
 * 2. count the total number of special kmer         *
 * 3. count the total number of unipath and u_offset */
int kmer_gen(const char *prefix, hash_idx *h_idx)
{
    gzFile fp; kseq_t *fa_seq;
    fp = xzopen(prefix, "r"); fa_seq = kseq_init(fp);

    kmer_int_t kmerInt; ref_offset_t seq_i; uint8_t pre_nt, next_nt, cur_spe_flag, next_spe_flag;
    hash_int_t pre_hashKey, cur_hashKey; kmer_num_t pre_k_i, cur_k_i;

    while (kseq_read(fa_seq) >= 0) {
        seq_i = 0;
        cur_spe_flag = 0;
        err_printf("[%s] Generating and sorting kmer for %s ...\n",__func__, fa_seq->name.s);
        if (get_first_kmer(fa_seq, &seq_i, h_idx->hp.k, &kmerInt, &next_nt)) {        // first kmer
#ifdef __DEBUG__
                err_printf("seq_i: %d\n", seq_i);
#endif
            next_spe_flag = hash_check(h_idx, kmerInt, nt_N, pre_hashKey, pre_k_i, &cur_hashKey, &cur_k_i, next_nt, cur_spe_flag);
            cur_spe_flag = next_spe_flag;
            pre_hashKey = cur_hashKey; pre_k_i = cur_k_i;
            while (get_shift_kmer(fa_seq, &seq_i, h_idx->hp.k, h_idx->hp.k_m, &kmerInt, &pre_nt, &next_nt)) { // generate kmer based on last kmer
#ifdef __DEBUG__
                err_printf("seq_i: %d\n", seq_i);
#endif
                next_spe_flag = hash_check(h_idx, kmerInt, pre_nt, pre_hashKey, pre_k_i, &cur_hashKey, &cur_k_i, next_nt, cur_spe_flag);
                cur_spe_flag = next_spe_flag;
                pre_hashKey = cur_hashKey; pre_k_i = cur_k_i;
            }
        }
    }
    
    err_printf("[%s] Generating and sorting done!\n",__func__);
    debwt_count_t k_count=0, uo_count=0, sk_count=0, k_ns;
    uint32_t i, j;
    for (i = 0; i < h_idx->hp.hash_size; ++i) {
        k_count += h_idx->kmer_num[i];
        k_ns = h_idx->kmer_c[i];
        for (j = 0; j < h_idx->kmer_num[i]; ++j) {
#ifdef __DEBUG__
            cons_kmer(h_idx->kmer_node[k_ns+j], i, *h_idx); // output kmer
#endif
        }
    }
    h_idx->kmer_real_count = k_count; // real count of normal kmer. kmer_tol_count is the collision-free count number

    for (i = 1; i < h_idx->uni_n+1; ++i) uo_count += h_idx->uni_offset_c[i];
    h_idx->off_tol_count = uo_count;

    err_printf("[%s] Total normal kmer real count: %lld\n",__func__, (long long)k_count);
    err_printf("[%s] Total unipath offset count: %lld\n",__func__, (long long)uo_count);

    for (i = 1; i < h_idx->hp.hash_size+1; ++i) {
        sk_count += h_idx->skmer_c[i];
        h_idx->skmer_c[i] = sk_count;
    }
    h_idx->skmer_tol_count = sk_count; 
    h_idx->uni_tol_count = sk_count / (h_idx->hp.k-1);

    err_printf("[%s] Total unipath count: %lld\n",__func__, (long long)h_idx->uni_tol_count);
    hash_init_snode(h_idx); 
    kseq_destroy(fa_seq); err_gzclose(fp);
    return 0;
}

/* Thrid-run :                     *
 * 1. generate special normal kmer *
 * 2. generate UID                 *
 * 3. generate u_offset_c          *
 * 4. generate u_offset            */
int spe_kmer_gen(const char *prefix, hash_idx *h_idx, de_bwt_t *db_idx)
{
    kmer_int_t kmerInt;
    gzFile fp; kseq_t *fa_seq;
    fp = xzopen(prefix, "r");
    fa_seq = kseq_init(fp);
    ref_offset_t chr_off, chr_len, seq_i; uint8_t pre_nt, next_nt;
    kmer_num_t k_i, pre_k_i;
    hash_int_t hashKey, pre_hashKey;

    chr_off = 0;
    while (kseq_read(fa_seq) >= 0) {
        chr_len = fa_seq->seq.l;
        seq_i = 0;
        err_printf("[%s] Re-generating and sorting special kmer for %s ...\n",__func__, fa_seq->name.s);
        pre_k_i = -1; pre_hashKey = 0;
        if (get_first_kmer(fa_seq, &seq_i, h_idx->hp.k, &kmerInt, &next_nt)) {        // first kmer
            k_i = hash_spe_check(h_idx, db_idx, seq_i+chr_off+1-h_idx->hp.k, kmerInt, nt_N, pre_hashKey, pre_k_i, &hashKey);
            pre_k_i = k_i; pre_hashKey = hashKey;
            while (get_shift_kmer(fa_seq, &seq_i, h_idx->hp.k, h_idx->hp.k_m, &kmerInt, &pre_nt, &next_nt)) { // generate kmer based on last kmer
                k_i = hash_spe_check(h_idx, db_idx, seq_i+chr_off+1-h_idx->hp.k, kmerInt, pre_nt, pre_hashKey, pre_k_i, &hashKey);
                pre_k_i = k_i; pre_hashKey = hashKey;
            }
        }
        chr_off += chr_len;
    }
 
    sort_spe_kmer(h_idx);
    err_printf("[%s] Re-generating and sorting done!\n",__func__);

    debwt_count_t sk_ns;
    uint32_t i, j;
    for (i = 0; i < h_idx->hp.hash_size; ++i) {
        sk_ns = h_idx->skmer_c[i];
        for (j = 0; j < h_idx->skmer_num[i]; ++j) {
#ifdef __DEBUG__
            cons_kmer(h_idx->skmer_node[sk_ns+j], i, *h_idx); // output kmer
#endif
        }
    }
    err_printf("[%s] Total special kmer count: %lld\n",__func__, (long long)h_idx->skmer_tol_count);

    kseq_destroy(fa_seq); err_gzclose(fp);
    return 0;
}

/* Final func of kmer_process *
 * 1. bwt_str[]               *
 * 2. uid[] XXX               *
 * 3. u_offset[]              */
int kmer_merge(hash_idx h_idx, de_bwt_t *de_idx)
{
    debwt_count_t i, j, k, k_ns, sk_ns;
    remn_int_t j_key, k_key;
    //debwt_count_t uid_i = 0;
    for (i = 0; i < h_idx.hp.hash_size; ++i) {
        k_ns = h_idx.kmer_c[i], sk_ns = h_idx.skmer_c[i];
        j = 0; k = 0;
        while (j < h_idx.kmer_num[i] && k < h_idx.skmer_num[i]) {
            j_key = get_kmer_key(h_idx.kmer_node[k_ns+j], h_idx);
            k_key = get_kmer_key(h_idx.skmer_node[sk_ns+k], h_idx);
            if (j_key >= k_key) { // '#' is always smaller
                kputc(BWT_STR[get_bwt_char(h_idx.skmer_node[sk_ns+k], h_idx)], de_idx->bwt_str);
                k++;
            } else {
                kputc(BWT_STR[get_bwt_char(h_idx.kmer_node[k_ns+j], h_idx)], de_idx->bwt_str);
                j++;
            }
        }
        while (j < h_idx.kmer_num[i]) {
            kputc(BWT_STR[get_bwt_char(h_idx.kmer_node[k_ns+j], h_idx)], de_idx->bwt_str);
            j++;
        }
        while (k < h_idx.skmer_num[i]) {
            kputc(BWT_STR[get_bwt_char(h_idx.skmer_node[sk_ns+k], h_idx)], de_idx->bwt_str);
            k++;
        }
    }
    // uni_offset_c, uni_offset
    for (i = 0; i <= de_idx->n_unipath; ++i)  de_idx->uni_offset_c[i] = h_idx.uni_offset_c[i];
    for (i = 0; i < de_idx->n_offset; ++i) de_idx->uni_offset[i] = h_idx.uni_offset[i];
#ifdef __DEBUG__
    for (i = 0; i < h_idx.hp.hash_size; ++i) {
        k_ns = h_idx.kmer_c[i];
        for (j = 0; j < h_idx.kmer_num[i]; ++j) {
            cons_kmer(h_idx.kmer_node[k_ns+j], i, h_idx);
            stdout_printf("\t\t%c\n", BWT_STR[get_bwt_char(h_idx.kmer_node[k_ns+j], h_idx)]);
        }
    }
    stdout_printf("\n");
    for (i = 0; i < h_idx.hp.hash_size; ++i) {
        sk_ns = h_idx.skmer_c[i];
        for (j = 0; j < h_idx.skmer_num[i]; ++j) {
            cons_kmer(h_idx.skmer_node[sk_ns+j], i, h_idx);
            stdout_printf("\t\t%c\n", BWT_STR[get_bwt_char(h_idx.skmer_node[sk_ns+j], h_idx)]);
        }
    }
#endif
    return 0;
}
