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

KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];

void hash_init(hash_idx *h_idx)
{
    h_idx->kmer_num = (kmer_num_t*)_err_calloc(h_idx->hp.hash_size, sizeof(kmer_num_t));
    h_idx->kmer_m = (kmer_num_t*)_err_malloc(h_idx->hp.hash_size * sizeof(kmer_num_t));
    h_idx->kmer_node = (kmer_node_t**)_err_calloc(h_idx->hp.hash_size, sizeof(kmer_node_t*));
    uint32_t i;
    for (i = 0; i < h_idx->hp.hash_size; ++i) {
        h_idx->kmer_m[i] = 8;
        h_idx->kmer_node[i] = (kmer_node_t*)_err_calloc(8, sizeof(kmer_node_t));
    }
}

void hash_free(hash_idx *h_idx)
{
    uint32_t i;
    for (i = 0; i < h_idx->hp.hash_size; ++i) {
        free(h_idx->kmer_node[i]);
    }
    free(h_idx->kmer_node); free(h_idx->kmer_num);
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
    int i = *seq_i;
    while (i < (int)seq->seq.l) {
        i++;
        if (*next_nt != nt_N) {
            *pre_nt = (*kmer) >> ((k-1)<<1);
            (*kmer) = ((*kmer) << 2 | (*next_nt)) & k_m; // 44 
            *next_nt = nst_nt4_table[(int)seq->seq.s[i]];
            *seq_i = i;
            return 1;
        } else {
            *seq_i = i;
            if (!get_first_kmer(seq, seq_i, k, kmer, next_nt)) return 0;
            *pre_nt = nt_N;
            return 1;
        }
    }
    // XXX printf("%c\n", seq->seq.s[i]);
    return 0;
}

#define first_n_bits(type, x, l, n)  (type)(x >> ((l)-(n)))
#define _first_n_bits(type, x, l, n) (type)((x >> ((l)-(n))) << ((l)-(n)))
#define last_n_bits(type, x, l, n)   (type)((x << ((l)-(n))) >> ((l)-(n)))
#define _last_n_bits(type, x, l, n)  (type)(x << ((l)-(n)))

uint32_t bin_search_kmer(uint8_t *ef, hash_idx h, kmer_num_t hashKey, uint32_t kmerKey)
{
    uint32_t tmpKey;
    int left=0, right = h.kmer_num[hashKey]-1, middle;
    if (right == -1) { *ef = 0; return 0; }

    // binary search for hashInt(hashKey, kmerKey)
    while (left <= right) {
        middle = (left+right)>>1;
        tmpKey = first_n_bits(uint32_t, h.kmer_node[hashKey][middle], _node_size, h.hp.remn_n);
        if (tmpKey == kmerKey) { *ef = 1; return middle; }
        else if (tmpKey > kmerKey) {
            if ((middle == 0) || (first_n_bits(uint32_t, h.kmer_node[hashKey][middle-1], _node_size, h.hp.remn_n) < kmerKey)) {
                *ef = 0; return middle;
            } else right = middle-1;
        } else left = middle+1;
    }
    *ef = 0;
    return h.kmer_num[hashKey];
}

#define set_kmer_key(n, kmerKey, h) { (*(n)) = (kmerKey << (h->hp.remn_ni)); }
#define set_MinEdge(n, h) { (*(n)) |= (1 << (h->hp.in_ni)); }
#define set_MoutEdge(n, h) { (*(n)) |= (1 << (h->hp.out_ni)); }
#define set_bwt_char(n, last_nt, h) { (*(n)) |= (last_nt << (h->hp.bwt_char_ni)); }
#define set_next_char(n, next_nt, h) { (*(n)) |= (next_nt << (h->hp.offset_ni)); }

#define has_MinEdge(n, h) (((n) >> (h->hp.in_ni)) & (h->hp.inout_m))
#define has_MoutEdge(n, h) (((n) >> (h->hp.out_ni)) & (h->hp.inout_m))
#define get_bwt_char(n, h) (((n) >> (h->hp.bwt_char_ni)) & (h->hp.char_m))
#define get_next_char(n, h) (((n) >> (h->hp.offset_ni)) & (h->hp.char_m))

int update_hash(uint32_t k_i, uint8_t ef, hash_idx *h, kmer_num_t hashKey, kmer_node_t kmerKey, uint8_t pre_nt, uint8_t next_nt)
{
    if (ef) {
        if (pre_nt != get_bwt_char(h->kmer_node[hashKey][k_i], h)) {
            set_MinEdge(h->kmer_node[hashKey]+k_i, h);
            set_bwt_char(h->kmer_node[hashKey]+k_i, nt_N, h);
        }
        if (next_nt != get_next_char(h->kmer_node[hashKey][k_i], h)) {
            set_MoutEdge(h->kmer_node[hashKey]+k_i, h);
            set_next_char(h->kmer_node[hashKey]+k_i, nt_N, h);
        }
    } else {
        kmer_num_t k_n = ++(h->kmer_num[hashKey]);
        if (k_n > h->kmer_m[hashKey]) {
            h->kmer_m[hashKey] <<= 1;
            h->kmer_node[hashKey] = (kmer_node_t*)_err_realloc(h->kmer_node[hashKey], h->kmer_m[hashKey] * sizeof(kmer_node_t));
        }
        // copy node
        memmove(h->kmer_node[hashKey]+k_i+1, h->kmer_node[hashKey]+k_i, (k_n-1-k_i)*sizeof(kmer_node_t));
        set_kmer_key(h->kmer_node[hashKey]+k_i, kmerKey, h);
        set_bwt_char(h->kmer_node[hashKey]+k_i, pre_nt, h);
        set_next_char(h->kmer_node[hashKey]+k_i, next_nt, h);
    }
    return 1;
}
        
int hash_check(hash_idx *h, kmer_int_t kmerInt, uint8_t pre_nt, uint8_t next_nt)
{
    kmer_num_t hashKey; uint32_t kmerKey, k_i;
    uint8_t equal_f=0;

    hashKey = first_n_bits(kmer_num_t, kmerInt, (h->hp.k)<<1, (h->hp.hash_k)<<1);
    kmerKey = last_n_bits(uint32_t, kmerInt, _node_size, h->hp.remn_n);

    k_i = bin_search_kmer(&equal_f, *h, hashKey, kmerKey);
    if (!update_hash(k_i, equal_f, h, hashKey, kmerKey, pre_nt, next_nt)) err_fatal_simple("ERROR: hash-update fail!\n");
    return 0;
}

void cons_kmer(kmer_node_t n, kmer_num_t hashKey, char **seq, uint8_t k, uint8_t remn_n, uint8_t remn_ni)
{
    kmer_int_t all = hashKey;
    all = (all << remn_n) | (n>>(remn_ni));
    int i;
    for (i = 0; i < k; ++i)
    {
        (*seq)[i] = "ACGTN"[(all >> (2 * (k-i-1)))&0x3];
    }
    fprintf(stdout, "%s\n", *seq);
}
void con_kmer(kmer_int_t all, char **seq, uint8_t k)
{
    int i;
    for (i = 0; i < k; ++i)
    {
        (*seq)[i] = "ACGTN"[(all >> (2 * (k-i-1)))&0x3];
    }
    fprintf(stdout, "%s %lld\n", *seq, all);
}
int kmer_gen(const char *prefix, hash_idx *h_idx)
{
    kmer_int_t kmerInt;
    gzFile fp; kseq_t *fa_seq;
    fp = xzopen(prefix, "r");
    fa_seq = kseq_init(fp);
    int seq_i; uint8_t pre_nt, next_nt;
    int count = 0;
    while (kseq_read(fa_seq) >= 0) {
        seq_i = 0;
        fprintf(stderr, "\n[%s] Generating and sorting kmer for %s ...",__func__, fa_seq->name.s);
        if (get_first_kmer(fa_seq, &seq_i, h_idx->hp.k, &kmerInt, &next_nt)) {        // first kmer
#ifdef __DEBUG__
            fprintf(stderr, "ref_i %d\n", seq_i);
#endif
            hash_check(h_idx, kmerInt, nt_N, next_nt);
            count++;
            while (get_shift_kmer(fa_seq, &seq_i, h_idx->hp.k, h_idx->hp.k_m, &kmerInt, &pre_nt, &next_nt)) { // generate kmer based on last kmer
#ifdef __DEBUG__
                fprintf(stderr, "ref_i %d\n", seq_i);
#endif
                hash_check(h_idx, kmerInt, pre_nt, next_nt);
                count++;
            }
        }
    }
    fprintf(stderr, "\n[%s] Generating and sorting done!\n",__func__);
    // output kmer
    char *kmerseq = (char*)malloc((h_idx->hp.k+1) * sizeof(char));
    uint32_t i, j;
    for (i = 0; i < h_idx->hp.hash_size; ++i) {
        for (j = 0; j < h_idx->kmer_num[i]; ++j) {
            cons_kmer(h_idx->kmer_node[i][j], i, &kmerseq, h_idx->hp.k, h_idx->hp.remn_n, h_idx->hp.remn_ni);
        }
    }
    free(kmerseq);
#ifdef __DEBUG__
    fprintf(stderr, "kmer count: %d\n", count);
#endif
    kseq_destroy(fa_seq); err_gzclose(fp);
    return 0;
}
