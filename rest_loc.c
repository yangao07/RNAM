#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "rest_aln.h"
#include "kmer_hash.h"
#include "debwt.h"
#include "bntseq.h"
#include "kseq.h"
#include "utils.h"

#define MEM_LEN 19

int rest_read_seq(kseq_t *read_seq, int chunk_read_n)
{
    kseq_t *s = read_seq; int n = 0;
    while (kseq_read(s+n) >= 0) {
        if (chunk_read_n-1 == n++) break;
    }
    return n;
}

void aux_free(rest_aux_t *aux)
{
    int i;
    free(aux->ap);
    debwt_index_free(aux->db); free(aux->db);
    free(aux->pac); bns_destroy(aux->bns);

    for (i = 0; i < CHUNK_READ_N; ++i) {
        if (aux->w_seqs+i != NULL) {
            free((aux->w_seqs+i)->name.s);
            free((aux->w_seqs+i)->comment.s);
            free((aux->w_seqs+i)->seq.s);
            free((aux->w_seqs+i)->qual.s);
        }
    }
    ks_destroy(aux->w_seqs->f);
    free(aux->w_seqs);
    free(aux);
}

void set_seed_loc(loc_t *loc, int read_off, int read_loc_len, uni_sa_t uid, ref_off_t uni_off, ref_off_t uni_loc_len)
{
    loc->uid = uid;
    loc->uni_off = uni_off, loc->len1 = uni_loc_len;
    loc->read_off = read_off, loc->len2 = read_loc_len;
}

/*
 * uni_mem: directly match with pac XXX
 */
/*
 * read         [----- off1 ----- off2 -----]
 * unipath [ ... ----- off1 ----- off2 ---------- ... ]
 */
int uni_mem(uint8_t *read_seq, int read_len, uint32_t read_off1, uint32_t read_off2, ref_off_t uni_off1, ref_off_t uni_off2, debwt_t *db, uni_sa_t uid, bntseq_t *bns, uint8_t *pac, int *l1, int *l2)
{
    uni_sa_t uni_pos_i = db->uni_pos_c[uid];
    ref_off_t uni_pos = db->uni_pos[uni_pos_i], who_uni_len = db->uni_len[uid];
    ref_off_t uni_start, uni_len, uni_len1, uni_len2, off1, off2; // off1/2: for both read and unipath
    ref_off_t pac_coor;
    uint8_t *_read_seq = read_seq;

    if (uni_off1 >= read_off1) {
        off1 = read_off1;
        uni_len1 = read_off1; uni_start = uni_off1-read_off1;
        who_uni_len -= (uni_off1-read_off1), uni_off2 -= (uni_off1-read_off1);
    } else {
        off1 = uni_off1; 
        uni_len1 = uni_off1; uni_start = 0;
        _read_seq = read_seq+(read_off1-uni_off1);
        read_len -= (read_off1-uni_off1), read_off2 -= (read_off1-uni_off1);
    }
    off2 = read_off2; // read_off2 == uni_off2
    if ((who_uni_len-uni_off2) >= (read_len-read_off2)) {
        uni_len2 = read_len-read_off2-1;
    } else {
        uni_len2 = who_uni_len-uni_off2-1;
    }
    uni_len = uni_len1+off2-off1+1+uni_len2;

    uint8_t is_rev = _debwt_get_strand(db->uni_pos_strand, uni_pos_i);
    if (is_rev) pac_coor = uni_pos-uni_start-uni_len;
    else pac_coor = uni_pos+uni_start-1;
    uint8_t *uni_seq = _bns_get_seq(bns->l_pac, pac, pac_coor, uni_len, is_rev);
    int mem_l = bi_extend(_read_seq, uni_seq, off1, off2, uni_len, l1, l2);
    return mem_l;
}

int gen_mem_loc(int cur_i, uint8_t *q, int qlen, debwt_t *db, bntseq_t *bns, uint8_t *pac, rest_aln_para *ap, loc_t *loc) {
    int old_i;
    debwt_count_t k, l, il, i, uni_occ_thd = ap->debwt_uni_occ_thd;

    for (; cur_i >= 0; --cur_i) {
        old_i = cur_i;
        // debwt hash 
        int hk = get_hash_value(q+cur_i, _BWT_HASH_K);
        k = db->bwt_hash[hk];
        if ((il = db->bwt_hash_il[hk]) == 0) continue;
        l = k+il-1;
        // bwt backtrack
        while (il > uni_occ_thd && cur_i >= 1) {
            il = debwt_exact_match_alt(db, 1, q+cur_i-1, &k, &l);
            if (il == 0 || cur_i ==1) break;
            --cur_i;
        }
        if (il == 0) continue;
        // base extend (back/forward) // NOT span unipaths
        ref_off_t uni_off, max_uni_off, max_loc_len1; uni_sa_t uid, max_uid; int m_len, max_len=0, l1, l2, max_read_off, max_loc_len2;
        for (i = 0; i < il; ++i) {
            uid = debwt_sa(db, k+i, &uni_off);
            m_len = uni_mem(q, qlen, cur_i, old_i+_BWT_HASH_K-1, uni_off, uni_off+(old_i-cur_i)+_BWT_HASH_K-1, db, uid, bns, pac, &l1, &l2);
            if (m_len > max_len) {
                max_uid = uid;
                max_uni_off = uni_off-l1, max_loc_len1 = m_len;
                max_read_off = cur_i-l1, max_loc_len2 = m_len;
                max_len = m_len;
            }
        }
        cur_i = old_i;
        if (max_len >= MEM_LEN) {
            set_seed_loc(loc, max_read_off, max_loc_len2, max_uid, max_uni_off, max_loc_len1);
            return (max_read_off-_BWT_HASH_K);
        }
    }
    return 0;
}

loc_clu_t *init_loc_clu(void) {
    loc_clu_t *lc = (loc_clu_t*)_err_malloc(sizeof(loc_clu_t));
    lc->n = 0, lc->m = 8;
    lc->loc = (loc_t*)_err_malloc(8*sizeof(loc_t));
    lc->clu_w = (int*)_err_malloc(8*sizeof(int));
    return lc;
}

loc_clu_t *relloc_loc_clu(loc_clu_t *lc) {
    lc->m <<= 1;
    lc->loc = (loc_t*)_err_realloc(lc->loc, lc->m * sizeof(loc_t));
    lc->clu_w = (int*)_err_realloc(lc->clu_w, lc->m * sizeof(int));
    return lc;
}

void free_loc_clu(loc_clu_t *lc) { free(lc->loc); free(lc->clu_w); free(lc); }

int check_loc(loc_t l1, loc_t l2) {
    if (l1.uid == l2.uid) {
        if (abs(l1.read_off))
    }
}

int check_rough_loc(loc_t l, loc_clu_t *lc) {
    int i;
    for (i = 0; i < lc->n; ++i) {
        if (check_loc(l, lc->loc[i])) {
            if ((lc->clu_w[i]++) >= T) {
                return i+1;
            }
        }
    }
    return 0;
}

int rest_rough_loc(debwt_t *db, bntseq_t *bns, uint8_t *pac, rest_aln_para *ap, uint8_t *q, int qlen, int *sid, int *eid) {
    *sid = *eid = -1;
    loc_clu_t *lc = init_loc_clu();
    int cur_i, old_i;
    debwt_count_t k, l, il;
    // backwards: set eid
    for (cur_i = qlen-_BWT_HASH_K; cur_i >= 0; --cur_i) {
        // MEM seeding and locating
        loc_t loc;
        cur_i = gen_mem_loc(cur_i, q, qlen, db, bns, pac, ap, &loc);
        // push and check clu
        if (check_rough_loc(loc, lc)) { 
            // set sid
        }
    }
    // forwards: set sid
    for (j = qlen-1; j > i; --j) {
        // kmer
        // MEM seeding and locating
        // check with sid
        // if (check_is_OK)
        // set eid
    }
    return 0;
}

int rest_main_loc(rest_aux_t *aux)
{
    debwt_t *db = aux->db; bntseq_t *bns = aux->bns; uint8_t *pac = aux->pac;
    kseq_t *w_seqs = aux->w_seqs; int n_seqs = aux->n_seqs; 
    rest_aln_para *ap = aux->ap;
    int i_seq=0; uint64_t i;
    while (i_seq < n_seqs) {
        if (i_seq == n_seqs) break;
        kseq_t *seqs = w_seqs+i_seq;
        stdout_printf("%s\n%s\n", seqs->name.s, seqs->seq.s);
        uint8_t *bseq = (uint8_t*)_err_malloc(seqs->seq.l * sizeof(uint8_t));
        for (i = 0; i < seqs->seq.l; ++i) bseq[i] = nst_nt4_table[(int)(seqs->seq.s[i])];

        // rough seeding
        uint32_t sid, eid;
        rest_rough_loc(db, bns, pac, ap, bseq, seqs->seq.l, &sid, &eid);
        // push into buf and sort
        /// init read's buf
        /// push
        /// if(full) qsort and dump
        // dump to disk
        // multi-merge and dump to disk

        free(bseq);
        i_seq++;
    }
}

static void *rest_thread_main_loc(void *a)
{
    a=0;
    return a;
}

int THREAD_READ_I;
pthread_rwlock_t RWLOCK;

/*
 * function: roughly locate read
 * ret value: Start and End Rid
 * (Rid: id of specific ref-region)
 */
uint32_t rest_loc_core(const char *ref_fn, const char *read_fn, rest_aln_para *rest_ap)
{
    /* load index */
    err_printf("[rest_loc_core] Restoring ref-indices ... ");
    debwt_t *db_idx = debwt_restore_index(ref_fn);
    bntseq_t *bns = bns_restore(ref_fn);
    uint8_t *pac = (uint8_t*)_err_calloc(bns->l_pac/4+1, 1);
    fread(pac, 1, bns->l_pac/4+1, bns->fp_pac); err_printf("done!\n");
    
    /* read read.fa/fq */
    int n_seqs, i;
    gzFile readfp; kseq_t *read_seqs;
    readfp = err_xzopen_core(__func__, read_fn, "r");
    kstream_t *fs = ks_init(readfp);
    read_seqs = (kseq_t*)_err_calloc(CHUNK_READ_N, sizeof(kseq_t));
    for (i = 0; i < CHUNK_READ_N; ++i) read_seqs[i].f = fs;
    
    // alloc and init for auxiliary data
    if (rest_ap->n_thread < 1) rest_ap->n_thread = 1;
    rest_aux_t *aux = (rest_aux_t*)_err_calloc(rest_ap->n_thread, sizeof(rest_aux_t));
    for (i = 0; i < rest_ap->n_thread; ++i) {
        aux[i].tid = i; aux[i].ap = rest_ap;
        aux[i].db = db_idx; aux[i].bns = bns; aux[i].pac = pac;
    }

    if (rest_ap->n_thread <= 1) {
        while ((n_seqs = rest_read_seq(read_seqs, CHUNK_READ_N)) != 0) { 
            aux->n_seqs = n_seqs;
            aux->w_seqs = read_seqs;
            rest_main_loc(aux);
            // output
        }
    } else { // multi-threads
        pthread_rwlock_init(&RWLOCK, NULL);
        while ((n_seqs = rest_read_seq(read_seqs, CHUNK_READ_N)) != 0) { 
            THREAD_READ_I = 0;
            pthread_t *tid; pthread_attr_t attr;
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            tid = (pthread_t*)_err_calloc(rest_ap->n_thread, sizeof(pthread_t));
            int j;
            for (j = 0; j < rest_ap->n_thread; ++j) {
                aux[j].n_seqs = n_seqs;
                aux[j].w_seqs = read_seqs;
                pthread_create(&tid[j], &attr, rest_thread_main_loc, aux+j);
            }
            for (j = 0; j < rest_ap->n_thread; ++j) pthread_join(tid[j], 0);
            free(tid);
            // output
        }
        pthread_rwlock_destroy(&RWLOCK);
    }
    aux_free(aux);
    err_gzclose(readfp);
    return 0;
}
