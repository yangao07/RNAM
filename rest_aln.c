#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <pthread.h>
#include "rest_aln.h"
#include "debwt.h"
#include "utils.h"
#include "kseq.h"

int rest_aln_usage(void)
{
    fprintf(stderr, "rest aln [option] ref.fa read.fq/fa\n");
    return 0;
}

rest_aln_para *rest_init_ap(void)
{
    rest_aln_para *ap = (rest_aln_para*)calloc(1, sizeof(rest_aln_para));
    ap->n_thread = 1;
    ap->seed_len = REST_SEED_LEN;
    return ap;
}

void aux_free(rest_aux_t *aux, int t)
{
    int i;
    for (i = 0; i < t; ++i) {
    }
    free(aux->ap);
    debwt_index_free(aux->db); free(aux->db);
    free(aux->pac); bns_destroy(aux->bns);
    kseq_destroy(aux->w_seqs);
    free(aux);
}

int rest_read_seq(kseq_t *read_seq, int chunk_read_n)
{
    kseq_t *s = read_seq; int n = 0;
    while (kseq_read(s+n) >= 0) {
        if (chunk_read_n-1 == n++) break;
    }
    return n;
}
int THREAD_READ_I;
pthread_rwlock_t RWLOCK;

int rest_main_aln(rest_aux_t *aux)
{
    debwt_t *db = aux->db; bntseq_t *bns = aux->bns; uint8_t *pac = aux->pac;
    kseq_t *w_seqs = aux->w_seqs; int n_seqs = aux->n_seqs; 
    rest_aln_para *ap = aux->ap;
    int i_seq=0, i;
    while (i_seq < n_seqs) {
        if (i_seq == n_seqs) break;
        kseq_t *seqs = w_seqs+i_seq;
        stdout_printf("%s\n%s\n", seqs->name.s, seqs->seq.s);
        uint8_t *bseq = (uint8_t*)_err_malloc(seqs->seq.l * sizeof(uint8_t));
        for (i = 0; i < seqs->seq.l; ++i) bseq[i] = nst_nt4_table[(int)(seqs->seq.s[i])];

        // seeding and locating
        // stitching
        // local assembly


        debwt_count_t ok, ol, l;
        uint32_t uid, off, m;
        l = debwt_exact_match(db, seqs->seq.l, bseq, &ok, &ol);
        for (i = 0; i < l; ++i) {
            uid = debwt_sa(db, ok+i, &off);
            for (m = db->uni_pos_c[uid]; m < db->uni_pos_c[uid+1]; ++m)
                stdout_printf("%c%d\t%d\n", "+-"[_debwt_get_strand(db->uni_pos_strand, m)], (int)db->uni_pos[m], off);
        }
        free(bseq);
        i_seq++;
    }
    return 0;
}

static void *rest_thread_main_aln(void *a)
{
    rest_aux_t *aux = (rest_aux_t*)a;
    debwt_t *db = aux->db; bntseq_t *bns = aux->bns; uint8_t *pac = aux->pac;
    kseq_t *w_seqs = aux->w_seqs; int n_seqs = aux->n_seqs; 
    rest_aln_para *ap = aux->ap;
    int i;
    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i = THREAD_READ_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i >= n_seqs) break;
        kseq_t *seqs = w_seqs+i;
        stdout_printf("%s\n%s\n", seqs->name.s, seqs->seq.s);
        uint8_t *bseq = (uint8_t*)_err_malloc(seqs->seq.l * sizeof(uint8_t));
        int j;
        for (j = 0; j < seqs->seq.l; ++j) bseq[j] = nst_nt4_table[(int)(seqs->seq.s[j])];
        debwt_count_t ok, ol, l;
        uint32_t uid, off, m, n;
        l = debwt_exact_match(db, seqs->seq.l, bseq, &ok, &ol);
        for (j = 0; j < l; ++j) {
            uid = debwt_sa(db, ok+j, &off);
            for (m = db->uni_pos_c[uid]; m < db->uni_pos_c[uid+1]; ++m)
                stdout_printf("%c%d\t%d\n", "+-"[_debwt_get_strand(db->uni_pos_strand, m)], (int)db->uni_pos[m], off);
        }
        free(bseq);
    }
    return 0;
}

int rest_aln_core(const char *ref_fn, const char *read_fn, rest_aln_para *rest_ap)
{
    /* load index */
    err_printf("[rest_aln_core] Restoring ref-indices ... ");
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
            rest_main_aln(aux);
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
                pthread_create(&tid[j], &attr, rest_thread_main_aln, aux+j);
            }
            for (j = 0; j < rest_ap->n_thread; ++j) pthread_join(tid[j], 0);
            free(tid);
            // output
        }
        pthread_rwlock_destroy(&RWLOCK);
    }
    aux_free(aux, rest_ap->n_thread);
    for(i = 0; i < CHUNK_READ_N; ++i) {
        free((read_seqs+i)->name.s); free((read_seqs+i)->comment.s); free((read_seqs+i)->seq.s); free((read_seqs+i)->qual.s);
    } free(read_seqs); ks_destroy(fs); err_gzclose(readfp);
    return 0;
}

int rest_aln(int argc, char *argv[])
{
    int c;
    rest_aln_para *rest_ap = rest_init_ap();

    while ((c = getopt(argc, argv, "t:l:")) >= 0) {
        switch (c)
        {
            case 'l': rest_ap->seed_len = atoi(optarg); break;
            default: return rest_aln_usage();
        }
    }
    if (argc - optind != 2) return rest_aln_usage();

    rest_aln_core(argv[optind], argv[optind+1], rest_ap);
    return 0;
}
