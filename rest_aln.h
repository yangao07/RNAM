#ifndef _REST_ALN_H_
#define _REST_ALN_H_

#include "debwt.h"
#include "bntseq.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#define CHUNK_READ_N 100000
#define CHUNK_SIZE  10000000

#define REST_SEED_LEN 22



typedef struct {
    int n_thread;
    int seed_len;
} rest_aln_para;

typedef struct {

} rest_seq_t;

typedef struct {
    int tid;

    debwt_t *db;       // index
    uint8_t *pac;
    bntseq_t *bns;

    int n_seqs;         // read seqs
    kseq_t *w_seqs;

    rest_aln_para *ap; // aln parameters
    // aux data during alignment
    // alignment result
} rest_aux_t;


int rest_aln(int argc, char *argv[]);

#endif
