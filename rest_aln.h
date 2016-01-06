#ifndef _REST_ALN_H_
#define _REST_ALN_H_

#define REST_SEED_LEN 22

typedef struct {
    int seed_len;
} rest_aln_para;


int rest_aln(int argc, char *argv[]);

#endif
