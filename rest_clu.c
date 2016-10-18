#include <stdio.h>
#include <stdlib.h>
#include "rest_clu.h"
#include "rest_aln.h"
#include "utils.h"

int check_seed_clu(seed_clu_t *sc, loc_t loc) {
    int i;
    // check with max
    if (check) return sc->max;
    // check with others
    for (i = 0; i < sc->n; ++i) {
        if (i == sc->max) continue;
        if (check) return i;

    }
    return -1;
}

void update_seed_clu(seed_clu_t *sc, int ci, loc_t loc) {
    int i; int max = sc->max;

    ++sc->clu[ci].seed_num;
    if (ci != max) {
        if (sc->clu[ci].seed_num > sc->clu[max].seed_num)
            sc->max = ci;
    }
}

seed_clu_t *init_seed_clu(void) {
    seed_clu_t *sc;
    sc = (seed_clu_t*)_err_malloc(sizeof(seed_clu_t));
    sc->n = 0, sc->m = 10, sc->max = 0;
    sc->clu = (s_clu*)_err_calloc(sc->m, sizeof(s_clu));
}

void realloc_seed_clu(seed_clu_t *sc) {
    sc->m <<= 1;
    sc->clu = _err_realloc(sc->clu, sc->m * sizeof(s_clu));
}

void free_seed_clu(seed_clu_t *sc) {
    free(sc->clu); free(sc);
}

void add_seed_clu(seed_clu_t *sc, loc_t loc) {
    if (sc->n == sc->m) {
        realloc_seed_clu(sc);
    }
    // set new clu
    sc->clu[sc->n];
    ++sc->n;
}


/* @input:  seeds location of one read
 * @output: one seed cluster represents the read
 * @para:                                      */
int seed_clu(seed_loc_t *sl)
{
    int i, j, res;
    seed_clu_t *sc = init_seed_clu();

    for (i = 0; i < sl->n; ++i) {
        if ((res = check_seed_clu(sc, sl->loc[i])) >= 0) {
            update_seed_clu(sc, res, sl->loc[i]);
        } else add_seed_clu(sc, sl->loc[i]);
    }
    //
    sc->clu[sc->max];
    free_seed_clu(sc);
    return 0;
}

/* @input:   read's RPRST unipath/pos and whole seeds' locs
 * @process: set corresponding bin based on the RPRST
 * @para:                                                  */
int read_clu()
{
    return 0;
}
