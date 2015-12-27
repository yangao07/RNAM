#include <stdio.h> 
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include "build_debwt.h"
#include "rest_index.h"
#include "kmer_hash.h"
#include "kstring.h"
#include "utils.h"

void debwt_init(debwt_t *db, hash_idx h)
{
    // init debwt_bwt
    db->bwt_l = h.kmer_real_count + h.hp.k * h.uni_tol_count;
    db->n_occ = (db->bwt_l + _OCC_INV - 1) / _OCC_INV; // do NOT store last C[5]
    db->bwt_size = (db->bwt_l + _BWT_INV - 1) / _BWT_INV + db->n_occ * _OCC_C;
    db->bwt = (debwt_int_t*)_err_calloc(db->bwt_size, sizeof(debwt_int_t));
    db->bwt_unit = 0, db->bwt_i = 0, db->bwt_k = 0;


    // init debwt_sa
    // XXX db->bwt_l - db->n_unipath
    db->n_sa = db->bwt_l / _SA_INV;
    db->n_s_sa = h.uni_tol_count;

    db->sa_uid = (debwt_int_t*)calloc(db->n_sa, sizeof(debwt_int_t));
    db->s_sa_uid = (debwt_int_t*)calloc(db->n_s_sa, sizeof(debwt_int_t));

    // init C[]
    int i; for (i = 0; i < _OCC_C; ++i) db->C[i] = 0;

    // other things
    db->bwt_str = (kstring_t*)calloc(1, sizeof(kstring_t));

    db->n_unipath = h.uni_tol_count;
    db->n_offset = h.off_tol_count;

    db->uni_id = (debwt_count_t*)malloc(db->bwt_l * sizeof(debwt_count_t));
    db->uni_offset_c = (debwt_count_t*)malloc((db->n_unipath+1) * sizeof(debwt_count_t));
    db->uni_offset = (debwt_count_t*)malloc(db->n_offset * sizeof(debwt_count_t));
}

void debwt_free(debwt_t *db_idx)
{
    free(db_idx->bwt);
    free(db_idx->bwt_str->s); free(db_idx->bwt_str);

    free(db_idx->sa_uid); free(db_idx->s_sa_uid);

    free(db_idx->uni_id); free(db_idx->uni_offset); free(db_idx->uni_offset_c);
}

void print_debwt(debwt_t db_idx)
{
    stdout_printf("BWT_STR:\n%s\n", db_idx.bwt_str->s);
    stdout_printf("UID:\n");
    uint32_t i, j;
    /*
    for (i = 0; i < db_idx.bwt_l-db_idx.n_unipath; ++i) {
        fprintf(stdout, "# %d:\t%d\n", i+1, db_idx.uni_id[i]); 
    }*/
    stdout_printf("Uni_offset_c:\n");
    for (i = 0; i < db_idx.n_unipath; ++i) stdout_printf("U_c %d:\t%d\n", i+1, db_idx.uni_offset_c[i+1]);
    stdout_printf("Uni_offset:\n");
    for (i = 0; i < db_idx.n_unipath; ++i) {
        stdout_printf("U %d:\t", i+1);
        for (j = db_idx.uni_offset_c[i]; j < db_idx.uni_offset_c[i+1]; ++j)
            stdout_printf("%d ", db_idx.uni_offset[j]);
        stdout_printf("\n");
    }
}

// XXX (sa_i)>>_OCC_INV_B
#define debwt_bwt(db, sa_i) ((db).bwt[_OCC_C + ((sa_i)>>_OCC_INV_B<<_BWT_OCC_B) + ((sa_i)>>_OCC_INV_B) + (((sa_i)&(_OCC_INV_M))>>_BWT_INV_B)])
#define debwt_get_bwt_nt(db, sa_i) ((debwt_bwt(db, sa_i)>>((~(sa_i)&_BWT_INV_M)<<_BWT_NT_B))&_BWT_NT_M)
#define debwt_occ_a(db, sa_i) ((db).bwt + ((sa_i)>>_OCC_INV_B<<_BWT_OCC_B) + ((sa_i)>>_OCC_INV_B))
#define __occ_cnt4(table, b) (table[(b)&0xff] + table[(b)>>8&0xff] + table[(b)>>16&0xff] + table[(b)>>24])

static inline int __occ_aux(debwt_int_t b, debwt_t db, uint8_t c)
{
    b = ((c&4)?b:~b)>>2 & ((c&2)?b:~b)>>1 & ((c&1)?b:~b) & 0x1111111111111111ull;
    return __occ_cnt4(db.bit_table, b);
}

static inline debwt_count_t debwt_occ(const debwt_t db, debwt_count_t sa_uid_i, uint8_t bwt_nt)
{
    debwt_int_t *p, *end;
    debwt_count_t occ = (p = debwt_occ_a(db, sa_uid_i))[bwt_nt];
    p += _OCC_C;

    end = p + ((sa_uid_i & _OCC_INV_M) >> _BWT_INV_M);

    // OCC up to sa_uid_i/_OCC_INV
    for (; p < end; ++p) occ += __occ_aux(*p, db, bwt_nt);
    occ += __occ_aux(*p >> (((~sa_uid_i&_BWT_INV_M)+1)<<_BWT_NT_B), db, bwt_nt);
    if (bwt_nt == 0) occ -= ((~sa_uid_i&_BWT_INV_M)+1);
    return occ;
}

static inline debwt_int_t bwt_pre_sa_i(const debwt_t db, debwt_count_t sa_uid_i, uint8_t *bwt_nt)
{
    *bwt_nt = debwt_get_bwt_nt(db, sa_uid_i);
    if (*bwt_nt == nt_N) return sa_uid_i;
    else return db.C[*bwt_nt] + debwt_occ(db, sa_uid_i, *bwt_nt);
}

void debwt_cal_sa(debwt_t *db_idx)
{
    debwt_count_t sa_uid_i, sa_uid;
    int64_t i, j;
    //update C[4]
    for (i = 3; i >= 0; --i) {
        db_idx->C[i] = 0;
        for (j = i-1; j >= 0; --j)
            db_idx->C[i] += db_idx->C[j];
        db_idx->C[i] += db_idx->n_unipath;
    }
    // calculate cnt_table
    db_idx->bit_table[0] = 0;
    for (i = 0; i < 256; ++i) db_idx->bit_table[i] = (i&1) + db_idx->bit_table[i/2];
    // calculate uni_id
    uint8_t bwt_nt=0;
    for (i = 0; i < db_idx->n_unipath; ++i) {
        sa_uid_i = i/*0-base*/, sa_uid = i/*0-base*/;
        while (bwt_nt < nt_N) {
            if ((sa_uid_i+1) % _SA_INV == 0) db_idx->sa_uid[((sa_uid_i+1)>>_SA_INV_B)-1] = sa_uid;
            sa_uid_i = bwt_pre_sa_i(*db_idx, sa_uid_i, &bwt_nt);
        }
        // bwt_nt == nt_N, sa_uid_i
        db_idx->s_sa_uid[debwt_occ(*db_idx, sa_uid_i, nt_N)] = sa_uid;
    }
}

void push_debwt_bwt(uint8_t bwt_nt, debwt_t *db_idx)
{
    bwt_nt = bwt_nt > nt_N ? nt_N : bwt_nt;
    int i;
    if (db_idx->bwt_i % _OCC_INV == 0) {
        for (i = 0; i < _OCC_C; ++i) db_idx->bwt[db_idx->bwt_k++] = db_idx->C[i];
    }
    db_idx->bwt_unit <<= _BWT_NT_K;
    db_idx->bwt_unit |= bwt_nt;
    db_idx->C[bwt_nt]++;
    db_idx->bwt_i++;
                                                // last bwt
    if (db_idx->bwt_i % _BWT_INV == 0 || db_idx->bwt_i == db_idx->bwt_l) {
        db_idx->bwt[db_idx->bwt_k++] = db_idx->bwt_unit;
        db_idx->bwt_unit = 0;
    }
}

int build_debwt(char *prefix, hash_idx *h_idx, debwt_t *db_idx)
{
    err_printf("[build_debwt] Building de Bruijn-BWT index for genome ...\n");

    hash_init_num(h_idx); 
    // FIRST run for counting the total number of kmer
    kmer_tol_count(prefix, h_idx);
    // SECOND run for generating k-mer and countint spe-kmer
    //   update in/out flag
    //   update bwt_char
    //   count the total number of spe-kmer
    kmer_gen(prefix, h_idx);
    // THIRD run for generating unipath
    //   extra hash_num&node for $-contained kmer
    //   pos of unipath
    //   update bwt_char of head and tail kmer 
    //   uid and offset of kmer
    debwt_init(db_idx, *h_idx);
    spe_kmer_gen(prefix, h_idx, db_idx);
    // merge normal kmer and $-contained kmer
    kmer_merge(*h_idx, db_idx);
    // calculate SA
    debwt_cal_sa(db_idx);
    
    hash_free(h_idx);
    print_debwt(*db_idx);
    
    err_printf("[build_debwt] Building de Bruijn-BWT index done!\n");
    return 0;
}
