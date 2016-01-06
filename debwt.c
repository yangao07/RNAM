#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <math.h>
#include "debwt.h"
#include "rest_index.h"
#include "kmer_hash.h"
#include "kstring.h"
#include "utils.h"

#ifdef __LIT__
#define _BWT_HASH_K 2
#else
#define _BWT_HASH_K 13
#endif

void debwt_init(debwt_t *db, hash_idx *h)
{
    // init debwt_bwt
    db->bwt_l = h->kmer_real_count + h->hp.k * h->uni_tol_count;
    debwt_count_t n_occ = (db->bwt_l + _OCC_INV - 1) / _OCC_INV; // do NOT store last C[5]
    db->bwt_size = (db->bwt_l + _BWT_INV - 1) / _BWT_INV + n_occ * _OCC_C;
    db->bwt = (debwt_int_t*)_err_calloc(db->bwt_size, sizeof(debwt_int_t));
    db->bwt_unit = 0, db->bwt_i = 0, db->bwt_k = 0;

    // init C[]
    int i; for (i = 0; i < _OCC_C; ++i) db->C[i] = 0;

    // init debwt_sa
    // XXX db->bwt_l - db->n_unipath
    db->n_sa = db->bwt_l / _SA_INV;
    db->n_s_sa = h->uni_tol_count;

    db->sa_uid = (debwt_count_t*)calloc(db->n_sa, sizeof(debwt_count_t));
    db->sa_u_off = (debwt_count_t*)calloc(db->n_sa, sizeof(debwt_count_t));
    db->s_sa_uid = (debwt_count_t*)calloc(db->n_s_sa, sizeof(debwt_count_t));

    // other things
    db->n_unipath = h->uni_tol_count;
    db->n_offset = h->off_tol_count;

    db->uni_offset_c = (debwt_count_t*)malloc((db->n_unipath+1) * sizeof(debwt_count_t));
    db->uni_offset = (ref_offset_t*)malloc(db->n_offset * sizeof(ref_offset_t));
}

void debwt_init0(debwt_t *db, hash_idx *h)
{
    // init debwt_bwt
    db->bwt_l = h->kmer_real_count + h->hp.k * h->uni_tol_count;
    debwt_count_t n_occ = (db->bwt_l + _OCC_INV - 1) / _OCC_INV; // do NOT store last C[5]
    db->bwt_size = (db->bwt_l + _BWT_INV - 1) / _BWT_INV + n_occ * _OCC_C;
    db->bwt = (debwt_int_t*)_err_calloc(db->bwt_size, sizeof(debwt_int_t));
    db->bwt_unit = 0, db->bwt_i = 0, db->bwt_k = 0;

    // init C[]
    int i; for (i = 0; i < _OCC_C; ++i) db->C[i] = 0;
}

void debwt_init1(debwt_t *db, hash_idx *h)
{
    // init debwt_sa
    // XXX db->bwt_l - db->n_unipath
    db->n_sa = db->bwt_l / _SA_INV;
    db->n_s_sa = h->uni_tol_count;

    db->sa_uid = (debwt_count_t*)calloc(db->n_sa, sizeof(debwt_count_t));
    db->sa_u_off = (debwt_count_t*)calloc(db->n_sa, sizeof(debwt_count_t));
    db->s_sa_uid = (debwt_count_t*)calloc(db->n_s_sa, sizeof(debwt_count_t));

    // offsets
    db->n_unipath = h->uni_tol_count;
    db->n_offset = h->off_tol_count;

    db->uni_offset_c = (debwt_count_t*)malloc((db->n_unipath+1) * sizeof(debwt_count_t));
    db->uni_offset = (ref_offset_t*)malloc(db->n_offset * sizeof(ref_offset_t));
}

void debwt_free(debwt_t *db_idx)
{
    free(db_idx->bwt);

    free(db_idx->sa_uid); free(db_idx->sa_u_off); free(db_idx->s_sa_uid);

    free(db_idx->uni_offset); free(db_idx->uni_offset_c);

    free(db_idx->bwt_hash); free(db_idx->bwt_hash_il);
}

// argv:
// sa_uid_i: 0-base sa index
// off: 0-base offset from start of the unipath
debwt_count_t debwt_sa(const debwt_t *db, debwt_count_t sa_uid_i, ref_offset_t *off);

void print_debwt(debwt_t *db_idx)
{
    uint32_t i, j;
    stdout_printf("UID:\n");
    for (i = db_idx->n_unipath; i < db_idx->bwt_l; ++i) {
        stdout_printf("%d: %d\t", i, debwt_sa(db_idx, i, &j));
        stdout_printf("%d\n", j);
    }
    /*
    for (i = 0; i < db_idx->bwt_l-db_idx->n_unipath; ++i) {
        fprintf(stdout, "# %d:\t%d\n", i+1, db_idx->uni_id[i]); 
    }*/
    stdout_printf("Uni_offset_c:\n");
    for (i = 0; i < db_idx->n_unipath; ++i) stdout_printf("U_c %d:\t%d\n", i+1, db_idx->uni_offset_c[i+1]);
    stdout_printf("Uni_offset:\n");
    for (i = 0; i < db_idx->n_unipath; ++i) {
        stdout_printf("U %d:\t", i+1);
        for (j = db_idx->uni_offset_c[i]; j < db_idx->uni_offset_c[i+1]; ++j)
            stdout_printf("%d ", db_idx->uni_offset[j]);
        stdout_printf("\n");
    }

    stdout_printf("BWT-HASH:\n");
    for (i = 0; i < db_idx->bwt_hash_size; ++i) {
        stdout_printf("%d: %d %d\n", i, db_idx->bwt_hash[i], db_idx->bwt_hash_il[i]);
    }
}

// XXX (sa_i)>>_OCC_INV_B
#define debwt_bwt(db, sa_i) ((db)->bwt[_OCC_C + ((sa_i)>>_OCC_INV_B<<_BWT_OCC_B) + ((sa_i)>>_OCC_INV_B) + (((sa_i)&(_OCC_INV_M))>>_BWT_INV_B)])
#define debwt_get_bwt_nt(db, sa_i) ((debwt_bwt(db, sa_i)>>((~(sa_i)&_BWT_INV_M)<<_BWT_NT_B))&_BWT_NT_M)
#define debwt_occ_a(db, sa_i) ((db)->bwt + ((sa_i)>>_OCC_INV_B<<_BWT_OCC_B) + ((sa_i)>>_OCC_INV_B))
#define __occ_cnt4(table, b) (table[(b)&0x1111] + table[(b)>>16&0x1111] + table[(b)>>32&0x1111] + table[(b)>>48&0x1111])
#define __occ_cnt2(table, b) (table[(b)&0x11111111] + table[(b)>>32&0x11111111])

static inline int __occ_aux(debwt_int_t b, const debwt_t *db, uint8_t c)
{
    b = ((c&4)?b:~b)>>2 & ((c&2)?b:~b)>>1 & ((c&1)?b:~b) & 0x1111111111111111ull;
    return __occ_cnt4(db->bit_table16, b);
}

static inline debwt_count_t debwt_occ(const debwt_t *db, debwt_count_t sa_uid_i, uint8_t nt)
{
    debwt_int_t *p, *end;
    debwt_count_t occ = (p = debwt_occ_a(db, sa_uid_i))[nt];
    p += _OCC_C;

    // OCC up to sa_uid_i/_BWT_INV
    end = p + ((sa_uid_i & _OCC_INV_M) >> _BWT_INV_B);
    for (; p < end; ++p) occ += __occ_aux(*p, db, nt);

    occ += __occ_aux((*p >> ((~sa_uid_i&_BWT_INV_M)<<_BWT_NT_B))>>(1<<_BWT_NT_B), db, nt);
    if (nt == 0) occ -= ((~sa_uid_i&_BWT_INV_M)+1);
    //occ += __occ_aux(*p >> ((~(sa_uid_i)&_BWT_INV_M)<<_BWT_NT_B), db, nt);
    //occ--;
    //if (nt == 0) occ -= (~(sa_uid_i)&_BWT_INV_M);
    return occ;
}

static inline debwt_int_t bwt_pre_sa_i(const debwt_t *db, debwt_count_t sa_uid_i, uint8_t bwt_nt)
{
    if (bwt_nt == nt_N) return sa_uid_i;
    else return db->C[bwt_nt] + debwt_occ(db, sa_uid_i, bwt_nt);
}

void debwt_cnt_table16(debwt_t *db_idx)
{
    int i;
    db_idx->bit_table16[0] = 0;
    for (i = 0; i < 65536; ++i) 
        db_idx->bit_table16[i] = (i&1) + db_idx->bit_table16[i/2];
}

void debwt_cal_sa(debwt_t *db_idx)
{
    err_printf("[%s] Calculating SA ...\n", __func__);
    debwt_count_t sa_uid_i, sa_uid;
    int64_t i, j;
    //update C[5]
    for (i = 4; i >= 0; --i) {
        db_idx->C[i] = 0;
        for (j = i-1; j >= 0; --j)
            db_idx->C[i] += db_idx->C[j];
        db_idx->C[i] += db_idx->n_unipath;
    }
    // calculate cnt_table
    debwt_cnt_table16(db_idx);
    // calculate uni_id
    uint8_t bwt_nt=0; debwt_count_t off, off_n, off_m, *off_a;
    off_m = 100; off_a = (debwt_count_t*)malloc(off_m * sizeof(debwt_count_t));
    for (i = 0; i < db_idx->n_unipath; ++i) {
        sa_uid_i = i/*0-base*/, sa_uid = i/*0-base*/;
        off = 0; off_n = 0;
        while ((bwt_nt = debwt_get_bwt_nt(db_idx, sa_uid_i)) < nt_N) {
            if ((sa_uid_i+1) % _SA_INV == 0) {
                db_idx->sa_uid[((sa_uid_i+1)>>_SA_INV_B)-1] = sa_uid;
                db_idx->sa_u_off[((sa_uid_i+1)>>_SA_INV_B)-1] = off;
                if (off_n >= off_m) {
                    off_m <<= 1;
                    off_a = (debwt_count_t*)_err_realloc(off_a, off_m * sizeof(debwt_count_t));
                }
                off_a[off_n++] = ((sa_uid_i+1)>>_SA_INV_B)-1;
            }
            sa_uid_i = bwt_pre_sa_i(db_idx, sa_uid_i, bwt_nt);
            off++;
#ifdef __DEBUG__
            stdout_printf("sa_i: %d\t%c\n", sa_uid_i, "ACGTN"[bwt_nt]);
#endif
        }
        // bwt_nt == nt_N, sa_uid_i
        db_idx->s_sa_uid[debwt_occ(db_idx, sa_uid_i, nt_N)] = sa_uid;
        // reset all the sa_u_off: p -> P-p
        for (j = 0; j < off_n; ++j) 
            db_idx->sa_u_off[off_a[j]] = off-db_idx->sa_u_off[off_a[j]];
    }
    free(off_a);
    err_printf("[%s] Calculating SA done!\n", __func__);
}

debwt_count_t debwt_sa(const debwt_t *db, debwt_count_t sa_uid_i, ref_offset_t *off)
{
    *off = 0;
    uint8_t nt = debwt_get_bwt_nt(db, sa_uid_i);
    while (nt != nt_N && ((sa_uid_i+1) & _SA_INV_M)) {
        (*off)++;
        sa_uid_i = bwt_pre_sa_i(db, sa_uid_i, nt);

        nt = debwt_get_bwt_nt(db, sa_uid_i);
    }
    if (nt == nt_N) {
        return db->s_sa_uid[debwt_occ(db, sa_uid_i, nt_N)];
    } else { 
        *off += db->sa_u_off[((sa_uid_i+1)>>_SA_INV_B)-1];
        return db->sa_uid[((sa_uid_i+1)>>_SA_INV_B)-1];
    }
}

void push_debwt_bwt(uint8_t bwt_nt, debwt_t *db_idx)
{
    bwt_nt = bwt_nt > nt_N ? nt_N : bwt_nt;
    //stdout_printf("BWT: %c\n", "ACGTN"[bwt_nt]);
    int i;
    if (db_idx->bwt_i % _OCC_INV == 0) {
        for (i = 0; i < _OCC_C; ++i) db_idx->bwt[db_idx->bwt_k++] = db_idx->C[i];
    }
    db_idx->bwt_unit |= (debwt_int_t)(bwt_nt) << ((~db_idx->bwt_i&_BWT_INV_M)<<_BWT_NT_B);
    db_idx->C[bwt_nt]++;
    db_idx->bwt_i++;
                                                // last bwt
    if (db_idx->bwt_i % _BWT_INV == 0 || db_idx->bwt_i == db_idx->bwt_l) {
        db_idx->bwt[db_idx->bwt_k++] = db_idx->bwt_unit;
        db_idx->bwt_unit = 0;
    }
}

#define debwt_set_intv(d, nt, i) {(i).u = (d)->C[(nt)], (i).l = (d)->C[(nt)+1]-1;}

int debwt_backtrack(debwt_t *db_idx, debwt_intv_t old_intv, uint8_t nt, debwt_intv_t *new_intv)
{
    debwt_count_t c1 = debwt_occ(db_idx, old_intv.u, nt);
    debwt_count_t c2 = debwt_occ(db_idx, old_intv.l+1, nt);
    if (c2 > c1) {
        new_intv->u = db_idx->C[nt] + c1;
        new_intv->l = db_idx->C[nt] + c2 - 1;
        return 1;
    }
    else return 0;
}

// before calculate hash, bwt and sa have been calculated
void debwt_cal_hash(debwt_t *db_idx)
{
    err_printf("[%s] Calculating bwt-hash ...\n", __func__);
    db_idx->bwt_hash_k = _BWT_HASH_K; // 13
    db_idx->bwt_hash_size = (uint32_t)pow(4, _BWT_HASH_K);
    db_idx->bwt_hash = (debwt_count_t*)_err_calloc(db_idx->bwt_hash_size, sizeof(debwt_count_t));
    db_idx->bwt_hash_il = (debwt_count_t*)_err_calloc(db_idx->bwt_hash_size, sizeof(debwt_count_t));

    int i, j; uint8_t *hash_nt = (uint8_t*)malloc(_BWT_HASH_K * sizeof(uint8_t));
    debwt_intv_t intv, new_intv;

    for (i = 0; i < db_idx->bwt_hash_size; ++i) {
        get_hash_nt(hash_nt, i, _BWT_HASH_K);
        debwt_set_intv(db_idx, hash_nt[_BWT_HASH_K-1], intv);
        for (j = _BWT_HASH_K-2; j >= 0; --j) {
            if (!debwt_backtrack(db_idx, intv, hash_nt[j], &new_intv)) {
                new_intv.u = new_intv.l = 0; // not exist in debwt
                break;
            }
            intv = new_intv;
        }
        db_idx->bwt_hash[i] = new_intv.u;
        db_idx->bwt_hash_il[i] = new_intv.l - new_intv.u + 1;
        //stdout_printf("%d:\t %d %d\n", i, db_idx->bwt_hash[i], db_idx->bwt_hash_il[i]);
    }
    free(hash_nt);
    err_printf("[%s] Calculating bwt-hash done!\n", __func__);
}

void debwt_dump(const char *prefix, debwt_t *db)
{
    char *fn = (char*)calloc(strlen(prefix)+15, sizeof(char));
    // dump bwt
    err_printf("[%s] Writing bwt to file ...\n", __func__);
    strcpy(fn, prefix); strcat(fn, ".rest.debwt");
    FILE *fp = xopen(fn, "wb");
    err_fwrite(&db->bwt_l, sizeof(debwt_count_t), 1, fp);
    err_fwrite(&db->bwt_size, sizeof(debwt_count_t), 1, fp);
    err_fwrite(db->C, sizeof(debwt_count_t), 5, fp);
    err_fwrite(db->bwt, sizeof(debwt_int_t), db->bwt_size, fp);
    err_fflush(fp);
    err_fclose(fp);

    // dump sa
    err_printf("[%s] Writing sa to file ...\n", __func__);
    strcpy(fn, prefix); strcat(fn, ".rest.usa");
    fp = xopen(fn, "wb");
    err_fwrite(&db->bwt_l, sizeof(debwt_count_t), 1, fp);
    err_fwrite(&db->n_sa, sizeof(debwt_count_t), 1, fp);
    err_fwrite(&db->n_s_sa, sizeof(debwt_count_t), 1, fp); // n_unipath
    err_fwrite(&db->n_offset, sizeof(debwt_count_t), 1, fp);

    err_fwrite(db->sa_uid, sizeof(debwt_count_t), db->n_sa, fp);
    err_fwrite(db->sa_u_off, sizeof(debwt_count_t), db->n_sa, fp);
    err_fwrite(db->s_sa_uid, sizeof(debwt_count_t), db->n_s_sa, fp);

    err_fwrite(db->uni_offset_c+1, sizeof(debwt_count_t), db->n_unipath, fp);
    err_fwrite(db->uni_offset, sizeof(ref_offset_t), db->n_offset, fp);
    err_fflush(fp);
    err_fclose(fp);

    // dump bwt-hash
    err_printf("[%s] Writing bwt-hash to file ...\n", __func__);
    strcpy(fn, prefix); strcat(fn, ".rest.dhs");
    fp = xopen(fn, "wb");
    err_fwrite(&db->bwt_hash_size, sizeof(uint32_t), 1, fp);
    err_fwrite(db->bwt_hash, sizeof(debwt_count_t), db->bwt_hash_size, fp);
    err_fwrite(db->bwt_hash_il, sizeof(debwt_count_t), db->bwt_hash_size, fp);
    err_fflush(fp);
    err_fclose(fp);

    free(fn);
    err_printf("[%s] Writing index done.\n", __func__);
}

static void fread_fix(FILE *fp, debwt_count_t size, void *a)
{
	const int bufsize = 0x1000000; // 16M block
	debwt_count_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	//return offset;
}

debwt_t *debwt_restore_index(const char *prefix)
{
    char *fn = (char*)calloc(strlen(prefix)+15, sizeof(char));
    debwt_t *db;
    // restore bwt
    strcpy(fn, prefix); strcat(fn, ".rest.debwt");
    FILE *fp = xopen(fn, "rb");

    db = (debwt_t*)_err_calloc(1, sizeof(debwt_t));
    err_fread_noeof(&db->bwt_l, sizeof(debwt_count_t), 1, fp);
    err_fread_noeof(&db->bwt_size, sizeof(debwt_count_t), 1, fp);
    err_fread_noeof(db->C, sizeof(debwt_count_t), 5, fp);
    db->bwt = (debwt_int_t*)_err_calloc(db->bwt_size, sizeof(debwt_int_t));
    fread_fix(fp, db->bwt_size * sizeof(debwt_int_t), db->bwt);
    err_fclose(fp); 

    // restore sa
    strcpy(fn, prefix); strcat(fn, ".rest.usa");
    fp = xopen(fn, "rb");
    debwt_count_t bwt_l;
    err_fread_noeof(&bwt_l, sizeof(debwt_count_t), 1, fp);
    xassert(bwt_l == db->bwt_l, "SA-BWT inconsistency: bwt length is not the same.");
    err_fread_noeof(&db->n_sa, sizeof(debwt_count_t), 1, fp);
    err_fread_noeof(&db->n_s_sa, sizeof(debwt_count_t), 1, fp);
    db->n_unipath = db->n_s_sa;
    err_fread_noeof(&db->n_offset, sizeof(debwt_count_t), 1, fp);

    db->sa_uid = (debwt_count_t*)_err_calloc(db->n_sa, sizeof(debwt_count_t));
    fread_fix(fp, db->n_sa * sizeof(debwt_count_t), db->sa_uid);
    db->sa_u_off = (debwt_count_t*)_err_calloc(db->n_sa, sizeof(debwt_count_t));
    fread_fix(fp, db->n_sa * sizeof(debwt_count_t), db->sa_u_off);
    db->s_sa_uid = (debwt_count_t*)_err_calloc(db->n_s_sa, sizeof(debwt_count_t));
    fread_fix(fp, db->n_s_sa * sizeof(debwt_count_t), db->s_sa_uid);

    db->uni_offset_c = (debwt_count_t*)_err_calloc(db->n_unipath+1, sizeof(debwt_count_t));
    fread_fix(fp, db->n_unipath * sizeof(debwt_count_t), db->uni_offset_c+1);
    db->uni_offset = (ref_offset_t*)_err_calloc(db->n_offset, sizeof(ref_offset_t));
    fread_fix(fp, db->n_offset * sizeof(ref_offset_t), db->uni_offset);
    err_fclose(fp); 

    // restore bwt-hash
    strcpy(fn, prefix); strcat(fn, ".rest.dhs");
    fp = xopen(fn, "rb");
    err_fread_noeof(&db->bwt_hash_size, sizeof(uint32_t), 1, fp);
    db->bwt_hash = (debwt_count_t*)_err_calloc(db->bwt_hash_size, sizeof(debwt_count_t));
    db->bwt_hash_il = (debwt_count_t*)_err_calloc(db->bwt_hash_size, sizeof(debwt_count_t));
    fread_fix(fp, db->bwt_hash_size * sizeof(debwt_count_t), db->bwt_hash);
    fread_fix(fp, db->bwt_hash_size * sizeof(debwt_count_t), db->bwt_hash_il);
    err_fclose(fp);

    // count table16
    debwt_cnt_table16(db);

    free(fn);
    return db;
}

int pac_build_debwt(const char *prefix, debwt_pac_t *db_pac, debwt_count_t l_pac, hash_idx *h_idx, debwt_t *db_idx)
{
    err_printf("[pac_build_debwt] Building de Bruijn-BWT index for genome ...\n");
    debwt_kmer_hash_init_num(h_idx); 

    pac_count_kmer(db_pac, l_pac, h_idx); // count tolal number of kmer
    pac_gen_kmer(db_pac, l_pac, h_idx);   // generate normal kmer

    debwt_init0(db_idx, h_idx);
    pac_gen_skmer(db_pac, l_pac, h_idx, db_idx); // generate special kmer and partial bwt 
    free(db_pac);
    debwt_kmer_hash_free0(h_idx);

    debwt_init1(db_idx, h_idx);
    kmer_merge(h_idx, db_idx); // merge normal kmer and $-contained kmer, generate whole bwt
    debwt_kmer_hash_free1(h_idx);

    debwt_cal_sa(db_idx);   // calculate SA
    debwt_cal_hash(db_idx); // calculate bwt_hash range
    debwt_dump(prefix, db_idx); // dump index files to disk
#ifdef __DEBUG__
    print_debwt(db_idx);
    debwt_t *d = debwt_restore_index(prefix);
    print_debwt(d);
    debwt_free(d); free(d);
#endif
    debwt_free(db_idx);
    err_printf("[pac_build_debwt] Building de Bruijn-BWT index done!\n");
    return 0;
}
