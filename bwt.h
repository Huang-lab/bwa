/* The MIT License

   Copyright (c) 2018-     Dana-Farber Cancer Institute
                 2009-2018 Broad Institute, Inc.
                 2008-2009 Genome Research Ltd. (GRL)

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <hli@jimmy.harvard.edu> */

#ifndef BWA_BWT_H
#define BWA_BWT_H

#include <stdint.h>
#include <stddef.h>

// Disk format (128-position interval, used by bwa index / bwt_bwtupdate_core)
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)
// Runtime one-hot format (64-position interval with popcount, converted at load time)
#define OCC_ONEHOT_SHIFT 6

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef unsigned char ubyte_t;
#endif

typedef uint64_t bwtint_t;

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	// occurance array, separated to two parts
	uint32_t cnt_table[256];
	// suffix array
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
	int is_onehot; // 1 if converted to one-hot format at runtime
} bwt_t;

typedef struct {
	bwtint_t x[3], info;
} bwtintv_t;

typedef struct { size_t n, m; bwtintv_t *a; } bwtintv_v;

// Dual-mode: supports both old (128-pos, 2-bit packed) and one-hot (64-pos, bitvectors)
static inline uint32_t *bwt_occ_intv(const bwt_t *bwt, bwtint_t k) {
	return bwt->is_onehot
		? bwt->bwt + ((k >> OCC_ONEHOT_SHIFT) << 4)   // one-hot: 64-pos blocks
		: bwt->bwt + ((k >> OCC_INTV_SHIFT) << 4);     // old: 128-pos blocks
}

// Old format helpers (for bwa index path)
#define bwt_bwt_old(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])

static inline int bwt_B0(const bwt_t *bwt, bwtint_t k) {
	if (bwt->is_onehot) {
		const uint64_t *oh = (const uint64_t*)bwt_occ_intv(bwt, k) + 4;
		int bit = k & 63;
		return (int)(((oh[1] | oh[3]) >> bit) & 1) | ((int)(((oh[2] | oh[3]) >> bit) & 1) << 1);
	} else {
		return (bwt_bwt_old(bwt, k) >> ((~(k) & 0xf) << 1)) & 3;
	}
}

#define bwt_set_intv(bwt, c, ik) ((ik).x[0] = (bwt)->L2[(int)(c)]+1, (ik).x[2] = (bwt)->L2[(int)(c)+1]-(bwt)->L2[(int)(c)], (ik).x[1] = (bwt)->L2[3-(c)]+1, (ik).info = 0)

#ifdef __cplusplus
extern "C" {
#endif

	void bwt_dump_bwt(const char *fn, const bwt_t *bwt);
	void bwt_dump_sa(const char *fn, const bwt_t *bwt);

	bwt_t *bwt_restore_bwt(const char *fn);
	void bwt_restore_sa(const char *fn, bwt_t *bwt);

	void bwt_destroy(bwt_t *bwt);

	void bwt_bwtgen(const char *fn_pac, const char *fn_bwt); // from BWT-SW
	void bwt_bwtgen2(const char *fn_pac, const char *fn_bwt, int block_size); // from BWT-SW
	void bwt_cal_sa(bwt_t *bwt, int intv);

	void bwt_bwtupdate_core(bwt_t *bwt);

	bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c);
	void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4]);
	bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k);

	// more efficient version of bwt_occ/bwt_occ4 for retrieving two close Occ values
	void bwt_gen_cnt_table(bwt_t *bwt);
	void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
	void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4]);

	int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end);
	int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0);

	/**
	 * Extend bi-SA-interval _ik_
	 */
	void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back);

	/**
	 * Given a query _q_, collect potential SMEMs covering position _x_ and store them in _mem_.
	 * Return the end of the longest exact match starting from _x_.
	 */
	int bwt_smem1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]);
	int bwt_smem1a(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]);

	int bwt_seed_strategy1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem);

#ifdef __cplusplus
}
#endif

#endif
