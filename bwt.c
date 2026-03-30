/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

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

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include "utils.h"
#include "bwt.h"
#include "kvec.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k) // compute inverse CSA
{
	bwtint_t x = k - (k > bwt->primary);
	x = bwt_B0(bwt, x);
	x = bwt->L2[x] + bwt_occ(bwt, k, x);
	return k == bwt->primary? 0 : x;
}

// bwt->bwt and bwt->occ must be precalculated
void bwt_convert_to_onehot(bwt_t *bwt); // forward declaration

void bwt_cal_sa(bwt_t *bwt, int intv)
{
	bwtint_t isa, sa, i; // S(isa) = sa
	int intv_round = intv;

	kv_roundup32(intv_round);
	xassert(intv_round == intv, "SA sample interval is not a power of 2.");
	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");
	if (!bwt->is_onehot) bwt_convert_to_onehot(bwt); // convert before SA computation

	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	// calculate SA value
	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % intv == 0) bwt->sa[isa/intv] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;
	bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t sa = 0, mask = bwt->sa_intv - 1;
	// prefetch the SA entry — k will change but the final SA block may be nearby
	__builtin_prefetch(&bwt->sa[k / bwt->sa_intv], 0, 0);
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	return sa + bwt->sa[k/bwt->sa_intv];
}

// Old-format occurrence counting (used during bwa index before one-hot conversion)
static inline int __occ_aux(uint64_t y, int c)
{
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

static bwtint_t bwt_occ_old(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
	bwtint_t n;
	uint32_t *p, *end;
	k -= (k >= bwt->primary);
	n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
	p += sizeof(bwtint_t);
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
	n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31;
	return n;
}

// One-hot occurrence counting — single popcount per base
static inline uint64_t occ_mask(bwtint_t k)
{
	int y = k & 63;
	return ~0ULL >> (63 ^ y);
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (bwtint_t)(-1)) return 0;
	if (!bwt->is_onehot) return bwt_occ_old(bwt, k, c);
	k -= (k >= bwt->primary);
	const uint64_t *p = (const uint64_t*)bwt_occ_intv(bwt, k);
	return p[c] + __builtin_popcountll(p[4 + c] & occ_mask(k));
}

void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
{
	bwtint_t _k, _l;
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if ((bwt->is_onehot ? (_l>>OCC_ONEHOT_SHIFT != _k>>OCC_ONEHOT_SHIFT) : (_l>>OCC_INTV_SHIFT != _k>>OCC_INTV_SHIFT)) || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		*ok = bwt_occ(bwt, k, c);
		*ol = bwt_occ(bwt, l, c);
	} else {
		if (k >= bwt->primary) --k;
		if (l >= bwt->primary) --l;
		const uint64_t *p = (const uint64_t*)bwt_occ_intv(bwt, k);
		uint64_t bits = p[4 + c];
		*ok = p[c] + __builtin_popcountll(bits & occ_mask(k));
		*ol = p[c] + __builtin_popcountll(bits & occ_mask(l));
	}
}

void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	k -= (k >= bwt->primary);
	const uint64_t *p = (const uint64_t*)bwt_occ_intv(bwt, k);
	uint64_t mask = occ_mask(k);
	cnt[0] = p[0] + __builtin_popcountll(p[4] & mask);
	cnt[1] = p[1] + __builtin_popcountll(p[5] & mask);
	cnt[2] = p[2] + __builtin_popcountll(p[6] & mask);
	cnt[3] = p[3] + __builtin_popcountll(p[7] & mask);
}

void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	bwtint_t _k, _l;
	_k = k - (k >= bwt->primary);
	_l = l - (l >= bwt->primary);
	if ((bwt->is_onehot ? (_l>>OCC_ONEHOT_SHIFT != _k>>OCC_ONEHOT_SHIFT) : (_l>>OCC_INTV_SHIFT != _k>>OCC_INTV_SHIFT)) || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		k -= (k >= bwt->primary);
		l -= (l >= bwt->primary);
		const uint64_t *p = (const uint64_t*)bwt_occ_intv(bwt, k);
		uint64_t mk = occ_mask(k), ml = occ_mask(l);
		cntk[0] = p[0] + __builtin_popcountll(p[4] & mk);
		cntk[1] = p[1] + __builtin_popcountll(p[5] & mk);
		cntk[2] = p[2] + __builtin_popcountll(p[6] & mk);
		cntk[3] = p[3] + __builtin_popcountll(p[7] & mk);
		cntl[0] = p[0] + __builtin_popcountll(p[4] & ml);
		cntl[1] = p[1] + __builtin_popcountll(p[5] & ml);
		cntl[2] = p[2] + __builtin_popcountll(p[6] & ml);
		cntl[3] = p[3] + __builtin_popcountll(p[7] & ml);
	}
}

int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t k, l, ok, ol;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) break; // no match
	}
	if (k > l) return 0; // no match
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
	int i;
	bwtint_t k, l, ok, ol;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}

/*********************
 * Bidirectional BWT *
 *********************/

void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back)
{
	bwtint_t tk[4], tl[4];
	int i;
	bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
	ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
	ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
}

static void bwt_reverse_intvs(bwtintv_v *p)
{
	if (p->n > 1) {
		int j;
		for (j = 0; j < p->n>>1; ++j) {
			bwtintv_t tmp = p->a[p->n - 1 - j];
			p->a[p->n - 1 - j] = p->a[j];
			p->a[j] = tmp;
		}
	}
}
void bwt_extend_fwd(const bwt_t *restrict bwt, const bwtintv_t *restrict ik, bwtintv_t ok[restrict 4]);
void bwt_extend_back(const bwt_t *restrict bwt, const bwtintv_t *restrict ik, bwtintv_t ok[restrict 4]);

// NOTE: $max_intv is not currently used in BWA-MEM
int bwt_smem1a(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	int i, j, c, ret;
	bwtintv_t ik, ok[4];
	bwtintv_v a[2], *prev, *curr, *swap;

	mem->n = 0;
	if (q[x] > 3) return x + 1;
	if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
	kv_init(a[0]); kv_init(a[1]);
	prev = tmpvec && tmpvec[0]? tmpvec[0] : &a[0]; // use the temporary vector if provided
	curr = tmpvec && tmpvec[1]? tmpvec[1] : &a[1];
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (ik.x[2] < max_intv) { // an interval small enough
			kv_push(bwtintv_t, *curr, ik);
			break;
		} else if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend_fwd(bwt, &ik, ok);
			if (ok[c].x[2] != ik.x[2]) { // change of the interval size
				kv_push(bwtintv_t, *curr, ik);
				if (ok[c].x[2] < min_intv) break; // the interval size is too small to be extended further
			}
			ik = ok[c]; ik.info = i + 1;
		} else { // an ambiguous base
			kv_push(bwtintv_t, *curr, ik);
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}
	if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
	bwt_reverse_intvs(curr); // s.t. smaller intervals (i.e. longer matches) visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? -1 : q[i] < 4? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			bwtintv_t *p = &prev->a[j];
			if (c >= 0 && ik.x[2] >= max_intv) bwt_extend_back(bwt, p, ok);
			if (c < 0 || ik.x[2] < max_intv || ok[c].x[2] < min_intv) { // keep the hit if reaching the beginning or an ambiguous base or the intv is small enough
				if (curr->n == 0) { // test curr->n>0 to make sure there are no longer matches
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
						kv_push(bwtintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			} else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]) {
				ok[c].info = p->info;
				kv_push(bwtintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

	if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
	if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
	return ret;
}

int bwt_smem1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	return bwt_smem1a(bwt, len, q, x, min_intv, 0, mem, tmpvec);
}

int bwt_seed_strategy1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem)
{
	int i, c;
	bwtintv_t ik, ok[4];

	memset(mem, 0, sizeof(bwtintv_t));
	if (q[x] > 3) return x + 1;
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	for (i = x + 1; i < len; ++i) { // forward search
		if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend_fwd(bwt, &ik, ok);
			if (ok[c].x[2] < max_intv && i - x >= min_len) {
				*mem = ok[c];
				mem->info = (uint64_t)x<<32 | (i + 1);
				return i + 1;
			}
			ik = ok[c];
		} else return i + 1;
	}
	return len;
}

/*************************
 * Read/write BWT and SA *
 *************************/

void bwt_dump_bwt(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	err_fwrite(bwt->bwt, 4, bwt->bwt_size, fp);
	err_fflush(fp);
	err_fclose(fp);
}

void bwt_dump_sa(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	err_fwrite(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fwrite(&bwt->seq_len, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	err_fflush(fp);
	err_fclose(fp);
}

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = xopen(fn, "rb");
	err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	err_fread_noeof(skipped, sizeof(bwtint_t), 4, fp); // skip
	err_fread_noeof(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	fread_fix(fp, sizeof(bwtint_t) * (bwt->n_sa - 1), bwt->sa + 1);
	err_fclose(fp);
}

// Convert old interleaved BWT format (OCC_INTERVAL=128, 2-bit packed) to
// one-hot format (OCC_INTERVAL=64, per-base bitvectors) at load time.
// Old block layout: [4×u64 counts][8×u32 2-bit packed BWT data] = 64 bytes per 128 positions
// New block layout: [4×u64 counts][4×u64 one-hot bitvectors] = 64 bytes per 64 positions
void bwt_convert_to_onehot(bwt_t *bwt)
{
	bwtint_t i;
	bwtint_t old_intv = 128; // the file was written with OCC_INTERVAL=128
	bwtint_t n_occ_new = (bwt->seq_len + (1LL<<OCC_ONEHOT_SHIFT) - 1) / (1LL<<OCC_ONEHOT_SHIFT) + 1;
	bwtint_t new_size = n_occ_new * 16; // 16 uint32_t per block
	uint32_t *new_bwt = (uint32_t*)calloc(new_size, 4);

	// First, extract the raw BWT string (2-bit packed, no interleaving)
	// from the old interleaved format
	bwtint_t n_raw_words = (bwt->seq_len + 15) / 16;
	uint32_t *raw = (uint32_t*)calloc(n_raw_words, 4);
	{
		// Old block at position i*128: bwt->bwt[(i<<4)]
		// Counts at offset 0..3 (as uint64_t[4] = uint32_t[0..7])
		// BWT data at offset 8..15 (uint32_t[8..15], 8 words = 128 positions)
		bwtint_t n_old_blocks = (bwt->seq_len + old_intv - 1) / old_intv;
		bwtint_t raw_idx = 0;
		for (i = 0; i < n_old_blocks; ++i) {
			uint32_t *old_block = bwt->bwt + (i << 4);
			// BWT data starts after 4 × uint64_t = 8 × uint32_t
			uint32_t *bwt_data = old_block + 8;
			int n_words = 8; // 128/16
			bwtint_t remaining = bwt->seq_len - i * old_intv;
			if (remaining < old_intv) n_words = (remaining + 15) / 16;
			int j;
			for (j = 0; j < n_words && raw_idx < n_raw_words; ++j)
				raw[raw_idx++] = bwt_data[j];
		}
	}

	// Now build the one-hot format from the raw BWT
	{
		bwtint_t c[4] = {0, 0, 0, 0};
		bwtint_t block_idx = 0;
		for (i = 0; i < bwt->seq_len; ) {
			// Start of a new block
			uint64_t *out = (uint64_t*)(new_bwt + (block_idx << 4));
			// Store cumulative counts
			out[0] = c[0]; out[1] = c[1]; out[2] = c[2]; out[3] = c[3];
			// Build one-hot bitvectors for positions [i, i+64)
			uint64_t oh[4] = {0, 0, 0, 0};
			bwtint_t end = i + 64;
			if (end > bwt->seq_len) end = bwt->seq_len;
			bwtint_t j;
			for (j = i; j < end; ++j) {
				int base = (raw[j >> 4] >> ((~j & 15) << 1)) & 3;
				oh[base] |= 1ULL << (j - i);
				++c[base];
			}
			out[4] = oh[0]; out[5] = oh[1]; out[6] = oh[2]; out[7] = oh[3];
			++block_idx;
			i = end;
		}
		// Final block with just counts (boundary)
		if (block_idx < n_occ_new) {
			uint64_t *out = (uint64_t*)(new_bwt + (block_idx << 4));
			out[0] = c[0]; out[1] = c[1]; out[2] = c[2]; out[3] = c[3];
			out[4] = out[5] = out[6] = out[7] = 0;
		}
	}

	free(raw);
	free(bwt->bwt);
	bwt->bwt = new_bwt;
	bwt->bwt_size = new_size;
	bwt->is_onehot = 1;
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = xopen(fn, "rb");
	err_fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	err_fseek(fp, 0, SEEK_SET);
	err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	fread_fix(fp, bwt->bwt_size<<2, bwt->bwt);
	bwt->seq_len = bwt->L2[4];
	err_fclose(fp);
	bwt_gen_cnt_table(bwt); // needed for old format; one-hot conversion deferred

	return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}

// Specialized bwt_extend for backward search (is_back=1) — eliminates indexing overhead
void bwt_extend_back(const bwt_t *restrict bwt, const bwtintv_t *restrict ik, bwtintv_t ok[restrict 4])
{
	bwtint_t tk[4], tl[4];
	int i;
	bwt_2occ4(bwt, ik->x[0] - 1, ik->x[0] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}
	ok[3].x[1] = ik->x[1] + (ik->x[0] <= bwt->primary && ik->x[0] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
	ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
	ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
}

// Specialized bwt_extend for forward search (is_back=0)
void bwt_extend_fwd(const bwt_t *restrict bwt, const bwtintv_t *restrict ik, bwtintv_t ok[restrict 4])
{
	bwtint_t tk[4], tl[4];
	int i;
	bwt_2occ4(bwt, ik->x[1] - 1, ik->x[1] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}
	ok[3].x[0] = ik->x[0] + (ik->x[1] <= bwt->primary && ik->x[1] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
	ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
	ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
}
