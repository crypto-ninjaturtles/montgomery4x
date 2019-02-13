/*
 * This file is in public domain.
 *
 * Authors: Huseyin Hisil, Berkan Egrice, Mert Yassi
 * Title:   Fast 4 way vectorized ladder for the complete set of Montgomery curves
 * Version: 2019-02-13
 *
 */

#include "thirdparty/pmp-inv-master/SL-DCC/1/include/p25519.h"
#include "thirdparty/hp-ecc-vec-master/include/faz_fp_avx2.h"
#include "thirdparty/hp-ecc-vec-master/include/simd_avx2.h"
#include "thirdparty/hp-ecc-vec-master/src/eltfp25519_4w_redradix.c"

#include "datatype.h"
#include "montgomery.h"

void intmul_sb(vec *C, vec *A, vec *B);
void intsqr(vec *C);
void compress(vec *C);

const ALIGN vec vecp25519[5] = {
	{p25519_05, p25519_05, p25519_05, p25519_05},
	{p25519_16, p25519_16, p25519_16, p25519_16},
	{p25519_27, p25519_27, p25519_27, p25519_27},
	{p25519_38, p25519_38, p25519_38, p25519_38},
	{p25519_49, p25519_49, p25519_49, p25519_49}
};

const ALIGN vec vecmask01 = { mask100, mask302, mask504, mask706 };
const ALIGN vec vecmask10 = { mask504, mask706, mask100, mask302 };
const ALIGN vec vecmask32 = { mask32, mask32, mask32, mask32 };
const ALIGN vec vecmask26 = { mask26, mask26, mask26, mask26 };
const ALIGN vec vecmask25 = { mask25, mask25, mask25, mask25 };
const ALIGN vec vecmask2625 = { mask2625, mask2625, mask2625, mask2625 };
const ALIGN vec vecmask2526 = { mask2526, mask2526, mask2526, mask2526 };
vec vecshft2625 = { shift2625, shift2625, shift2625, shift2625 };
vec vecshft2526 = { shift2526, shift2526, shift2526, shift2526 };

ALIGN INLINE void mul19fast(vec *A){
	vec vecmultfast19 = {
		1UL | (19UL << 32), 1UL | (19UL << 32),
		1UL | (19UL << 32), 1UL | (19UL << 32)
	};
	*A=VSHF(MULW(*A,vecmultfast19), _MM_SHUFFLE(2,3,0,1));
}

ALIGN INLINE void gfe4x_reducefast(vec *c){
	vec h05, h16, h27, h38, h49;
	h05=SHRW(c[0],vecshft2625); c[0]&=vecmask2625; c[1]=VADD(c[1],h05);
	h16=SHRW(c[1],vecshft2526); c[1]&=vecmask2526; c[2]=VADD(c[2],h16);
	h27=SHRW(c[2],vecshft2625); c[2]&=vecmask2625; c[3]=VADD(c[3],h27);
	h38=SHRW(c[3],vecshft2526); c[3]&=vecmask2526; c[4]=VADD(c[4],h38);
	h49=SHRW(c[4],vecshft2625); c[4]&=vecmask2625;
	mul19fast(&h49);                               c[0]=VADD(c[0],h49);
	h05=SHRW(c[0],vecshft2625); c[0]&=vecmask2625; c[1]=VADD(c[1],h05);
}

ALIGN INLINE void gfe4x_reduce_sequential(vec *c){
	const ALIGN vec nineteen = { 19UL, 19UL, 19UL, 19UL };
	vec h;
	h=VSHR(c[0],26); c[0]&=vecmask26; c[1]=VADD(c[1],h);
	h=VSHR(c[1],25); c[1]&=vecmask25; c[2]=VADD(c[2],h);
	h=VSHR(c[2],26); c[2]&=vecmask26; c[3]=VADD(c[3],h);
	h=VSHR(c[3],25); c[3]&=vecmask25; c[4]=VADD(c[4],h);
	h=VSHR(c[4],26); c[4]&=vecmask26; c[5]=VADD(c[5],h);
	h=VSHR(c[5],25); c[5]&=vecmask25; c[6]=VADD(c[6],h);
	h=VSHR(c[6],26); c[6]&=vecmask26; c[7]=VADD(c[7],h);
	h=VSHR(c[7],25); c[7]&=vecmask25; c[8]=VADD(c[8],h);
	h=VSHR(c[8],26); c[8]&=vecmask26; c[9]=VADD(c[9],h);
	h=VSHR(c[9],25); c[9]&=vecmask25;
	h=VMUL(h,nineteen);               c[0]=VADD(c[0],h);
}

ALIGN INLINE void gfe4x_mul(vec *C, vec *A, vec *B){
	vec a[NLIMB], b[NLIMB], c[NLIMB];

	a[0] = A[0];   b[0] = B[0];
	a[2] = A[1];   b[2] = B[1];
    a[4] = A[2];   b[4] = B[2];
    a[6] = A[3];   b[6] = B[3];
	a[8] = A[4];   b[8] = B[4];
	a[1] = A[5];   b[1] = B[5];
    a[3] = A[6];   b[3] = B[6];
	a[5] = A[7];   b[5] = B[7];
	a[7] = A[8];   b[7] = B[8];
	a[9] = A[9];   b[9] = B[9];

	intmul_sb(c,a,b);
	compress(c);

	C[0] = c[0];
	C[1] = c[2];
	C[2] = c[4];
	C[3] = c[6];
	C[4] = c[8];
	C[5] = c[1];
	C[6] = c[3];
	C[7] = c[5];
	C[8] = c[7];
	C[9] = c[9];
}

ALIGN INLINE void gfe4x_squ(vec *C, vec *A){
	vec a[NLIMB];

	a[0] = A[0];
	a[2] = A[1];
    a[4] = A[2];
    a[6] = A[3];
	a[8] = A[4];
	a[1] = A[5];
    a[3] = A[6];
	a[5] = A[7];
	a[7] = A[8];
	a[9] = A[9];

	intsqr(a);
	compress(a);

	C[0] = a[0];
	C[1] = a[2];
	C[2] = a[4];
	C[3] = a[6];
	C[4] = a[8];
	C[5] = a[1];
	C[6] = a[3];
	C[7] = a[5];
	C[8] = a[7];
	C[9] = a[9];
}

ALIGN INLINE void gfe4x_squeeze(vec *ip) {
	int i;
	for(i = 0; i < 5; i++){
		ip[i] = VSHL(ip[i+5], 32) ^ ip[i];
	}
}

ALIGN INLINE void gfe4x_unsqueeze(vec *ip) {
	int i;
	for(i = 0; i < 5; i++){
		ip[i+5] = VSHR(ip[i],32);
	}
}

ALIGN INLINE void gfe4x_mask32(vec *ip) {
	int i;
	for(i = 0; i < 5; i++){
		ip[i] = ip[i] & vecmask32;
	}
}

ALIGN INLINE void gfe4x_swap(vec *ip, int swap) {
	swap=-swap; int i;
	vec vc = _mm256_set_epi32(swap, swap, swap, swap, swap, swap, swap, swap);
	vec perm_mask = _mm256_blendv_epi8(vecmask01, vecmask10, vc);
	for(i = 0; i < 5; i++){
		ip[i] = _mm256_permutevar8x32_epi32(ip[i], perm_mask);
	}
}

ALIGN INLINE void gfe4x_align(vec *op, vec *ip) {
	int i;
	for(i = 0; i < 5; i++){
		op[i] = _mm256_permute4x64_epi64(ip[i], _MM_SHUFFLE(3,2,2,3));
	}
}

ALIGN INLINE void gfe4x_hadamard(vec *ip) {
	vec temp1, temp2, temp3;
	int i;
	for(i = 0; i < 5; i++){
		temp1 = VSHF(ip[i],_MM_SHUFFLE(1,0,3,2));
		temp3 = VSUB(vecp25519[i], ip[i]);
		temp2 = VBLD(ip[i], temp3, _MM_SHUFFLE(3,0,3,0));
		ip[i] = VADD(temp1, temp2);
	}
}

ALIGN INLINE void gfe4x_premix(vec *z1x1x2A, vec *x5z5z2z4, vec *z1x1_A, vec *x3z3x2z2, vec *x5z5x4z4){
	int i;
	for(i = 0; i < 5; i++){
		z1x1x2A[i] = VBLD(z1x1_A[i], x3z3x2z2[i], 48);
		x3z3x2z2[i] = VSHF(x3z3x2z2[i],_MM_SHUFFLE(1,0,3,2));
		x5z5z2z4[i] = VBLD(x5z5x4z4[i], x3z3x2z2[i], 48);
	}
}

ALIGN INLINE void postmix(vec *x3z3x2z2, vec *x7z7x6z6, vec *x5z5x4z4){
	vec x9z9x8z8[NLIMB], x1z1x0z0[NLIMB];
	int i;
	for (i = 0; i < 5; i++) {
		x9z9x8z8[i] = VADD(x7z7x6z6[i], x7z7x6z6[i]);
		x9z9x8z8[i] = VSHF(x9z9x8z8[i],_MM_SHUFFLE(1,0,3,2));
		x1z1x0z0[i] = VSHF(x5z5x4z4[i],_MM_SHUFFLE(1,0,3,2));
		x9z9x8z8[i] = VADD(x9z9x8z8[i], x7z7x6z6[i]);
		x1z1x0z0[i] = VSUB(VADD(x5z5x4z4[i], vecp25519[i]), x1z1x0z0[i]);
		x3z3x2z2[i] = VBLD(x7z7x6z6[i], x1z1x0z0[i], 48);
		x3z3x2z2[i] = VBLD(x3z3x2z2[i], x9z9x8z8[i], 192);
	}
}

ALIGN INLINE void ladder_4x1(vec *x3z3x2z2, vec *z1x1_A, const unsigned char *n) {
	int i, j, bit, swap, prevbit = 0;
	vec x5z5x4z4[NLIMB], z2x2x2z2[NLIMB], z1x1x2A[NLIMB], x5z5z2z4[NLIMB], x7z7x6z6[NLIMB];
	gfe4x_squeeze(x3z3x2z2);
	gfe4x_squeeze(z1x1_A);
	j = 6;
	for(i = 31; i >= 0; i--){
		for(; j >= 0; j--){
			bit = (n[i] >> j) & 1;
			swap = bit ^ prevbit;
			prevbit = bit;
			gfe4x_swap(x3z3x2z2, swap);
			gfe4x_hadamard(x3z3x2z2);
			gfe4x_align(z2x2x2z2, x3z3x2z2);
			gfe4x_unsqueeze(z2x2x2z2);
			gfe4x_unsqueeze(x3z3x2z2);
			gfe4x_mul(x3z3x2z2, x3z3x2z2, z2x2x2z2);
			gfe4x_squeeze(x3z3x2z2);
			gfe4x_hadamard(x3z3x2z2);
			gfe4x_unsqueeze(x3z3x2z2);
			gfe4x_squ(x5z5x4z4, x3z3x2z2);
			gfe4x_squeeze(x5z5x4z4);
			gfe4x_premix(z1x1x2A, x5z5z2z4, z1x1_A, x3z3x2z2, x5z5x4z4);
			gfe4x_unsqueeze(x5z5z2z4);
			gfe4x_unsqueeze(z1x1x2A);
			gfe4x_mul(x7z7x6z6, z1x1x2A, x5z5z2z4);
			gfe4x_squeeze(x7z7x6z6);
			postmix(x3z3x2z2, x7z7x6z6, x5z5x4z4);
			gfe4x_reducefast(x3z3x2z2);
		}
		j = 7;
	}
	gfe4x_swap(x3z3x2z2, bit);
	gfe4x_unsqueeze(x3z3x2z2);
	gfe4x_mask32(x3z3x2z2);
	gfe4x_reduce_sequential(x3z3x2z2);
	gfe4x_reduce_sequential(x3z3x2z2); //Do not remove this 2nd reduction.
}

ALIGN void scalar_mult_var_base(unsigned char *q, const unsigned char *n, const unsigned char *p, const unsigned char *A){
	vec z1x1_A[NLIMB], x3z3x2z2[NLIMB];
	gfe25519 x64[1], z64[1], t[1], inv[1];
	gfe t0[1], t1[1];
	unsigned char k[32], x1[32];
	int i;

	//Change scalar. Do clamping as in curve25519.
	for(i = 0; i < 32; i++){ k[i] = n[i]; }
	k[0] &= 248; k[31] &= 127; k[31] |= 64;
	for(i = 0; i < 32; i++){ x1[i] = p[i]; }
	x1[31] &= 127;

	char_to_gfe(t0, x1);
	char_to_gfe(t1, A);
	set_base(z1x1_A, t0, t1);
	set_vector(x3z3x2z2, t0);

	ladder_4x1(x3z3x2z2, z1x1_A, k);

	get_channel(t0, x3z3x2z2, 2);
	get_channel(t1, x3z3x2z2, 3);
	reduce_unique(t0);
	reduce_unique(t1);
	gfe64_pack((gfe64 *)x64, t0);
	gfe64_pack((gfe64 *)z64, t1);

	gfp25519inv(inv, z64);
	gfp25519mul(t,inv,x64);
	gfp25519reduce(t);
	gfp25519makeunique(t);

	for (i = 0; i < 4; i++) {
		((u64 *)q)[i] = t->l[i];
	}
}
