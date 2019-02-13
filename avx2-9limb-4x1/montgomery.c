/*
 * This file is in public domain.
 *
 * Authors: Huseyin Hisil, Berkan Egrice, Mert Yassi
 * Title:   Fast 4 way vectorized ladder for the complete set of Montgomery curves
 * Version: 2019-02-13
 *
 */

#include "datatype.h"
#include "montgomery.h"
#include "thirdparty/pmp-inv-master/SL-DCC/1/include/p25519.h"

const ALIGN vec ptimes2[5] = {
	{ p25519x2_40, p25519x2_40, p25519x2_40, p25519x2_40 },
	{ p25519x2_51, p25519x2_51, p25519x2_51, p25519x2_51 },
	{ p25519x2_62, p25519x2_62, p25519x2_62, p25519x2_62 },
	{ p25519x2_73, p25519x2_73, p25519x2_73, p25519x2_73 },
	{ p25519x2_8_, p25519x2_8_, p25519x2_8_, p25519x2_8_ }
};
const ALIGN vec ptimes3[5] = {
	{ p25519x3_40, p25519x3_40, p25519x3_40, p25519x3_40 },
	{ p25519x3_51, p25519x3_51, p25519x3_51, p25519x3_51 },
	{ p25519x3_62, p25519x3_62, p25519x3_62, p25519x3_62 },
	{ p25519x3_73, p25519x3_73, p25519x3_73, p25519x3_73 },
	{ p25519x3_8_, p25519x3_8_, p25519x3_8_, p25519x3_8_ }
};

const ALIGN vec nineteen = { 19UL, 19UL, 19UL, 19UL };
const ALIGN vec vecmask01 = { mask100, mask302, mask504, mask706 };
const ALIGN vec vecmask10 = { mask504, mask706, mask100, mask302 };
const ALIGN vec vecmask32 = { mask32, mask32, mask32, mask32 };
const ALIGN vec vecmask29 = { mask29, mask29, mask29, mask29 };
const ALIGN vec vecmask28 = { mask28, mask28, mask28, mask28 };
const ALIGN vec vecmask2829 = { mask2829, mask2829, mask2829, mask2829 };
const ALIGN vec vecmask2928 = { mask2928, mask2928, mask2928, mask2928 };
const ALIGN vec vecmask2828 = { mask2828, mask2828, mask2828, mask2828 };
const ALIGN vec vecmask2832 = { mask2832, mask2832, mask2832, mask2832 };
ALIGN vec vecshft2800 = { shift28, shift28, shift28, shift28 };
ALIGN vec vecshft2829 = { shift2829, shift2829, shift2829, shift2829 };
ALIGN vec vecshft2928 = { shift2928, shift2928, shift2928, shift2928 };
ALIGN vec vecshft2828 = { shift2828, shift2828, shift2828, shift2828 };
ALIGN vec vecshft2832 = { shift2832, shift2832, shift2832, shift2832 };

ALIGN INLINE void mul19fast(vec *A){
	const vec vecmultfast19 = {
		19UL | (1UL << 32), 19UL | (1UL << 32),
		19UL | (1UL << 32), 19UL | (1UL << 32)
	};
	*A=VSHF(MULW(*A,vecmultfast19), _MM_SHUFFLE(2,3,0,1));
}

ALIGN INLINE void gfe4x_reducefast(vec *W){ /*W must have 9 digits. Works in squeezed form.*/
	vec vtemp;
	W[8]&=vecmask32;
	vtemp=SHRW(W[4],vecshft2829); W[4]&=vecmask2829; W[5]=VADD(W[5],vtemp);
	vtemp=SHRW(W[5],vecshft2828); W[5]&=vecmask2828; W[6]=VADD(W[6],vtemp);
	vtemp=SHRW(W[6],vecshft2928); W[6]&=vecmask2928; W[7]=VADD(W[7],vtemp);
	vtemp=SHRW(W[7],vecshft2829); W[7]&=vecmask2829; W[8]=VADD(W[8],vtemp);
	vtemp=SHRW(W[8],vecshft2800); W[8]&=vecmask2832;
	mul19fast(&vtemp);                               W[4]=VADD(W[4],vtemp);
	vtemp=SHRW(W[4],vecshft2829); W[4]&=vecmask2829; W[5]=VADD(W[5],vtemp);
}

ALIGN INLINE void gfe4x_reduce(vec *W){ /*W must have 11 digits.*/
	vec vtemp0, vtemp1;
	// [9] -> [10]              // [4] -> [5]
	vtemp1=VSHR(W[9],29);       vtemp0=VSHR(W[4],28);
	W[9]&=vecmask29;            W[4]&=vecmask28;
	W[10]=VADD(W[10],vtemp1);   W[5]=VADD(W[5],vtemp0);
	// [1] += [10]*19           // [0] += [9]*19
	VMUL19(W[10]);              W[9]=VMUL(nineteen,W[9]);
	W[1]=VADD(W[1],W[10]);      W[0]=VADD(W[0],W[9]);
	// [0] -> [1]               // [5] -> [6]
	vtemp1=VSHR(W[0],29);       vtemp0=VSHR(W[5],28);
	W[0]&=vecmask29;            W[5]&=vecmask28;
	W[1]=VADD(W[1],vtemp1);     W[6]=VADD(W[6],vtemp0);
	// [1] -> [2]               // [6] -> [7]
	vtemp1=VSHR(W[1],28);       vtemp0=VSHR(W[6],29);
	W[1]&=vecmask28;            W[6]&=vecmask29;
	W[2]=VADD(W[2],vtemp1);     W[7]=VADD(W[7],vtemp0);
	// [2] -> [3]               // [7] -> [8]
	vtemp1=VSHR(W[2],28);       vtemp0=VSHR(W[7],28);
	W[2]&=vecmask28;            W[7]&=vecmask28;
	W[3]=VADD(W[3],vtemp1);     W[8]=VADD(W[8],vtemp0);
	// [3] -> [4]               // [8] -> [0]
	vtemp1=VSHR(W[3],29);        vtemp0=VSHR(W[8],28);
	W[3]&=vecmask29;            W[8]&=vecmask28;
	                            VMUL19(vtemp0);
	W[4]=VADD(W[4],vtemp1);     W[0]=VADD(W[0],vtemp0);
	// [4] -> [5]               // [0] -> [1]
	vtemp1=VSHR(W[4],28);       vtemp0=VSHR(W[0],29);
	W[4]&=vecmask28;            W[0]&=vecmask29;
	W[5]=VADD(W[5],vtemp1);     W[1]=VADD(W[1],vtemp0);
}

ALIGN INLINE void gfe4x_reduce_sequential(vec *W){ /*W must have 9 digits.*/
	vec t;
	t=VSHR(W[0],29);   W[0]&=vecmask29;   W[1]=VADD(W[1],t);
	t=VSHR(W[1],28);   W[1]&=vecmask28;   W[2]=VADD(W[2],t);
	t=VSHR(W[2],28);   W[2]&=vecmask28;   W[3]=VADD(W[3],t);
	t=VSHR(W[3],29);   W[3]&=vecmask29;   W[4]=VADD(W[4],t);
	t=VSHR(W[4],28);   W[4]&=vecmask28;   W[5]=VADD(W[5],t);
	t=VSHR(W[5],28);   W[5]&=vecmask28;   W[6]=VADD(W[6],t);
	t=VSHR(W[6],29);   W[6]&=vecmask29;   W[7]=VADD(W[7],t);
	t=VSHR(W[7],28);   W[7]&=vecmask28;   W[8]=VADD(W[8],t);
	t=VSHR(W[8],28);   W[8]&=vecmask28;
	t=VMUL(t,nineteen);                   W[0]=VADD(W[0],t);
}

ALIGN INLINE void gfe4x_mul(vec *r, vec *m, vec *n){
	vec AA[3], BB[3], CC[3], DD[3], EE[3], FF[3];
	vec GG[5], HH[5], II[5], JJ[5], KK[5], LL[5];
	vec MM[5], NN[5], PP[5], RR[5], TT[5], UU[5];
	vec VV[5], XX[5], ZZ[5], SS[5], QQ[5], WW[5];
	vec YY[5], AB[5], CD[5], EF[5];
	vec *u0 = &m[0], *u1 = &m[3], *u2 = &m[6];
	vec *v0 = &n[0], *v1 = &n[3], *v2 = &n[6];
	vec w[11];

	mul_3_(GG, u0, v0);
	add_n_(AA, u0, u1, 3);
	add_n_(BB, u0, u2, 3);
	add_n_(CC, u1, u2, 3);
	add_n_(DD, v0, v1, 3);
	add_n_(EE, v0, v2, 3);
	add_n_(FF, v1, v2, 3);
	mul_3_(HH, u1, v1);
	mul_3_(II, u2, v2);
	mul_3_(JJ, AA, DD);
	mul_3_(KK, BB, EE);
	add_n_(MM, HH, II, 5);
	add_n_(NN, GG, II, 5);
	sub_n_(PP, JJ, NN, 5);
	mul_3_(LL, CC, FF);
	sub_n_(RR, KK, NN, 5);
	sub_n_(SS, LL, MM, 5);
	sub_n_(TT, PP, HH, 5);
	add_n_(UU, RR, HH, 5);
	shl_n_(QQ, SS, 4, 5);
	add_n_(VV, GG, SS, 5);
	add_n_(WW, SS, SS, 5);
	add_n_(ZZ, VV, WW, 5);
	add_n_(YY, ZZ, QQ, 5);
	shl_n_(XX, II, 2, 5);
	add_n_(AB, XX, II, 5);
	shl_n_(CD, AB, 2, 5);
	add_n_(EF, TT, CD, 5);

	w[0]  = YY[0];
	w[1]  = YY[1];
	w[2]  = YY[2];
	w[3]  = VADD(YY[3], EF[0]);
	w[4]  = VADD(YY[4], EF[1]);
	w[5]  = EF[2];
	w[6]  = VADD(EF[3], UU[0]);
	w[7]  = VADD(EF[4], UU[1]);
	w[8]  = UU[2];
	w[9]  = UU[3];
	w[10] = UU[4];

	gfe4x_reduce(w); /*Do not eliminate w with r! (w[11],r[9])*/

	r[0] = w[0];
	r[1] = w[1];
	r[2] = w[2];
	r[3] = w[3];
	r[4] = w[4];
	r[5] = w[5];
	r[6] = w[6];
	r[7] = w[7];
	r[8] = w[8];
}

ALIGN INLINE void gfe4x_squ(vec *r, vec *m){
	vec AA[5], BB[5], CC[5], DD[5], EE[5], FF[5], GG[3], HH[3];
	vec II[3], KK[5], LL[5], MM[5], NN[5], PP[5], RR[5], TT[5];
	vec UU[5], VV[5], XX[5], ZZ[5], SS[5], QQ[5], WW[5], YY[5];
	vec *u0 = &m[0], *u1 = &m[3], *u2 = &m[6];
	vec w[11];

	squ_3_(AA, u0);
	squ_3_(BB, u1);
	squ_3_(CC, u2);
	add_n_(GG, u1, u2, 3);
	add_n_(HH, u0, u1, 3);
	squ_3_(DD, GG);
	add_n_(II, u0, u2, 3);
	add_n_(KK, AA, CC, 5);
	shl_n_(PP, CC, 2, 5);
	shl_n_(RR, CC, 4, 5);
	add_n_(LL, RR, PP, 5);
	sub_n_(MM, KK, LL, 5);
	sub_n_(NN, DD, BB, 5);
	shl_n_(QQ, NN, 4, 5);
	shl_n_(SS, NN, 1, 5);
	add_n_(TT, QQ, SS, 5);
	add_n_(UU, TT, NN, 5);
	squ_3_(EE, HH);
	add_n_(VV, UU, MM, 5);
	sub_n_(WW, EE, BB, 5);
	squ_3_(FF, II);
	add_n_(YY, FF, BB, 5);
	sub_n_(XX, WW, MM, 5);
	sub_n_(ZZ, YY, KK, 5);

	w[0]  = VV[0];
	w[1]  = VADD(VV[1],VV[1]);
	w[2]  = VADD(VV[2],VV[2]);
	w[3]  = VADD(VSHL(VV[3],2),XX[0]);
	w[4]  = VADD(VV[4],VADD(XX[1],XX[1]));
	w[5]  = VADD(XX[2],XX[2]);
	w[6]  = VADD(ZZ[0],VSHL(XX[3],2));
	w[7]  = VADD(VADD(ZZ[1],ZZ[1]),XX[4]);
	w[8]  = VADD(ZZ[2],ZZ[2]);
	w[9]  = VSHL(ZZ[3],2);
	w[10] = ZZ[4];

	gfe4x_reduce(w); /*Do not eliminate w with r! (w[11],r[9])*/

	r[0] = w[0];
	r[1] = w[1];
	r[2] = w[2];
	r[3] = w[3];
	r[4] = w[4];
	r[5] = w[5];
	r[6] = w[6];
	r[7] = w[7];
	r[8] = w[8];
}

ALIGN INLINE void gfe4x_squeeze(vec *ip) {
	ip[4]^=VSHL(ip[0],32);
	ip[5]^=VSHL(ip[1],32);
	ip[6]^=VSHL(ip[2],32);
	ip[7]^=VSHL(ip[3],32);
}

ALIGN INLINE void gfe4x_unsqueeze(vec *ip) {
	ip[0]=VSHR(ip[4],32);
	ip[1]=VSHR(ip[5],32);
	ip[2]=VSHR(ip[6],32);
	ip[3]=VSHR(ip[7],32);
}

ALIGN INLINE void gfe4x_mask32(vec *ip) {
	ip[4]&=vecmask32;
	ip[5]&=vecmask32;
	ip[6]&=vecmask32;
	ip[7]&=vecmask32;
}

ALIGN INLINE void gfe4x_swap(vec *ip, int swap) {
	int i;
	swap=-swap;
	vec vc = _mm256_set_epi32(swap, swap, swap, swap, swap, swap, swap, swap);
	vec perm_mask = _mm256_blendv_epi8(vecmask01, vecmask10, vc);
	for(i = 4; i < NLIMB; i++){
		ip[i] = _mm256_permutevar8x32_epi32(ip[i], perm_mask);
	}
}

ALIGN INLINE void gfe4x_align(vec *op, vec *ip) {
	int i;
	for(i = 0; i < NLIMB; i++){
		op[i] = _mm256_permute4x64_epi64(ip[i], _MM_SHUFFLE(3,2,2,3));
	}
}

ALIGN INLINE void gfe4x_hadamard(vec *ip) {
	vec temp1[NLIMB], temp2[NLIMB], temp3[NLIMB];
	int i;
	for(i = 4; i < NLIMB; i++){
		temp1[i] = VSHF(ip[i],_MM_SHUFFLE(1,0,3,2));
		temp3[i] = VSUB(ptimes3[i-4], ip[i]);
		temp2[i] = VBLD(ip[i], temp3[i], _MM_SHUFFLE(3,0,3,0));
		ip[i] = VADD(temp1[i], temp2[i]);
	}
}

ALIGN INLINE void gfe4x_premix(vec *z1x1x2A, vec *x5z5z2z4, vec *z1x1_A, vec *x3z3x2z2, vec *x5z5x4z4){
	vec x3z3z2x2[NLIMB];
	int i;
	for(i = 4; i < NLIMB; i++){
		z1x1x2A[i] = VBLD(z1x1_A[i], x3z3x2z2[i], 48);
		x3z3z2x2[i] = VSHF(x3z3x2z2[i],_MM_SHUFFLE(1,0,3,2));
		x5z5z2z4[i] = VBLD(x5z5x4z4[i], x3z3z2x2[i], 48);
	}
}

ALIGN INLINE void gfe4x_postmix(vec *x3z3x2z2, vec *x7z7x6z6, vec *x5z5x4z4){
	vec x9z9x8z8[NLIMB], x1z1x0z0[NLIMB];
	int i;
	for(i = 4; i < NLIMB; i++){
		x9z9x8z8[i] = VADD(x7z7x6z6[i], x7z7x6z6[i]);
		x9z9x8z8[i] = VSHF(x9z9x8z8[i],_MM_SHUFFLE(1,0,3,2));
		x1z1x0z0[i] = VSHF(x5z5x4z4[i],_MM_SHUFFLE(1,0,3,2));
		x9z9x8z8[i] = VADD(x9z9x8z8[i], x7z7x6z6[i]);
		x1z1x0z0[i] = VSUB(VADD(x5z5x4z4[i], ptimes2[i-4]), x1z1x0z0[i]);
		x3z3x2z2[i] = VBLD(x7z7x6z6[i], x1z1x0z0[i], 48);
		x3z3x2z2[i] = VBLD(x3z3x2z2[i], x9z9x8z8[i], 192);
	}
}

ALIGN INLINE void ladder_4x1(vec *x3z3x2z2, vec *z1x1_A, const unsigned char *n) {
	vec z1x1x2A[NLIMB], z2x2x2z2[NLIMB], x5z5x4z4[NLIMB], x5z5z2z4[NLIMB], x7z7x6z6[NLIMB];
	int i, j, bit, swap, prevbit = 0;
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
			gfe4x_reducefast(x3z3x2z2);
			gfe4x_unsqueeze(x3z3x2z2);
			gfe4x_align(z2x2x2z2, x3z3x2z2);
			gfe4x_mul(x3z3x2z2, x3z3x2z2, z2x2x2z2);
			gfe4x_squeeze(x3z3x2z2);
			gfe4x_hadamard(x3z3x2z2);
			gfe4x_reducefast(x3z3x2z2);
			gfe4x_unsqueeze(x3z3x2z2);
			gfe4x_squ(x5z5x4z4, x3z3x2z2);
			gfe4x_squeeze(x5z5x4z4);
			gfe4x_premix(z1x1x2A, x5z5z2z4, z1x1_A, x3z3x2z2, x5z5x4z4);
			gfe4x_unsqueeze(x5z5z2z4);
			gfe4x_unsqueeze(z1x1x2A);
			gfe4x_mul(x7z7x6z6, z1x1x2A, x5z5z2z4);
			gfe4x_squeeze(x7z7x6z6);
			gfe4x_postmix(x3z3x2z2, x7z7x6z6, x5z5x4z4);
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
	gfe25519 x64[1], z64[1], t[1], temp[1];
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

	ladder_4x1(x3z3x2z2, z1x1_A, n);

	get_channel(t0, x3z3x2z2, 2);
	get_channel(t1, x3z3x2z2, 3);
	reduce_unique(t0);
	reduce_unique(t1);
	gfe64_pack((gfe64 *)x64, t0);
	gfe64_pack((gfe64 *)z64, t1);

	gfp25519inv(temp, z64);
	gfp25519mul(t,temp,x64);
	gfp25519reduce(t);
	gfp25519makeunique(t);

	for (i = 0; i < 4; i++) {
		((u64 *)q)[i] = t->l[i];
	}
}
