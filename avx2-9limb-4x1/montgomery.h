/*
 * This file is in public domain.
 *
 * Authors: Huseyin Hisil, Berkan Egrice, Mert Yassi
 * Title:   Fast 4 way vectorized ladder for the complete set of Montgomery curves
 * Version: 2019-02-13
 *
 */

#define VADD _mm256_add_epi64
#define VSUB _mm256_sub_epi64
#define VMUL _mm256_mul_epu32
#define VSHL _mm256_slli_epi64
#define VSHR _mm256_srli_epi64
#define VBLD _mm256_blend_epi32
#define VSHF _mm256_shuffle_epi32
#define SHRW _mm256_srlv_epi32
#define MULW _mm256_mullo_epi32

#define mask100 0x100000000UL
#define mask302 0x300000002UL
#define mask504 0x500000004UL
#define mask706 0x700000006UL
#define mask32 0x0ffffffffUL
#define mask29 0x01fffffffUL
#define mask28 0x00fffffffUL
#define mask2829 0x0fffffffUL | (0x1fffffffUL << 32)
#define mask2928 0x1fffffffUL | (0x0fffffffUL << 32)
#define mask2828 0x0fffffffUL | (0x0fffffffUL << 32)
#define mask2832 0x0fffffffUL | (0xffffffffUL << 32)

#define p25519x2_40  p25519_4(2) | (p25519_0(2) << 32)
#define p25519x2_51  p25519_5(2) | (p25519_1(2) << 32)
#define p25519x2_62  p25519_6(2) | (p25519_2(2) << 32)
#define p25519x2_73  p25519_7(2) | (p25519_3(2) << 32)
#define p25519x2_8_  p25519_8(2)

#define p25519x3_40  p25519_4(3) | (p25519_0(3) << 32)
#define p25519x3_51  p25519_5(3) | (p25519_1(3) << 32)
#define p25519x3_62  p25519_6(3) | (p25519_2(3) << 32)
#define p25519x3_73  p25519_7(3) | (p25519_3(3) << 32)
#define p25519x3_8_  p25519_8(3)

#define shift28 28UL
#define shift2829 28UL | (29UL << 32)
#define shift2832 28UL | (32UL << 32)
#define shift2928 29UL | (28UL << 32)
#define shift2828 28UL | (28UL << 32)

#define VMUL19(A){ A=VADD(VADD(VSHL(A,4),VSHL(A,1)),A); }

#define add_n_(Z,A,B,n){               \
	do{                                \
		int i;                         \
		for(i=0; i < n; i++){          \
			Z[i]=VADD(B[i],A[i]);      \
		}                              \
	}while(0);                         \
}

#define sub_n_(Z,A,B,n){               \
	do{                                \
		int i;                         \
		for(i=0; i < n; i++){          \
			Z[i]=VSUB(A[i],B[i]);      \
		}                              \
	}while(0);                         \
}

#define shl_n_(Z,A,c,n){               \
	do{                                \
		int i;                         \
		for(i=0; i < n; i++){          \
			Z[i]=VSHL(A[i],c);         \
		}                              \
	}while(0);                         \
}

#define mul_3_(Z,A,B){                 \
	do{                                \
		vec t00, t01, t03, t04, t05;   \
		vec t06, t07, t08, t09;        \
		Z[0]=VMUL(A[0],B[0]);          \
		Z[4]=VMUL(A[2],B[2]);          \
		t04=VMUL(B[0],A[1]);           \
		t03=VMUL(B[1],A[0]);           \
		Z[1]=VADD(t04,t03);            \
		t00=VADD(A[1],A[1]);           \
		t01=VADD(A[2],A[2]);           \
		t05=VMUL(A[0],B[2]);           \
		t06=VMUL(B[0],A[2]);           \
		t07=VMUL(t00,B[1]);            \
		t08=VMUL(t00,B[2]);            \
		t09=VMUL(B[1],t01);            \
		Z[2]=VADD(t07,VADD(t05,t06));  \
		Z[3]=VADD(t08,t09);            \
	}while(0);                         \
}

#define squ_3_(Z,A){                                 \
	do{                                              \
		Z[0]=VMUL(A[0],A[0]);                        \
		Z[1]=VMUL(A[0],A[1]);                        \
		Z[2]=VADD(VMUL(A[0],A[2]),VMUL(A[1],A[1]));  \
		Z[3]=VMUL(A[1],A[2]);                        \
		Z[4]=VMUL(A[2],A[2]);                        \
	}while(0);                                       \
}

ALIGN void scalar_mult_var_base(unsigned char *q, const unsigned char *n, const unsigned char *p, const unsigned char *A);
