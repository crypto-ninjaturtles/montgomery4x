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
#define mask32 0xffffffffUL
#define mask26 0x03ffffffUL
#define mask25 0x01ffffffUL
#define mask2625 mask26 | (mask25 << 32)
#define mask2526 mask25 | (mask26 << 32)

#define p25519_05 p25519_0(2) | (p25519_5(2) << 32)
#define p25519_16 p25519_1(2) | (p25519_6(2) << 32)
#define p25519_27 p25519_2(2) | (p25519_7(2) << 32)
#define p25519_38 p25519_3(2) | (p25519_8(2) << 32)
#define p25519_49 p25519_4(2) | (p25519_9(2) << 32)

#define shift2625 26UL | (25UL << 32)
#define shift2526 25UL | (26UL << 32)

ALIGN void scalar_mult_var_base(unsigned char *q, const unsigned char *n, const unsigned char *p, const unsigned char *A);
