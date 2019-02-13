/**
 * Copyright 2017 Armando Faz Hernández
 * This file is part of faz_ecc_avx2.
 *
 * faz_ecc_avx2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * faz_ecc_avx2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with faz_ecc_avx2.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * Modified by HEY.
 * Date 07/02/2019.
 * Notes:
 * 1.Unused parts are eliminated.
 * 2.The body of functions that are in use are intact.
 **/


#define MUL19(X) ADD(X,ADD(ADD(X,X),SHL(X,4)))

inline void intmul_sb(argElement_4w C, argElement_4w A, argElement_4w B) {
    const __m256i times19 = SET164(19);
    __m256i c0,c1,c2,c3,c4,c5,c6,c7,c8,c9;
    __m256i    d1,d2,d3,d4,d5,d6,d7,d8,d9;
    __m256i b0,b1,b2,b3,b4,b5,b6,b7,b8,b9;
    __m256i ai,a2i;

    b0 = B[0];	b5 = B[1];
    b1 = B[2];	b6 = B[3];
    b2 = B[4];	b7 = B[5];
    b3 = B[6];	b8 = B[7];
    b4 = B[8];	b9 = B[9];

	ai = A[0];
	c0 = MUL(ai,b0);			c5 = MUL(ai,b5);
	c1 = MUL(ai,b1);			c6 = MUL(ai,b6);
	c2 = MUL(ai,b2);			c7 = MUL(ai,b7);
	c3 = MUL(ai,b3);			c8 = MUL(ai,b8);
	c4 = MUL(ai,b4);			c9 = MUL(ai,b9);

    ai = A[2];
    a2i = ADD(A[2],A[2]);
    d9 = MUL(b9,times19);
    c0 = ADD(c0,MUL(a2i,d9));	c5 = ADD(c5,MUL(ai, b4));
    c1 = ADD(c1,MUL(ai, b0));	c6 = ADD(c6,MUL(a2i,b5));
    c2 = ADD(c2,MUL(a2i,b1));	c7 = ADD(c7,MUL(ai, b6));
    c3 = ADD(c3,MUL(ai, b2));	c8 = ADD(c8,MUL(a2i,b7));
    c4 = ADD(c4,MUL(a2i,b3));	c9 = ADD(c9,MUL(ai, b8));

    ai = A[4];
	d8 = MUL(b8,times19);
	c0 = ADD(c0,MUL(ai,d8));	c5 = ADD(c5,MUL(ai,b3));
	c1 = ADD(c1,MUL(ai,d9));	c6 = ADD(c6,MUL(ai,b4));
	c2 = ADD(c2,MUL(ai,b0));	c7 = ADD(c7,MUL(ai,b5));
	c3 = ADD(c3,MUL(ai,b1));	c8 = ADD(c8,MUL(ai,b6));
	c4 = ADD(c4,MUL(ai,b2));	c9 = ADD(c9,MUL(ai,b7));

	ai = A[6];
	a2i = ADD(A[6],A[6]);
	d7 = MUL(b7,times19);
	c0 = ADD(c0,MUL(a2i,d7));	c5 = ADD(c5,MUL(ai,b2));
	c1 = ADD(c1,MUL(ai,d8));	c6 = ADD(c6,MUL(a2i,b3));
	c2 = ADD(c2,MUL(a2i,d9));	c7 = ADD(c7,MUL(ai,b4));
	c3 = ADD(c3,MUL(ai,b0));	c8 = ADD(c8,MUL(a2i,b5));
	c4 = ADD(c4,MUL(a2i,b1));	c9 = ADD(c9,MUL(ai,b6));

	ai = A[8];
	d6 = MUL(b6,times19);
	c0 = ADD(c0,MUL(ai,d6));	c5 = ADD(c5,MUL(ai,b1));
	c1 = ADD(c1,MUL(ai,d7));	c6 = ADD(c6,MUL(ai,b2));
	c2 = ADD(c2,MUL(ai,d8));	c7 = ADD(c7,MUL(ai,b3));
	c3 = ADD(c3,MUL(ai,d9));	c8 = ADD(c8,MUL(ai,b4));
	c4 = ADD(c4,MUL(ai,b0));	c9 = ADD(c9,MUL(ai,b5));

	ai = A[1];
	a2i = ADD(A[1],A[1]);
	d5 = MUL(b5,times19);
	c0 = ADD(c0,MUL(a2i,d5));	c5 = ADD(c5,MUL(ai,b0));
	c1 = ADD(c1,MUL(ai,d6));	c6 = ADD(c6,MUL(a2i,b1));
	c2 = ADD(c2,MUL(a2i,d7));	c7 = ADD(c7,MUL(ai,b2));
	c3 = ADD(c3,MUL(ai,d8));	c8 = ADD(c8,MUL(a2i,b3));
	c4 = ADD(c4,MUL(a2i,d9));	c9 = ADD(c9,MUL(ai,b4));

	ai = A[3];
	d4 = MUL(b4,times19);
	c0 = ADD(c0,MUL(ai,d4));	c5 = ADD(c5,MUL(ai,d9));
	c1 = ADD(c1,MUL(ai,d5));	c6 = ADD(c6,MUL(ai,b0));
	c2 = ADD(c2,MUL(ai,d6));	c7 = ADD(c7,MUL(ai,b1));
	c3 = ADD(c3,MUL(ai,d7));	c8 = ADD(c8,MUL(ai,b2));
	c4 = ADD(c4,MUL(ai,d8));	c9 = ADD(c9,MUL(ai,b3));

	ai = A[5];
	a2i = ADD(A[5],A[5]);
	d3 = MUL(b3,times19);
	c0 = ADD(c0,MUL(a2i,d3));	c5 = ADD(c5,MUL(ai,d8));
	c1 = ADD(c1,MUL(ai,d4));	c6 = ADD(c6,MUL(a2i,d9));
	c2 = ADD(c2,MUL(a2i,d5));	c7 = ADD(c7,MUL(ai,b0));
	c3 = ADD(c3,MUL(ai,d6));	c8 = ADD(c8,MUL(a2i,b1));
	c4 = ADD(c4,MUL(a2i,d7));	c9 = ADD(c9,MUL(ai,b2));

	ai = A[7];
	d2 = MUL(b2,times19);
	c0 = ADD(c0,MUL(ai,d2));	c5 = ADD(c5,MUL(ai,d7));
	c1 = ADD(c1,MUL(ai,d3));	c6 = ADD(c6,MUL(ai,d8));
	c2 = ADD(c2,MUL(ai,d4));	c7 = ADD(c7,MUL(ai,d9));
	c3 = ADD(c3,MUL(ai,d5));	c8 = ADD(c8,MUL(ai,b0));
	c4 = ADD(c4,MUL(ai,d6));	c9 = ADD(c9,MUL(ai,b1));

	ai = A[9];
	a2i = ADD(A[9],A[9]);
	d1 = MUL(b1,times19);
	c0 = ADD(c0,MUL(a2i,d1));	c5 = ADD(c5,MUL(ai,d6));
	c1 = ADD(c1,MUL(ai,d2));	c6 = ADD(c6,MUL(a2i,d7));
	c2 = ADD(c2,MUL(a2i,d3));	c7 = ADD(c7,MUL(ai,d8));
	c3 = ADD(c3,MUL(ai,d4));	c8 = ADD(c8,MUL(a2i,d9));
	c4 = ADD(c4,MUL(a2i,d5));	c9 = ADD(c9,MUL(ai,b0));

    C[0] = c0;	C[1] = c5;
    C[2] = c1;	C[3] = c6;
    C[4] = c2;	C[5] = c7;
    C[6] = c3;	C[7] = c8;
    C[8] = c4;	C[9] = c9;
}

inline void intsqr(argElement_4w c) {
	const __m256i times19 = SET164(19);
	__m256i c0,c1,c2,c3,c4,c5,c6,c7,c8,c9;
	__m256i          d3,d4,d5,d6,d7,d8,d9;
	__m256i b0,b1,b2,b3,b4,b5,b6,b7,b8,b9;
	__m256i ai, aj, ak;

	b0 = c[0];	b5 = c[1];
	b1 = c[2];	b6 = c[3];
	b2 = c[4];	b7 = c[5];
	b3 = c[6];	b8 = c[7];
	b4 = c[8];	b9 = c[9];

	ai = c[0];
	aj = ADD(ai,ai);
	c0 = MUL(ai,b0);				c5 = MUL(ai,b5);
	c1 = MUL(aj,b1);				c6 = MUL(aj,b6);
	c2 = MUL(aj,b2);				c7 = MUL(aj,b7);
	c3 = MUL(aj,b3);				c8 = MUL(aj,b8);
	c4 = MUL(aj,b4);				c9 = MUL(aj,b9);

	ai = c[2];
	aj = ADD(ai,ai);
	ak = ADD(aj, aj);
	d9 = MUL(b9,times19);
	c0 = ADD(c0,MUL(ak,d9));		c5 = ADD(c5,MUL(aj,b4));
	/*c1 = ADD(c1,MUL(ai,b0));		c6 = ADD(c6,MUL(aj,b5));*/
	c2 = ADD(c2,MUL(aj,b1));		c7 = ADD(c7,MUL(ai,b6));
	c3 = ADD(c3,MUL(aj,b2));		c8 = ADD(c8,MUL(ak,b7));
	c4 = ADD(c4,MUL(ak,b3));		c9 = ADD(c9,MUL(aj,b8));

	ai = c[4];
	aj = ADD(ai,ai);
	ak = ADD(aj,aj);
	d8 = MUL(b8,times19);
	c0 = ADD(c0,MUL(aj,d8));		c5 = ADD(c5,MUL(aj,b3));
	c1 = ADD(c1,MUL(aj,d9));		c6 = ADD(c6,MUL(aj,b4));
	/*c2 = ADD(c2,MUL(ai,b0));		c7 = ADD(c7,MUL(ai,b5));*/
	/*c3 = ADD(c3,MUL(ai,b1));		c8 = ADD(c8,MUL(ai,b6));*/
	c4 = ADD(c4,MUL(ai,b2));		c9 = ADD(c9,MUL(ai,b7));

	ai = c[6];
	aj = ADD(ai,ai);
	ak = ADD(aj,aj);
	d7 = MUL(b7,times19);
	/*c0 = ADD(c0,MUL(aj,d7));		c5 = ADD(c5,MUL(ai,b2));*/
	c1 = ADD(c1,MUL(ai,d8));		c6 = ADD(c6,MUL(aj,b3));
	c2 = ADD(c2,MUL(ak,d9));		c7 = ADD(c7,MUL(aj,b4));
	/*c3 = ADD(c3,MUL(ai,b0));		c8 = ADD(c8,MUL(aj,b5));*/
	/*c4 = ADD(c4,MUL(aj,b1));		c9 = ADD(c9,MUL(ai,b6));*/

	ai = c[8];
	aj = ADD(ai,ai);
	ak = ADD(aj,aj);
	d6 = MUL(b6,times19);
	/*c0 = ADD(c0,MUL(ai,d6));		c5 = ADD(c5,MUL(ai,b1));*/
	/*c1 = ADD(c1,MUL(ai,d7));		c6 = ADD(c6,MUL(ai,b2));*/
	/*c2 = ADD(c2,MUL(ai,d8));		c7 = ADD(c7,MUL(ai,b3));*/
	c3 = ADD(c3,MUL(ai,d9));		c8 = ADD(c8,MUL(ai,b4));
	/*c4 = ADD(c4,MUL(ai,b0));		c9 = ADD(c9,MUL(ai,b5));*/

	ai = c[1];
	aj = ADD(ai,ai);
	ak = ADD(aj,aj);
	d5 = MUL(b5,times19);
	c0 = ADD(c0,MUL(aj,d5));		c5 = ADD(c5,MUL(ai,b0));
	c1 = ADD(c1,MUL(aj,d6));		c6 = ADD(c6,MUL(ak,b1));
	c2 = ADD(c2,MUL(ak,d7));		c7 = ADD(c7,MUL(aj,b2));
	c3 = ADD(c3,MUL(aj,d8));		c8 = ADD(c8,MUL(ak,b3));
	c4 = ADD(c4,MUL(ak,d9));		c9 = ADD(c9,MUL(aj,b4));

	ai = c[3];
	aj = ADD(ai,ai);
	ak = ADD(aj,aj);
	d4 = MUL(b4,times19);
	c0 = ADD(c0,MUL(aj,d4));		c5 = ADD(c5,MUL(aj,d9));
	/*c1 = ADD(c1,MUL(ai,d5));		c6 = ADD(c6,MUL(ai,b0));*/
	c2 = ADD(c2,MUL(ai,d6));		c7 = ADD(c7,MUL(ai,b1));
	c3 = ADD(c3,MUL(aj,d7));		c8 = ADD(c8,MUL(aj,b2));
	c4 = ADD(c4,MUL(aj,d8));		c9 = ADD(c9,MUL(aj,b3));

	ai = c[5];
	aj = ADD(ai,ai);
	ak = ADD(aj,aj);
	d3 = MUL(b3,times19);
	c0 = ADD(c0,MUL(ak,d3));		c5 = ADD(c5,MUL(aj,d8));
	c1 = ADD(c1,MUL(aj,d4));		c6 = ADD(c6,MUL(ak,d9));
	/*c2 = ADD(c2,MUL(aj,d5));		c7 = ADD(c7,MUL(ai,b0));*/
	/*c3 = ADD(c3,MUL(ai,d6));		c8 = ADD(c8,MUL(aj,b1));*/
	c4 = ADD(c4,MUL(aj,d7));		c9 = ADD(c9,MUL(ai,b2));

	ai = c[7];
	aj = ADD(ai,ai);
	/*ak = ADD(aj,aj);*/
	/*c0 = ADD(c0,MUL(ai,d2));		c5 = ADD(c5,MUL(ai,d7));*/
	c1 = ADD(c1,MUL(ai,d3));		c6 = ADD(c6,MUL(ai,d8));
	c2 = ADD(c2,MUL(aj,d4));		c7 = ADD(c7,MUL(aj,d9));
	/*c3 = ADD(c3,MUL(ai,d5));		c8 = ADD(c8,MUL(ai,b0));*/
	/*c4 = ADD(c4,MUL(ai,d6));		c9 = ADD(c9,MUL(ai,b1));*/

	ai = c[9];
	aj = ADD(ai,ai);
	/*ak = ADD(aj,aj);*/
	/*c0 = ADD(c0,MUL(aj,d1));		c5 = ADD(c5,MUL(ai, d6));*/
	/*c1 = ADD(c1,MUL(ai,d2));		c6 = ADD(c6,MUL(aj, d7));*/
	/*c2 = ADD(c2,MUL(aj,d3));		c7 = ADD(c7,MUL(ai, d8));*/
	c3 = ADD(c3,MUL(ai,d4));		c8 = ADD(c8,MUL(aj, d9));
	/*c4 = ADD(c4,MUL(aj,d5));		c9 = ADD(c9,MUL(ai, b0));*/

	c[0] = c0;	c[1] = c5;
	c[2] = c1;	c[3] = c6;
	c[4] = c2;	c[5] = c7;
	c[6] = c3;	c[7] = c8;
	c[8] = c4;	c[9] = c9;
}

inline void compress(argElement_4w C) {
	const uint64_t ones26 = ((uint64_t) 1 << BASE0_FP25519) - 1;
	const uint64_t ones25 = ((uint64_t) 1 << BASE1_FP25519) - 1;
	const __m256i mask26  = SET164(ones26);
	const __m256i mask25  = SET164(ones25);

	__m256i c0 = LOAD(C+0);__m256i c5 = LOAD(C+1);
	__m256i c1 = LOAD(C+2);__m256i c6 = LOAD(C+3);
	__m256i c2 = LOAD(C+4);__m256i c7 = LOAD(C+5);
	__m256i c3 = LOAD(C+6);__m256i c8 = LOAD(C+7);
	__m256i c4 = LOAD(C+8);__m256i c9 = LOAD(C+9);

	__m256i h0,h1,h2,h3,h4,h5,h6,h7,h8,h9;

	h0 = SHR(c0, BASE0_FP25519);
	c0 = AND(c0, mask26);
	c1 = ADD(c1, h0);

	h1 = SHR(c1, BASE1_FP25519);
	c1 = AND(c1, mask25);
	c2 = ADD(c2, h1);

	h2 = SHR(c2, BASE0_FP25519);
	c2 = AND(c2, mask26);
	c3 = ADD(c3, h2);

	h3 = SHR(c3, BASE1_FP25519);
	c3 = AND(c3, mask25);
	c4 = ADD(c4, h3);

	h4 = SHR(c4, BASE0_FP25519);
	c4 = AND(c4, mask26);
	c5 = ADD(c5, h4);

	h5 = SHR(c5, BASE1_FP25519);
	c5 = AND(c5, mask25);
	c6 = ADD(c6, h5);

	h6 = SHR(c6, BASE0_FP25519);
	c6 = AND(c6, mask26);
	c7 = ADD(c7, h6);

	h7 = SHR(c7, BASE1_FP25519);
	c7 = AND(c7, mask25);
	c8 = ADD(c8, h7);

	h8 = SHR(c8, BASE0_FP25519);
	c8 = AND(c8, mask26);
	c9 = ADD(c9, h8);

	h9 = SHR(c9, BASE1_FP25519);
	c9 = AND(c9, mask25);
	c0 = ADD(c0, MUL19(h9));

	h0 = SHR(c0, BASE0_FP25519);
	c0 = AND(c0, mask26);
	c1 = ADD(c1, h0);

	STORE(C+0,c0);	STORE(C+1,c5);
	STORE(C+2,c1);	STORE(C+3,c6);
	STORE(C+4,c2);	STORE(C+5,c7);
	STORE(C+6,c3);	STORE(C+7,c8);
	STORE(C+8,c4);	STORE(C+9,c9);
}
