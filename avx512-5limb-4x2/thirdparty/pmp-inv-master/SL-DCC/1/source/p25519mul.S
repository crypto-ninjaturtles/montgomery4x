/*
+-----------------------------------------------------------------------------+
| This code corresponds to the Galois field F(2^255-19) from the paper        |
| "Efficient Inversion In (Pseudo-)Mersenne Prime Order Fields" authored by   |
| Kaushik Nath,  Indian Statistical Institute, Kolkata, India, and            |
| Palash Sarkar, Indian Statistical Institute, Kolkata, India.	              |
+-----------------------------------------------------------------------------+
| Copyright (c) 2018, Kaushik Nath and Palash Sarkar.                         |
|                                                                             |
| Permission to use this code is granted.                          	      |
|                                                                             |
| Redistribution and use in source and binary forms, with or without          |
| modification, are permitted provided that the following conditions are      |
| met:                                                                        |
|                                                                             |
| * Redistributions of source code must retain the above copyright notice,    |
|   this list of conditions and the following disclaimer.                     |
|                                                                             |
| * Redistributions in binary form must reproduce the above copyright         |
|   notice, this list of conditions and the following disclaimer in the       |
|   documentation and/or other materials provided with the distribution.      |
|                                                                             |
| * The names of the contributors may not be used to endorse or promote       |
|   products derived from this software without specific prior written        |
|   permission.                                                               |
+-----------------------------------------------------------------------------+
| THIS SOFTWARE IS PROVIDED BY THE AUTHORS ""AS IS"" AND ANY EXPRESS OR       |
| IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES   |
| OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.     |
| IN NO EVENT SHALL THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,      |
| INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT    |
| NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,   |
| DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY       |
| THEORY LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING |
| NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,| 
| EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                          |
+-----------------------------------------------------------------------------+
*/

.p2align 5
.globl gfp25519mul
gfp25519mul:

    push    %rbp
    push    %rbx
    push    %r12
    push    %r13
    push    %r14
    push    %r15
    
    movq    %rdx, %rbx
    
    xorq    %r13, %r13
    movq    0(%rbx), %rdx    

    mulx    0(%rsi), %r8, %r9
    mulx    8(%rsi), %rcx, %r10
    adcx    %rcx, %r9     

    mulx    16(%rsi), %rcx, %r11
    adcx    %rcx, %r10    

    mulx    24(%rsi), %rcx, %r12
    adcx    %rcx, %r11
    adcx    %r13, %r12

    xorq    %r14, %r14
    movq    8(%rbx), %rdx
   
    mulx    0(%rsi), %rcx, %rbp
    adcx    %rcx, %r9
    adox    %rbp, %r10
    
    mulx    8(%rsi), %rcx, %rbp
    adcx    %rcx, %r10
    adox    %rbp, %r11
    
    mulx    16(%rsi), %rcx, %rbp
    adcx    %rcx, %r11
    adox    %rbp, %r12
    
    mulx    24(%rsi), %rcx, %rbp
    adcx    %rcx, %r12
    adox    %rbp, %r13
    adcx    %r14, %r13

    xorq    %r15, %r15

    movq    16(%rbx), %rdx
    
    mulx    0(%rsi), %rcx, %rbp
    adcx    %rcx, %r10
    adox    %rbp, %r11
    
    mulx    8(%rsi), %rcx, %rbp
    adcx    %rcx, %r11
    adox    %rbp, %r12
    
    mulx    16(%rsi), %rcx, %rbp
    adcx    %rcx, %r12
    adox    %rbp, %r13
    
    mulx    24(%rsi), %rcx, %rbp
    adcx    %rcx, %r13
    adox    %rbp, %r14
    adcx    %r15, %r14

    xorq    %rax, %rax

    movq    24(%rbx), %rdx
    
    mulx    0(%rsi), %rcx, %rbp
    adcx    %rcx, %r11
    adox    %rbp, %r12
    
    mulx    8(%rsi), %rcx, %rbp
    adcx    %rcx, %r12
    adox    %rbp, %r13
    
    mulx    16(%rsi), %rcx, %rbp
    adcx    %rcx, %r13
    adox    %rbp, %r14
    
    mulx    24(%rsi), %rcx, %rbp
    adcx    %rcx, %r14
    adox    %rbp, %r15			
    adcx    %rax, %r15

    xorq    %rbp, %rbp
    movq    $38, %rdx

    mulx    %r12, %rax, %r12 
    adcx    %rax, %r8
    adox    %r12, %r9

    mulx    %r13, %rcx, %r13
    adcx    %rcx, %r9
    adox    %r13, %r10

    mulx    %r14, %rcx, %r14
    adcx    %rcx, %r10
    adox    %r14, %r11

    mulx    %r15, %rcx, %r15
    adcx    %rcx, %r11
    adox    zero, %r15
    adcx    zero, %r15

    shld    $1, %r11, %r15
    andq    mask63, %r11

    imul    $19, %r15, %r15
    addq    %r15, %r8
    adcq    $0, %r9
    adcq    $0, %r10
    adcq    $0, %r11

    movq    %r8,   0(%rdi)
    movq    %r9,   8(%rdi)
    movq    %r10, 16(%rdi)
    movq    %r11, 24(%rdi)

    pop     %r15
    pop     %r14
    pop     %r13
    pop     %r12
    pop     %rbx
    pop     %rbp

    ret


.p2align 5
.globl gfp25519reduce
gfp25519reduce:

    movq    0(%rdi),   %r8
    movq    8(%rdi),   %r9
    movq    24(%rdi), %r10

    movq    %r10, %r11
    shrq    $63, %r11
    andq    mask63, %r10

    imul    $19, %r11, %r11
    addq    %r11, %r8
    adcq    $0, %r9

    movq    %r8,   0(%rdi)
    movq    %r9,   8(%rdi)
    movq    %r10, 24(%rdi)

    ret
