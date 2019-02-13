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
.global gfp25519nsqr
gfp25519nsqr:

    push    %rbp
    push    %rbx
    push    %r12
    push    %r13
    push    %r14
    push    %r15
    push    %rdi

    movq    0(%rsi),  %rbx  
    movq    8(%rsi),  %rbp  
    movq    16(%rsi), %rax
    movq    24(%rsi), %rsi

    movq    %rdx, %rdi

    .START:

    decq    %rdi

    xorq    %r13, %r13
    movq    %rbx, %rdx
    
    mulx    %rbp, %r9, %r10

    mulx    %rax, %rcx, %r11
    adcx    %rcx, %r10
    
    mulx    %rsi, %rcx, %r12
    adcx    %rcx, %r11
    adcx    %r13, %r12

    movq    %rbp, %rdx
    xorq    %r14, %r14
    
    mulx    %rax, %rcx, %rdx
    adcx    %rcx, %r11
    adox    %rdx, %r12
    
    movq    %rbp, %rdx
    mulx    %rsi, %rcx, %rdx
    adcx    %rcx, %r12
    adox    %rdx, %r13
    adcx    %r14, %r13

    xorq    %r15, %r15
    movq    %rax, %rdx
    
    mulx    %rsi, %rcx, %r14
    adcx    %rcx, %r13
    adcx    %r15, %r14

    shld    $1, %r14, %r15
    shld    $1, %r13, %r14
    shld    $1, %r12, %r13
    shld    $1, %r11, %r12
    shld    $1, %r10, %r11
    shld    $1, %r9, %r10
    shlq    $1, %r9
     
    xorq    %rdx, %rdx
    movq    %rbx, %rdx
    mulx    %rdx, %r8, %rdx
    adcx    %rdx, %r9

    movq    %rbp, %rdx
    mulx    %rdx, %rcx, %rdx
    adcx    %rcx, %r10
    adcx    %rdx, %r11

    movq    %rax, %rdx
    mulx    %rdx, %rcx, %rdx
    adcx    %rcx, %r12
    adcx    %rdx, %r13

    movq    %rsi, %rdx
    mulx    %rdx, %rcx, %rdx
    adcx    %rcx, %r14
    adcx    %rdx, %r15

    xorq    %rbp, %rbp
    movq    $38, %rdx    		

    mulx    %r12, %rbx, %rbp
    adcx    %r8, %rbx
    adox    %r9, %rbp

    mulx    %r13, %rcx, %rax
    adcx    %rcx, %rbp
    adox    %r10, %rax

    mulx    %r14, %rcx, %rsi
    adcx    %rcx, %rax
    adox    %r11, %rsi

    mulx    %r15, %rcx, %r15
    adcx    %rcx, %rsi
    adox    zero, %r15
    adcx    zero, %r15

    shld    $1, %rsi, %r15
    andq    mask63, %rsi

    imul   $19, %r15, %r15
    addq   %r15, %rbx
    adcq   $0, %rbp
    adcq   $0, %rax
    adcq   $0, %rsi

    cmpq   $0, %rdi

    jne    .START

    pop     %rdi

    movq    %rbx,  0(%rdi)
    movq    %rbp,  8(%rdi)
    movq    %rax, 16(%rdi)
    movq    %rsi, 24(%rdi)

    pop     %r15
    pop     %r14
    pop     %r13
    pop     %r12
    pop     %rbx
    pop     %rbp

    ret
