  .section	.rdata,"dr"
	.p2align	4               # @_ZZL7rsqrt64yiE11exp_adj_tab
_rsqrt64_tables:
  # rsqrt_exp_adj_tab
	.long	-1
	.long	3037000500
  # rsqrt_tab
  .byte 252, 245,  238,  232,  226,  221,  216,  211
  .byte 207, 203,  199,  195,  192,  189,  185,  182
	.addrsig

	.text
	.def	 @feat.00;
	.scl	3;
	.type	0;
	.endef
	.globl	@feat.00
.set @feat.00, 0
	.file	"uu.cpp"
	.def	 asm_sqrt448;
	.scl	2;
	.type	32;
	.endef
	.globl	asm_sqrt448             # -- Begin function asm_sqrt448

	.p2align	4, 0x90
.smul2:
# signed multiply 128-bit number a, by 64-bit number b, return upper 128 bits of the product
# inputs:
# rdx:rax - a
# rbx     - b
# outputs:
# rdx:rax - result
# rbx     - b (preserved)
# other registers affected:
# rdi, rbp
  mov   %rdx, %rdi
  sar   $63,  %rdx                # rdx = sign(err)
  xor   %rdx, %rax
  xor   %rdx, %rdi                # rdi:rax = abs(err)
  mov   %rdx, %rbp                # rbp = sign(err)
  # Sqrt_adj = abs(err)*invS/2, scaled by 2**(127-exp1)
  mulq  %rbx                      # rdx:rax = invS[3]*abs_err[0]
  mov   %rdi, %rax                # rax = abs_err[1]
  mov   %rdx, %rdi                # rdi = mulh(invS[3],abs_err[0])
  mulq  %rbx                      # rdx:rax = invS[3]*abs_err[1]
  add   %rdi, %rax
  adc   $0,   %rdx                # rdx:rax = invS[3]*abs_err[1] + mulh(invS[3],abs_err[0])
  #  = abs(err)*invS/2, scaled by 2**(128-exp1)
  not   %rbp
  xor   %rbp, %rax
  xor   %rbp, %rdx                # rdx:rax = -err*invS/2, scaled by 2**(128-exp1)
  ret

.smul4:
# signed multiply 192-bit number a, by 128-bit number b, return upper 192 bits of the product
# inputs:
# rdx:rax:r8 - a
# rbx:r14    - b
# outputs:
# rdx:rax:r8 - result
# rbx:r14    - b (preserved)
# rbp        - ~sign(a)
# other registers affected:
# rdi, r10, r11
  mov   %rdx, %rdi
  sar   $63,  %rdx                # rdx = sign(a)
  xor   %rdx, %rax
  xor   %rdx, %rdi
  xor   %rdx, %r8
  # rdi:rax:r8 = abs(a)
  mov   %rdx, %rbp                # rbp = sign(a)
  mov   %rax, %r10                # r10 = abs_a[1]
  mulq  %r14                      # rdx:rax = abs_a[1]*b[0]
  mov   %r8,  %rax
  mov   %rdx, %r8                 # r8 = prod[0] = (abs_a[1]*b[0]).HI
  mulq  %rbx                      # rdx:rax = abs_a[0]*b[1]
  xor   %r11d, %r11d
  add   %rdx, %r8                 # r8  = prod[0] += (abs_a[1]*b[0]).HI
  adc   %r11, %r11                # r11 = prod[1] = carry
  mov   %r10, %rax                # rax = abs_a[1]
  mulq  %rbx                      # rdx:rax = abs_a[1]*b[1]
  xor   %r10d, %r10d
  add   %rax, %r8                 # r8  = prod[0] += (abs_a[1]*b[1]).LO
  adc   %rdx, %r11                # r11 = prod[1] += (abs_a[1]*b[1]).LO + carry
  adc   %r10, %r10                # r10 = prod[2] = carry
  mov   %rdi, %rax                # rax = abs_a[2]
  mulq  %r14                      # rdx:rax = abs_a[2]*b[0]
  add   %rax, %r8                 # r8  = prod[0] += (abs_a[2]*b[0]).LO
  adc   %rdx, %r11                # r11 = prod[1] += (abs_a[2]*b[0]).LO + carry
  adc   $0,   %r10                # r10 = prod[2] += carry
  mov   %rdi, %rax                # rax = abs_a[2]
  mulq  %rbx                      # rdx:rax = abs_a[2]*b[1]
  add   %r11, %rax                # rax = prod[1] += (abs_a[2]*b[1]).LO
  adc   %r10, %rdx                # rdx = prod[2] += (abs_a[2]*b[0]).HI + carry
  # rdx:rax:r8 = abs(a)*b / 2**128
  not   %rbp
  xor   %rbp, %r8
  xor   %rbp, %rax
  xor   %rbp, %rdx
  # rdx:rax:r8 = -a*b / 2**128
  ret

	.p2align	4, 0x90
asm_sqrt448:                            # @asm_sqrt448
.seh_proc asm_sqrt448
# %bb.0:
	pushq	%r15
	.seh_pushreg 15
	pushq	%r14
	.seh_pushreg 14
	pushq	%r13
	.seh_pushreg 13
	pushq	%r12
	.seh_pushreg 12
	pushq	%rsi
	.seh_pushreg 6
	pushq	%rdi
	.seh_pushreg 7
	pushq	%rbp
	.seh_pushreg 5
	pushq	%rbx
	.seh_pushreg 3
	subq	$96, %rsp
	.seh_stackalloc 96
	.seh_endprologue

  movq  %rcx, 56(%rsp)      # save dst pointer
	movq	48(%rdx), %rax      # r10 = x = src[6]
  mov   %rdx, %rsi          # rsi = src
	movq	%rax, %rbx          # rbx = x
	movq	%rax, %r10          # r10 = x
	leaq	_rsqrt64_tables(%rip), %rdi
	shr   $59, %rax
	and   $15, %eax             # rax = idx = (x >> 59) & 15;
	movzbl 8(%rax,%rdi), %eax   # rax = y = rsqrt_tab[idx]; y scaled by 2^8

  # 1st NR step - y = y*(3+eps1 - y*y*x)/2
  mov  %eax, %edx             # rdx = y
  imul %eax, %eax             # rax = yy = y*y;    scaled by 2^16
  shr  $49,  %rbx             # rbx = x >> (64-15)
  imul %ebx, %eax             # rax = yyx = yy*(x >> (64-15)) = y*y*x scaled by 2^(63+16+15-64)=2^30
  NR1_3 = ((3)<<30)+((3)<<17);
  mov  $(NR1_3), %ebx
  sub  %eax, %ebx
  imul %rbx, %rdx
  shr  $7,   %rdx             # rdx = y = (y * ((3u << 30)+(3u<<17)-yyx))>>(6+1); scaled by 2^(8+30-6)=2^32

  # 2nd NR step - y = y*(3+eps2 - y*y*x)/2
  mov  %rdx, %rax             # rax = y
  imul %rdx, %rdx             # rdx = yy=y*y       scaled by 2^64
  mov  %r10, %rbx
  shr  $32,  %rdx             # rdx = yy=(y*y)>>32 scaled by 2^32
  shr  $32,  %rbx             # rbx = x >> 32;     scaled by 2^31
  imul %rbx, %rdx
  shr  $33,  %rdx             # rdx = yyx = (yy*(x >> 32)) >> 33 = y*y*x scaled by 2^(32+63-32-33)=2^30
  NR2_3 = ((3)<<30)+((3)<<4);
  mov  $(NR2_3), %ebx
  sub  %edx, %ebx             # rbx = (3u << 30)+(3u<<4)-yyx
  imul %rbx, %rax             # rax = y * ((3u << 30)+(3u<<4)-yyx)
  mov  %r8d, %ecx             # rcx = exp1
  shr  $31,  %rax             # rax = y = (y * ((3u << 30)+(3u<<4)-yyx))>>(30+1); scaled by 2^(32+30-30)=2^32

  # handle exponent
  mull (%rdi,%rcx,4)          # rdx = y = (y * exp_adj_tabb[exp1]) >> 32

  # 3rd NR step
  # Use 2nd order polynomial: y =  y - y*(err/2 - 3/8*err**2) = y - y/2*(err - 3/4*err**2)
  mov  %rdx, %rbp             # rbp = y
  imul %rdx, %rdx             # rdx = y*y
  lea  (%r10,%r10), %rax      # rax = x1 = x*2 = x - 1           scaled by 2**64
  shl  %cl,  %rdx             # rdx = yy = (y * y) << exp1; scaled by 2**64
  mov  %rdx, %rbx             # rbx = yy
  mul  %rdx                   # rdx:rax = x1*yy
  add  %rbx, %rdx             # rdx = err = int64_t(mulh(yy, x1)+yy) = y*y*x-1 scaled by 2**64
  mov  %rdx, %rax             # rax = err
  mov  %rdx, %rbx             # rbx = err
  imul %rax                   # rdx = err2 = imulh(err, err) = err*err scaled by 2**64
  lea  (%rdx,%rdx,2),%rax     # rax = err2*3
  shl  $2,   %rbx             # rbx = err*4
  sub  %rbx, %rax             # rax = m2 = err2*3 - err*4 = 3*err**2 - err*4 scaled by 2**64
  mov  %rbp, %rbx             # rbx = y
  shl  $29,  %rbp             # rbp = y << 29
  imul %rbp                   # rdx = adj = imulh(m2, int64_t(y<<29)) = m2*(y/8)
  shl  $32,  %rbx             # rbx = y << 32
  add  %rdx, %rbx             # rbx = invSx1 = (y<<32) + adj;

  # Calculate sqrt() as src*invSqrt
  movq  %rbx, %rax                # rax = invSx1
  mulq  %r10                      # rdx = mulh(invSx1, src[6])
  shl   %cl,  %rdx                # rdx = mulh(invSx1, src[6]) << expl
  add   %rcx, %rdx                # rdx = (mulh(invSx1, src[6]) << expl) + expl
  xor   $1,   %ecx                # ecx = exp1^1
  movq  %rdx, %r12                # r12 = Sqrt[7]

  # improve precision to 127-eps bits
  mov   %rdx, %rax
  mulq  %rdx                      # rdx:rax = sqr(sqrt(src)) scaled by 2**126
	shldq	%cl, %rax, %rdx
  shl   %cl, %rax                 # rdx:rax = rdx:rax << (exp1^1) = sqr(sqrt(src)) scaled by 2**(127-exp1)
  # Sqrt_sqr = err = Sqrt**2 - src
  sub   40(%rsi), %rax
  sbb   %r10,     %rdx
  # rdx:rax = err = Sqrt**2 - src
  call  .smul2
  # rdx:rax = -err*invS/2, scaled by 2**(128-exp1)
	shrdq	%cl, %rdx, %rax
  sar   %cl, %rdx                 # rdx:rax = abs(err)*invS/2, scaled by 2**127
  # Sqrt -= err*invS/2
  add   %rdx, %r12                # r12 = Sqrt[7]
  mov   %rax, %r13                # r13 = Sqrt[6]

  # adjust invS
  # Rsqrt_err = sqrt*rsqrt, scaled by 2**127
  mulq  %rbx                      # rdx:rax = invS[3]*Sqrt[6]
  mov   %r12, %rax                # rax = Sqrt[7]
  mov   %rdx, %rdi                # rdi = mulh(invS[3],Sqrt[6])
  mulq  %rbx                      # rdx:rax = invS[3]*Sqrt[7]
  add   %rdi, %rax
  adc   $0,   %rdx                # rdx:rax = invS[3]*Sqrt[7] + mulh(invS[3],Sqrt[6])
  # rdx:rax = Rsqrt_err = sqrt*rsqrt, scaled by 2**127
	shldq	$1, %rax, %rdx
  shl   $1, %rax                  # rdx:rax = rdx:rax << 1 = Rsqrt_err = sqrt*rsqrt, scaled by 2**128
  call  .smul2
  # rdx:rax = -Rsqrt_err*invS, scaled by 2**127
  # invS -= Rsqrt_err*invS
  add   %rdx, %rbx
  jnz .no_ovf2
    movq $-1, %rbx
    movq %rbx,%rax
  .no_ovf2:
  mov   %rax, %r14                # r14 = invS[2]

  # Registers usage
  # rax - free (copy of invS[2])
  # rdx - free
  # rcx - (exp1^1)
  # rbx - invS[3]
  # rsi - src
  # rdi - free
  # rbp - free
  # r8  - free
  # r9  - free
  # r10 - free (copy of src[6])
  # r11 - free
  # r12 - Sqrt[7]
  # r13 - Sqrt[6]
  # r14 - invS[2]
  # r15 - free

  # improve precision to 255-eps bits
  # Sqrt_sqr = sqr(sqrt(src)) scaled by 2**254
  mov  %r13, %rax                  # rax = Sqrt[6]
  mulq %rax                        # rdx:rax = Sqrt[6]*Sqrt[6]
  mov  %rax, %r8                   # r8  = Sqrt_sqr[0]
  mov  %rdx, %r10                  # r10 = Sqrt_sqr[1]
  mov  %r13, %rax                  # rax = Sqrt[6]
  mulq %r12                        # rdx:rax = Sqrt[6]*Sqrt[7]
  add  %rax, %rax
  adc  %rdx, %rdx                  # rdx:rax = Sqrt[6]*Sqrt[7]*2
  mov  %r12, %r11
  imul %r11, %r11                  # r11 = Sqrt[7]*Sqrt[7]
  add  %r10, %rax                  # rax = Sqrt_sqr[1]
  adc  %r11, %rdx                  # rdx = Sqrt_sqr[2]
  # rdx:rax:r8 = sqr(sqrt(src))[2:0] scaled by 2**254
	shldq	%cl, %rax, %rdx
	shldq	%cl, %r8 , %rax
  shl   %cl, %r8
  # rdx:rax:r8 = sqr(sqrt(src))[2:0] scaled by 2**(255-exp1)
  sub   24(%rsi),%r8               # Sqrt_sqr[0] -= src[3]
  sbb   32(%rsi),%rax              # Sqrt_sqr[1] -= src[4]
  sbb   40(%rsi),%rdx              # Sqrt_sqr[2] -= src[5]
  call  .smul4
  # rdx:rax:r8 = Sqrt_adj = -err*invS/2, scaled by 2**(255-exp1)
	shrdq	%cl, %rax, %r8
	shrdq	%cl, %rdx, %rax
  sar   %cl, %rdx
  # rdx:rax:r8 = Sqrt_adj = -err*invS/2, scaled by 2**255
  # Sqrt += Sqrt_adj
  add  %rdx, %r13
  adc  %rbp, %r12
  movq %r8,  %r9                   # r9 = Sqrt[4]

  # adjust invS

  # Rsqrt_err = sqrt*rsqrt, scaled by 2**255
  # Multiply a=Sqrt[7:4] (4 qwords) by b=invS[3:2] (2 qwords),
  # Rsqrt_err[2:0] = words[4:2] of the product
  mov  %rax, %r15                  # r15 = Sqrt[5]
  mulq %r14                        # rdx:rax = a[1]*b[0]
  mov  %r8,  %rax                  # rax = a[0]==Sqrt[4]
  mov  %rdx, %r8                   # r8  = prod2 = (a[1]*b[0]).HI
  mulq %rbx                        # rdx:rax = a[0]*b[1]
  xor  %edi, %edi                  # rdi = prod3 = 0
  add  %rdx, %r8                   # r8  = prod2+= (a[1]*b[0]).HI
  adc  %rdi, %rdi                  # rdi = prod3+= carry
  mov  %r15, %rax                  # rax = a[1]==Sqrt[5]
  mulq %rbx                        # rdx:rax = a[1]*b[1]
  add  %rax, %r8                   # r8  = prod2+= (a[1]*b[1]).LO
  adc  %rdx, %rdi                  # rdi = prod3+= (a[1]*b[1]).HI + carry. Here carry out can't happen
  mov  %r13, %rax                  # rax = a[2]==Sqrt[6]
  mulq %r14                        # rdx:rax = a[2]*b[0]
  xor  %ebp, %ebp                  # rbp = prod4 = 0
  add  %rax, %r8                   # r8  = prod2+= (a[1]*b[1]).LO
  adc  %rdx, %rdi                  # rdi = prod3+= (a[1]*b[1]).HI + carry.
  adc  %rbp, %rbp                  # rbp = prod4+= carry
  mov  %r13, %rax                  # rax = a[2]==Sqrt[6]
  mulq %rbx                        # rdx:rax = a[2]*b[1]
  add  %rax, %rdi                  # rdi = prod3+= (a[2]*b[1]).LO
  adc  %rdx, %rbp                  # rbp = prod4+= (a[2]*b[1]).HI + carry.
  mov  %r12, %rax                  # rax = a[3]==Sqrt[7]
  mulq %r14                        # rdx:rax = a[3]*b[0]
  add  %rdi, %rax                  # rax = prod3+= (a[3]*b[0]).LO
  adc  %rdx, %rbp                  # rbp = prod4+= (a[3]*b[0]).HI + carry.
  mov  %r12, %rdx                  # rdx = a[3]==Sqrt[7]
  imul %rbx, %rdx                  # rdx = (a[3]*b[1]).LO
  add  %rbp, %rdx                  # rdx = prod4+= (a[3]*b[1]).LO
  # rdx:rax:r8 = Rsqrt_err = Rsqrt_err << 1
	shldq	$1, %rax, %rdx
	shldq	$1, %r8,  %rax
  shl   $1, %r8
  # Rsqrt_err = sqrt*rsqrt-1, scaled by 2**256
  call  .smul4
  # rdx:rax:r8 = -Rsqrt_err*invS, scaled by 2**255
  # invS -= Rsqrt_err*invS
  add  %rdx, %r14
  adc  %rbp, %rbx
  jnz .no_ovf4
    # overflow in temporary result - set invS to maximum
    movq $-1, %rbx
    movq %rbx,%r14
    movq %rbx,%r8
    movq %rbx,%rax
  .no_ovf4:
  movq %r8,    (%rsp)              # wrkbuf[0]  = invS[0] = r8
  movq %rax,  8(%rsp)              # wrkbuf[1]  = invS[1] = rax
  movq %r14, 16(%rsp)              # wrkbuf[2]  = invS[2] = r14
  movq %rbx, 24(%rsp)              # wrkbuf[3]  = invS[3] = rbx

  # Registers usage
  # rax - invS[1]
  # rdx - free
  # rcx - (exp1^1)
  # rbx - invS[3]
  # rsi - src
  # rdi - free
  # rbp - free
  # r8  - invS[0]
  # r9  - Sqrt[4]
  # r10 - free
  # r11 - free
  # r12 - Sqrt[7]
  # r13 - Sqrt[6]
  # r14 - invS[2]
  # r15 - Sqrt[5]

  # improve precision of Sqrt to 511-eps bits
  # Sqrt_sqr = sqr(Sqrt) scaled by 2**510
  # Calculate 5 LS qwords of the product
  mov %r9,  %rax                   # rax = a[0] = Sqrt[4]
  mul %rax                         # rdx:rax = a[0]*a[0]
  mov %rax, %r10                   # r10 = Sqrt_sqr[0] = (a[0]*a[0]).LO
  mov %rdx, %rdi                   # rdi = Sqrt_sqr[1] = (a[0]*a[0]).HI

  mov %r15, %rax                   # rax = a[1] = Sqrt[5]
  mul %r9                          # rdx:rax = a[0]*a[1]
  mov %rdx, %rbx                   # rbx = Sqrt_sqr[2] = (a[0]*a[1]).HI
  add %rax, %rdi                   # rdi = Sqrt_sqr[1]+= (a[0]*a[1]).LO
  adc $0,   %rbx                   # rbx = Sqrt_sqr[2]+= carry
  xor %r11d, %r11d                 # r11 = Sqrt_sqr[3] = 0
  add %rax, %rdi                   # rdi = Sqrt_sqr[1]+= (a[0]*a[1]).LO
  adc %rdx, %rbx                   # rbx = Sqrt_sqr[2]+= (a[0]*a[1]).HI + carry
  adc %r11, %r11                   # r11 = Sqrt_sqr[3]+= carry

  mov %r13, %rax                   # rax = a[2] = Sqrt[6]
  mul %r9                          # rdx:rax = a[0]*a[2]
  xor %r8d, %r8d                   # r8  = Sqrt_sqr[4] = 0
  add %rax, %rax
  adc %rdx, %rdx
  adc %r8,  %r8                    # r8:rdx:rax = a[0]*a[2]*2
  add %rax, %rbx                   # rbx = Sqrt_sqr[2]+= (a[0]*a[2]*2).LO
  adc %rdx, %r11                   # r11 = Sqrt_sqr[3]+= (a[0]*a[2]*2).HI + carry
  adc $0,   %r8                    # r8  = Sqrt_sqr[4]+= carry
  mov %r15, %rax                   # rax = a[1] = Sqrt[5]
  mul %rax                         # rdx:rax = a[1]*a[1]
  add %rax, %rbx                   # rbx = Sqrt_sqr[2]+= (a[1]*a[1]).LO
  adc %rdx, %r11                   # r11 = Sqrt_sqr[3]+= (a[1]*a[1]).HI + carry
  adc $0,   %r8                    # r8  = Sqrt_sqr[4]+= carry

  mov %r12, %rax                   # rax = a[3] = Sqrt[7]
  mul %r9                          # rdx:rax = a[0]*a[3]
  mov %rax, %rbp                   # rbp = (a[0]*a[3]).LO
  mov %rdx, %r14                   # r14 = (a[0]*a[3]).HI
  mov %r13, %rax                   # rax = a[2] = Sqrt[6]
  mul %r15                         # rdx:rax = a[1]*a[2]
  add %rbp, %rax                   #
  adc %r14, %rdx                   # rdx:rax = a[0]*a[3]+a[1]*a[2]
  add %rax, %rax
  adc %rdx, %rdx                   # rdx:rax = (a[0]*a[3]+a[1]*a[2])*2
  add %rax, %r11                   # r11 = Sqrt_sqr[3]+= ((a[0]*a[3]+a[1]*a[2])*2).LO + carry
  adc %rdx, %r8                    # r8  = Sqrt_sqr[4]+= ((a[0]*a[3]+a[1]*a[2])*2).HI + carry

  mov %r12, %rax                   # rax = a[3] = Sqrt[7]
  imul %r15,%rax                   # rax = (a[1]*a[3]).LO
  add %rax, %rax                   # rax = (a[1]*a[3]*2).LO
  add %rax, %r8                    # r8  = Sqrt_sqr[4]+= (a[1]*a[3]*2).LO
  mov %r13, %rax                   # rax = a[2] = Sqrt[6]
  imul %rax,%rax                   # rax = (a[2]*a[3])
  add %rax, %r8                    # r8  = Sqrt_sqr[4]+= (a[2]*a[2]).LO
  # r8:r11:rbx:rdi:r10 = Sqrt_sqr = sqr(Sqrt) scaled by 2**510

  # Sqrt_sqr <<= (exp1 ^ 1)
	shldq	%cl, %r11, %r8
	shldq	%cl, %rbx, %r11
	shldq	%cl, %rdi, %rbx
	shldq	%cl, %r10, %rdi
  shl   %cl, %r10
  # r8:r11:rbx:rdi:r10 = Sqrt_sqr = sqr(Sqrt) scaled by 2**(511-exp1)

  # Sqrt_sqr_err[4:1} <<= Sqrt_sqr[4:1] - src[3:0]
  sub   (%rsi), %rdi
  sbb  8(%rsi), %rbx
  sbb 16(%rsi), %r11
  sbb 24(%rsi), %r8

  # at this point src[] is used no longer, so register rsi is free

  mov %r8,  %rbp                   # rbp = Sqrt_sqr_err[4]
  sar $63,  %r8                    # r8  = sign(Sqrt_sqr_err)

  # Registers usage
  # rax - free
  # rdx - free
  # rcx - (exp1^1)
  # rbx - Sqrt_sqr_err[2]
  # rsi - free
  # rdi - Sqrt_sqr_err[1]
  # rbp - Sqrt_sqr_err[4]
  # r8  - sign(Sqrt_sqr_err)
  # r9  - Sqrt[4]
  # r10 - Sqrt_sqr_err[0]
  # r11 - Sqrt_sqr_err[3]
  # r12 - Sqrt[7]
  # r13 - Sqrt[6]
  # r14 - free
  # r15 - Sqrt[5]

  # Sqrt_sqr_err = abs(Sqrt_sqr_err)
  xor  %r8, %r10
  xor  %r8, %rdi
  xor  %r8, %rbx
  xor  %r8, %r11
  xor  %r8, %rbp

  # Sqrt_adj = abs(err)*invS/2, scaled by 2**(511-exp1)
  # multiply 5 qwords by 4 qwords, Sqrt_adj = 5 MS qwords of the product
  mov 24(%rsp),%rax                # rax         = b[3] = invS[3]
  mul %r10                         # rdx:rax     = a[0]*b[3]
  mov %rdx,    %r10                # r10 = prod4 = (a[0]*b[3]).HI
  mov 16(%rsp),%rax                # rax         = b[2] = invS[2]
  mul %rdi                         # rdx:rax     = a[1]*b[2]
  xor %esi,    %esi                # rsi = prod5 = 0
  add %rdx,    %r10                # r10 = prod4+=(a[1]*b[2]).HI
  adc %rsi,    %rsi                # rsi = prod5+= carry
  mov 8(%rsp), %rax                # rax         = b[1] = invS[1]
  mul %rbx                         # rdx:rax     = a[2]*b[1]
  add %rdx,    %r10                # r10 = prod4+=(a[2]*b[1]).HI
  adc $0,      %rsi                # rsi = prod5+= carry
  mov (%rsp),  %rax                # rax         = b[0] = invS[0]
  mul %r11                         # rdx:rax     = a[3]*b[0]
  add %rdx,    %r10                # r10 = prod4+=(a[3]*b[0]).HI
  adc $0,      %rsi                # rsi = prod5+= carry

  mov 24(%rsp),%rax                # rax         = b[3] = invS[3]
  mul %rdi                         # rdx:rax     = a[1]*b[3]
  xor %edi,    %edi                # rdi = prod6 = 0
  add %rax,    %r10                # r10 = prod4+=(a[1]*b[3]).LO
  adc %rdx,    %rsi                # rsi = prod5+=(a[1]*b[3]).HI + carry
  adc %rdi,    %rdi                # rdi = prod6+= carry
  mov 16(%rsp),%rax                # rax         = b[2] = invS[2]
  mul %rbx                         # rdx:rax     = a[2]*b[2]
  add %rax,    %r10                # r10 = prod4+=(a[2]*b[2]).LO
  adc %rdx,    %rsi                # rsi = prod5+=(a[2]*b[2]).HI + carry
  adc $0,      %rdi                # rdi = prod6+= carry
  mov 8(%rsp), %rax                # rax         = b[1] = invS[1]
  mul %r11                         # rdx:rax     = a[3]*b[1]
  add %rax,    %r10                # r10 = prod4+=(a[3]*b[1]).LO
  adc %rdx,    %rsi                # rsi = prod5+=(a[3]*b[1]).HI + carry
  adc $0,      %rdi                # rdi = prod6+= carry
  mov (%rsp),  %rax                # rax         = b[0] = invS[0]
  mul %rbp                         # rdx:rax     = a[4]*b[0]
  add %rax,    %r10                # r10 = prod4+=(a[4]*b[0]).LO
  adc %rdx,    %rsi                # rsi = prod5+=(a[4]*b[0]).HI + carry
  adc $0,      %rdi                # rdi = prod6+= carry

  mov 24(%rsp),%rax                # rax         = b[3] = invS[3]
  mul %rbx                         # rdx:rax     = a[2]*b[3]
  xor %ebx,    %ebx                # rbx = prod7 = 0
  add %rax,    %rsi                # rsi = prod5+=(a[2]*b[3]).LO
  adc %rdx,    %rdi                # rdi = prod6+=(a[2]*b[3]).HI + carry
  adc %rbx,    %rbx                # rbx = prod7+= carry
  mov 16(%rsp),%rax                # rax         = b[2] = invS[2]
  mul %r11                         # rdx:rax     = a[3]*b[2]
  add %rax,    %rsi                # rsi = prod5+=(a[3]*b[2]).LO
  adc %rdx,    %rdi                # rdi = prod6+=(a[3]*b[2]).HI + carry
  adc $0,      %rbx                # rbx = prod7+= carry
  mov 8(%rsp), %rax                # rax         = b[1] = invS[1]
  mul %rbp                         # rdx:rax     = a[4]*b[1]
  add %rax,    %rsi                # rsi = prod5+=(a[1]*b[4]).LO
  adc %rdx,    %rdi                # rdi = prod6+=(a[1]*b[4]).HI + carry
  adc $0,      %rbx                # rbx = prod7+= carry

  mov 24(%rsp),%rax                # rax         = b[3] = invS[3]
  mul %r11                         # rdx:rax     = a[3]*b[3]
  xor %r11d,   %r11d               # r11 = prod8 = 0
  add %rax,    %rdi                # rdi = prod6+=(a[3]*b[3]).LO
  adc %rdx,    %rbx                # rbx = prod7+=(a[3]*b[3]).HI + carry
  adc %r11,    %r11                # r11 = prod8+= carry
  mov 16(%rsp),%rax                # rax         = b[2] = invS[2]
  mul %rbp                         # rdx:rax     = a[4]*b[2]
  add %rax,    %rdi                # rdi = prod6+=(a[4]*b[2]).LO
  adc %rdx,    %rbx                # rbx = prod7+=(a[4]*b[2]).HI + carry
  adc $0,      %r11                # r11 = prod8+= carry

  mov 24(%rsp),%rax                # rax         = b[3] = invS[3]
  mul %rbp                         # rdx:rax     = a[4]*b[3]
  add %rax,    %rbx                # rbx = prod7+=(a[4]*b[3]).LO
  adc %rdx,    %r11                # r11 = prod8+=(a[4]*b[3]).HI + carry
  # r11:rbx:rdi:rsi:r10 = abs_sqrt_adj = abs(err)*invS/2, scaled by 2**(511-exp1)

  # change sign to the opposite of original sign of err
  not %r8
  xor %r8, %r10
  xor %r8, %rsi
  xor %r8, %rdi
  xor %r8, %rbx
  xor %r8, %r11
  # r11:rbx:rdi:rsi:r10 = Sqrt_adj = -err*invS/2, scaled by 2**(511-exp1)

  # r11:rbx:rdi:rsi:r10 = Sqrt_adj = -err*invS/2, scaled by 2**(511-exp1)
	shrdq	%cl, %rsi, %r10
	shrdq	%cl, %rdi, %rsi
	shrdq	%cl, %rbx, %rdi
	shrdq	%cl, %r11, %rbx
  sar   %cl, %r11
  # r11:rbx:rdi:rsi:r10 = Sqrt_adj = -err*invS/2, scaled by 2**511

  # Sqrt += Sqrt_adj
  add  %r11, %r9
  adc  %r8,  %r15
  adc  %r8,  %r13
  adc  %r8,  %r12

  # Registers usage
  # rax - free
  # rdx - free
  # rcx - free
  # rbx - Sqrt[3]
  # rsi - Sqrt[1]
  # rdi - Sqrt[2]
  # rbp - free
  # r8  - free
  # r9  - Sqrt[4]
  # r10 - Sqrt[0] = lsw
  # r11 - free
  # r12 - Sqrt[7]
  # r13 - Sqrt[6]
  # r14 - free
  # r15 - Sqrt[5]

  # check if lsw is close to the middle of qword range
  UINT64_MID = (1) << 63;
  MAX_ERR    = (1) << 20; # probably less, but it does not cost much to be on the safe side
	movabsq	$(UINT64_MID - MAX_ERR),%rax
  mov %r10, %rdx
  sub %rax, %rdx                   # rdx = lsw - (UINT64_MID-MAX_ERR)
  cmp $(MAX_ERR*2), %rdx
  jae .round_and_copy
    # We are very close to tip point
    # so more precise calculations are needed to decide if we are above or below
    # Calculate TP_sqr = (sqrt(7:1)*2**64 + 2**63)**2 == (sqrt(7:1)**2 + sqrt(7:1))*2**128 + 2**126
    # Only TP_sqr[8] is of interest, more specifically only bit 62 of TP_sqr[8]
    # When TP_sqr[8].b62 = 0 it means that TP_sqr > src, so the result = sqrt(7:1)
    # When TP_sqr[8].b62 = 1 it means that TP_sqr < src, so the result = sqrt(7:1)+1

    # The routine is supposed to reach this point very rarely, so calculation of TP_sqr[8] optimized for size rather than for speed

    # store sqrt(7:1) in work buffer
    mov %rsi, 0(%rsp)
    mov %rdi, 8(%rsp)
    mov %rbx, 16(%rsp)
    mov %r9,  24(%rsp)
    mov %r15, 32(%rsp)
    mov %r13, 40(%rsp)
    mov %r12, 48(%rsp)

    mov %rsi, %rcx                # rcx  = acc0 = a[0]
    xor %ebp, %ebp                # rbp  = acc1
    mov $1,   %r8d                # r8d  = K = 1
    .sqr_tip_point_k:
      xor %r11d, %r11d            # r11 = acc2
      add (%rsp,%r8,8), %rbp      # rbp = acc1 += a[K];
      adc %r11, %r11              # r11 = acc2 += carry
      xor %r10d, %r10d            # r10d = i = 0
      lea -1(%r8), %r14           # r14d = j = K-1
      .sqr_tip_point_j:
        mov  (%rsp,%r10,8),%rax   # rax = a[i]
        mulq (%rsp,%r14,8)        # rdx:rax = a[i]*a[j]
        add %rax, %rcx            # rcx = acc0 += (a[i]*a[j]).LO
        adc %rdx, %rbp            # rbp = acc1 += (a[i]*a[j]).HI + carry
        adc $0,   %r11            # r11 = acc2 += carry
        inc %r10d                 # ++i
      sub $1, %r14d               # --j
      jge .sqr_tip_point_j        # j >= 0

      mov %rbp, %rcx              # rcx = acc0 = acc1
      mov %r11, %rbp              # rbp = acc1 = acc2
    inc %r8d                      # ++K
    cmp $7, %r8d                  # while (K < 7)
    jb  .sqr_tip_point_k

    xor %r8d, %r8d                # r8d  = i = 0
    mov $6,   %r14d               # r14d = j = 6
    .sqr_tip_point_j2:
      mov  (%rsp,%r8,8),%rax      # rax = a[i]
      imul (%rsp,%r14,8),%rax     # rax = (a[i]*a[j]).LO
      add   %rax, %rcx            # rcx = acc0 += (a[i]*a[j]).LO
      inc %r8d                    # ++i
    sub $1, %r14d                 # --j
    jge .sqr_tip_point_j2         # j >= 0

    # rcx = TP_sqr[8]

    lea (%rcx,%rcx), %r10         # r10= lsw = acc0 << 1, that moves bit 62 in position 63, where it inspected by the rest of the code
  .round_and_copy:
  mov 56(%rsp), %rcx              # rcx = dst
  xor %eax, %eax                  # rax = 0
  add %r10, %r10                  # carry = bit 63 of lsw
  # add carry
  adc %rax, %rsi
  adc %rax, %rdi
  adc %rax, %rbx
  adc %rax, %r9
  adc %rax, %r15
  adc %rax, %r13
  adc %rax, %r12
  # strore result at dst
  mov %rsi,  0(%rcx)
  mov %rdi,  8(%rcx)
  mov %rbx, 16(%rcx)
  mov %r9,  24(%rcx)
  mov %r15, 32(%rcx)
  mov %r13, 40(%rcx)
  mov %r12, 48(%rcx)

	addq	$96, %rsp
	popq	%rbx
	popq	%rbp
	popq	%rdi
	popq	%rsi
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	retq
	.seh_handlerdata
	.text
	.seh_endproc
                                        # -- End function
	.addrsig
