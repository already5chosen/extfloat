  .section	.rdata,"dr"
	.p2align	4               # @_ZZL7rsqrt64yiE11exp_adj_tab
_rsqrt64_tab:
  .byte 252, 178, 245, 173, 238, 168, 232, 164
  .byte 226, 160, 221, 156, 216, 153, 211, 149
  .byte 207, 146, 203, 143, 199, 141, 195, 138
  .byte 192, 136, 189, 133, 185, 131, 182, 129
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
asm_sqrt448:                      # @asm_sqrt448
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
	.seh_endprologue

  # rcx = dst, rdx = src, r8d = exp1
	movq	48(%rdx), %rax       # rax = src[6]
  movq  %rcx, 72(%rsp)       # save dst pointer at RCX home
	movq	%rax, %rbx           # rbx = src[6]
	movq	%rax, %rsi           # rsi = src[6]
	leaq	_rsqrt64_tab(%rip), %rdi
	shr   $58, %rax
	and   $30, %eax            # rax = idx = (x >> 58) & 30;
  add   %r8d,%eax            # rax = idx = ((x >> 58) & 30) + exp1;
	movzbl (%rax,%rdi), %eax   # rax = y = rsqrt_tab[idx]; y scaled by 2^8
  mov   %r8d, %ecx           # rcx = exp1
  xor   $1,   %ecx           # ecx = exp0 = exp1^1
  shr  %cl,  %rbx            # rbx = x = ssrc[7] = src[6] >> exp0

  # 1st NR step - y = y*(3+eps1 - y*y*x)/2
  mov  %eax, %r9d             # r9  = y
  imul %eax, %eax             # rax = yy = y*y;    scaled by 2^16
  mov  %rbx, %r8              # r8  = x
	movq	  (%rdx), %rbp        # rbp = src[0]
  shr  $48,  %rbx             # rbx = x >> (64-16)
  imul %ebx, %eax             # rax = yyx = yy*(x >> (64-16)) = y*y*x scaled by 2^(62+16+16-64)=2^30
  NR1_3 = ((3)<<30)+((3)<<17);
  mov  $(NR1_3), %ebx
	movq	 8(%rdx), %rdi        # rdi = src[1]
  sub  %eax, %ebx             # rbx = (3u << 30) + (3u<<17) - yyx
  imul %rbx, %r9
  shr  $7,   %r9              # r9  = y = (y * ((3u << 30)+(3u<<17)-yyx))>>(6+1); scaled by 2^(8+30-6)=2^32
	movq	16(%rdx), %r10        # r10 = src[2]

  # 2nd NR step - y = y*(3+eps2 - y*y*x)/2
  mov  %r9,  %rax             # rax = y
  imul %r9,  %r9              # r9  = yy=y*y       scaled by 2^64
  mov  %r8,  %rbx             # rbx = x
	movq	24(%rdx), %r11        # r11 = src[3]
  shr  $32,  %r9              # r9  = yy=(y*y)>>32 scaled by 2^32
  shr  $32,  %r8              # r8  = x >> 32;     scaled by 2^30
  imul %r8,  %r9
	movq	32(%rdx), %r12        # r12 = src[4]
  shr  $32,  %r9              # r9  = yyx = (yy*(x >> 32)) >> 33 = y*y*x scaled by 2^(32+62-32-32)=2^30
  NR2_3 = ((3)<<30)+((3)<<4);
  mov  $(NR2_3), %r8d
  sub  %r9d, %r8d             # r8  = (3u << 30)+(3u<<4)-yyx
  imul %r8,  %rax             # rax = y = y * ((3u << 30)+(3u<<4)-yyx)
  shr  $31,  %rax             # rax = y = (y * ((3u << 30)+(3u<<4)-yyx))>>(30+1); scaled by 2^(32+30-30)=2^32
	movq	40(%rdx), %r13        # r13 = src[5]

  # 3rd NR step
  # Use 2nd order polynomial: y =  y - y*(err/2 - 3/8*err**2) = y - y/2*(err - 3/4*err**2)
  mov  %rax, %r8              # r8  = y
  imul %rax, %rax             # rax = y*y, scaled by 2^64
  mul  %rbx                   # rdx:rax - y*y*x, scaled by 2^(64+62)=2^126
  shrd $62,  %rdx, %rax       # rax - err = y*y*x >> 62 = y*y*x-1 scaled by 2**64
  mov  %rax, %r9              # r9  - err
  imul %rax                   # rdx = err2 = imulh(err, err) = err*err scaled by 2**64
  lea  (%rdx,%rdx,2),%rax     # rax = err2*3
  shl  $2,   %r9              # r9  = err*4
  sub  %r9,  %rax             # rax = m2 = err2*3 - err*4 = 3*err**2 - err*4 scaled by 2**64
  mov  %r8,  %r9              # r9  = y
  shl  $29,  %r8              # r8  = y << 29
  imul %r8                    # rdx = adj = imulh(m2, int64_t(y<<29)) = m2*(y/8)
  shl  $32,  %r9              # r9  = y << 32
  add  %rdx, %r9              # r9  = Rsqrt[1] = (y<<32) + adj;

  # Calculate sqrt() as src*invSqrt
  movq  %r9,  %rax            # rax = Rsqrt[1]
  mulq  %rbx                  # rdx = mulh(Rsqrt[1], ssrc[7])
  lea   2(%rdx,%rdx),%rdx     # rdx = Sqrt[7] = mulh(Rsqrt[1], ssrc[7])*2 + 2
  movq  %rdx, %r15            # r15 = Sqrt[7]

  # Registers usage
  # rax - free
  # rdx - free (copy of Sqrt[7])
  # rcx - exp0 = exp1^1
  # rbx - ssrc[7]
  # rsi - src[6]
  # rdi - src[1]
  # rbp - src[0]
  # r8  - free
  # r9  - Rsqrt[1]
  # r10 - src[2]
  # r11 - src[3]
  # r12 - src[4]
  # r13 - src[5]
  # r14 - free
  # r15 - Sqrt[7]

  xor   %eax, %eax           # rax = 0
	shrdq	%cl, %rbp, %rax      # rax = ssrc[0] = ((src[0]:0) >> exp0).LO
	shrdq	%cl, %rdi, %rbp      # rbp = ssrc[1] = ((src[1:0]) >> exp0).LO
	shrdq	%cl, %r10, %rdi      # rdi = ssrc[2] = ((src[2:1]) >> exp0).LO
	shrdq	%cl, %r11, %r10      # r10 = ssrc[3] = ((src[3:2]) >> exp0).LO
	shrdq	%cl, %r12, %r11      # r11 = ssrc[4] = ((src[4:3]) >> exp0).LO
	shrdq	%cl, %r13, %r12      # r12 = ssrc[5] = ((src[5:4]) >> exp0).LO
	shrdq	%cl, %rsi, %r13      # r13 = ssrc[6] = ((src[6:5]) >> exp0).LO
  mov   %rbx,%rsi            # rsi = ssrc[7] = src[6] >> exp0
  mov   %r9, %rbx            # rbx = Rsqrt[1]

  mov  %rax, 80(%rsp)        # save ssrc[0] at RDX home
  mov  %rbp, 88(%rsp)        # save ssrc[1] at R8  home

  # Registers usage
  # rax - free
  # rdx - free (copy of Sqrt[7])
  # rcx - free
  # rbx - Rsqrt[1]
  # rsi - ssrc[7]
  # rdi - ssrc[2]
  # rbp - free
  # r8  - free
  # r9  - free (copy of Rsqrt[1])
  # r10 - ssrc[3]
  # r11 - ssrc[4]
  # r12 - ssrc[5]
  # r13 - ssrc[6]
  # r14 - free
  # r15 - Sqrt[7]

  # improve precision of Sqrt to 64*2-eps bits
  mov   %rdx, %rax            # rax = Sqrt[7]
  mulq  %rax                  # rdx:rax = sqr(Sqrt[7]) scaled by 2**126
  # Sqrt_sqr = err = Sqrt**2 - ssrc[7:6]
  sub   %r13, %rax
  sbb   %rsi, %rdx
  # rdx:rax = err = Sqrt**2 - src
  call  .smul2
  # rdx:rax = err*invS/2, scaled by 2**127
  # Sqrt -= err*invS/2
  neg   %rax                  # rax = Sqrt[6] = -adj[0]
  sbb   %rdx, %r15            # r15 = Sqrt[7] -= adj[1] - borrow
  mov   %rax, %r14            # r14 = Sqrt[6]

  # improve precision of Rsqrt to 64*2-eps bits
  # Multiply sqrt[7:6] by Rsqrt[1], take upper 128 bits of the product
  mul   %rbx                  # rdx:rax = Rsqrt[1]*Sqrt[6]
  mov   %r15, %rax            # rax     = Sqrt[7]
  mov   %rdx, %rcx            # rcx     = mulh(invS[3],Sqrt[6])
  mulq  %rbx                  # rdx:rax = Rsqrt[1]*Sqrt[7]
  add   %rcx, %rax
  adc   $0,   %rdx            # rdx:rax = Rsqrt[1]*Sqrt[7] + mulh(Rsqrt[1],Sqrt[6])
  # rdx:rax = Rsqrt_err = sqrt*rsqrt, scaled by 2**127
	add  	%rax, %rax
  adc   %rdx, %rdx            # rdx:rax = rdx:rax << 1 = Rsqrt_err = sqrt*rsqrt, scaled by 2**128
  call  .smul2
  # rdx:rax = adj = Rsqrt_err*Rsqrt[1], scaled by 2**127
  # Rsqrt[1:0] -= adj
  neg %rax                    # rax = Rsqrt[0] = -adj[0]
  sbb %rdx, %rbx              # rbx = Rsqrt[1] -= adj[1] - borrow
  mov %rax, %rsi              # rsi = Rsqrt[0]
  mov %rbx, %rcx              # rcx = Rsqrt[1]

  # Registers usage
  # rax - free
  # rdx - free
  # rdi - ssrc[2]
  # rbx - free
  # rsi - Rsqrt[0]
  # rcx - Rsqrt[1]
  # rbp - free
  # r8  - free
  # r9  - free
  # r10 - ssrc[3]
  # r11 - ssrc[4]
  # r12 - ssrc[5]
  # r13 - ssrc[6]
  # r14 - Sqrt[6]
  # r15 - Sqrt[7]

  # improve precision to 64*4-eps bits
  # Sqrt_sqr[2:0] = sqr(Sqrt[7:6])[2:0], scaled by 2**(64*3-2)
  mov   %r14, %rax            # rax = Sqrt[6]
  mul   %rax                  # rdx:rax = Sqrt[6]*Sqrt[6]
  mov   %rax, %r8             # r8  = Sqrt_sqr[0] = (Sqrt[6]*Sqrt[6]).LO
  mov   %rdx, %r9             # r9  = Sqrt_sqr[1] = (Sqrt[6]*Sqrt[6]).HI
  mov   %r15, %rbx            # rbx = Sqrt[7]
  imul  %rbx, %rbx            # rbx = (Sqrt[7]*Sqrt[7]).LO
  mov   %r15, %rax            # rax = Sqrt[7]
  mul   %r14                  # rdx:rax = Sqrt[7]*Sqrt[6]
  add   %rax, %rax
  adc   %rdx, %rdx            # rdx:rax = Sqrt[7]*Sqrt[6]*2
  add   %rax, %r9             # r9  = Sqrt_sqr[1] += (Sqrt[7]*Sqrt[6]*2).LO
  adc   %rbx, %rdx            # rdx = Sqrt_sqr[2] += (Sqrt[7]*Sqrt[7]*2).LO + (Sqrt[7]*Sqrt[7]).LO + carry
  # rdx:r9:r8   - sqr[2:0],   scaled by 2**(64*3-2)
  mov   %r11, %rax            # rax = ssrc[4]
  sub   %r8,  %rax            # rax = err[0] = ssrc[4] - sqr[0]
  sbb   %r9,  %r12            # rbp = err[1] = ssrc[5] - sqr[1] - borrow
  sbb   %rdx, %r13            # rdi = err[2] = ssrc[5] - sqr[2] - borrow
  # r13:r12:rax - err[2:0], scaled by 2**(64*3-2)
  shrd  $6,   %r12, %rax
  shrd  $6,   %r13, %r12
  # r12:rax - err[1:0], scaled by 2**(64*3-2-6)

  # mul2sX2u - multiply signed 128-bit number by unsigned 128-bit number.
  # return upper 128 bits of the product
  # inputs:
  # r12:rax - a[1:0], signed multiplicand
  # rcx:rsi - b[1:0], unsigned multiplicand
  # outputs:
  # r13:r8   - res[1:0], signed result
  # rcx:rsi - b[1:0], (preserved)
  # r12     - sign(adj)
  # Other registers affected:
  # rdx
  # start unsigned multiplication
  mul   %rcx                  # rdx:rax = a[0]*b[1]
  mov   %rdx, %r8             # r8      = res[0] = (a[0]*b[1]).HI
  mov   %r12, %rax            # rax     = a[1]
  mul   %rsi                  # rdx:rax = a[1]*b[0]
  xor   %r13d,%r13d           # r13     = res[1] = 0
  add   %rdx, %r8             # r8      = res[0] += (a[1]*b[0]).HI
  adc   %r13, %r13            # r13     = res[1] += carry
  mov   %r12, %rax            # rax     = a[1]
  mul   %rcx                  # rdx:rax = a[1]*b[1]
  add   %rax, %r8             # r8      = res[0] += (a[1]*b[1]).LO
  adc   %rdx, %r13            # r13     = res[1] += (a[1]*b[1]).HI + carry
  # unsigned multiplication done

  # convert unsigned product to signed
  sar   $63,  %r12            # r12 = sign(adj) = (err[1] < 0) ? -1 : 0
  mov   %rsi, %rax
  and   %r12, %rax            # rax = b[0] & sign(adj)
  mov   %rcx, %rdx
  and   %r12, %rdx            # rdx = b[1] & sign(adj)
  sub   %rax, %r8             # r8  = res[0] -= (b[0] & sign(adj))
  sbb   %rdx, %r13            # r13 = res[1] -= (b[1] & sign(adj)) - borrow

  # r13:r8 = res[1:0] = adj[1:0] = err*Rsqrt[1:0], scaled by 2**(64*3-1-6)
  mov   %r13, %rdx
  shld  $6,   %r8, %r13
  sar   $58,  %rdx
  shl   $6,   %r8
  # rdx:r13:r8   - adj[2:0], scaled by 2**(64*3-1)
  add   %rdx,%r14             # r14 = Sqrt[6] += adj[2]
  adc   %r12,%r15             # r15 = Sqrt[7] += sign(adj) + carry
  mov   %r8, %r12             # r12 = Sqrt[4]  = adj[0]
                              # r13 = Sqrt[5]  = adj[1]

  # Registers usage
  # rax - free
  # rdx - free
  # rdi - ssrc[2]
  # rbx - free
  # rsi - Rsqrt[0]
  # rcx - Rsqrt[1]
  # rbp - free
  # r8  - free (copy of Sqrt[4])
  # r9  - free
  # r10 - ssrc[3]
  # r11 - ssrc[4]
  # r12 - Sqrt[4]
  # r13 - Sqrt[5]
  # r14 - Sqrt[6]
  # r15 - Sqrt[7]

  # improve precision to 64*6-eps bits
  # Sqrt_sqr[2:0] = sqr(Sqrt[7:4])[4:2], scaled by 2**(64*3-2)
  # non-diagonal
  mov   %r8,  %rax            # rax     = Sqrt[4]
  mul   %r13                  # rdx:rax = Sqrt[4]*Sqrt[5]
  mov   %rdx, %r8             # r8      = Sqrt_sqr[0] = (Sqrt[4]*Sqrt[5]).HI
  mov   %r12, %rax            # rax     = Sqrt[4]
  mul   %r14                  # rdx:rax = Sqrt[4]*Sqrt[6]
  xor   %r9d, %r9d            # r9      = Sqrt_sqr[1] = 0
  add   %rax, %r8             # r8      = Sqrt_sqr[0] += (Sqrt[4]*Sqrt[6]).LO
  adc   %rdx, %r9             # r9      = Sqrt_sqr[1] += (Sqrt[4]*Sqrt[6]).HI + carry
  mov   %r12, %rax            # rax     = Sqrt[4]
  mul   %r15                  # rdx:rax = Sqrt[4]*Sqrt[7]
  xor   %ebx, %ebx            # rbx     = Sqrt_sqr[2] = 0
  add   %rax, %r9             # r9      = Sqrt_sqr[1] += (Sqrt[4]*Sqrt[7]).LO
  adc   %rdx, %rbx            # rbx     = Sqrt_sqr[2] += (Sqrt[4]*Sqrt[7]).HI + carry
  mov   %r13, %rax            # rax     = Sqrt[5]
  mul   %r14                  # rdx:rax = Sqrt[5]*Sqrt[6]
  add   %rax, %r9             # r9      = Sqrt_sqr[1] += (Sqrt[5]*Sqrt[6]).LO
  adc   %rdx, %rbx            # rbx     = Sqrt_sqr[2] += (Sqrt[5]*Sqrt[6]).HI + carry
  mov   %r13, %rax            # rax     = Sqrt[5]
  imul  %r15, %rax            # rax     = (Sqrt[5]*Sqrt[7]).LO
  add   %rax, %rbx            # rbx     = Sqrt_sqr[2] += (Sqrt[5]*Sqrt[7]).LO

  add   %r8,  %r8
  adc   %r9,  %r9
  adc   %rbx, %rbx            # non-diagonal *= 2

  mov   %r13, %rax            # rax     = Sqrt[5]
  mul   %rax                  # rdx:rax = Sqrt[5]*Sqrt[5]
  mov   %r14, %rbp            # rbp     = Sqrt[6]
  imul  %rbp, %rbp            # rbp     = (Sqrt[6]*Sqrt[6]).LO

  add   %rax, %r8
  adc   %rdx, %r9
  adc   %rbp, %rbx            # non-diagonal += diagonal
  # rbx:r9:r8   - sqr[2:0],   scaled by 2**(64*3-2)

  mov %rdi,     %rax          # rax = ssrc[2]
  sub %r8,      %rax          # rax = err[0] = ssrc[2] - sqr[0]
  sbb %r9,      %r10          # rbp = err[1] = ssrc[3] - sqr[1] - borrow
  sbb %rbx,     %r11          # rdi = err[2] = ssrc[4] - sqr[2] - borrow
  # r11:r10:rax - err[2:0], scaled by 2**(64*3-2)
  shrd $10, %r10, %rax
  shrd $10, %r11, %r10
  # r10:rax - err[1:0], scaled by 2**(64*3-2-10)

  # mul2sX2u - multiply signed 128-bit number by unsigned 128-bit number.
  # return upper 128 bits of the product
  # inputs:
  # r10:rax - a[1:0], signed multiplicand
  # rcx:rsi - b[1:0], unsigned multiplicand
  # outputs:
  # r11:r8   - res[1:0], signed result
  # rcx:rsi - b[1:0], (preserved)
  # r10     - sign(adj)
  # Other registers affected:
  # rdx
  # start unsigned multiplication
  mul   %rcx                  # rdx:rax = a[0]*b[1]
  mov   %rdx, %r8             # r8      = res[0] = (a[0]*b[1]).HI
  mov   %r10, %rax            # rax     = a[1]
  mul   %rsi                  # rdx:rax = a[1]*b[0]
  xor   %r11d,%r11d           # r11      = res[1] = 0
  add   %rdx, %r8             # r8      = res[0] += (a[1]*b[0]).HI
  adc   %r11, %r11            # r11      = res[1] += carry
  mov   %r10, %rax            # rax     = a[1]
  mul   %rcx                  # rdx:rax = a[1]*b[1]
  add   %rax, %r8             # r8      = res[0] += (a[1]*b[1]).LO
  adc   %rdx, %r11            # r11      = res[1] += (a[1]*b[1]).HI + carry
  # unsigned multiplication done

  # convert unsigned product to signed
  sar   $63,  %r10            # r10 = sign(adj) = (err[1] < 0) ? -1 : 0
  mov   %rsi, %rax
  and   %r10, %rax            # rax = b[0] & sign(adj)
  mov   %rcx, %rdx
  and   %r10, %rdx            # rdx = b[1] & sign(adj)
  sub   %rax, %r8             # r8  = res[0] -= (b[0] & sign(adj))
  sbb   %rdx, %r11            # r11  = res[1] -= (b[1] & sign(adj)) - borrow
  # r11:r8   - adj[1:0] = err*Rsqrt[1:0], scaled by 2**(64*3-1-10)

  mov   %r11, %rdx
  shld  $10,  %r8, %r11
  shl   $10,  %r8
  sar   $54,  %rdx
  # rdx:r11:r8 - adj[2:0], scaled by 2**(64*3-1)
  add   %rdx, %r12            # r12 = Sqrt[4] += adj[2]
  adc   $0,   %r10            # r10 = incdec = sign(adj) + carry
  jnz   .incdec_sqrt_5to7
  .incdec_sqrt_5to7_done:
  mov   %r8, %r10             # r10 = Sqrt[2]  = adj[0]
                              # r11 = Sqrt[3]  = adj[1]

  # Registers usage
  # rax - free
  # rdx - free
  # rdi - ssrc[2]
  # rbx - free
  # rsi - Rsqrt[0]
  # rcx - Rsqrt[1]
  # rbp - free
  # r8  - free (copy of Sqrt[2])
  # r9  - free
  # r10 - Sqrt[2]
  # r11 - Sqrt[3]
  # r12 - Sqrt[4]
  # r13 - Sqrt[5]
  # r14 - Sqrt[6]
  # r15 - Sqrt[7]

  # improve precision to 64*6-eps bits
  # Sqrt_sqr[2:0] = sqr(Sqrt[7:2])[6:4], scaled by 2**(64*3-2)
  # non-diagonal
  mov   %r8,  %rax            # rax     = Sqrt[2]
  mul   %r13                  # rdx:rax = Sqrt[2]*Sqrt[5]
  mov   %rdx, %r8             # r8      = Sqrt_sqr[0] = (Sqrt[2]*Sqrt[5]).HI
  mov   %r10, %rax            # rax     = Sqrt[2]
  mul   %r14                  # rdx:rax = Sqrt[2]*Sqrt[6]
  xor   %r9d, %r9d            # r9      = Sqrt_sqr[1] = 0
  add   %rax, %r8             # r8      = Sqrt_sqr[0] += (Sqrt[2]*Sqrt[6]).LO
  adc   %rdx, %r9             # r9      = Sqrt_sqr[1] += (Sqrt[2]*Sqrt[6]).HI + carry
  mov   %r11, %rax            # rax     = Sqrt[3]
  mul   %r12                  # rdx:rax = Sqrt[3]*Sqrt[4]
  add   %rdx, %r8             # r8      = Sqrt_sqr[0] += (Sqrt[3]*Sqrt[4]).HI
  adc   $0,   %r9             # r9      = Sqrt_sqr[1] += carry
  mov   %r11, %rax            # rax     = Sqrt[3]
  mul   %r13                  # rdx:rax = Sqrt[3]*Sqrt[5]
  xor   %ebx, %ebx            # rbx     = Sqrt_sqr[2] = 0
  add   %rax, %r8             # r8      = Sqrt_sqr[0] += (Sqrt[3]*Sqrt[5]).LO
  adc   %rdx, %r9             # r9      = Sqrt_sqr[1] += (Sqrt[3]*Sqrt[5]).HI + carry
  adc   $0,   %rbx            # rbx     = Sqrt_sqr[2] += carry
  mov   %r10, %rax            # rax     = Sqrt[2]
  mul   %r15                  # rdx:rax = Sqrt[2]*Sqrt[7]
  add   %rax, %r9             # r9      = Sqrt_sqr[1] += (Sqrt[2]*Sqrt[7]).LO
  adc   %rdx, %rbx            # rbx     = Sqrt_sqr[2] += (Sqrt[2]*Sqrt[7]).HI + carry
  mov   %r11, %rax            # rax     = Sqrt[3]
  mul   %r14                  # rdx:rax = Sqrt[3]*Sqrt[6]
  add   %rax, %r9             # r9      = Sqrt_sqr[1] += (Sqrt[3]*Sqrt[6]).LO
  adc   %rdx, %rbx            # rbx     = Sqrt_sqr[2] += (Sqrt[3]*Sqrt[6]).HI + carry
  mov   %r12, %rax            # rax     = Sqrt[4]
  mul   %r13                  # rdx:rax = Sqrt[4]*Sqrt[5]
  add   %rax, %r9             # r9      = Sqrt_sqr[1] += (Sqrt[4]*Sqrt[5]).LO
  adc   %rdx, %rbx            # rbx     = Sqrt_sqr[2] += (Sqrt[4]*Sqrt[5]).HI + carry
  mov   %r11, %rax            # rax     = Sqrt[3]
  imul  %r15, %rax            # rax     = (Sqrt[3]*Sqrt[7]).LO
  add   %rax, %rbx            # rbx     = Sqrt_sqr[2] += (Sqrt[3]*Sqrt[7]).LO
  mov   %r12, %rax            # rax     = Sqrt[4]
  imul  %r14, %rax            # rax     = (Sqrt[4]*Sqrt[6]).LO
  add   %rax, %rbx            # rbx     = Sqrt_sqr[2] += (Sqrt[4]*Sqrt[6]).LO

  add   %r8,  %r8
  adc   %r9,  %r9
  adc   %rbx, %rbx            # non-diagonal *= 2

  # diagonal elements
  mov   %r12, %rax            # rax     = Sqrt[4]
  mul   %rax                  # rdx:rax = Sqrt[4]*Sqrt[4]
  mov   %r13, %rbp            # rbp     = Sqrt[5]
  imul  %rbp, %rbp            # rbp     = (Sqrt[5]*Sqrt[5]).LO

  add   %rax, %r8
  adc   %rdx, %r9
  adc   %rbp, %rbx            # non-diagonal += diagonal
  # rbx:r9:r8   - sqr[2:0],   scaled by 2**(64*3-2)

  mov 80(%rsp), %rax          # rax = ssrc[0]
  mov 88(%rsp), %rbp          # rbp = ssrc[1]
  sub %r8,      %rax          # rax = err[0] = ssrc[0] - sqr[0]
  sbb %r9,      %rbp          # rbp = err[1] = ssrc[1] - sqr[1] - borrow
  sbb %rbx,     %rdi          # rdi = err[2] = ssrc[2] - sqr[2] - borrow
  # rdi:rbp:rax - err[2:0], scaled by 2**(64*3-2)
  shrd $14, %rbp, %rax
  shrd $14, %rdi, %rbp
  # rbp:rax - err[1:0], scaled by 2**(64*3-2-14)

  # mul2sX2u - multiply signed 128-bit number by unsigned 128-bit number.
  # return upper 128 bits of the product
  # inputs:
  # rbp:rax - a[1:0], signed multiplicand
  # rcx:rsi - b[1:0], unsigned multiplicand
  # outputs:
  # r9:r8   - res[1:0], signed result
  # rcx:rsi - b[1:0], (preserved)
  # rbp     - sign(adj)
  # Other registers affected:
  # rdx
  # start unsigned multiplication
  mul   %rcx                  # rdx:rax = a[0]*b[1]
  mov   %rdx, %r8             # r8      = res[0] = (a[0]*b[1]).HI
  mov   %rbp, %rax            # rax     = a[1]
  mul   %rsi                  # rdx:rax = a[1]*b[0]
  xor   %r9d, %r9d            # r9      = res[1] = 0
  add   %rdx, %r8             # r8      = res[0] += (a[1]*b[0]).HI
  adc   %r9,  %r9             # r9      = res[1] += carry
  mov   %rbp, %rax            # rax     = a[1]
  mul   %rcx                  # rdx:rax = a[1]*b[1]
  add   %rax, %r8             # r8      = res[0] += (a[1]*b[1]).LO
  adc   %rdx, %r9             # r9      = res[1] += (a[1]*b[1]).HI + carry
  # unsigned multiplication done

  # convert unsigned product to signed
  sar   $63,  %rbp            # rbp = sign(adj) = (err[1] < 0) ? -1 : 0
  mov   %rsi, %rax
  and   %rbp, %rax            # rax = b[0] & sign(adj)
  mov   %rcx, %rdx
  and   %rbp, %rdx            # rdx = b[1] & sign(adj)
  sub   %rax, %r8             # r8  = res[0] -= (b[0] & sign(adj))
  sbb   %rdx, %r9             # r9  = res[1] -= (b[1] & sign(adj)) - borrow

  # r9:r8   - adj[1:0] = err*Rsqrt[1:0], scaled by 2**(64*3-1-14)
  # rbp     - sign(adj)
  mov  %r9, %rdx
  shld $14, %r8, %r9
  shl  $14, %r8
  sar  $50, %rdx
  # rdx:r9:r8   - adj[2:0] = err*Rsqrt[1:0], scaled by 2**(64*3-1)
                              # r8  = Sqrt[0]  = adj[0]
                              # r9  = Sqrt[1]  = adj[1]
  add   %rdx,%r10             # r10 = Sqrt[2] += adj[2]
  adc   $0,  %rbp             # rbp = incdec = sign(adj) + carry
  jnz   .incdec_sqrt_3to7
  .incdec_sqrt_3to7_done:

  # Registers usage
  # rax - free
  # rdx - free
  # rcx - free
  # rbx - free
  # rsi - free
  # rdi - free
  # rbp - free
  # r8  - Sqrt[0] = lsw
  # r9  - Sqrt[1]
  # r10 - Sqrt[2]
  # r11 - Sqrt[3]
  # r12 - Sqrt[4]
  # r13 - Sqrt[5]
  # r14 - Sqrt[6]
  # r15 - Sqrt[7]

  mov 72(%rsp), %rcx          # rcx = dst, restore from RCX home
  # check if lsw is close to the middle of qword range
  UINT64_MID = (1) << 63;
  MAX_ERR    = (1) << 17;     # probably less, but it does not cost much to be on the safe side
	movabsq	$(UINT64_MID - MAX_ERR),%rax
  mov %r8, %rdx
  sub %rax, %rdx              # rdx = lsw - (UINT64_MID-MAX_ERR)
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
    mov %r9,   0(%rcx)
    mov %r10,  8(%rcx)
    mov %r11, 16(%rcx)
    mov %r12, 24(%rcx)
    mov %r13, 32(%rcx)
    mov %r14, 40(%rcx)
    mov %r15, 48(%rcx)

    # r9   = acc0 = a[0]
    xor %ebp, %ebp                # rbp  = acc1
    mov $1,   %ebx                # ebx  = K = 1
    .sqr_tip_point_k:
      xor %esi, %esi              # rsi = acc2
      add (%rcx,%rbx,8), %rbp     # rbp = acc1 += a[K];
      adc %rsi, %rsi              # rsi = acc2 += carry
      xor %r8d, %r8d              # r8d = i = 0
      lea -1(%rbx), %rdi          # edi = j = K-1
      .sqr_tip_point_j:
        mov  (%rcx,%r8, 8),%rax   # rax = a[i]
        mulq (%rcx,%rdi,8)        # rdx:rax = a[i]*a[j]
        add %rax, %r9             # r9  = acc0 += (a[i]*a[j]).LO
        adc %rdx, %rbp            # rbp = acc1 += (a[i]*a[j]).HI + carry
        adc $0,   %rsi            # rsi = acc2 += carry
        inc %r8d                  # ++i
      sub $1, %edi                # --j
      jge .sqr_tip_point_j        # j >= 0

      mov %rbp, %r9               # r9  = acc0 = acc1
      mov %rsi, %rbp              # rbp = acc1 = acc2
    inc %ebx                      # ++K
    cmp $7, %ebx                  # while (K < 7)
    jb  .sqr_tip_point_k

    xor %ebx, %ebx                # ebx = i = 0
    mov $6,   %edi                # edi = j = 6
    .sqr_tip_point_j2:
      mov  (%rcx,%rbx,8),%rax     # rax = a[i]
      imul (%rcx,%rdi,8),%rax     # rax = (a[i]*a[j]).LO
      add   %rax, %r9             # r9  = acc0 += (a[i]*a[j]).LO
      inc %ebx                    # ++i
    sub $1, %edi                  # --j
    jge .sqr_tip_point_j2         # j >= 0

    # r9  = TP_sqr[8]

    lea (%r9 ,%r9 ), %r8          # r8= lsw = acc0 << 1, that moves bit 62 in position 63, where it inspected by the rest of the code
    mov 0(%rcx), %r9              # restore r9=Sqrt[1]
  .round_and_copy:
  xor %eax,%eax                   # rax = 0
  add %r8, %r8                    # carry = bit 63 of lsw
  # add carry
  adc %rax, %r9
  jnc .store_results
  adc %rax, %r10
  adc %rax, %r11
  adc %rax, %r12
  adc %rax, %r13
  adc %rax, %r14
  adc %rax, %r15
  # store result at dst
  .store_results:
  mov %r9,   0(%rcx)
  mov %r10,  8(%rcx)
  mov %r11, 16(%rcx)
  mov %r12, 24(%rcx)
  mov %r13, 32(%rcx)
  mov %r14, 40(%rcx)
  mov %r15, 48(%rcx)

	popq	%rbx
	popq	%rbp
	popq	%rdi
	popq	%rsi
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	retq

.p2align	4, 0x90
.smul2:
# signed multiply 128-bit number a, by 64-bit number b, return upper 128 bits of the product
# inputs:
# rdx:rax - a[1:0]
# rbx     - b, b >= 2**63
# outputs:
# rdx:rax - result
# rbx     - b (preserved)
# other registers affected:
# rcx, rbp
  mov   %rdx, %rbp                # rbp     = a[1]
  mul   %rbx                      # rdx:rax = a[0]*b
  mov   %rdx, %rcx                # rcx     = (a[0]*b).HI
  mov   %rbp, %rax                # rax     = a[1]
  imul  %rbx                      # rdx:rax = a[1]*(b-2**64)
  add   %rcx, %rax                # rdx:rax += (a[0]*b).HI
  adc   %rbp, %rdx                # rdx:rax += a[1] + carry (emulate signed*unsigned multiplication)
  ret

.p2align	4, 0x90
  .incdec_sqrt_5to7:
  mov %r10, %rax             # rax - incdec
  sar $1,   %r10             # r10 - sign(incdec) = incdec < 0 ? -1 : 0
  add %rax, %r13             # r11 = Sqrt[5] += incdec
  adc %r10, %r14             # r14 = Sqrt[6] += sign(incdec)+carry
  adc %r10, %r15             # r15 = Sqrt[7] += sign(incdec)+carry
  jmp .incdec_sqrt_5to7_done

  .incdec_sqrt_3to7:
  mov %rbp, %rax             # rax - incdec
  sar $1,   %rbp             # rbp - sign(incdec) = incdec < 0 ? -1 : 0
  add %rax, %r11             # r11 = Sqrt[3] += incdec
  adc %rbp, %r12             # r12 = Sqrt[4] += sign(incdec)+carry
  adc %rbp, %r13             # r13 = Sqrt[5] += sign(incdec)+carry
  adc %rbp, %r14             # r14 = Sqrt[6] += sign(incdec)+carry
  adc %rbp, %r15             # r15 = Sqrt[7] += sign(incdec)+carry
  jmp .incdec_sqrt_3to7_done

	.seh_handlerdata
	.text
	.seh_endproc
                                        # -- End function
	.addrsig
