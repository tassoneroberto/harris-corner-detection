%include "init.nasm"
%include "sseutils64.nasm"

extern printf

section .data		
f	equ	8	; matrice f ->rdi
g	equ	16	; matrice g ->rsi
n	equ	24	; rows->rdx
m	equ	32	; cols->rcx
ris	equ 	40	; matrice ris -> r8

dim 	equ	8	; dimensione ogni locazione

p 	equ	4	; parallel

pdim	equ 	32
pdim2 	equ	64
pdim3	equ	96

section .bss

section .text			


global prodotto	

prodotto: go


	mov rax,rdx
	imul rax,rcx	; n*m
	mov rbx,32
	xor rdx,rdx
	idiv rbx		; rax=quoziente , rdx=resto	
	mov rcx,rax		; rcx=quoziente (for quoziente)

	;add rsi, rbp		; indirizzo di g
	;add rdi, rbp		; indirizzo di f
	;add r8, rbp

	cmp rcx, 1
	jl .cicloR

.cicloQ:	
	vmovups ymm1, [rdi]		; prelevo f[i:i+3]
	vmovups ymm2, [rsi]		; prelevo g[i:i+3]	
	vmulps ymm1,ymm2
	vmovups [r8],ymm1

	vmovups ymm3, [rdi+pdim]	; prelevo f[i:i+3]
	vmovups ymm4, [rsi+pdim]	; prelevo g[i:i+3]	
	vmulps ymm3,ymm4
	vmovups [r8+pdim],ymm3

	vmovups ymm5, [rdi+pdim2]	; prelevo f[i:i+3]
	vmovups ymm6, [rsi+pdim2]	; prelevo g[i:i+3]	
	vmulps ymm5,ymm6
	vmovups [r8+pdim2],ymm5

	vmovups ymm0, [rdi+pdim3]	; prelevo f[i:i+3]
	vmovups ymm7, [rsi+pdim3]	; prelevo g[i:i+3]	
	vmulps ymm0,ymm7
	vmovups [r8+pdim3],ymm0

	add rsi, 128
	add rdi, 128
	add r8, 128
	sub rcx,1
	cmp rcx,0
	jg .cicloQ
	
	cmp rdx, 0
	jle .fine
	
.cicloR:
	vmovss xmm1, [rdi]	; prelevo f[i]
	vmovss xmm2, [rsi]	; prelevo g[i]	
	vmulss xmm1,xmm2
	vmovss [r8],xmm1

	add rsi, 4
	add rdi, 4
	add r8, 4
	sub rdx,1
	jg .cicloR




.fine:

end
