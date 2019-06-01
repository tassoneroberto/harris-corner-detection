%include "init.nasm"
%include "sseutils64.nasm"

section .data		
f	equ	8	; matrice f -> rdi
n	equ	16	; rows -> rsi
m	equ	24	; cols -> rdx
s	equ 	32	; soglia -> xmm0

dim 	equ	4	; dimensione ogni locazione

section .bss

section .text			

global sogliatura	

sogliatura: go

	imul rsi,rdx ;rsi=n*m
	xor rbx,rbx
	;movss xmm0,[rbp+s]
.for:	
	movss xmm1, [rdi]	; prelevo f[i]
	cmpss xmm1,xmm0,5
	movq rbx,xmm1
	cmp rbx,0
	jg .nonazzero
	mov [rdi], dword 0
.nonazzero:
	add rdi, dim
	dec rsi
	jnz .for

end


	;movss xmm1, [rdi]	; prelevo f[i]
	;movss xmm2, [rsi]	; prelevo g[i]	
	;mulss xmm1,xmm2
	;movss [r8],xmm1

	;add rsi, 4
	;add rdi, 4
	;add r8, 4
	;sub rdx,1
	;jg .cicloR

