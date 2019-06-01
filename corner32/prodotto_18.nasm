%include "init.nasm"

section .data		
f	equ	8	; matrice f
g	equ	12	; matrice g
n	equ	16	; rows
m	equ	20	; cols
ris	equ 	24	; matrice ris
pdim	equ 	16
pdim2 	equ	32
pdim3	equ	48

section .bss

section .text			

global prodotto	

prodotto: go

	mov eax,[ebp+n]
	imul eax,[ebp+m]	; n*m
	mov ebx,16
	mov edx,0
	div ebx			; eax=quoziente , edx=resto	
	mov ecx,eax		; ecx=quoziente (for quoziente)

	mov esi, [ebp+f]	; indirizzo di f
	mov edi, [ebp+g]	; indirizzo di g
	mov eax, [ebp+ris]

	cmp ecx, 1
	jl .cicloR

.cicloQ:	
	movups xmm1, [esi]		; prelevo f[i:i+3]
	movups xmm2, [edi]		; prelevo g[i:i+3]	
	mulps xmm1,xmm2
	movups [eax],xmm1

	movups xmm3, [esi+pdim]	; prelevo f[i:i+3]
	movups xmm4, [edi+pdim]	; prelevo g[i:i+3]	
	mulps xmm3,xmm4
	movups [eax+pdim],xmm3

	movups xmm5, [esi+pdim2]	; prelevo f[i:i+3]
	movups xmm6, [edi+pdim2]	; prelevo g[i:i+3]	
	mulps xmm5,xmm6
	movups [eax+pdim2],xmm5

	movups xmm0, [esi+pdim3]	; prelevo f[i:i+3]
	movups xmm7, [edi+pdim3]	; prelevo g[i:i+3]	
	mulps xmm0,xmm7
	movups [eax+pdim3],xmm0

	add esi, 64
	add edi, 64
	add eax, 64
	dec ecx
	jnz .cicloQ
	
	cmp edx, 0
	jz .fine
	
.cicloR:
	movss xmm1, [esi]	; prelevo f[i]
	movss xmm2, [edi]	; prelevo g[i]	
	mulss xmm1,xmm2
	movss [eax],xmm1

	add esi, 4
	add edi, 4
	add eax, 4
	dec edx
	jnz .cicloR


.fine:

end
