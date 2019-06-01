%include "init.nasm"

section .data		
f	equ	8	; matrice f
n	equ	12	; rows
m	equ	16	; cols
s	equ 	20	; soglia
dim 	equ	4	; dimensione ogni locazione

section .bss

section .text			

global sogliatura	

sogliatura: go
	mov ecx,[ebp+n];
	imul ecx,[ebp+m];	
	mov edx,[ebp+s];
	mov esi, [ebp+f]	; indirizzo di f
	mov ebx, 0
.for:	
	mov eax, [esi+ebx]	; prelevo f[i]
	sub eax, edx;
	jge .nonazzero
	mov [esi+ebx], dword 0
	;movss xmm4, qword[0]
.nonazzero:
	add ebx, dim
	dec ecx
	jnz .for
end
