;
; macro inizializzazione codice NASM 32
;

%macro	go	0
	push ebp	; salva il Base Pointer
	mov ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push ebx	; salva i registri da preservare
	push esi
	push edi
%endmacro

%macro	end		0
	pop edi		; ripristina i registri da preservare
	pop esi
	pop ebx
	mov esp, ebp	; ripristina lo Stack Pointer
	pop ebp		; ripristina il Base Pointer
	ret		; ritorna alla funzione chiamante
%endmacro
