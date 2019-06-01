;
; macro inizializzazione codice NASM a 64 bit
;

%macro	go		0
		push	rbp				; salva il Base Pointer
		mov		rbp, rsp		; il Base Pointer punta al Record di Attivazione corrente
		pushaq					; salva i registri generali
%endmacro

%macro	end		0
		popaq					; ripristina i registri generali
		mov		rsp, rbp		; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
%endmacro
