/**************************************************************************************
 *
 * CdL Magistrale in Ingegneria Informatica
 * Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2016/17
 *
 * Progetto di un algoritmo di Corner Detection
 * in linguaggio assembly x86-64 + AVX
 *
 * Fabrizio Angiulli, 20 aprile 2016
 *
 **************************************************************************************/

/*

 Software necessario per l'esecuzione:

     NASM (www.nasm.us)
     GCC (gcc.gnu.org)

 entrambi sono disponibili come pacchetti software
 installabili mediante il packaging tool del sistema
 operativo; per esempio, su Ubuntu, mediante i comandi:

     sudo apt-get install nasm
     sudo apt-get install gcc

 potrebbe essere necessario installare le seguenti librerie:

     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
     sudo apt-get install libc6-dev-i386

 Per generare il file eseguibile:

 nasm -f elf32 corner32.nasm && gcc -O0 -m64 -msse corner64.o corner64c.c -o corner64c && ./corner64c

 oppure

 ./runcorner64

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>
#include <omp.h>

#define	CORNER	0
#define	FILTERING 1

#define	MATRIX		float*
#define	VECTOR		float*
#define	IMAGE		float*

typedef struct {
    int x;
    int y;
} location;

typedef struct {
    int method; // +0
    char* img_name; // +4
    IMAGE img; // +8
    int N; // +12
    int M; // +16
// corner detection parameters
    float sigma; // +20
    float theta; // +24
    int t; // +28
// filtering parameters
    char* filt_name; // +32
    IMAGE filt; // +36
    int N_filt; // +40
    int M_filt; // +44
// generic parameters
    int silent; // +48
    int display; // +52
// filtering output
    IMAGE filtered_img; // +56
// corner detection output
    IMAGE R; // +60
    int n_corners; // +64
    location* corners; // +68
} params;

params* input;

float* temp1,temp2,temp3;
int corners=0;

/* Changelog */
// Prima versione: Execution time = 3.479 seconds
// Seconda versione: Execution time = 0.203 seconds (filtro separabile + altre ottimizzazioni)
// Terza versione: Execution time = 0.188 seconds (sse)
// Quarta versione: Execution time = 0.170 seconds (avx)

/* Inizio metodi */
float* create(int,int);
float* filterGauss(float*, int, int, float* , int ,float*,float* );
float* filter(float*, int, int, float* , int , int );
float filtN(float*,int, int,float*, int , int,int,int );
float filtS(float*,int, int,float*, int , int,int,int );
float filtE(float*,int, int,float*, int , int,int,int );
float filtO(float*,int, int,float*, int , int,int,int );
float filtV(float*,int, int,float*, int , int,int,int );
float filtInner(float*, int,float*, int , int,int,int );
void cancella(float*,int,int);
void print(float*,int,int);
float tr(float*,int,int);
float* smoothing(float,int);
float* smooth2(float,int);
float smooth(int,int,int,float,float);
void derive(float*,int,int,float*,float*);
void sogliatura(float*,int,int,float);
float* soppressione(float*,int , int , int );
int soppr(float* ,int , int , int , int , int );
void intorno(float*,int, int,int,int,int);
void findCorners(float*,int, int, location*);

typedef enum {false,true} bool;



// Funzione NASM x86-64 + AVX
extern void sogliatura(float*, int, int, float);
extern void prodotto(float*,float*, int, int, float*);

int r1,c1,nc,mc;

float* create(int r,int c) {
    float *m;
    m=(float*)malloc(r*c*sizeof(float));
    return m;
}

float* smoothing(float theta, int n) {
    float smooth1=1/(sqrt(2*M_PI*theta*theta));
    float smooth2=2*theta*theta;
    int np=n>>1; //centro filtro
    float somma=0; // serve per normalizzare il filtro
    float* ris=create(n,n);

    // valori semi-diagonale e semi-bordo
    for(int i=0; i<np; i++) {
	int ini=i*n+i;
	int innp=i*n+np;
        ris[ini]=smooth(i,i,np,smooth1,smooth2);
        ris[innp]=smooth(i,np,np,smooth1,smooth2);
        somma+=ris[ini]*4;
        somma+=ris[innp]*4;
    }
    // valori "triangoli" interni
    for(int j=1; j<np; j++) {

        for(int i=0; i<j; i++) {
	    int inj=i*n+j;
            ris[inj]=smooth(i,j,np,smooth1,smooth2);
            somma+=ris[inj]*8;
        }
    }
    // centro filtro
    int npnnp=np*n+np;
    ris[npnnp]=smooth(np,np,np,smooth1,smooth2);
    somma+= ris[npnnp];


    // normalizzazione
    for(int i=0; i<np; i++) {
        int ni1=n-i-1;
        int ini=i*n+i;
        ris[ini]/=somma;// diagonale in alto a sinistra
        ris[i*n+np]/=somma;// bordo in alto

        ris[(ni1)*n+(i)]=ris[ini];//diagonale in basso a sinistra
        ris[(ni1)*n+(ni1)]=ris[ini];//diagonale in basso a destra
        ris[(i)*n+(ni1)]=ris[ini];//diagonale in alto a destra
        ris[(np)*n+(i)]=ris[i*n+np];//bordo sinistro
        ris[(np)*n+(ni1)]=ris[i*n+np];//bordo destro
        ris[(ni1)*n+(np)]= ris[i*n+np];//bordo in basso

    }


    for(int j=1; j<np; j++) {

        for(int i=0; i<j; i++) {
            int nj1=n-j-1;
            int ni1=n-i-1;
            int in=i*n;
            int jn=j*n;
            int nj1n=(nj1)*n;
            int ni1n=(ni1)*n;
            int inj=in+j;

            ris[in+j]/=somma;
            ris[jn+i]=ris[inj];
            ris[nj1n+i]=ris[inj];
            ris[ni1n+j]=ris[inj];
            ris[ni1n+nj1]=ris[inj];
            ris[nj1n+ni1]=ris[inj];
            ris[jn+ni1]=ris[inj];
            ris[in+nj1]=ris[inj];
        }
    }

    ris[np*n+np]/=somma;
    return ris;
}

float smooth(int x,int y,int np,float smooth1,float smooth2) {
    return (smooth1*exp(-((x-np)*(x-np)+(y-np)*(y-np))/(smooth2)));
}

float* smooth2(float theta,int n) {
    int np=n>>1;
    float* ris=(float*)malloc(n*sizeof(float));
    float smooth1=1/((sqrt(2*M_PI)*theta));
    float smooth2=-2*theta*theta;
    int nc=n>>1;
    int nm1=n-1;
    float t, sum=0;
    int inp;
    for(int i=0; i<nc; i++) {
        inp=i-np;
        t=smooth1*exp((inp*inp)/smooth2);
        ris[i]=t;
        ris[nm1-i]=t;
        sum+=2*t;
    }
    ris[nc]=smooth1*exp(((nc-np)*(nc-np))/smooth2);
    sum+=ris[nc];
    for(int i=0; i<nc; i++) {
        ris[i]/=sum;
        ris[nm1-i]=ris[i];
    }
    ris[nc]/=sum;
    return ris;
}

void derive(float* f, int n, int m, float* fx, float* fy) {
    int i,j;
    int imj,im=0;
    int nm=n*m;
    int nmmm=nm-m;
    int nmpm=nm+m;
    int m2=m<<1;
    int nm2m=nm-m2;
    int nm1=n-1;
    int mm1=m-1;
    int mm2=mm1-1;
    int mp1=m+1;
    int mmm1=-m-1;
    int mmp1=-m+1;
    int imjpmmm1,imjpmmp1,imjpmp1;
    // valori interni
    for(i=1; i<nm1; i++) {
            im+=m;
            for(j=1; j<mm1; j++) {
            imj=im+j;
	    imjpmmm1=imj+mmm1;
	    imjpmmp1=imj+mmp1;
	    imjpmp1=imj+mp1;
            fx[imj]=-f[imjpmmm1]-f[imj-m]-f[imjpmmp1]+f[imj+mm1]+f[imj+m]+f[imjpmp1];
            fy[imj]=-f[imjpmmm1]+f[imjpmmp1]-f[imj-1]+f[imj+1]-f[imj+mm1]+f[imjpmp1];
        }
    }
    // bordi verticali
    int immm,impm,imp2m,imm1,imm2,impmp1;
    im=0;
    for(i=1; i<nm1; i++) {
        im+=m;
        imm1=im-1;
        imm2=imm1-1;
        immm=im-m;
        impm=im+m;
        impmp1=impm+1;
        imp2m=impm+m;
        fx[im]=-f[immm]-f[immm]-f[immm+1]+f[impm]+f[impm]+f[impmp1];
        fx[impm-1]=-f[imm2]-f[imm1]-f[imm1]+f[imp2m-2]+f[imp2m-1]+f[imp2m-1];
        fy[im]=-f[immm]+f[immm+1]-f[im]+f[im+1]-f[impm]+f[impmp1];
        fy[impm-1]=-f[imm2]+f[imm1]-f[impm-2]+f[impm-1]-f[imp2m-2]+f[imp2m-1];
    }
    // bordi orizzontali
    int jm1,jp1;
    for(j=1; j<m-1; j++) {
	jm1=j-1; jp1=j+1;
        fx[j]=-f[jm1]-f[j]-f[jp1]+f[m+jm1]+f[m+j]+f[m+jp1];
        fx[nmmm+j]=-f[nmmm-m+jm1]-f[nmmm-m+j]-f[nmmm-m+jp1]+f[nmmm+jm1]+f[nmmm+j]+f[nmmm+jp1];
        fy[j]=-f[jm1]+f[jp1]-f[jm1]+f[jp1]-f[m+jm1]+f[m+jp1];
        fy[nmmm+j]=-f[nmmm-m+jm1]+f[nmmm-m+jp1]-f[nmmm+jm1]+f[nmmm+jp1]-f[nmmm+jm1]+f[nmmm+jp1];
    }
    //angoli
    fx[0]=-f[0]-f[0]-f[1]+f[m]+f[m]+f[mp1];
    fx[mm1]=-f[m2]-f[mm1]-f[mm1]+f[m2-2]+f[m2-1]+f[m2-1];
    fx[nmmm]=-f[nm2m]-f[nm2m]-f[nm2m+1]+f[nmmm]+f[nmmm]+f[nmmm+1];
    fx[nm-1]=-f[nmpm-2]-f[nmpm-1]-f[nmpm-1]+f[nm-2]+f[nm-1]+f[nm-1];

    fy[0]=-f[0]+f[1]-f[0]+f[1]-f[m]+f[mp1];
    fy[mm1]=-f[mm2]+f[mm1]-f[mm2]+f[mm1]-f[m2-2]+f[m2-1];
    fy[nmmm]=-f[nm2m]+f[nm2m+1]-f[nmmm]+f[nmmm+1]-f[nmmm]+f[nmmm+1];
    fy[nm-1]=-f[nmpm-2]+f[nmpm-1]-f[nm-2]+f[nm-1]-f[nm-2]+f[nm-1];


}


float* filterGauss(float* f, int r, int c, float* w, int n,float* ris,float* ris2) {
    nc=n>>1;
    int i,j,k,icj,ic,nck,icc1,kc,ncckc;
    int c1=c-1;
    int r1=r-1;
    int r1c=r1*c;
    int ncc=nc*c;
    int ncp1=nc+1;
    int ncm1=nc-1;
    float s,s2;
    // FILTRO ORIZZONTALE
    // valori interni
    ic=-c;
    for(i=0; i<r; i++) {
        ic+=c;
	icj=ic+ncm1;
        for(j=nc; j<c-nc; j++) {
            icj++;
            s=0;
	    nck=ncp1;
            for(k=0; k<nc; k++) {
                nck--;
                s+=w[k]*(f[icj-nck]+f[icj+nck]);
            }
            ris[icj]=s+(w[nc]*f[icj]);
        }
    }
    // valori esterni
    ic=-c;
    for(i=0; i<r; i++) {
    	ic+=c;
	icj=ic;
	icc1=ic+c1;
	s=0;
	s2=0;
	// prima e ultima colonna
	nck=ncp1;
	for(k=0;k<nc;k++){
		nck--;
		s+=w[k]*(f[ic]+f[ic+nck]);
		s2+=w[k]*(f[icc1-nck]+f[icc1]);
	}
	ris[ic]=s+(w[nc]*f[ic]);
    	ris[icc1]=s2+(w[nc]*f[icc1]);
	// sezioni ovest-est
    	for(j=1; j<nc; j++) {
		icj++;
		s=0;
		s2=0;
		//sfora ovest-est
		nck=ncp1;
		for(k=0;k<nc-j;k++){
			nck--;
			s+=w[k]*(f[ic]+f[icj+nck]);// sezione ovest
			s2+=w[k]*(f[icc1-j-nck]+f[icc1]);// sezione est
		}
		//valori rimanenti che non sforano
		for(k=nc-j; k<nc; k++) {
			nck--;
			s+=w[k]*(f[icj-nck]+f[icj+nck]);
			s2+=w[k]*(f[icc1-j-nck]+f[icc1-j+nck]);
		}
		ris[icj]=s+(w[nc]*f[icj]);
            	ris[icc1-j]=s2+(w[nc]*f[icc1-j]);
	}
    }


    // FILTRO VERTICALE
    // valori interni
    ic=ncc-c;
    for(i=nc; i<r-nc; i++) {
        ic+=c;
	icj=ic-1;
        for(j=0; j<c; j++) {
            icj++;
            s=0;
	    kc=-c;
            for(k=0; k<nc; k++) {
		kc+=c;
		ncckc=ncc-kc;
                s+=w[k]*(ris[ic-ncckc+j]+ris[ic+ncckc+j]);
            }
            ris2[icj]=s+(w[nc]*ris[icj]);
        }
    }
   // valori esterni

// prima e ultima riga
    	ic=0;
	icj=ic-1;
    	for(j=0; j<c; j++) {
		icj++;
		s=0;
		s2=0;
		//sfora nord-sud
		kc=-c;
		for(k=0;k<nc;k++){
			kc+=c;
			ncckc=ncc-kc;
			s+=w[k]*(ris[j]+ris[ncckc+j]);// sezione nord
			s2+=w[k]*(ris[r1c+j]+ris[r1c-ncckc+j]);// sezione sud
		}
		ris2[icj]=s+(w[nc]*ris[icj]);
            	ris2[r1c+j]=s2+(w[nc]*ris[r1c+j]);
	}





    for(i=1; i<nc; i++) {
    	ic+=c;
	icj=ic-1;
    	for(j=0; j<c; j++) {
		icj++;
		s=0;
		s2=0;
		//sfora nord-sud
		kc=-c;
		for(k=0;k<nc-i;k++){
			kc+=c;
			ncckc=ncc-kc;
			s+=w[k]*(ris[j]+ris[ic+ncckc+j]);// sezione nord
			s2+=w[k]*(ris[r1c+j]+ris[r1c-ic-ncckc+j]);// sezione sud
		}
		//valori rimanenti che non sforano
		for(k=nc-i; k<nc; k++) {
			kc+=c;
			ncckc=ncc-kc;
			s+=w[k]*(ris[ic-ncckc+j]+ris[ic+ncckc+j]);
			s2+=w[k]*(ris[r1c-ic-ncckc+j]+ris[r1c-ic+ncckc+j]);
		}
		ris2[icj]=s+(w[nc]*ris[icj]);
            	ris2[r1c-ic+j]=s2+(w[nc]*ris[r1c-ic+j]);
	}
    }
    
    return ris2;
}


float* filter(float* f, int r, int c, float* filter, int n, int m) {
    float* ris=create(r,c);
    nc=n>>1;
    mc=m>>1;
    r1=r-1;
    c1=c-1;
    int i,j,ic;
    int r1c=r1*c;
    int nccc=nc*c-c;
    
    // valori interni (non sfora mai)
    ic=nccc;
    for(i=nc;i<r-nc;i++){
	ic+=c;
	for(j=nc;j<c-nc;j++){
		ris[ic+j]=filtInner(f,c,filter,n,m,i,j);
	}
    }
    // valori righe nord-sud (sfora solo a nord-sud)
    ic=-c;
    for(i=0; i<nc; i++){
	ic+=c;
        for(j=nc; j<c-nc; j++){
            ris[ic+j]=filtN(f,r,c,filter,n,m,i,j);
	    ris[r1c-ic+j]=filtS(f,r,c,filter,n,m,r1-i,j);
	}
    }
    // valori colonne est-ovest (sfora solo a est-ovest)
    ic=nccc;
    for(i=nc; i<r-nc; i++){
	ic+=c;
        for(j=0; j<nc; j++){
            ris[ic+j]=filtE(f,r,c,filter,n,m,i,j);
	    ris[ic+c1-j]=filtO(f,r,c,filter,n,m,i,c1-j);
	}
    }
    // valori nei "vertici" (sfora in qualunque direzione)
    ic=-c;
    for(i=0; i<nc; i++){
	ic+=c;
        for(j=0; j<nc; j++){
            ris[ic+j]=filtV(f,r,c,filter,n,m,i,j);
	    ris[r1c-ic+j]=filtV(f,r,c,filter,n,m,r1-i,j);
	}
    }
    return ris;
}

float filtN(float* f,int r, int c,float* w, int n, int m, int x, int y) {
    int i,j;
    float s=0;
    int exti,extj,im;
    int xnc=x-nc;
    int ymc=y-mc;
    for(i=0; i<n; i++){
	im=i*m;
	exti=xnc+i;
	for(j=0; j<m; j++) {
		extj=ymc+j;
		if(exti<0) //nord
		    exti=0;
		s+=w[im+j]*f[exti*c+extj];
	    }
    }
    return s;
}

float filtS(float* f,int r, int c,float* w, int n, int m, int x, int y) {
    int i,j;
    float s=0;
    int exti,extj,im;
    int xnc=x-nc;
    int ymc=y-mc;
    for(i=0; i<n; i++){
	im=i*m;
	exti=xnc+i;
	for(j=0; j<m; j++) {
		extj=ymc+j;
		if(exti>r1) //sud
		    exti=r1;
		s+=w[im+j]*f[exti*c+extj];
	    }
    }
    return s;
}

float filtO(float* f,int r, int c,float* w, int n, int m, int x, int y) {
    int i,j;
    float s=0;
    int exti,extj,im;
    int xnc=x-nc;
    int ymc=y-mc;
    for(i=0; i<n; i++){
	im=i*m;
	exti=xnc+i;
	for(j=0; j<m; j++) {
		extj=ymc+j;
		if(extj<0) //ovest
		    extj=0;
		s+=w[im+j]*f[exti*c+extj];
	    }
    }
    return s;
}

float filtE(float* f,int r, int c,float* w, int n, int m, int x, int y) {
    int i,j;
    float s=0;
    int exti,extj,im;
    int xnc=x-nc;
    int ymc=y-mc;
    for(i=0; i<n; i++){
	im=i*m;
	exti=xnc+i;
	for(j=0; j<m; j++) {
		extj=ymc+j;
		if(extj>c1) //est
		    extj=c1;
		s+=w[im+j]*f[exti*c+extj];
	    }
    }
    return s;
}


float filtV(float* f,int r, int c,float* w, int n, int m, int x, int y) {
    int i,j;
    float s=0;
    int exti,extj,im;
    int xnc=x-nc;
    int ymc=y-mc;
    for(i=0; i<n; i++){
	im=i*m;
	exti=xnc+i;
	for(j=0; j<m; j++) {
		extj=ymc+j;
		if(exti<0) //nord
		    exti=0;
		else if(exti>r1) //sud
		    exti=r1;
		if(extj<0) //ovest
		    extj=0;
		else if(extj>c1) //est
		    extj=c1;
		s+=w[im+j]*f[exti*c+extj];
	    }
    }
    return s;
}


float filtInner(float* f, int c,float* w, int n, int m, int x, int y) {
    int i,j,im,ic;
    float s=0;
    int xnccymc=(x-nc)*c+y-mc;
    im=-m;
    ic=-c;
    for(i=0; i<n; i++){
	im+=m;
	ic+=c;
	for(j=0; j<m; j++) {
		s+=w[im+j]*f[xnccymc+ic+j];
	    }
    }
    return s;
}


float* soppressione(float* f,int r, int c, int t) {
    float* ris=create(r,c);
    int i,j,icj,ic;
    for(i=0; i<r; i++){
	ic=i*c;
        for(j=0; j<c; j++) {
            icj=ic+j;
            ris[icj]=soppr(f,r,c,t,i,j);
            if(ris[icj]==1)
                corners++;
        }
    }
    input->n_corners=corners;
    return ris;
}

int soppr(float* f,int r, int c, int t, int x, int y) {
    int exti,extj;
    int i,j;
    int cm1=c-1;
    int rm1=r-1;
    int t2p1=(t<<1)+1;
    int xmt=x-t;
    int ymt=y-t;
    int xcy=x*c+y;
    for(i=0; i<t2p1; i++){
	exti=xmt+i;
	// nord
        if(exti<0)
            continue;
        // sud
        else if(exti>rm1)
            continue;
        for(j=0; j<t2p1; j++) {
            if(i==t && j==t)
                continue;
            extj=ymt+j;
            // ovest
            if(extj<0)
                continue;
            // est
            else if(extj>cm1)
                continue;
            if(f[xcy]<=f[exti*c+extj])
                return 0;
        }
    }
    return 1;
}

void findCorners(float* r,int n, int m, location* c_locations){
	int i,j,index=0;
	for(i=0;i<n;i++){
		for(j=0;j<m;j++){
			if(r[i*m+j]){
				c_locations[index].x=i;
				c_locations[index].y=j;
				index++;
			}
		}
	}

    FILE* fp;
    char fpath[256];
    sprintf(fpath, "%s_corners.txt", input->img_name);
    fp = fopen(fpath, "w");
    for (i = 0; i < input->n_corners; i++)
        fprintf(fp, "%d %d\n", c_locations[i].x, c_locations[i].y);
    fclose(fp);

}

void print(float*m,int r,int c) {
    int i,j;
    printf("\t");
    for(i=0; i<c; i++) {
        printf("%d\t\t",i);
    }
    printf("\n");
    for(i=0; i<r; i++) {
        printf("%d\t",i);
        for(j=0; j<c; j++) {
            printf("%f\t",m[i*c+j]);
        }
        printf("\n");
    }
}


void intorno(float*m,int r, int c,int x,int y,int t) {
    int i,j;
    printf("\t");
    for(i=y-t; i<y+t; i++) {
        printf("%d\t\t",i);
    }
    printf("\n");
    for(i=x-t; i<x+t; i++) {
        printf("%d\t",i);
        for(j=y-t; j<y+t; j++) {
            printf("%f\t",m[i*c+j]);
        }
        printf("\n");
    }
}

/*Fine metodi*/


/*
 *
 *	Le funzioni sono state scritte assumento che le matrici siano memorizzate
 * 	mediante un array (float*), in modo da occupare un unico blocco
 * 	di memoria, ma a scelta del candidato possono essere
 * 	memorizzate mediante array di array (float**).
 *
 * 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
 * 	matrici per righe (row-major order) o per colonne (column major-order).
 *
 * 	L'assunzione corrente è che le matrici siano in row-major order.
 *
 */


void* get_block(int size, int elements) {
    return _mm_malloc(elements*size,16);
}


void free_block(void* p) {
    _mm_free(p);
}


MATRIX alloc_matrix(int rows, int cols) {
    return (MATRIX) get_block(sizeof(float),rows*cols);
}


void dealloc_matrix(MATRIX mat) {
    free_block(mat);
}


/*
 *
 * 	load_input
 * 	===========
 *
 *	Legge da file l'immagine codificata come una matrice di N righe
 * 	e M colonne e la memorizza in un array lineare in row-major order
 *
 * 	Codifica del file:
 * 	primi 4 byte: numero di colonne (M) --> numero intero
 * 	4 byte successivi: numero di righe (N) --> numero intero
 * 	N*M*4 byte successivi: imag data in row-major order --> numeri floating-point a precisione singola
 *
 *****************************************************************************
 *	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
 * 	dell'immagine.
 *****************************************************************************
 *
 */
IMAGE load_img(char* filename, int *n, int *m) {
    FILE* fp;
    int rows, cols, status, i;
    char fpath[256];

    sprintf(fpath, "%s.img", filename);
    fp = fopen(fpath, "rb");

    if (fp == NULL) {
        printf("'%s' : bad image file name!\n", fpath);
        exit(0);
    }

    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);
    IMAGE data = alloc_matrix(rows,cols);
    status = fread(data, sizeof(float), rows*cols, fp);
    fclose(fp);

    *n = rows;
    *m = cols;

    return data;
}


void save_img(char* filename, IMAGE img, int n, int m) {
    FILE* fp;
    char fpath[256];

    sprintf(fpath, "%s_output.img", filename);
    fp = fopen(fpath, "wb");
    fwrite(&m, sizeof(int), 1, fp);
    fwrite(&n, sizeof(int), 1, fp);
    if (img != NULL)
        fwrite(img, sizeof(float), n*m, fp);
    fclose(fp);
}


void save_corners(char* filename, int n_corners, location* corners) {
    FILE* fp;
    int i;
    char fpath[256];

    sprintf(fpath, "%s_corners.txt", filename);
    fp = fopen(fpath, "w");
    for (i = 0; i < n_corners; i++)
        fprintf(fp, "%d %d\n", corners[i].x, corners[i].y);
    fclose(fp);
}


/*
 *	corner
 * 	====
 *
 *	img contiene l'immagine codificato come una matrice di N righe
 * 	ed M colonne memorizzata in un array lineare in row-major order
 *
 *	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
 * 	dell'immagine.
 *
 *
 */
void corner(params* input) {

    if(input->method == CORNER) {

        float* f=input->img;
        int n=input->N;
        int m=input->M;
        int nm=n*m;
	int dim=nm*sizeof(float);
        int i,j;
        float* temp1=(float*)malloc(dim);
        float* r=(float*)malloc(dim);
        // 1) Smoothing immagine con filtro Gauss (sigma)
	int sigma6=input->sigma*6;
        int nf= ceil(sigma6);
        if((nf & 1)==0)
            nf++;
        float* t2=smooth2(input->sigma,nf);
        float* fg=filterGauss(f,n,m,t2,nf,temp1,r);
        // 2) Derivate parziali Fx e Fy
        float* fx=(float*)malloc(dim);
        float* fxx=(float*)malloc(dim);
        float* fxy=(float*)malloc(dim);
        float* fyy=(float*)malloc(dim);
        float* fy=(float*)malloc(dim);
        derive(fg,n,m,fx,fy);
        // 3) Calcolo Fxx, Fxy e Fyy
        prodotto(fx,fx,n,m,fxx);
        prodotto(fx,fy,n,m,fxy);
        prodotto(fy,fy,n,m,fyy);
        // 4) Smoothing Fxx, Fxy e Fyy con filtro Gauss (2*sigma)
        int nf2= ceil(2*sigma6);
        if((nf2 & 1)==0)
            nf2++;
        float* t3=smooth2(input->sigma*2, nf2);
        float* sxx=filterGauss(fxx,n,m,t3,nf2,input->img,fg);
        float* sxy=filterGauss(fxy,n,m,t3,nf2,input->img,fx);
        float* syy=filterGauss(fyy,n,m,t3,nf2,input->img,fy);
        // 5-6) Calcolo matrice R
        float Rmax=-1;
        int index;
        float tr;
	int im;
        for (i=0; i<n; i++) {
	    im=i*m;
            for(j=0; j<m; j++) {
                index=im+j;
                tr=sxx[index]+syy[index];
                r[index]=(sxx[index]*syy[index]-sxy[index]*sxy[index])-0.05*tr*tr;
                if(r[index]>Rmax)
                    Rmax=r[index];
            }
        }
        // 7) Sogliatura: azzerare valori minori di Rmax*theta
        sogliatura(r,n,m,Rmax*input->theta);
        // 8) Soppressione: matrice binaria con valori strettamente maggiori nel loro intorno di raggio t
        float* final=soppressione(r,n,m,input->t);
	location c_locations[corners];
	findCorners(final,n,m,c_locations);
    } else if(input->method == FILTERING) {
        input->filtered_img=filter(input->img,input->N,input->M,input->filt,input->N_filt,input->M_filt);
    }
}

int main(int argc, char** argv) {
    input = malloc(sizeof(params));

    input->method = CORNER;
    input->img_name = "";
    input->img = NULL;
    input->N = 0;
    input->M = 0;
    input->sigma = 1.5f;
    input->theta = 0.2f;
    input->t = 1;
    input->filt_name = "";
    input->filt = NULL;
    input->N_filt = 0;
    input->M_filt = 0;
    input->silent = 0;
    input->display = 0;
    input->filtered_img = NULL;
    input->R = NULL;
    input->n_corners = 0;
    input->corners = NULL;

    int i, j;

    int par = 1;
    while (par < argc) {
        if (par == 1) {
            input->img_name = argv[par];
            par++;
        } else if (strcmp(argv[par],"-s") == 0) {
            input->silent = 1;
            par++;
        } else if (strcmp(argv[par],"-d") == 0) {
            input->display = 1;
            par++;
        } else if (strcmp(argv[par],"-sigma") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing sigma value!\n");
                exit(1);
            }
            input->sigma = atof(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-theta") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing theta value!\n");
                exit(1);
            }
            input->theta = atof(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-t") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing t value!\n");
                exit(1);
            }
            input->t = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-corner") == 0) {
            input->method = CORNER;
            par++;
        } else if (strcmp(argv[par],"-filtering") == 0) {
            input->method = FILTERING;
            par++;
            if (par >= argc) {
                printf("Missing filter file name!\n");
                exit(1);
            }
            input->filt_name = argv[par];
            par++;
        } else
            par++;
    }

    if (!input->silent) {
        printf("Usage: %s <img_file_name> [-d][-s] [[-corner][-sigma <value>][-theta <value>][-t <value>]] [-filtering <filter_file_name>]\n", argv[0]);
        printf("\nParameters:\n");
        printf("\t-corner: perform corner detection (default)\n");
        printf("\t-d : display input and output\n");
        printf("\t-s : silent\n");
        printf("\t-sigma <value> : st.dev. of Gaussian filter (default 1.5)\n");
        printf("\t-theta <value> : thresholding parameter (default 0.2)\n");
        printf("\t-t <value> : non-maxima suppression parameter (default 1)\n");
        printf("\t-filtering <filter_file_name> : perform image filtering\n");
        printf("\n");
    }

    if (strlen(input->img_name) == 0) {
        printf("Missing image file name!\n");
        exit(1);
    }

    input->img = load_img(input->img_name, &input->N, &input->M);

    if (input->method == FILTERING)
        input->filt = load_img(input->filt_name, &input->N_filt, &input->M_filt);

    if (!input->silent) {
        printf("Image file name: '%s'\n", input->img_name);
        printf("Image rows: %d\n", input->N);
        printf("Image columns: %d\n", input->M);
        if (input->method == CORNER) {
            printf("Parameter sigma: %f\n", input->sigma);
            printf("Parameter theta: %f\n", input->theta);
            printf("Parameter t: %d\n", input->t);
        } else if (input->method == FILTERING) {
            printf("Filter file name: '%s'\n", input->filt_name);
            printf("Filter rows: %d\n", input->N_filt);
            printf("Filter columns: %d\n", input->M_filt);
        }
    }

    clock_t t = clock();
    corner(input);
    t = clock() - t;

    if (!input->silent)
        printf("\nExecution time = %.3f seconds\n", ((float)t)/CLOCKS_PER_SEC);
    else
        printf("%.3f\n", ((float)t)/CLOCKS_PER_SEC);

    if (!input->silent && input->display) {
        printf("\nCorner locations:\n");
        for (i = 0; i < input->n_corners; i++) {
            printf("(%d,%d)\n", input->corners[i].x, input->corners[i].y);
        }
    }

    if (input->method == CORNER) {
        //save_corners(input->img_name, input->n_corners, input->corners);
        save_img(input->img_name, input->R, input->N, input->M);
    } else if (input->method == FILTERING) {
        save_img(input->img_name, input->filtered_img, input->N, input->M);
    }

    return 0;
}
