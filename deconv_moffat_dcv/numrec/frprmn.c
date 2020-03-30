#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn(FLOAT p[], int n, FLOAT ftol, int *iter, FLOAT *fret,
	FLOAT (*func)(FLOAT []), void (*dfunc)(FLOAT [], FLOAT []))
{
/* JLP2002:
	void linmin(FLOAT p[], FLOAT xi[], int n, FLOAT *fret,
		FLOAT (*func)(FLOAT []));
* I use dlinmin, with gradient information instead: */
	void dlinmin(FLOAT p[], FLOAT xi[], int n, FLOAT *fret,
		FLOAT (*func)(FLOAT []), void (*dfunc)(FLOAT [], FLOAT []));
	int j,its;
	FLOAT gg,gam,fp,dgg;
	FLOAT *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
/* JLP2002: I use dlinmin, with gradient information instead: 
		linmin(p,xi,n,fret,func);
*/
printf("calling dlinmin with n=%d\n",n);
		dlinmin(p,xi,n,fret,func,dfunc);
printf("exit from dlinmin\n");
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *$!$!6)$. */
