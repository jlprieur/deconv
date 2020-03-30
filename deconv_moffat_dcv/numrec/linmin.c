#define NRANSI
#include "nrutil.h"
#define TOL 2.0e-4

int ncom;
FLOAT *pcom,*xicom,(*nrfunc)(FLOAT []);

void linmin(FLOAT p[], FLOAT xi[], int n, FLOAT *fret, FLOAT (*func)(FLOAT []))
{
	FLOAT brent(FLOAT ax, FLOAT bx, FLOAT cx,
		FLOAT (*f)(FLOAT), FLOAT tol, FLOAT *xmin);
	FLOAT f1dim(FLOAT x);
	void mnbrak(FLOAT *ax, FLOAT *bx, FLOAT *cx, FLOAT *fa, FLOAT *fb,
		FLOAT *fc, FLOAT (*func)(FLOAT));
	int j;
	FLOAT xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *$!$!6)$. */
