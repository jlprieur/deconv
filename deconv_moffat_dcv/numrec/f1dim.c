#define NRANSI
#include "nrutil.h"

extern int ncom;
extern FLOAT *pcom,*xicom,(*nrfunc)(FLOAT []);

FLOAT f1dim(FLOAT x)
{
	int j;
	FLOAT f,*xt;

printf("f1dim: ncom=%d\n",ncom);

	xt=vector(1,ncom);
printf("f1dim: OK\n");
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_vector(xt,1,ncom);
	return f;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software *$!$!6)$. */
