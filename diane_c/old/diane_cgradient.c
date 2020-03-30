/***************************************************************/
/*                                                             */
/* PROGRAM: diane_cgradient                                    */
/*                                                             */
/* PURPOSE: Deconvolution with multiresolution constraint      */
/*          using a mean square minimization involving         */
/*          conjugate gradient method                          */
/*                                                             */
/* INPUT:  argv[1] = generic_name                              */
/*         argv[2] = threshold for residual                    */
/*         argv[3] = error1                                    */
/*         argv[4] = error2                                    */
/*         argv[5] = error3                                    */
/*                                                             */
/* NB: *.PHIT  first approximation of the object       (input) */
/*     *.PCR   real part of FT of PHIT                 (input) */
/*     *.PCI   imaginary part of FT of PHIT            (input) */
/*     *.V     multiresolution support                 (input) */
/*     *.G     interpolation function                  (input) */
/*     *.PHIN  final solution                         (output) */
/*                                                             */
/* VERSION: March 1993                                         */
/*                                                             */
/* AUTHOR: JLP, SR, JV                                         */
/*         translated to C by Karim BOUYOUCEF                  */
/*                                                             */
/***************************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

extern   void      fourn();
extern   float **  alloc_square_float();
extern   void      free_square_float();
extern   int       decomposition_mallat_ncs();
extern   int       recomposition_mallat_ncs();

/***************************************************************/
main(argc,argv)
int argc;
char *argv[];
/***************************************************************/
{
/*
      DECLARATIONS
*/
auto     int      dim, size[4];
auto     char     file[100], type[16], mode[8], nature[16], comments[1024];
auto     int      nbr_lin, nbr_col, err, d_order;
auto     float    eps, eps1, eps2, eps3;
register int      i, j, k, I;
auto     int      isign, nn[2], ndim;
auto     float    major, xabs[111], xrn[111], xdn[111], g2, xnorrn, psdz;
auto     float    **phit, **phin, **pcr, **pci, **g, **v, **phitr, **phiti;
auto     float    *phitc, **rn, **zn, **dnr, **dni, **dntemp;
auto     float    xnorfin, xnorrnp, x, y, omegan, isig, tmp, xx, xnun, pas, rac[45];
auto     float    test, omega_n, gamma_n, deltan, xnfinc, ron, tete, tetap, xmum;
auto     int      m, itest, num, nomb, nrac, mult[45], icn[110];
/**********************************************************************************/

/*      VERIFICATION DE LA LIGNE DE COMMANDE
*/
if (argc != 5)
  {
  (void) printf("\nUnexpected number of arguments\n");
  (void) printf("\nUSAGE:\n");
  (void) printf("diane_gradient  generic_name  d_order  eps\n\n");
  (void) printf("\t eps1 = error1\n");
  (void) printf("\t eps2 = error2\n");
  (void) printf("\t eps3 = error3\n");
  (void) printf("\t \n\n");
  exit(-1);
  }

/*       LECTURE DES DIFFERENTS FICHIERS ET DES DATAs
*/
dim = 2;
(void) strcpy(type,"float");
(void) strcpy(nature,"real");

sprintf(file,"%s.PHIT",argv[1]);
phit = (float **) read_data(file,dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];
/*
split_image(phit,nbr_lin,nbr_col);
*/

sprintf(file,"%s.PCR",argv[1]);
pcr = (float **) read_data(file,dim,size,type,mode,nature,comments);
if ((nbr_lin != size[0]) || (nbr_col != size[1]))
 {
 fprintf(stderr,"\a\nFATAL ERROR : [%s.PCR & %s.PHIT] have different sizes\n",argv[1],argv[1]);
 exit(-1);
 }
split_image(pcr,nbr_lin,nbr_col);

sprintf(file,"%s.PCI",argv[1]);
pci = (float **) read_data(file,dim,size,type,mode,nature,comments);
if ((nbr_lin != size[0]) || (nbr_col != size[1]))
 {
 fprintf(stderr,"\a\nFATAL ERROR : [%s.PCI & %s.PHIT] have different sizes\n",argv[1],argv[1]);
 exit(-1);
 }
split_image(pci,nbr_lin,nbr_col);

sprintf(file,"%s.G",argv[1]);
g = (float **) read_data(file,dim,size,type,mode,nature,comments);
if ((nbr_lin != size[0]) || (nbr_col != size[1]))
  {
  fprintf(stderr,"\a\nFATAL ERROR : [%s.G & %s.PHIT] have different sizes\n",argv[1],argv[1]);
  exit(-1);
  }
split_image(g,nbr_lin,nbr_col);

sprintf(file,"%s.V",argv[1]);
v = (float **) read_data(file,dim,size,type,mode,nature,comments);
if ((nbr_lin != size[0]) || (nbr_col != size[1]))
  {
  fprintf(stderr,"\a\nFATAL ERROR : [%s.V & %s.PHIT] have different sizes\n",argv[1],argv[1]);
  exit(-1);
  }

err = sscanf(comments,"J=%d",&d_order);
if ((err != 1) || (d_order < 0) || (d_order > 10))
   d_order = 0;

err = sscanf(argv[2],"%f",&eps1);
if ((err != 1) || (eps1 < 0.0))
  fprintf(stderr,"\nWARNING : incorrect error parameter [%.10f]\n",eps1);

err = sscanf(argv[3],"%f",&eps2);
if ((err != 1) || (eps2 < 0.0))
  fprintf(stderr,"\nWARNING : incorrect error parameter [%.10f]\n",eps2);

err = sscanf(argv[4],"%f",&eps3);
if ((err != 1) || (eps3 < 0.0))
  fprintf(stderr,"\nWARNING : incorrect error parameter [%.10f]\n",eps3);

/*
      MEMORY ALLOCATION
*/
phin = (float**)alloc_square_float(nbr_lin, nbr_col);
phitr = (float**)alloc_square_float(nbr_lin, nbr_col);
phiti = (float**)alloc_square_float(nbr_lin, nbr_col);
phitc = (float *)malloc(2*nbr_lin*nbr_col*sizeof(float));
rn = (float**)alloc_square_float(nbr_lin, nbr_col);
zn = (float**)alloc_square_float(nbr_lin, nbr_col);
dnr = (float**)alloc_square_float(nbr_lin, nbr_col);
dni = (float**)alloc_square_float(nbr_lin, nbr_col);
dntemp = (float**)alloc_square_float(nbr_lin, nbr_col);

fprintf(stderr,"\n*****************************");
fprintf(stderr,"\nPROGRAM :      diane_cgradient");
fprintf(stderr,"\n*****************************");
fprintf(stderr,"\ndecomposition order of the support constraint = %d",d_order);
fprintf(stderr,"\nthreshold for residual error = %.10f",eps);
fprintf(stderr,"\nerror eps1 = %.10f",eps1);
fprintf(stderr,"\nerror eps2 = %.10f",eps2);
fprintf(stderr,"\nerror eps3 = %.10f",eps3);

/*************************************************************
* 
*      CENTRAL RECONSTRUCTION
*
*     compute eigenvalues of A*A
*     conjugate gradients method
*     search polynom whose roots converge towards eigenvalues
*
**************************************************************/
/*  upper limit for 1/mu  */
major = 10.0;

/*********************************************************
   initialisation of polynoms
   xabs = absices , xrn = polynom rn  , xdn = polynom dn
   icn = polynom cn of sign changes
*********************************************************/
for (k=1; k<=101; k++)
  {
  xabs[k] = (float)(k-1.0)*0.01;
  xrn[k] = 1.0;
  xdn[k] = 1.0;
  icn[k] = 0;
  }

/*   phi0 = V * phit    */
/*
err = (int) decomposition_mallat_ncs(phit,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the decomposition\n\n");
  exit(-1);
  }
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  phin[i][j] = v[i][j] * phit[i][j];
err = (int) recomposition_mallat_ncs(phin,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  exit(-1);
  }
err = (int) recomposition_mallat_ncs(phit,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  exit(-1);
  }
*/
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  phin[i][j] = phit[i][j] * v[i][j];

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phitr[i][j] = phin[i][j];
  phiti[i][j] = 0.0;
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phitc[2*(j+i*nbr_col)] = phitr[i][j];
  phitc[2*(j+i*nbr_col)+1] = phiti[i][j];
  }
ndim = 2;
isign = -1;
nn[0] = nbr_lin;
nn[1] = nbr_col;
(void) fourn(phitc-1,nn-1,ndim,isign);

for (k=0; k<nbr_lin*nbr_col; k++)
  {
  i = k/nbr_col;
  j = k-i*nbr_col;
  phitr[i][j] = phitc[2*k];
  phiti[i][j] = phitc[2*k+1];
  }

/*  r0 = v U* g2   (phitc-phi0 chap)  */
/*  compute d0 and ||rn||**2 for n=0  */
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  g2 = g[i][j] * g[i][j];
  phitr[i][j] = (pcr[i][j] - phitr[i][j]) * g2;
  phiti[i][j] = (pci[i][j] - phiti[i][j]) * g2;
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phitc[2*(j+i*nbr_col)] = phitr[i][j];
  phitc[2*(j+i*nbr_col)+1] = phiti[i][j];
  }
isign = 1;
(void) fourn(phitc-1,nn-1,ndim,isign);

for (k=0; k<nbr_lin*nbr_col; k++)
  {
  i = k/nbr_col;
  j = k-i*nbr_col;
  phitr[i][j] = phitc[2*k]/((float)(nbr_lin*nbr_col));
  phiti[i][j] = phitc[2*k+1]/((float)(nbr_lin*nbr_col));
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  phitr[i][j] = (float)sqrt(phitr[i][j]*phitr[i][j]+phiti[i][j]*phiti[i][j]);

xnorrn = 0.0;
/*
err = (int) decomposition_mallat_ncs(phitr,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the decomposition\n\n");
  exit(-1);
  }
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phiti[i][j] = 0.0;
  rn[i][j] = phitr[i][j] * v[i][j];
  }
err = (int) recomposition_mallat_ncs(rn,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  exit(-1);
  }
err = (int) recomposition_mallat_ncs(phitr,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  exit(-1);
  }
*/
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phiti[i][j] = 0.0;
  rn[i][j] = phitr[i][j] * v[i][j];
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  dnr[i][j] = rn[i][j];
  xnorrn += (dnr[i][j]*dnr[i][j]);
  }

/*****************************
        ITERATION
******************************/
fprintf(stderr,"\n   N    Test\n");
I = 0;
do
{
I++;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  dntemp[i][j] = dnr[i][j];
  dni[i][j] = 0.0;
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phitc[2*(j+i*nbr_col)] = dnr[i][j];
  phitc[2*(j+i*nbr_col)+1] = dni[i][j];
  }
isign = -1;
(void) fourn(phitc-1,nn-1,ndim,isign);

for (k=0; k<nbr_lin*nbr_col; k++)
  {
  i = k/nbr_col;
  j = k-i*nbr_col;
  dnr[i][j] = phitc[2*k];
  dni[i][j] = phitc[2*k+1];
  }

/*  compute zn and omega n  */
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  g2 = g[i][j] * g[i][j];
  dnr[i][j] *= g2;
  dni[i][j] *= g2;
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phitc[2*(j+i*nbr_col)] = dnr[i][j];
  phitc[2*(j+i*nbr_col)+1] = dni[i][j];
  }
isign = 1;
(void) fourn(phitc-1,nn-1,ndim,isign);

for (k=0; k<nbr_lin*nbr_col; k++)
  {
  i = k/nbr_col;
  j = k-i*nbr_col;
  dnr[i][j] = phitc[2*k]/((float)(nbr_lin*nbr_col));
  dni[i][j] = phitc[2*k+1]/((float)(nbr_lin*nbr_col));
  }


/*
err = (int) decomposition_mallat_ncs(dnr,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the decomposition\n\n");
  exit(-1);
  }
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  zn[i][j] = dnr[i][j] * v[i][j];
err = (int) recomposition_mallat_ncs(zn,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  exit(-1);
  }
err = (int) recomposition_mallat_ncs(dnr,nbr_lin,nbr_col,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  exit(-1);
  }
*/
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  zn[i][j] = dnr[i][j] * v[i][j];

psdz = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  psdz = psdz + (dntemp[i][j]*zn[i][j]);
if (psdz != 0.0)
   omegan = xnorrn / psdz;
else
   {
   fprintf(stderr,"FATAL ERROR: emergency exit psdz=0\n\n");
   exit(-1);
   }

/*
    compute phi n+1 , ||phi n+1||**2
    r n+1 ,  ||r n+1||**2         
*/
xnorfin = 0.0;
xnorrnp = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  phin[i][j] = phin[i][j] + (omegan * dntemp[i][j]);
  x = phin[i][j] * phin[i][j];
  xnorfin = xnorfin + x;
  rn[i][j] = rn[i][j] - (omegan * zn[i][j]);
  y = rn[i][j] * rn[i][j];
  xnorrnp = xnorrnp + y;
  }

/*
    compute polynoms rn and cn
    number of degrees of freedom : icn[94]
    isig correspond to polynom Sn(x)
*/
itest = 0;
for (k=1; k<=101; k++)
  {
  if (xrn[k] >= 0.0)        isig = 1.;
  else if (xrn[k] < 0.0)   isig = -1.;
  xx = omegan * xabs[k] * xdn[k];
  xrn[k] = xrn[k] - xx;
  if (xrn[k] >= 0.0)       tmp = 1.;
  else if (xrn[k] < 0.0)  tmp = -1.;
  if (k <= 95)
    {
    if (isig != tmp)
      {
      itest = 1;
      icn[k] = icn[k] + 1;
      }
    }
  }

test = xnorrnp / xnorfin;
test = major * (float)sqrt(test);
fprintf(stderr,"  %d       %.10f\n",I,test);
  
/*  compute nu n  */
xnun = xnorrnp / xnorrn;

/*  compute d n+1  */
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  dnr[i][j] = rn[i][j] + xnun * dntemp[i][j];

/*  compute polynom d n+1  */
for (k=1; k<=101; k++)
  xdn[k] = xrn[k] + xnun * xdn[k];

xnorrn = xnorrnp;
}
while((test > 0.000000001) && (I < 30) && (itest != 0));
fprintf(stderr,"\ntest=%.10f I=%d itest=%d\n",test,I,itest);

/*****************************
    END OF ITERATION
*****************************/

fprintf(stderr,"\nNumber of degrees of freedom : %d\n",icn[95]);
/*******************************
  Compute eigenvalues of A*A
  rac is the root of rn
  num is the order number
  mult is the multiplicity
  nomb is the number of roots (taking into account multiplicity)
********************************/
pas = 1./(2.*100.);
m = 94;
num = 0;
nomb = 0;
fprintf(stderr,"\nN     Eigenvalues   Multiplicity\n");
for (k=1; k<=m; k++)
  {
  nrac = icn[k+1] - icn[k];
  if (nrac != 0)
    {
    num = num + 1;
    rac[num] = xabs[k+1] - pas;
    mult[num] = nrac;
    nomb = nomb + nrac;
    fprintf(stderr,"%d    %.10f      %d\n",num,rac[num],mult[num]);
    }
  }

/********************************
  Compute upper limit for error
********************************/
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  phitr[i][j] = g[i][j] * (phitr[i][j] - phin[i][j]);

/*  compute L2 norm   */
deltan = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  deltan += (phitr[i][j]*phitr[i][j] + phiti[i][j]*phiti[i][j]);
deltan = (float)sqrt(deltan/((float)(nbr_lin*nbr_col)));

xnfinc = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  xnfinc += (phin[i][j]*phin[i][j]);
xnfinc = (float)sqrt(xnfinc);
ron = (float)sqrt(fabs(eps1-deltan*deltan));

/*    upper limit for error  */
tete = (float)sqrt(rac[1] * xnfinc);
tete = ron / tete;

/**** compute 1/mu mean and corresponding errors  *****/
xmum = 0.0;
for (k=1; k<=num; k++)
  xmum = xmum + ((1.0/rac[k])*(float)mult[k]);
xmum = xmum / (float)nomb;
xmum = (xmum+1.)/2.;
xmum = (float)sqrt(xmum);
tetap = (xmum*ron)/xnfinc;

fprintf(stderr,"\nNorm of the solution ||Phi n|| = %f",xnfinc);
fprintf(stderr,"\nWeighted mean of the inverses of the eigenvalues (1/mu) = %f",xmum);
fprintf(stderr,"\nError of image reconstruction (sqrt(eps1)) = %f",sqrt(eps1));
fprintf(stderr,"\nError due to the choice of Hr (epsilon_o) = %f",sqrt(eps3));
fprintf(stderr,"\nError related to the SNR ratio (epsilon_i) = %f",sqrt(eps2));
fprintf(stderr,"\nTerm to correct the error (delta_n) = %f",deltan);
fprintf(stderr,"\nRelative error for object reconstruction = %f",tetap);
fprintf(stderr,"\n\n");

/*
     STORE FINAL SOLUTION
*/
dim = 2;
size[0] = nbr_lin;
size[1] = nbr_col;
strcpy(type,"float");
strcpy(mode,"bin");
strcpy(nature,"real");
sprintf(file,"%s.PHIN",argv[1]);
sprintf(comments,"diane");
write_data(file,phin,dim,size,type,mode,nature,comments);

/*
     FREE MEMORY
*/
free_square_float(phit,nbr_lin);
free_square_float(phin,nbr_lin);
free_square_float(pcr,nbr_lin);
free_square_float(pci,nbr_lin);
free_square_float(v,nbr_lin);
free_square_float(g,nbr_lin);
free_square_float(phitr,nbr_lin);
free_square_float(phiti,nbr_lin);
free_square_float(rn,nbr_lin);
free_square_float(zn,nbr_lin);
free_square_float(dnr,nbr_lin);
free_square_float(dni,nbr_lin);
free_square_float(dntemp,nbr_lin);
}
