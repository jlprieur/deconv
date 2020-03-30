/***************************************************************
*                                   
* PROGRAM: dc_cgradient (previously: diane_cgradient)
*                                 
* PURPOSE: Deconvolution with multiresolution constraint 
*          using a mean square minimization involving   
*          conjugate gradient method                   
*                                                     
* INPUT:  argv[1] = generic_name                     
*         argv[2] = threshold for residual          
*         argv[3] = error1,error2,error3                         
*                                               
* INPUT:
*     *.PHIT  first approximation of the object 
*     *.PCR   real part of FT of PHIT          
*     *.PCI   imaginary part of FT of PHIT    
*     *.V     multiresolution support        
*     *.G     interpolation function        
* OUTPUT:
*     *.PHIN  final solution
*                                          
* From Karim's version of March 1993                     
*
* JLP: change to fftw Fourier transform
*                                        
* AUTHOR: JLP, SR, JV                   
*         translated to C by Karim BOUYOUCEF  
* 
* JLP
* Version 10/03/2008
***************************************************************/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <jlp_ftoc.h>

/***************************************************************/
int main(int argc, char *argv[])
{
/*
      DECLARATIONS
*/
char     filename[61], generic_name[61], comments[81];
INT_PNTR pntr_ima;
INT4     nx, ny, nx1, ny1;
float    eps, eps1, eps2, eps3;
register int  i, k, iter;
float    major, xabs[111], xrn[111], xdn[111], g2, xnorrn, psdz;
float    *phit, *phin, *pcr, *pci, *g, *v, *phitr, *phiti;
float    *phitc, *rn, *zn, *dnr, *dni, *dntemp;
float    xnorfin, xnorrnp, x, y, omegan, isig, tmp, xx, xnun, pas, rac[45];
float    test, deltan, xnfinc, ron, tete, tetap, xmum;
int      m, itest, num, nomb, nrac, mult[45], icn[110];
/***************************************************************************/

/*  Syntax checking: */ 
/* Carefull: 7 parameters always, using JLP "runs" */
if(argc == 7 && *argv[2]) argc = 3;
if(argc != 3) {
   printf("\nUnexpected number of arguments\n");
   printf("\nUSAGE:\n");
   printf("diane_gradient  generic_name  eps1,eps2,eps3\n\n");
   printf("\t eps1 = error1\n");
   printf("\t eps2 = error2\n");
   printf("\t eps3 = error3\n");
   printf("\t \n\n");
  return(-1);
  }

/*   Input of the parameters and files 
*/
strcpy(generic_name,argv[1]);

if(sscanf(argv[2],"%f,%f,%f",&eps1,&eps2,&eps3) != 3)
  {
  fprintf(stderr,"\nWARNING : incorrect error parameters [%s]\n", argv[2]);
  return(-1);
  }
if((eps1 < 0.0) || (eps2 < 0.0) || (eps3 < 0.0)) {
  fprintf(stderr,"\nWARNING : incorrect error parameters [%f, %f, %f]\n",
          eps1, eps2, eps3);
  return(-1);
  }

JLP_INQUIFMT();

sprintf(filename,"%sPHIT",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, filename, comments);
phit = (float *) pntr_ima; 
/*
RECENT_FFT(phit,phit,&nx,&ny,&nx);
*/

sprintf(filename,"%sPCR",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
pcr = (float *) pntr_ima; 
if ((nx1 != nx) || (ny1 != ny))
 {
 fprintf(stderr,"\a\nFATAL ERROR : [%s.PCR & %s.PHIT] have different sizes\n",
  generic_name, generic_name);
 return(-1);
 }
RECENT_FFT(pcr,pcr,&nx,&ny,&nx);

sprintf(filename,"%sPCI",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
pci = (float *) pntr_ima; 
if ((nx1 != nx) || (ny1 != ny))
 {
 fprintf(stderr,"\a\nFATAL ERROR : [%s.PCI & %s.PHIT] have different sizes\n",
  generic_name, generic_name);
 return(-1);
 }
RECENT_FFT(pci,pci,&nx,&ny,&nx);

sprintf(filename,"%sG",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
g = (float *) pntr_ima; 
if ((nx1 != nx) || (ny1 != ny))
 {
 fprintf(stderr,"\a\nFATAL ERROR : [%s.G & %s.PHIT] have different sizes\n",
  generic_name, generic_name);
 return(-1);
 }
RECENT_FFT(g,g,&nx,&ny,&nx);

sprintf(filename,"%sV",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
v = (float *) pntr_ima; 
if ((nx1 != nx) || (ny1 != ny))
 {
 fprintf(stderr,"\a\nFATAL ERROR : [%s.V & %s.PHIT] have different sizes\n",
  generic_name, generic_name);
 return(-1);
 }

/*
      MEMORY ALLOCATION
*/
phin = (float *) malloc(nx * ny * sizeof(float));
phitr = (float *) malloc(nx * ny * sizeof(float));
phiti = (float *) malloc(nx * ny * sizeof(float));
phitc = (float *) malloc(2*ny*nx*sizeof(float));
rn = (float *) malloc(nx * ny * sizeof(float));
zn = (float *) malloc(nx * ny * sizeof(float));
dnr = (float *) malloc(nx * ny * sizeof(float));
dni = (float *) malloc(nx * ny * sizeof(float));
dntemp = (float *) malloc(nx * ny * sizeof(float));

printf("\n*****************************");
printf("\nPROGRAM :      diane_cgradient");
printf("\n*****************************");
printf("\nthreshold for residual error = %.10f",eps);
printf("\nerror eps1 = %.10f",eps1);
printf("\nerror eps2 = %.10f",eps2);
printf("\nerror eps3 = %.10f",eps3);

/*************************************************************
*  CENTRAL RECONSTRUCTION
*
*  Compute eigenvalues of A*A
*  with conjugate gradients method
*  and search for polynom whose roots converge towards eigenvalues
*
**************************************************************/
/*  upper limit for 1/mu  */
major = 10.0;

/*********************************************************
  Initialisation of polynoms
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

for (i = 0; i < nx * ny; i++) phin[i] = v[i] * phit[i];

for (i = 0; i < nx * ny; i++) {
  phitr[i] = phin[i];
  phiti[i] = 0.0;
  }

fftw_float(phitr,phiti,nx,ny,-1);

/*  r0 = v U* g2   (phitc-phi0 chap)  */
/*  compute d0 and ||rn||**2 for n=0  */
for (i = 0; i < nx * ny; i++) 
  {
  g2 = g[i] * g[i];
  phitr[i] = (pcr[i] - phitr[i]) * g2;
  phiti[i] = (pci[i] - phiti[i]) * g2;
  }

fftw_float(phitr,phiti,nx,ny,1);

for (i = 0; i < nx * ny; i++) 
  phitr[i] = (float)sqrt(phitr[i]*phitr[i]+phiti[i]*phiti[i]);

for (i = 0; i < nx * ny; i++) 
  {
  phiti[i] = 0.0;
  rn[i] = phitr[i] * v[i];
  }

xnorrn = 0.0;
for (i = 0; i < nx * ny; i++) 
  {
  dnr[i] = rn[i];
  xnorrn += (dnr[i]*dnr[i]);
  }

/*****************************
        ITERATION
******************************/
printf("\n   N    Test\n");
iter = 0;
/* Beginning of do ... while */
do {
iter++;
for (i = 0; i < nx * ny; i++) {
  dntemp[i] = dnr[i];
  dni[i] = 0.0;
  }

fftw_float(dnr,dni,nx,ny,-1);

/* Compute zn and omega n  */
for (i = 0; i < nx * ny; i++) {
  g2 = g[i] * g[i];
  dnr[i] *= g2;
  dni[i] *= g2;
  }

fftw_float(dnr,dni,nx,ny,1);

for (i = 0; i < nx * ny; i++) zn[i] = dnr[i] * v[i];

psdz = 0.0;
for (i = 0; i < nx * ny; i++) psdz += dntemp[i]*zn[i];

if (psdz != 0.0) {
  omegan = xnorrn / psdz;
} else {
  fprintf(stderr,"FATAL ERROR: emergency return psdz=0\n\n");
  return(-1);
}

/* Compute phi n+1 , ||phi n+1||**2
    r n+1 ,  ||r n+1||**2         
*/
xnorfin = 0.0;
xnorrnp = 0.0;
for (i = 0; i < nx * ny; i++) 
  {
  phin[i] = phin[i] + (omegan * dntemp[i]);
  x = phin[i] * phin[i];
  xnorfin = xnorfin + x;
  rn[i] = rn[i] - (omegan * zn[i]);
  y = rn[i] * rn[i];
  xnorrnp = xnorrnp + y;
  }

/* Compute polynoms rn and cn
   Number of degrees of freedom : icn[94]
   isig corresponds to polynom Sn(x)
*/
itest = 0;
for (k=1; k<=101; k++) {
  isig =  (xrn[k] >= 0.0) ? 1. : -1.;
  xx = omegan * xabs[k] * xdn[k];
  xrn[k] -= xx;
  tmp = (xrn[k] >= 0.0) ? 1. : -1;  
  if (k <= 95) {
    if (isig != tmp) {
      itest = 1;
      icn[k] += 1;
      }
  }
}

test = xnorrnp / xnorfin;
test = major * (float)sqrt(test);
printf("  %d       %.10f\n",iter,test);
  
/*  compute nu n  */
xnun = xnorrnp / xnorrn;

/*  compute d n+1  */
for (i = 0; i < nx * ny; i++) 
  dnr[i] = rn[i] + xnun * dntemp[i];

/*  compute polynom d n+1  */
for (k=1; k<=101; k++)
  xdn[k] = xrn[k] + xnun * xdn[k];

xnorrn = xnorrnp;
} while((test > 0.000000001) && (iter < 30) && (itest != 0));
/* End of do while loop */

printf("\n End of main loop: test=%.10f iter=%d itest=%d\n",test,iter,itest);
/*****************************
    END OF ITERATIONS
*****************************/

printf("\nNumber of degrees of freedom : %d\n",icn[95]);
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
printf("\nN     Eigenvalues   Multiplicity\n");
for (k=1; k<=m; k++)
  {
  nrac = icn[k+1] - icn[k];
  if (nrac != 0)
    {
    num = num + 1;
    rac[num] = xabs[k+1] - pas;
    mult[num] = nrac;
    nomb = nomb + nrac;
    printf("%d    %.10f      %d\n",num,rac[num],mult[num]);
    }
  }

/********************************
  Compute upper limit for error
********************************/
for (i = 0; i < nx * ny; i++) 
  phitr[i] = g[i] * (phitr[i] - phin[i]);

/*  compute L2 norm   */
deltan = 0.0;
for (i = 0; i < nx * ny; i++) 
  deltan += (phitr[i]*phitr[i] + phiti[i]*phiti[i]);
deltan = (float)sqrt(deltan/((float)(ny*nx)));

xnfinc = 0.0;
for (i = 0; i < nx * ny; i++) 
  xnfinc += (phin[i]*phin[i]);
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

printf("Norm of the solution ||Phi n|| = %f\n",xnfinc);
printf("Weighted mean of the inverses of the eigenvalues (1/mu) = %f\n",xmum);
printf("Error of image reconstruction: sqrt(eps1)=%f\n",sqrt(eps1));
printf("Error due to the choice of Hr: epsilon_o=%f\n",sqrt(eps3));
printf("Error related to the SNR ratio: epsilon_i=%f\n",sqrt(eps2));
printf("Term to correct the error: delta_n=%f\n",deltan);
printf("Relative error for object reconstruction: %f\n",tetap);

/*
     STORE FINAL SOLUTION
*/
sprintf(filename,"%sPHIN",generic_name);
sprintf(comments,"diane");
JLP_WRITEIMAG(phin, &nx, &ny, &nx, filename, comments);

/*
     FREE MEMORY
*/
free(phit);
free(phin);
free(pcr);
free(pci);
free(v);
free(g);
free(phitr);
free(phiti);
free(rn);
free(zn);
free(dnr);
free(dni);
free(dntemp);

return(0);
}
