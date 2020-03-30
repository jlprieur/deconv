/***************************************************************
*                                   
* PROGRAM: dc_cgradient_mallat (previously: diane_cgradient)
*          same as dc_cgradient but with Mallat decomposition of the support
*          => another input parameter (= order)
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

static int recomposition_mallat_ncs(float *mat, int nx, int ny, int order);
static int decomposition_mallat_ncs(float *mat, int nx, int ny, int order);

/***************************************************************/
int main(int argc, char *argv[])
{
/*
      DECLARATIONS
*/
char     filename[61], generic_name[61], comments[81];
INT_PNTR pntr_ima;
INT4     nx, ny, nx1, ny1;
int      err, d_order;
float    eps, eps1, eps2, eps3;
register int  i, k, ii;
float    major, xabs[111], xrn[111], xdn[111], g2, xnorrn, psdz;
float    *phit, *phin, *pcr, *pci, *g, *v, *phitr, *phiti;
float    *phitc, *rn, *zn, *dnr, *dni, *dntemp;
float    xnorfin, xnorrnp, x, y, omegan, isig, tmp, xx, xnun, pas, rac[45];
float    test, deltan, xnfinc, ron, tete, tetap, xmum;
int      m, itest, num, nomb, nrac, mult[45], icn[110];
/***************************************************************************/

/*  Syntax checking: */ 
if (argc != 4 && argc != 7)
  {
   printf("\nUnexpected number of arguments\n");
   printf("\nUSAGE:\n");
   printf("diane_gradient  generic_name  d_order  eps1,eps2,eps3\n\n");
   printf("\t eps1 = error1\n");
   printf("\t eps2 = error2\n");
   printf("\t eps3 = error3\n");
   printf("\t \n\n");
  return(-1);
  }

/*   Input of the parameters and files 
*/
strcpy(generic_name,argv[1]);

err = sscanf(argv[2],"%d",&d_order);
if ((err != 1) || (d_order < 0) || (d_order > 10))
   d_order = 0;

err = sscanf(argv[3],"%f,%f,%f",&eps1,&eps2,&eps3);
if ((err != 3) || (eps1 < 0.0))
  fprintf(stderr,"\nWARNING : incorrect error parameter [%.10f]\n",eps1);

if ((err != 3) || (eps2 < 0.0))
  fprintf(stderr,"\nWARNING : incorrect error parameter [%.10f]\n",eps2);

if ((err != 3) || (eps3 < 0.0))
  fprintf(stderr,"\nWARNING : incorrect error parameter [%.10f]\n",eps3);


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
printf("\ndecomposition order of the support constraint = %d",d_order);
printf("\nthreshold for residual error = %.10f",eps);
printf("\nerror eps1 = %.10f",eps1);
printf("\nerror eps2 = %.10f",eps2);
printf("\nerror eps3 = %.10f",eps3);

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
err = (int) decomposition_mallat_ncs(phit,nx,ny,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the decomposition\n\n");
  return(-1);
  }
for (i = 0; i < nx * ny; i++) phin[i] = v[i] * phit[i];
err = (int) recomposition_mallat_ncs(phin,nx,ny,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  return(-1);
  }
err = (int) recomposition_mallat_ncs(phit,nx,ny,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  return(-1);
  }
*/
for (i = 0; i < nx * ny; i++) phin[i] = v[i] * phit[i];

for (i = 0; i < nx * ny; i++) 
  {
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

/*
err = (int) decomposition_mallat_ncs(phitr,nx,ny,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the decomposition\n\n");
  return(-1);
  }
for (i = 0; i < nx * ny; i++) 
  {
  phiti[i] = 0.0;
  rn[i] = phitr[i] * v[i];
  }
err = (int) recomposition_mallat_ncs(rn,ny,nx,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  return(-1);
  }
err = (int) recomposition_mallat_ncs(phitr,nx,ny,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  return(-1);
  }
*/
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
ii = 0;
do
{
ii++;
for (i = 0; i < nx * ny; i++) 
  {
  dntemp[i] = dnr[i];
  dni[i] = 0.0;
  }

fftw_float(dnr,dni,nx,ny,-1);

/*  compute zn and omega n  */
for (i = 0; i < nx * ny; i++) 
  {
  g2 = g[i] * g[i];
  dnr[i] *= g2;
  dni[i] *= g2;
  }

fftw_float(dnr,dni,nx,ny,1);

/*
err = (int) decomposition_mallat_ncs(dnr,nx,ny,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the decomposition\n\n");
  return(-1);
  }
for (i = 0; i < nx * ny; i++) 
   zn[i] = dnr[i] * v[i];
err = (int) recomposition_mallat_ncs(zn,nx,ny,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  return(-1);
  }
err = (int) recomposition_mallat_ncs(dnr,nx,ny,d_order);
if (err != 0)
  {
  fprintf(stderr,"\a\nFATAL ERROR : problem in the recomposition\n\n");
  return(-1);
  }
*/
for (i = 0; i < nx * ny; i++) 
   zn[i] = dnr[i] * v[i];

psdz = 0.0;
for (i = 0; i < nx * ny; i++) 
  psdz = psdz + (dntemp[i]*zn[i]);
if (psdz != 0.0)
   omegan = xnorrn / psdz;
else
   {
   fprintf(stderr,"FATAL ERROR: emergency return psdz=0\n\n");
   return(-1);
   }

/*
    compute phi n+1 , ||phi n+1||**2
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
printf("  %d       %.10f\n",ii,test);
  
/*  compute nu n  */
xnun = xnorrnp / xnorrn;

/*  compute d n+1  */
for (i = 0; i < nx * ny; i++) 
  dnr[i] = rn[i] + xnun * dntemp[i];

/*  compute polynom d n+1  */
for (k=1; k<=101; k++)
  xdn[k] = xrn[k] + xnun * xdn[k];

xnorrn = xnorrnp;
}
while((test > 0.000000001) && (ii < 30) && (itest != 0));
printf("\ntest=%.10f ii=%d itest=%d\n",test,ii,itest);

/*****************************
    END OF ITERATION
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

printf("\nNorm of the solution ||Phi n|| = %f",xnfinc);
printf("\nWeighted mean of the inverses of the eigenvalues (1/mu) = %f",xmum);
printf("\nError of image reconstruction: sqrt(eps1)=%f",sqrt(eps1));
printf("\nError due to the choice of Hr: epsilon_o=%f",sqrt(eps3));
printf("\nError related to the SNR ratio: epsilon_i=%f",sqrt(eps2));
printf("\nTerm to correct the error: delta_n=%f",deltan);
printf("\nRelative error for object reconstruction: %f",tetap);

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

/****************************************************************************
*                                                                        
* FUNCTION: decomposition_mallat_ncs                                    
*                                                                      
* PURPOSE: approximates an image with Mallat algorithm using cubic splines.
*                                                                 
* INPUT:  mat[0..ny-1][0..nx-1] = image to be decomposed
*	   ny,nx = dimensions of the image             
*	   order = order of decomposition                       
*                                                              
* OUTPUT: mat[0..ny-1][0..nx-1] = image decomposed and normalized
*                                                      
* RETURN: 0 if every thing OK                         
*         -1 if memory allocation failure            
*                                                   
* VERSION: November 1992                           
*                                                 
* AUTHORS: Karim BOUYOUCEF and Eric ANTERRIEU    
*                                               
*****************************************************************************/
#define  NPT    41
#define  NPTS2  20   /*  (NPT-1)/2  */
/*****************************************************************************/
static int  decomposition_mallat_ncs(float *mat, int nx, int ny, int order)
/*****************************************************************************/
{
register  int    i, j, ii, jj, f;
auto      int    dim, iter, fmin;
auto      int    iy, iy2, imin, imax, i2;
auto      int    ix, ix2, jmin, jmax, j2;
auto      float  *vect, H[NPT], G[NPT], xx, yy, zz;

dim = ny;  if (nx > ny) dim = nx;
vect = (float *) malloc(dim*sizeof(float));
if (vect == NULL) return(-1);

/* h filter
*/
H[0]  =  0.000000; H[1]  = -0.000164; H[2]  = -0.000202; H[3]  =  0.000327; 
H[4]  =  0.000396; H[5]  = -0.000656; H[6]  = -0.000780; H[7]  =  0.001331; 
H[8]  =  0.001546; H[9]  = -0.002745; H[10] = -0.003079; H[11] =  0.005799; 
H[12] =  0.006141; H[13] = -0.012715; H[14] = -0.012145; H[15] =  0.029747; 
H[16] =  0.022685; H[17] = -0.077808; H[18] = -0.035498; H[19] =  0.306830;
H[20] =  0.541736; H[21] =  0.306830; H[22] = -0.035498; H[23] = -0.077808; 
H[24] =  0.022685; H[25] =  0.029747; H[26] = -0.012145; H[27] = -0.012715; 
H[28] =  0.006141; H[29] =  0.005799; H[30] = -0.003079; H[31] = -0.002745; 
H[32] =  0.001546; H[33] =  0.001331; H[34] = -0.000780; H[35] = -0.000656; 
H[36] =  0.000396; H[37] =  0.000327; H[38] = -0.000202; H[39] = -0.000164;
H[40] =  0.000000;

/* g filter
*/
G[0]  =  0.000000; G[1]  =  0.000000; G[2]  =  0.000164; G[3]  = -0.000202; 
G[4]  = -0.000327; G[5]  =  0.000396; G[6]  =  0.000656; G[7]  = -0.000780; 
G[8]  = -0.001331; G[9]  =  0.001546; G[10] =  0.002745; G[11] = -0.003079; 
G[12] = -0.005799; G[13] =  0.006141; G[14] =  0.012715; G[15] = -0.012145; 
G[16] = -0.029747; G[17] =  0.022685; G[18] =  0.077808; G[19] = -0.035498;
G[20] = -0.306830; G[21] =  0.541736; G[22] = -0.306830; G[23] = -0.035498; 
G[24] =  0.077808; G[25] =  0.022685; G[26] = -0.029747; G[27] = -0.012145; 
G[28] =  0.012715; G[29] =  0.006141; G[30] = -0.005799; G[31] = -0.003079; 
G[32] =  0.002745; G[33] =  0.001546; G[34] = -0.001331; G[35] = -0.000780; 
G[36] =  0.000656; G[37] =  0.000396; G[38] = -0.000327; G[39] = -0.000202;
G[40] =  0.000164;

for (iter=0; iter<order; iter++)
  {
  iy = ny >> iter;
  ix = nx >> iter;
  iy2 = iy/2;
  ix2 = ix/2;

  for (i=0; i<iy; i++) /* transformation of lines */
    {
    for (j = 0; j < ix; j++) vect[j] = mat[j + i * nx];
    for (j = 0; j < ix2; j++)
      {
      j2 = j+j;
      jmin = j2 - NPTS2; if (jmin < 0)  jmin = 0;
      jmax = j2 + NPTS2; if (jmax > (ix-1)) jmax = ix-1;
      yy = zz = 0.;
      f = fmin = NPTS2 + (jmin - j2);
      for (jj=jmin; jj<=jmax; jj++)
        {
        xx = vect[jj];
        yy += (xx*H[f]);
        zz += (xx*G[f]);
        f++;
        }
      mat[j + i * nx] = yy;
      mat[ix2+j + i * nx] = zz;
      }
    }

  for (j=0; j<ix; j++) /* transformation of columns */
    {
    for (i=0; i<iy; i++) vect[i] = mat[j + i * nx];
    for (i=0; i<iy2; i++)
      {
      i2 = i+i;
      imin = i2 - NPTS2; if (imin < 0)  imin = 0;
      imax = i2 + NPTS2; if (imax > (iy-1)) imax = iy-1;
      yy = zz = 0.;
      f = fmin = NPTS2 + (imin - i2);
      for (ii=imin; ii<=imax; ii++)
        {
        xx = vect[ii];
        yy += (xx*H[f]);
        zz += (xx*G[f]);
        f++;
        }
      mat[j + i * nx] = 2.0*yy;
      mat[j + (iy2+i)*nx] = 2.0*zz;
      }
    }
  }

free((char *) vect);
return(0);
}

/*****************************************************************************
*                                                                   
* FUNCTION: recomposition_mallat_ncs                                 
*                                                                     
* PURPOSE: rebuilts an image with Mallat algorithm using cubic splines.
*                                                                  
* INPUT:  mat[0..ny-1][0..nx-1] = image to be recomposed  
*         ny,nx = dimensions of the image          
*         order = order of decomposition                      
*                                                              
* OUTPUT: mat[0..ny-1][0..nx-1] = image recomposed    
*                                        
* RETURN: 0 if every thing OK             
*         -1 if memory allocation failure  
*                                           
* VERSION: November 1992                     
*                                             
* AUTHORS: Karim BOUYOUCEF and Eric ANTERRIEU  
*                                               
*****************************************************************************/
#define  NPT    41
#define  NPTS2  20   /*  (NPT-1)/2  */
/*****************************************************************************/
static int  recomposition_mallat_ncs(float *mat, int nx, int ny, int order)
/*****************************************************************************/
{
register  int    i, j, ii, jj, f;
auto      int    dim, iter, fmin;
auto      int    iy, iy2, imin, imax, i2;
auto      int    ix, ix2, jmin, jmax, j2;
auto      float  *xvect, *yvect, H[NPT], G[NPT], yy, zz;

dim = ny;  if (nx > ny) dim = nx;
xvect = (float *) malloc(dim*sizeof(float));
if (xvect == NULL) return(-1);
yvect = (float *) malloc(dim*sizeof(float));
if (yvect == NULL) return(-1);

/* h filter
*/
H[0]  =  0.000000; H[1]  = -0.000164; H[2]  = -0.000202; H[3]  =  0.000327;
H[4]  =  0.000396; H[5]  = -0.000656; H[6]  = -0.000780; H[7]  =  0.001331;
H[8]  =  0.001546; H[9]  = -0.002745; H[10] = -0.003079; H[11] =  0.005799;
H[12] =  0.006141; H[13] = -0.012715; H[14] = -0.012145; H[15] =  0.029747;
H[16] =  0.022685; H[17] = -0.077808; H[18] = -0.035498; H[19] =  0.306830;
H[20] =  0.541736; H[21] =  0.306830; H[22] = -0.035498; H[23] = -0.077808;
H[24] =  0.022685; H[25] =  0.029747; H[26] = -0.012145; H[27] = -0.012715;
H[28] =  0.006141; H[29] =  0.005799; H[30] = -0.003079; H[31] = -0.002745;
H[32] =  0.001546; H[33] =  0.001331; H[34] = -0.000780; H[35] = -0.000656;
H[36] =  0.000396; H[37] =  0.000327; H[38] = -0.000202; H[39] = -0.000164;
H[40] =  0.000000;

/* g filter
*/
G[0]  =  0.000000; G[1]  =  0.000000; G[2]  =  0.000164; G[3]  = -0.000202;
G[4]  = -0.000327; G[5]  =  0.000396; G[6]  =  0.000656; G[7]  = -0.000780;
G[8]  = -0.001331; G[9]  =  0.001546; G[10] =  0.002745; G[11] = -0.003079;
G[12] = -0.005799; G[13] =  0.006141; G[14] =  0.012715; G[15] = -0.012145;
G[16] = -0.029747; G[17] =  0.022685; G[18] =  0.077808; G[19] = -0.035498;
G[20] = -0.306830; G[21] =  0.541736; G[22] = -0.306830; G[23] = -0.035498;
G[24] =  0.077808; G[25] =  0.022685; G[26] = -0.029747; G[27] = -0.012145;
G[28] =  0.012715; G[29] =  0.006141; G[30] = -0.005799; G[31] = -0.003079;
G[32] =  0.002745; G[33] =  0.001546; G[34] = -0.001331; G[35] = -0.000780;
G[36] =  0.000656; G[37] =  0.000396; G[38] = -0.000327; G[39] = -0.000202;
G[40] =  0.000164;

if (order == 0)   return(0);

for (iter=1; iter<=order; iter++)
  {
  iy = ny >> (order-iter);
  ix = nx >> (order-iter);
  iy2 = iy/2;
  ix2 = ix/2;

  for (j=0; j<ix; j++) /* transformation of columns */
    {
    for (i=0; i<iy; i++) {xvect[i] = 0.0; yvect[i] = mat[j + i * nx];}
    for (i=0; i<iy2; i++)
      {
      i2 = i+i;
      imin = i2 - NPTS2; if (imin < 0)  imin = 0;
      imax = i2 + NPTS2; if (imax > (iy-1)) imax = iy-1;
      f = fmin = NPTS2 + (imin - i2);
      yy = yvect[i];
      zz = yvect[iy2+i];
      for (ii=imin; ii<=imax; ii++)
        {
        xvect[ii] += (yy*H[f] + zz*G[f]);
        f++;
        }
      }
    for (i=0; i<iy; i++)   mat[j + i * nx] = xvect[i];
    }

  for (i=0; i<iy; i++) /* transformation of lines */
    {
    for (j=0; j<ix; j++)  {xvect[j]=0.; yvect[j] = mat[j + i * nx];}
    for (j=0; j<ix2; j++)
      {
      j2 = j+j;
      jmin = j2 - NPTS2; if (jmin < 0)  jmin = 0;
      jmax = j2 + NPTS2; if (jmax > (ix-1)) jmax = ix-1;
      f = fmin = NPTS2 + (jmin - j2);
      yy = yvect[j];
      zz = yvect[ix2+j];
      for (jj=jmin; jj<=jmax; jj++)
        {
        xvect[jj] += (yy*H[f] + zz*G[f]);
        f++;
        }
      }
    for (j=0; j<ix; j++)   mat[j + i * nx] = 2.0*xvect[j];
    }	
  }

free((char *) xvect);

free((char *) yvect);
return(0);
}
