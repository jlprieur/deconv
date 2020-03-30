/******************************************************************************
*                                                                 
* PROGRAM: dc_lissage_1D (diane_lissage)
* 1D-Version of dc_lissage
*                                                                   
* PURPOSE: compute smoothing function and raw approximation          
*          of the object                                              
*                                                                      
* INPUT:  argv[1] = generic name                                        
*         argv[2] = lower thresold for SNR (alphat)                      
*         argv[3] = PSF radius                                            
*                                                                          
* INPUT:
*     *.SNR  signal to noise ratio in fourier domain
*     *.FTB  bounded modulus FT of .PSF            
*     *.RTF  real part of fourier transform of input image
*     *.ITF  imaginary part of fourier transform of input image
* OUTPUT:
*     *.HR   synthetic aperture
*     *.FLI  smoothing filter 
*     *.PCR  FT real part of first approximation of the object
*     *.PCR  FT imaginary part of first approximation of the object
*     *.PHIT first approximation of the object
*                                                                         
* AUTHOR: JLP, SR, JV                                                      
*         translated to C by Karim BOUYOUCEF                                
*         adapted to 1-D by JLP
*                                                                            
* JLP
* Version 01/02/00
******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

#define DEBUG
#define pi 3.14159265358979323846

void  main(argc,argv)
int   argc;
char  *argv[];
{
/*    DECLARATIONS
*/
char     filename[61], comments[81], generic_name[61], FTB_name[61];
register int i, j, ix, iy;
INT4     ny, nx, nx1, ny1, one;
INT_PNTR pntr_ima;
int      err, flag;
float    *SNR, *FTB, *RFT, *IFT, *sectr, *elipr, *smooth0, *smooth1;
float    *temp, *visib;
float    alphat, radius, cumul;
float    sum, error, eta0, epsilon, X2, Y2;
float    ratio, norm, rr, rp, test_real, p_axe;
int      itere, c2, d2, is, it_max;
int      ixc, iyc, irp, r2p, ic;

/*    TEST OF COMMAND LINE
*/
/* JLP99: to allow for batch commands... */
if(argc == 7) argc = 5;
if (argc != 5)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndc_lissage_1D lower_snr_thresold(alphat) generic_name(with dot)  PSF_radius transfer_function(.FTB)\n\n");
  fprintf(stderr,"lower_snr_thresold = alphat");
  fprintf(stderr,"\nPSF_radius = radius for synthetic aperture *.HR\n\n");
  exit(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
err = sscanf(argv[1],"%f",&alphat);
if (err !=1)
  {
  fprintf(stderr,"\nFATAL ERROR: Lower thresold [%s] incorrect\n\n",argv[2]);
  exit(-1);
  }
strcpy(generic_name,argv[2]);
strcpy(FTB_name,argv[4]);

JLP_BEGIN();
JLP_INQUIFMT();

/*    READ FILE SIZE
*/
JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, FTB_name, comments);
if ((ny < 1) || (ny > 1024) || (nx < 1) || (nx > 1024))
   {
   fprintf(stderr,"\nFATAL ERROR: invalid dimensions [%d,%d] from %s\n\n",
           ny,nx,FTB_name);
   exit(-1);
   }
FTB = (float *)pntr_ima;
if (FTB == NULL)
  {
  fprintf(stderr,"FATAL ERROR: when loading file [%s]\n\n",FTB_name);
  exit(-1);
  }
RECENT_FFT_1D_Y_FLOAT(FTB, FTB, &nx, &ny, &nx);

err = sscanf(argv[3],"%f",&radius);
if ((err !=1) || (radius < 1.) || (radius > (float)(ny/2-1)) )
  {
  fprintf(stderr,"\nFATAL ERROR: radius for synthetic aperture [%s] incorrect\n\n",argv[3]);
  exit(-1);
  }

/*   DISPLAY  alphat radius
*/
printf("PROGRAM :    dc_lissage_1D");
printf("**************************\n");
printf("alphat = %f\n",alphat);
printf("radius of Hr (fourier space) = %f\n",radius);

/*    MEMORY ALLOCATION
*/
temp = (float*) malloc(nx * ny * sizeof(float));
smooth1 = (float*) malloc(ny * sizeof(float));
sectr = (float*) malloc(ny * sizeof(float));

/* Computing 1-D smoothing function (to be applied to columns)
*/
eta0 = 1.55;
it_max = 50;
epsilon = 0.00001;
iyc = ny/2;
rr = eta0*ny/(radius*pi);
fprintf(stderr,"\nradius of smoothing function (direct space) = %f",rr);
/* Building "sectr" as the support in Fourier domain */
rp = radius;
irp = (int)rp+1;
r2p = irp*irp;
is = 0;
for (j=0; j<ny; j++)
   {
   c2 = (j-iyc)*(j-iyc);
     if (c2 <= r2p)
       {
       sectr[j] = 1.0;
/* Integral in Fourier domain: */
       is++;
       }
   }
one = 1;
RECENT_FFT_1D_Y_FLOAT(sectr, sectr, &one, &ny, &one);

/* ellipsoide of resolution: building "elipr" as the support in direct space*/
/* I take the radius computed for direct space: */
p_axe = rr;
ic = 0;
p_axe *= p_axe;
elipr = (float*) malloc(ny * sizeof(float));
for (j = 0; j < ny; j++)
   {
   c2 = (j-iyc)*(j-iyc);
   test_real = (float)c2/p_axe;
   if (test_real <= 1.0)
      {
      elipr[j] = 1.0;
/* Integral in direct space: */
      ic++;
      }
   smooth1[j] = elipr[j];
   }

eta0 = sqrt((float)(ic*is)/(float)ny);
#ifdef DEBUG
printf("is = %d ic = %d eta0=%e\n",is, ic, eta0);
#endif

/* compute smoothing function by power method */
smooth0 = (float*) malloc(ny * sizeof(float));
itere = 0;
error = epsilon + 1.0;

/******* do while loop: begin ***********************************************/
do
   {
   for (i = 0; i < ny; i++)
     {
     smooth0[i] = smooth1[i];
     temp[i] = 0.0;
     }
/* FFT and troncation in Fourier domain */
   fftw_1D_Y_float(smooth1, temp,(int)one,(int)ny,1);
   for (i = 0; i < ny; i++)
     {
     if (sectr[i] == 0.0)
       {
       smooth1[i] = 0.0;
       temp[i] = 0.0;
       }
     }
/* FFT-1 and troncation in direct space */
   fftw_1D_Y_float(smooth1, temp,(int)one,(int)ny,-1);
   for (i = 0; i < ny; i++)
     {
     if (elipr[i] == 0.0)
       {
       smooth1[i] = 0.0;
       temp[i] = 0.0;
       }
     }
/* Normalization and convergence test: */
   norm = 0.0;
   for (i = 0; i < ny; i++)
     {
     X2 = smooth1[i] * smooth1[i];
     Y2 = temp[i] * temp[i];
     norm += (X2+Y2);
     smooth1[i] = (float)sqrt(X2+Y2);
     }
   norm = (float)sqrt((double)norm);
   printf(" norm = %f iter=%d \n",norm,itere);
   error = 0.0;
   for (i = 0; i < ny; i++)
     {
/* normalization: */
     smooth1[i] /= norm;
/* (difference between present value (smooth1) and previous (smooth0)*/
     X2 = (smooth1[i] - smooth0[i]) * (smooth1[i] - smooth0[i]);
     error = error + X2;
     }
   error = (float)sqrt((double)error);
   itere++;
   }
while ((itere < it_max) && (error > epsilon));
/******* do while loop: end *************************************************/

fprintf(stderr,"\nnumber of iterations = %d",itere);
fprintf(stderr,"\nError = %f",error);
fprintf(stderr,"\nKhi = %f",norm);

/*   FREE MEMORY
*/
free(elipr);
free(smooth0);

/* computes smoothing function  */
sum = 0.0;
for (i = 0; i < ny; i++)
   {
   temp[i] = 0.0;
   sum += smooth1[i];
   if(i < 5) printf(" smooth1[%d] = %e \n",i,smooth1[i]);
   }
for (i = 0; i < ny; i++)
   smooth1[i] /= sum;

fftw_1D_Y_float(smooth1, temp,(int)one,(int)ny,1);
for (i = 0; i < ny; i++)
  smooth1[i] = (float)sqrt((double)(smooth1[i]*smooth1[i]+temp[i]*temp[i]));  

/*  STORE RESULT  *.HR
*/
RECENT_FFT_1D_Y_FLOAT(sectr, sectr, &one, &ny, &one);
sprintf(filename,"%sHR",generic_name);
sprintf(comments,"dc_lissage_1D");
JLP_WRITEIMAG(sectr, &ny, &one, &ny, filename, comments);

/*   COMPUTES PHIT
     first approximation of the object
*/
sprintf(filename,"%sSNR",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
SNR = (float *)pntr_ima;
if (SNR == NULL)
  {
  fprintf(stderr,"FATAL ERROR: when loading file [%s]\n\n",filename);
  exit(-1);
  }
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",
          filename);
  exit(-1);
  } 
RECENT_FFT_1D_Y_FLOAT(SNR, SNR, &nx, &ny, &nx);

sprintf(filename,"%sRFT",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
RFT = (float *)pntr_ima;
if (RFT == NULL)
  {
  fprintf(stderr,"FATAL ERROR: when loading file [%s]\n\n",filename);
  exit(-1);
  }
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",
          filename);
  exit(-1);
  } 
RECENT_FFT_1D_Y_FLOAT(RFT, RFT, &nx, &ny, &nx);

sprintf(filename,"%sIFT",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
IFT = (float *)pntr_ima;
if (IFT == NULL)
  {
  fprintf(stderr,"FATAL ERROR: when loading file [%s]\n\n",filename);
  exit(-1);
  }
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",
          filename);
  exit(-1);
  } 
RECENT_FFT_1D_Y_FLOAT(IFT, IFT, &nx, &ny, &nx);

/* Compute  .PCR  .PCI
*/
for (ix = 0; ix < nx; ix++)
{
for (iy = 0; iy < ny; iy++)
   {
   i = ix + iy * nx;
   if ((SNR[i] >= alphat) && (FTB[i] != 0.0))
      {
      ratio = smooth1[iy]/FTB[i];
      RFT[i] *= ratio;
      IFT[i] *= ratio;
      }
   else
      {
      RFT[i] = 0.0;
      IFT[i] = 0.0;
      }
    }
}

/* Compute pseudo visibility (modulus of the ratio of the moduli): */
visib = (float*) malloc(nx * ny * sizeof(float));
for (ix = 0; ix < nx; ix++)
{
for (iy = 0; iy < ny; iy++)
   {
   i = ix + iy * nx;
   if ((SNR[i] >= alphat) && (FTB[i] != 0.0))
   {
   visib[i] = (float)sqrt((double)(RFT[i]*RFT[i]+IFT[i]*IFT[i]))/smooth1[iy];
   }
   else
   visib[i] = 0.; 
   }
}
/*  Store results  .FLI  .PCR  .PCI
*/
RECENT_FFT_1D_Y_FLOAT(smooth1, smooth1, &one, &ny, &one);
sprintf(filename,"%sFLI",generic_name);
sprintf(comments,"dc_lissage_1D");
JLP_WRITEIMAG(smooth1, &ny, &one, &ny, filename, comments);

RECENT_FFT_1D_Y_FLOAT(RFT, temp, &nx, &ny, &nx);
sprintf(filename,"%sPCR",generic_name);
sprintf(comments,"dc_lissage_1D");
JLP_WRITEIMAG(temp, &nx, &ny, &nx, filename, comments);

RECENT_FFT_1D_Y_FLOAT(IFT, temp, &nx, &ny, &nx);
sprintf(filename,"%sPCI",generic_name);
sprintf(comments,"dc_lissage_1D");
JLP_WRITEIMAG(temp, &nx, &ny, &nx, filename, comments);

/* STORE RESULT  .VISI
*/
RECENT_FFT_1D_Y_FLOAT(visib, visib, &nx, &ny, &nx);
sprintf(filename,"%sVISI",generic_name);
sprintf(comments,"dc_lissage_1D");
JLP_WRITEIMAG(visib, &nx, &ny, &nx, filename, comments);

/*    FREE MEMORY
*/
free(SNR);
free(FTB);
free(smooth1);
free(temp);

/* COMPUTE PHIT
*/
fftw_1D_Y_float(RFT,IFT,(int)nx,(int)ny,-1);

for (i = 0; i < nx * ny; i++)
   RFT[i] = (float)sqrt((double)(RFT[i]*RFT[i]+IFT[i]*IFT[i]));

/* STORE RESULT  .PHIT
*/
sprintf(filename,"%sPHIT",generic_name);
sprintf(comments,"dc_lissage_1D");
RECENT_FFT_1D_Y_FLOAT(RFT, RFT, &nx, &ny, &nx);
JLP_WRITEIMAG(RFT, &nx, &ny, &nx, filename, comments);
free(RFT);
free(IFT);
fprintf(stderr,"\n\n");

JLP_END();
}
