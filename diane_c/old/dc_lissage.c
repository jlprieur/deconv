/******************************************************************************
*                                                                 
* PROGRAM: dc_lissage (diane_lissage)                                           
*                                                                   
* PURPOSE: compute smoothing function and raw approximation          
*          of the object                                              
*                                                                      
* INPUT:  argv[1] = generic name                                        
*         argv[2] = lower thresold for SNR (alphat)
*         argv[3] = transf. function radius
*         argv[4] = transf. function (.FTB)
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
* From Karims'version of December 1992 
*                                                                         
* AUTHOR: JLP, SR, JV                                                      
*         translated to C by Karim BOUYOUCEF                                
*                                                                            
* JLP
* Version 10/03/08
******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

int  main(int argc, char *argv[])
{
/*    DECLARATIONS
*/
char     filename[61], comments[81], generic_name[61], FTB_name[61];
INT4     ny, nx, nx1, ny1;
INT_PNTR pntr_ima;
float    *SNR, *FTB, *RFT, *IFT, *sectr, *elipr, *stock, *tempr, *tempi, *visib;
float    alphat, radius;
float    somme, error, eta0, epsilon, X2, Y2;
float    ratio, norm, norme, rr, rp, test_real, p_axe, gr_axe;
int      err, itere, c2, d2, is, it_max;
int      ixc, iyc, irp, r2p, ic;
register int i, j;

/*    TEST OF COMMAND LINE
*/
/* JLP99: to allow for batch commands... */
if(argc == 7) argc = 5;
if (argc != 5)
  {
  printf("Fatal error/unexpected number of arguments\n");
  printf("\nUSAGE:\n");
  printf("\ndc_lissage  lower_snr_threshold generic_name(with dot)  transf_function_radius transfer_function(.FTB)\n\n");
  printf("lower_snr_thresold = alphat");
  printf("\ntransf_function_radius = radius for synthetic aperture *.HR\n\n");
  printf("transfer_function(.FTB) or 1 if unity (double stars...)\n");
  return(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
err = sscanf(argv[1],"%f",&alphat);
if (err !=1) {
  printf("\nFATAL ERROR: Lower thresold [%s] incorrect\n\n",argv[2]);
  return(-1);
  }
strcpy(generic_name,argv[2]);
strcpy(FTB_name,argv[4]);

JLP_INQUIFMT();

/*  Read file size with *.SNR
*/
sprintf(filename,"%sSNR",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, filename, comments);
SNR = (float *)pntr_ima;
if (SNR == NULL) {
  printf("FATAL ERROR: when loading file [%s]\n\n",filename);
  return(-1);
  }
if ((ny < 1) || (ny > 1024) || (nx < 1) || (nx > 1024))
   {
   printf("\nFATAL ERROR: invalid dimensions [%d,%d] from %s\n\n",
           ny,nx,FTB_name);
   return(-1);
   }
RECENT_FFT(SNR, SNR, &nx, &ny, &nx);

sprintf(filename,"%sRFT",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
RFT = (float *)pntr_ima;
if (RFT == NULL) {
  printf("FATAL ERROR: when loading file [%s]\n\n",filename);
  return(-1);
  }
if ((nx != nx1) || (ny != ny1)) {
  printf("\nFATAL ERROR: file [%s] have incompatible size\n\n",
          filename);
  return(-1);
  } 
RECENT_FFT(RFT, RFT, &nx, &ny, &nx);

sprintf(filename,"%sIFT",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
IFT = (float *)pntr_ima;
if (IFT == NULL) {
  printf("FATAL ERROR: when loading file [%s]\n\n",filename);
  return(-1);
  }
if ((nx != nx1) || (ny != ny1)) {
  printf("\nFATAL ERROR: file [%s] have incompatible size\n\n",
          filename);
  return(-1);
  } 
RECENT_FFT(IFT, IFT, &nx, &ny, &nx);

/* Transfer function or unity field: */
if(FTB_name[0] == '1' && 
(FTB_name[1] == 0 || FTB_name[1] == '.' || FTB_name[1] == ' ')) {
  printf("OK: Unity field for transfer function\n");
  FTB = (float *)malloc(nx * ny * sizeof(float));
  for(i = 0; i < nx * ny; i++) FTB[i] = 1.;
  } else { 
  JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, FTB_name, comments);
  FTB = (float *)pntr_ima;
  if (FTB == NULL) {
    printf("FATAL ERROR: when loading file [%s]\n\n",FTB_name);
    return(-1);
    }
  if ((nx != nx1) || (ny != ny1)) {
    printf("\nFATAL ERROR: file [%s] have incompatible size\n\n",
            filename);
    return(-1);
    } 
  RECENT_FFT(FTB, FTB, &nx, &ny, &nx);
  }

err = sscanf(argv[3],"%f",&radius);
if ((err !=1) || (radius < 1.) || (radius > (float)(ny/2-1)) || (radius > (float)(nx/2-1)))
  {
  printf("\nFATAL ERROR: radius for synthetic aperture [%s] incorrect\n\n",argv[3]);
  return(-1);
  }

/*   DISPLAY  alphat radius
*/
printf("PROGRAM :    dc_lissage");
printf("**************************\n");
printf("alphat = %f\n",alphat);
printf("radius of Hr (fourier space) = %f\n",radius);

/*    MEMORY ALLOCATION
*/
tempr = (float*) malloc(ny * nx * sizeof(float));
sectr = (float*) malloc(ny * nx * sizeof(float));

/*    COMPUTING SMOOTHING FUNCTION
*/
eta0 = 1.55;
it_max = 50;
epsilon = 0.00001;
iyc = ny/2;
ixc = nx/2;
rr = eta0*ny/(radius*PI);
printf("\nradius of smoothing function (direct space) = %f",rr);
rp = radius;
irp = (int)rp+1;
r2p = irp*irp;
is = 0;
for (j=0; j<ny; j++)
   {
   c2 = (j-iyc)*(j-iyc);
   for (i=0; i<nx; i++)
     {
     if ((c2+(i-ixc)*(i-ixc)) <= r2p)
       {
       sectr[i + j * nx] = 1.0;
       is++;
       }
     }
   }
RECENT_FFT(sectr, sectr, &nx, &ny, &nx);

/* ellipsoide of resolution */
p_axe = rr;
gr_axe = p_axe;
ic = 0;
p_axe *= p_axe;
gr_axe *= gr_axe;
elipr = (float*) malloc(nx * ny * sizeof(float));
for (j = 0; j < ny; j++)
for (i = 0; i < nx; i++)
   {
   c2 = (j-iyc)*(j-iyc);
   d2 = (i-ixc)*(i-ixc);
   test_real = (float)c2/p_axe+(float)d2/gr_axe;
   if (test_real <= 1.0)
      {
      elipr[i + j * nx] = 1.0;
      ic++;
      }
   tempr[i + j * nx] = elipr[i + j * nx];
   }

eta0 = sqrt((float)(ic*is)/(float)ny);

/* Compute smoothing function with the power method */
stock = (float*) malloc(nx * ny * sizeof(float));
tempi = (float*) malloc(nx * ny * sizeof(float));
itere = 0;
error = epsilon + 1.0;

do
   {
   for (i = 0; i < nx * ny; i++)
     {
     stock[i] = tempr[i];
     tempi[i] = 0.0;
     }
   fftw_float(tempr, tempi,(int)nx,(int)ny,1);
   for (i = 0; i < nx * ny; i++)
     if (sectr[i] == 0.0)
       {
       tempr[i] = 0.0;
       tempi[i] = 0.0;
       }
   fftw_float(tempr, tempi,(int)nx,(int)ny,-1);
   for (i = 0; i < nx * ny; i++)
     if (elipr[i] == 0.0)
       {
       tempr[i] = 0.0;
       tempi[i] = 0.0;
       }
   /* normalization and convergence test */
   norme = 0.0;
   error = 0.0;
   for (i = 0; i < nx * ny; i++)
     {
     X2 = tempr[i] * tempr[i];
     Y2 = tempi[i] * tempi[i];
     norme += (X2+Y2);
     tempr[i] = (float)sqrt(X2+Y2);
     }
   norme = (float)sqrt((double)norme);
   for (i = 0; i < nx * ny; i++)
     {
     tempr[i] /= norme;
     X2 = (tempr[i] - stock[i]) * (tempr[i] - stock[i]);
     error = error + X2;
     }
   error = (float)sqrt((double)error);
   itere++;
   }
while ((itere < it_max) && (error > epsilon));

printf("\nnumber of iterations = %d",itere);
printf("\nError = %f",error);
printf("\nKhi = %f",norme);

/*   FREE MEMORY
*/
free(elipr);
free(stock);

/* computes smoothing function  */
somme = 0.0;
for (i = 0; i < nx * ny; i++)
   {
   tempi[i] = 0.0;
   somme += tempr[i];
   }
for (i = 0; i < nx * ny; i++)
   tempr[i] /= somme;

fftw_float(tempr, tempi,(int)nx,(int)ny,1);
for (i = 0; i < nx * ny; i++)
  tempr[i] = (float)sqrt((double)(tempr[i]*tempr[i]+tempi[i]*tempi[i]));  

/*   FREE MEMORY
*/
free(tempi);

/*   STORE RESULT  *.HR
*/
RECENT_FFT(sectr, sectr, &nx, &ny, &nx);
sprintf(filename,"%sHR",generic_name);
sprintf(comments,"dc_lissage");
JLP_WRITEIMAG(sectr, &nx, &ny, &nx, filename, comments);

/*    COMPUTES PHIT
      first approximation of the object
*/

/*    COMPUTE  .PCR  .PCI
*/
for (i = 0; i < nx * ny; i++)
   {
   if ((SNR[i] >= alphat) && (FTB[i] != 0.0))
      {
      ratio = tempr[i]/FTB[i];
      RFT[i] *= ratio;
      IFT[i] *= ratio;
      }
   else
      {
      RFT[i] = 0.0;
      IFT[i] = 0.0;
      }
    }

/* Compute pseudo visibility (modulus of the ratio of the moduli): */
visib = (float*) malloc(nx * ny * sizeof(float));
for (i = 0; i < nx * ny; i++)
   {
   if ((SNR[i] >= alphat) && (FTB[i] != 0.0))
   {
   visib[i] = (float)sqrt((double)(RFT[i]*RFT[i]+IFT[i]*IFT[i]))/tempr[i];
   }
   else
   visib[i] = 0.; 
   }
/*   STORE RESULTS  .FLI  .PCR  .PCI
*/
RECENT_FFT(tempr, tempr, &nx, &ny, &nx);
sprintf(filename,"%sFLI",generic_name);
sprintf(comments,"dc_lissage");
JLP_WRITEIMAG(tempr, &nx, &ny, &nx, filename, comments);

RECENT_FFT(RFT, RFT, &nx, &ny, &nx);
sprintf(filename,"%sPCR",generic_name);
sprintf(comments,"dc_lissage");
JLP_WRITEIMAG(RFT, &nx, &ny, &nx, filename, comments);

RECENT_FFT(IFT, IFT, &nx, &ny, &nx);
sprintf(filename,"%sPCI",generic_name);
sprintf(comments,"dc_lissage");
JLP_WRITEIMAG(IFT, &nx, &ny, &nx, filename, comments);

/* STORE RESULT  .VISI
*/
RECENT_FFT(visib, visib, &nx, &ny, &nx);
sprintf(filename,"%sVISI",generic_name);
sprintf(comments,"dc_lissage");
JLP_WRITEIMAG(visib, &nx, &ny, &nx, filename, comments);

/*    FREE MEMORY
*/
free(SNR);
free(FTB);
free(tempr);

/* COMPUTE PHIT
*/
norm = 0.0;
for (i = 0; i < nx * ny; i++)
   norm += (RFT[i]*RFT[i]+IFT[i]*IFT[i]);
norm = (float)sqrt((double)norm/(double)(ny*nx));
printf("\nnorm of smoothing function = %f",norm);

RECENT_FFT(RFT, RFT, &nx, &ny, &nx);
RECENT_FFT(IFT, IFT, &nx, &ny, &nx);
fftw_float(RFT,IFT,(int)nx,(int)ny,-1);

for (i = 0; i < nx * ny; i++)
   RFT[i] = (float)sqrt((double)(RFT[i]*RFT[i]+IFT[i]*IFT[i]));


/*   STORE RESULT  .PHIT
*/
sprintf(filename,"%sPHIT",generic_name);
sprintf(comments,"dc_lissage");
RECENT_FFT(RFT, RFT, &nx, &ny, &nx);
JLP_WRITEIMAG(RFT, &nx, &ny, &nx, filename, comments);
free(RFT);
free(IFT);
printf("\n\n");

return(0);
}
