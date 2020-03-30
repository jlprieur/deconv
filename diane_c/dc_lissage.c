/******************************************************************************
*                                                                 
* PROGRAM: dc_lissage (diane_lissage)                                           
*                                                                   
* PURPOSE: compute smoothing function and raw approximation          
*          of the object                                              
*                                                                      
* INPUT:  argv[1] = generic name                                        
*         argv[2] = lower thresold for SNR (alphat)
*         argv[3] = synthetic aperture: x and y radii
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
*     *.VISI visibility function (power spectrum)
*                                                                       
* From Karims'version of December 1992 (but very little has been left from
* his version...)
*                                                                         
* JLP
* Version 10/04/08
******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

int compute_smooth_function(float f_xradius, float f_yradius, float *SNR, 
                            float *FTB, float *FLI, 
                            float *HR, int nx, int ny, float alphat);
static int output_visib(float *SNR, float *RFT, float *IFT, float *FTB,
                        float alphat, int nx, int ny, char *generic_name);
static int output_phit_pcr_pci(float *SNR, float *RFT, float *IFT, float *FTB,
                               float *FLI, float alphat, INT4 nx, INT4 ny, 
                               char *generic_name);

int  main(int argc, char *argv[])
{
char     filename[61], comments[81], generic_name[61], FTB_name[61];
INT4     ny, nx, nx1, ny1;
INT_PNTR pntr_ima;
float    *SNR, *FTB, *RFT, *IFT, *FLI, *HR;
float    alphat, f_xradius, f_yradius;
int      ival;
register int i;

/*    TEST OF COMMAND LINE
*/
/* JLP99: to allow for batch commands... */
if(argc == 7) argc = 5;
if (argc != 5)
  {
  printf("Fatal error/unexpected number of arguments\n");
  printf("\nUSAGE:\n");
  printf("\ndc_lissage  snr_lower_threshold generic_name(with dot) f_xradius,f_yradius transfer_function(.FTB)\n\n");
  printf("snr_lower_thresold = alphat");
  printf("\nf_xradius, f_yradius = x (or) y radius for synthetic aperture *.HR\n\n");
  printf("transfer_function(.FTB) or 1 if unity (double stars...)\n");
  return(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
ival = sscanf(argv[1],"%f",&alphat);
if (ival !=1) {
  printf("\nFATAL ERROR: SNR lower thresold [%s] incorrect\n\n",argv[2]);
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
RECENT_FFT(SNR, SNR, &nx, &ny, &nx);

/* Reading real part of input experimental spectrum: RFT
*/
sprintf(filename,"%sRFT",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
RFT = (float *)pntr_ima;
if ((nx != nx1) || (ny != ny1)) {
  printf("\nFATAL ERROR: file [%s] have incompatible size\n\n",
          filename);
  return(-1);
  } 
RECENT_FFT(RFT, RFT, &nx, &ny, &nx);

/* Reading imaginary part of input experimental spectrum: IFT
*/
sprintf(filename,"%sIFT",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
IFT = (float *)pntr_ima;
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
  if ((nx != nx1) || (ny != ny1)) {
    printf("\nFATAL ERROR: file [%s] have incompatible size\n\n",
            filename);
    return(-1);
    } 
  RECENT_FFT(FTB, FTB, &nx, &ny, &nx);
  }

ival = sscanf(argv[3],"%f,%f",&f_xradius,&f_yradius);
if(ival == 1) f_yradius = f_xradius;
if ((ival != 1 && ival != 2) 
     || (f_xradius < 1.) || (f_xradius > (float)(nx/2-1)) 
     || (f_yradius < 1.) || (f_yradius > (float)(ny/2-1)))
  {
  printf("\ndc_lissage/FATAL ERROR: wrong radii for synthetic aperture [radii=%f,%f]!\n",
        f_xradius, f_yradius);
  if ((f_xradius > (float)(nx/2-1)) || (f_yradius > (float)(ny/2-1)))
    printf("Warning, maximum values: f_xradius=%.2f, f_yradius=%.2f\n",
           (float)(nx/2-1), (float)(ny/2-1));
  printf("\n");
  return(-1);
  }

/*  Display input parameters 
*/
printf("PROGRAM :    dc_lissage");
printf("**************************\n");
printf("alphat = %f\n",alphat);
printf("x and y radii of Hr (Fourier space) = %f,%f\n", f_xradius, f_yradius);

/*    MEMORY ALLOCATION
*/
FLI = (float*) malloc(ny * nx * sizeof(float));
HR = (float*) malloc(ny * nx * sizeof(float));

/* Compute the smoothing function FLI : */
compute_smooth_function(f_xradius, f_yradius, SNR, FTB, FLI, HR, nx, ny, alphat);

/*   STORE RESULT  *.HR
*/
RECENT_FFT(HR, HR, &nx, &ny, &nx);
sprintf(filename,"%sHR",generic_name);
sprintf(comments,"dc_lissage: f_xradius=%.3f f_yradius=%.3f", 
         f_xradius, f_yradius);
JLP_WRITEIMAG(HR, &nx, &ny, &nx, filename, comments);

printf("JLP OK1\n");
/* Computes the visibility of the first approximation image PHIT:
*/
 output_visib(SNR, RFT, IFT, FTB, alphat, nx, ny, generic_name);

printf("JLP OK1\n");
/* Computes the first approximation image PHIT and the
* corresponding spectrum PCR and PCI
*/
 output_phit_pcr_pci(SNR, RFT, IFT, FTB, FLI, alphat, nx, ny, generic_name);
printf("JLP OK1\n");


/* Save .FLI
*/
RECENT_FFT(FLI, FLI, &nx, &ny, &nx);
sprintf(filename,"%sFLI",generic_name);
sprintf(comments,"From dc_lissage (alphat=%.3f)", alphat);
JLP_WRITEIMAG(FLI, &nx, &ny, &nx, filename, comments);

/* Free memory
*/
free(SNR);
free(FTB);
free(FLI);
return(0);
}

/**********************************************************************
* Computes smoothing function "FLI" (in Fourier domain) 
* with successive FFTs and truncations both in direct and Fourier space
*
* INPUT:
*  xradius: Radius of HR in Fourier space along the X axis
*  yradius: Radius of HR in Fourier space along the Y axis
*           (*.HR=synthetic aperture)
*  axis_ratio (X/Y): axis ratio of the ellipse
*  SNR: signal to noise ratio of the input experimental spectrum
*  FTB: bounded transfer function
*  alphat: lower snr threshold to be applied on the experimental spectrum 
*
* OUTPUT:
*  HR: synthetic aperture
*  FLI: smoothing function (in Fourier domain)
**********************************************************************/
int compute_smooth_function(float f_xradius, float f_yradius, float *SNR, 
                            float *FTB, float *FLI, float *HR, int nx, int ny, 
                            float alphat)
{
float    *ellipr, *old_tempr, *tempi, *tempr;
float    sum, error, eta0, epsilon, w2;
float    rrx, rry, rrx2, rry2, rrfx2, rrfy2;
double   norm;
int      itere, is, it_max;
int      ixc, iyc, ic;
register int i, j;

old_tempr = (float*) malloc(nx * ny * sizeof(float));
tempi = (float*) malloc(nx * ny * sizeof(float));
tempr = (float*) malloc(nx * ny * sizeof(float));

/* Computing smoothing function
*/
/* If disk (in Fourier space) and disk (direct space) configuration:
* eta0 = 1.55 for Khi**2=0.95
*/
eta0 = 1.55;
/* Test value of the error, used to stop iterations for the power method */
epsilon = 1.e-5;
iyc = ny/2;
ixc = nx/2;
/* Remember: if radius in Fourier space,
*  f_xradius/nx spatial freq. thus radius = nx/f_xradius in direct space 
*/
rrx = eta0*nx/(f_xradius*PI);
rry = eta0*ny/(f_yradius*PI);
printf("The limit of resolution is set so that eta0=1.55 (disks)\n");
printf("X and Y radii of smoothing function (direct space) = %.3f %.3f\n", 
        rrx, rry);

/* Computation of the Fourier mask "HR",
* i.e. the characteristic function in Fourier space (to be used for truncation) 
*/
rrfx2 = SQUARE(f_xradius);
rrfy2 = SQUARE(f_yradius);
is = 0;
for (j=0; j<ny; j++)
   {
   for (i=0; i<nx; i++)
     {
     w2 = SQUARE(j-iyc) / rrfy2 + SQUARE(i-ixc) / rrfx2;
     if (w2 <= 1.0)
       {
       HR[i + j * nx] = 1.0;
       is++;
       }
      else HR[i + j * nx] = 0.;
     }
   }
/* Recenter HR to work in the frequency domain */
RECENT_FFT(HR, HR, &nx, &ny, &nx);

/* Computing the ellipsoid of resolution in direct space
* and the ellipsoid resolution function "ellipr" (to be used for truncation) */
ic = 0;
rrx2 = SQUARE(rrx);
rry2 = SQUARE(rry);
ellipr = (float*) malloc(nx * ny * sizeof(float));
for (j = 0; j < ny; j++)
for (i = 0; i < nx; i++)
   {
   w2 = SQUARE(j-iyc) / rry2 + SQUARE(i-ixc) / rrx2;
   if (w2 <= 1.0)
      {
      ellipr[i + j * nx] = 1.0;
      ic++;
      }
   else ellipr[i + j * nx] = 0.;
   tempr[i + j * nx] = ellipr[i + j * nx];
   }

/* Re-evaluation of eta0 */
eta0 = sqrt((float)(ic*is))/(float)ny;
  printf("Actual value of interpolation coefficient eta0: %.3f\n", eta0);
  printf("(Don't worry as long as eta0 is within [1.4 , 2.2])\n");

/* Compute smoothing function with the power method */
itere = 0;
error = epsilon + 1.0;

for (i = 0; i < nx * ny; i++) tempr[i] = 1.0;

/* Series of Fourier transforms and truncations */
it_max = 50;
do
   {
/* Load previous "tempr" in order to check the convergence */
   for (i = 0; i < nx * ny; i++)
     {
     old_tempr[i] = tempr[i];
     tempi[i] = 0.0;
     }

/* F.T. and truncation in Fourier space
*/
   fftw_float(tempr, tempi,(int)nx,(int)ny,1);
   for (i = 0; i < nx * ny; i++)
     if (HR[i] == 0.0)
       {
       tempr[i] = 0.0;
       tempi[i] = 0.0;
       }

/* Inverse F.T. and truncation in direct space
*/
   fftw_float(tempr, tempi,(int)nx,(int)ny,-1);
   for (i = 0; i < nx * ny; i++)
     if (ellipr[i] == 0.0)
       {
       tempr[i] = 0.0;
       tempi[i] = 0.0;
       }
/* Normalization of "tempr" and convergence test */
   norm = 0.0;
   error = 0.0;
   for (i = 0; i < nx * ny; i++)
     {
     w2 = SQUARE(tempr[i]) + SQUARE(tempi[i]);
     norm += w2;
     tempr[i] = (float)sqrt(w2);
     }
   norm = sqrt(norm);
   for (i = 0; i < nx * ny; i++)
     {
     tempr[i] /= norm;
     w2 = SQUARE(tempr[i] - old_tempr[i]);
     error = error + w2;
     }
   error = (float)sqrt((double)error);
   itere++;
   }
while ((itere < it_max) && (error > epsilon));

printf("\nEnd of computation of the smoothing function with the power method.\n");
printf("Number of iterations = %d, error=%.3e, Khi=%.3f\n",
        itere, error, norm);

/* Free memory
*/
free(ellipr);
free(old_tempr);

/* Computes smoothing function  */
/* Normalization of tempr (with L1 norm) in direct space: */
sum = 0.0;
for (i = 0; i < nx * ny; i++) sum += tempr[i];
for (i = 0; i < nx * ny; i++) tempr[i] /= sum;

/* 
* Final smoothing function is taken as the modulus of the FFT of that
* normalized functioni (i.e. assumes it has zero phase):
*/
for (i = 0; i < nx * ny; i++) tempi[i] = 0.0;
fftw_float(tempr, tempi,(int)nx,(int)ny,1);

for (i = 0; i < nx * ny; i++)
  FLI[i] = (float)sqrt((double)SQUARE(tempr[i])+(double)SQUARE(tempi[i]));  

/* Free memory
*/
free(tempr);
free(tempi);

return(0);
}
/*******************************************************************
* Computes visibility function and output first approximation image (PHIT)
* i.e., the ratio of the modulus of the experimental spectrum PCR,PCI by FTB.
*
* INPUT:
* SNR: snr of the input spectrum
* RFT, IFT: experimental input spectrum
* FTB: bounded tranfer function
* alphat: snr lower threshold to be applied on the input spectrum
* nx, ny: size of arrays
* generic_name: name to be used for the output file
*
* OUTPUT:
* none 
********************************************************************/
static int output_visib(float *SNR, float *RFT, float *IFT, float *FTB,
                        float alphat, INT4 nx, INT4 ny, char *generic_name)
{
register int i;
float *VISIB;
char filename[60], comments[80];

/* Compute the visibility of the first approximation image (PHIT) 
* (i.e., ratio of the modulus of the experimental spectrum PCR,PCI by FTB)
*/
VISIB = (float*) malloc(nx * ny * sizeof(float));
for (i = 0; i < nx * ny; i++)
   {
   if ((SNR[i] >= alphat) && (FTB[i] != 0.0))
   {
   VISIB[i] = (float)sqrt((double)(RFT[i]*RFT[i]+IFT[i]*IFT[i]))/FTB[i];
   }
   else
   VISIB[i] = 0.; 
   }

/* Store result in *.VISI
*/
RECENT_FFT(VISIB, VISIB, &nx, &ny, &nx);
sprintf(filename,"%sVISI",generic_name);
sprintf(comments,"dc_lissage");
JLP_WRITEIMAG(VISIB, &nx, &ny, &nx, filename, comments);

/* Free memory
*/
free(VISIB);
return(0);
}
/*******************************************************************
* Computes the first approximation image PHIT and the
* corresponding spectrum PCR and PCI
*
* INPUT:
* SNR: snr of the input spectrum
* RFT, IFT: experimental input spectrum
* FTB: bounded tranfer function
* alphat: snr lower threshold to be applied on the input spectrum
* nx, ny: size of arrays
* generic_name: name to be used for the output file
*
* OUTPUT:
* none 
********************************************************************/
static int output_phit_pcr_pci(float *SNR, float *RFT, float *IFT, float *FTB,
                               float *FLI, float alphat, INT4 nx, INT4 ny, 
                               char *generic_name)
{
float *tempr, *tempi, *PCR, *PCI, ratio;
double norm;
register int i;
char filename[60], comments[80];

tempr = (float *)malloc(nx * ny * sizeof(float));
tempi = (float *)malloc(nx * ny * sizeof(float));

/* Computes PHIT, the first approximation of the object (direct space)
* PCR and PCI: real and imaginary parts of the Fourier Transform of PHIT
*/

/* First compute PCR and PCI
*/
for (i = 0; i < nx * ny; i++)
   {
   if ((SNR[i] >= alphat) && (FTB[i] != 0.0))
      {
      ratio = FLI[i]/FTB[i];
      tempr[i] = RFT[i] * ratio;
      tempi[i] = IFT[i] * ratio;
      }
   else
      {
      tempr[i] = 0.0;
      tempi[i] = 0.0;
      }
    }

/* Save  PCR and  PCI to files
*/
PCR = (float *)malloc(nx * ny * sizeof(float));
RECENT_FFT(tempr, PCR, &nx, &ny, &nx);
sprintf(filename,"%sPCR",generic_name);
sprintf(comments,"From dc_lissage, alphat=%.3f", alphat);
JLP_WRITEIMAG(PCR, &nx, &ny, &nx, filename, comments);
free(PCR);

PCI = (float *)malloc(nx * ny * sizeof(float));
RECENT_FFT(tempi, PCI, &nx, &ny, &nx);
sprintf(filename,"%sPCI",generic_name);
sprintf(comments,"From dc_lissage, alphat=%.3f", alphat);
JLP_WRITEIMAG(PCI, &nx, &ny, &nx, filename, comments);
free(PCI);

/* Computes the norm of spectrum of PHIT: */
norm = 0.0;
for (i = 0; i < nx * ny; i++)
   norm += SQUARE(tempr[i]) + SQUARE(tempi[i]);
norm = sqrt(norm/(double)(ny*nx));
printf("norm of (PCR,PCI) (i.e. spectrum of PHIT) : %f\n",norm);

/* Compute PHIT
*Assume zero phase (real image), so I take the modulus of the inverse FFT:
*/
fftw_float(tempr, tempi, (int)nx, (int)ny, -1);

for (i = 0; i < nx * ny; i++)
   tempr[i] = (float)sqrt((double)SQUARE(tempr[i])+(double)SQUARE(tempi[i]));

/* Save PHIT to file:
*/
sprintf(filename,"%sPHIT",generic_name);
sprintf(comments,"From dc_lissage, alphat=%.3f", alphat);
RECENT_FFT(tempr, tempr, &nx, &ny, &nx);
JLP_WRITEIMAG(tempr, &nx, &ny, &nx, filename, comments);

/* Free memory
*/
free(tempr);
free(tempi);

return(0);
}
