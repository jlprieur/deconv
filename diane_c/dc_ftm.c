/***************************************************************
*                                                     
* PROGRAM: dc_ftm (diane_ftm)                                   
*                                                       
* PURPOSE: normalize PSF respecting L1 normalization     
*          performs fourier transform of PSF and bound it 
*                                                          
* INPUT:  argv[1] = PSF_name                            
*         argv[2] = FTB_name                            
*         argv[3] = Lower thresold for PSF in Fourier domain 
*                                                             
* NB:   generic_name.PSF point spread fuction (input)         
*       generic_name.FTB bounded imodulus FT of .PSF (output) 
*                                                             
* From Karim's version: December 1992                                      
*                                                             
* AUTHOR: JLP, SR, JV                                         
*         translated to C by Karim BOUYOUCEF                  
*                                                             
* JLP 
* Version 03/05/99
***************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

int main(int argc, char *argv[])
{
/*    DECLARATIONS
*/
INT4 nx, ny;
INT_PNTR pntr_ima;
int      status;
char     comments[81], ftb_name[60], psf_name[61];
float    *psf, max;
double   *tfr, *tfi;
double   alpha, cumul;
register int i;

 printf("NB: To determine the correct threshold, first use this program with 0.\n");
/*    TEST OF COMMAND LINE
*/
if (argc == 7) argc = 4;
if (argc != 4)
  {
  printf("\nUnexpected number of arguments \n");
  printf("\nUSAGE:\n");
  printf("\ndc_ftm  psf_file output_FTB  Lower_thresold\n\n");
  printf("\nLower_thresold = lower thresold of modulus fourier transform of *.PSF\n\n");
  exit(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
strcpy(psf_name,argv[1]);
strcpy(ftb_name,argv[2]);
status = sscanf(argv[3],"%lf",&alpha);
if ((status != 1) || (alpha > 1.))
  {
  printf("\nFATAL ERROR: Lower thresold [%s] incorrect\n\n",argv[2]);
  exit(-1);
  }

printf("\n***************************");
printf("\nPROGRAM :         dc_ftm");
printf("\n***************************");
printf("\nlower thresold for FT of PSF = %f \n",alpha);

/*    INPUT PSF
*/
JLP_INQUIFMT();
JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, psf_name, comments);
psf = (float *)pntr_ima;

/*    NORMALIZE L1 THE PSF
*/
cumul = 0.0;
for (i = 0; i < nx * ny; i++) cumul += psf[i];
printf(" cumul = %f\n", cumul);
for (i = 0; i < nx * ny; i++) psf[i] /= cumul;

tfr = (double *)malloc(ny * nx * sizeof(double));
for (i = 0; i < nx * ny; i++) tfr[i] = psf[i];
tfi = (double *)malloc(ny * nx * sizeof(double));
for (i = 0; i < nx * ny; i++) tfi[i] = 0.;
fftw_double(tfr,tfi,(int)nx,(int)ny,1);

/*    MODULUS OF PSF FOURIER TRANSFORM
*/
for (i = 0; i < nx * ny; i++) 
  tfr[i] = sqrt(tfr[i]*tfr[i]+tfi[i]*tfi[i]);

/*  Apply the threshold
*/
max = tfr[0];
for (i = 0; i < nx * ny; i++) max = MAXI(max, tfr[i]);

printf("Maximum of transfer function = %f\n", max);

if(alpha >= max) 
  {
  printf("Fatal error/Threshold is too high: max < alpha\n");
  exit(-1);
  }

for (i = 0; i < nx * ny; i++) 
    if (tfr[i] < alpha) tfr[i] = 0.0;

RECENT_FFT_DOUBLE(tfr, tfr, &ny, &nx, &nx);

/*    STORAGE OF THE RESULT   *.FTB
*/
strcpy(comments,"dc_ftm");
JLP_D_WRITEIMAG(tfr, &nx, &ny, &nx, ftb_name, comments);
printf("\n\n");

return(0);
}
