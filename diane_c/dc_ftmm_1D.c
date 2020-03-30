/***************************************************************
*                                                     
* PROGRAM: dc_ftmm_1D (diane_ftm)                                   
* 1D version of dc_ftmm
*                                                       
* PURPOSE: normalize Transfer Function with L1 normalization     
*          and bound it 
*                                                          
* INPUT:  argv[1] = modsq_name                            
*         argv[2] = FTB_name                            
*         argv[3] = Lower thresold for Transfer Function in Fourier domain 
*                                                             
* NB:     *.FTB bounded of Transfer Function (output) 
*
* From Karim's version: December 1992                                      
*                                                             
* AUTHOR: JLP, SR, JV                                         
*         translated to C by Karim BOUYOUCEF                  
*         adapted to 1-D by JLP
*                                                             
* JLP 
* Version 01/02/00
***************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

void  main(argc,argv)
int   argc;
char  *argv[];
{
/*    DECLARATIONS
*/
INT4 nx, ny;
INT_PNTR pntr_ima;
int      status;
char     comments[81], ftb_name[60], modsq_name[61];
register int i, j, ix, iy;
float    *modsq, min, max;
double   *ftr, *fti;
double   alpha, cumul, ww;

 printf(" dc_ftmm_1D, version 01/02/00 \n");
 printf("NB: To determine the correct threshold, first use this program with 0.\n");
/*    TEST OF COMMAND LINE
*/
if (argc == 7) argc = 4;
if (argc != 4)
  {
  fprintf(stderr,"\nUnexpected number of arguments \n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndc_ftm  modsq_file output_FTB  Lower_thresold\n\n");
  fprintf(stderr,"\nLower_thresold = lower thresold of Transfer Function\n\n");
  exit(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
strcpy(modsq_name,argv[1]);
strcpy(ftb_name,argv[2]);
status = sscanf(argv[3],"%lf",&alpha);
if ((status != 1) || (alpha > 1.))
  {
  fprintf(stderr,"\nFATAL ERROR: Lower thresold [%s] incorrect\n\n",argv[2]);
  exit(-1);
  }

fprintf(stderr,"\n***************************");
fprintf(stderr,"\nPROGRAM :   dc_ftmm_1D  Version 01/02/00");
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nlower thresold for Transfer Function = %f \n",alpha);

/*    INPUT MODSQ 
*/
JLP_BEGIN();
JLP_INQUIFMT();
JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, modsq_name, comments);
modsq = (float *)pntr_ima;

fti = (double *)malloc(nx * ny * sizeof(double));
ftr = (double *)malloc(nx * ny * sizeof(double));

/* Square root of modsq: */
for (i = 0; i < nx * ny; i++) 
   {
   ww = modsq[i];
   if(ww > 0)
#ifdef TOTO
      ftr[i] = sqrt(ww);
#endif
      ftr[i] = ww;
   else
      ftr[i] = 0.;
   }

/* First compute the PSF: */
RECENT_FFT_1D_Y(ftr, ftr, &nx, &ny, &nx);
for (i = 0; i < nx * ny; i++) fti[i] = 0.;
fftw_1D_Y(ftr, fti, (int)nx, (int)ny, -1);

/* Make sure that the PSF is zero on the edges: */
for(ix = 0; ix < nx; ix++)
{
  min = ftr[ix];
  for (iy = 1; iy < ny; iy++) 
     if(min > ftr[ix + iy*nx]) min = ftr[ix + iy*nx];
  for (iy = 0; iy < ny; iy++) ftr[ix + iy*nx] -= min;
}

/* Then normalize L1 the PSF for each column */
for(ix = 0; ix < nx; ix++)
{
  cumul = 0.0;
  for (iy = 0; iy < ny; iy++) cumul += ftr[ix + iy*nx];
#ifdef DEBUG
  if(ix < 1) printf(" cumul = %f (line#%d)\n", cumul,ix);
#endif
  for (iy = 0; iy < ny; iy++) ftr[ix + iy*nx] /= cumul;
}

/* Computes back the Transfer Function: */
for (i = 0; i < nx * ny; i++) fti[i] = 0.;
fftw_1D_Y(ftr, fti, (int)nx, (int)ny, 1);

/* Threshold according to ALPHAT:
*/
max = ftr[0];
for (i = 0; i < nx * ny; i++) 
   if(max < ftr[i]) max = ftr[i];
/* Just for an estimation (for a first go) */
printf("Maximum of Transfer Function is = %f\n", max);

for (i = 0; i < nx * ny; i++) 
    if (ftr[i] < alpha) ftr[i] = 0.0;
if(alpha >= max) 
  {
  printf("Fatal error/Threshold is too high: max < alpha\n");
  exit(-1);
  }

/*    STORAGE OF THE RESULT   *.FTB
*/
RECENT_FFT_1D_Y(ftr, ftr, &nx, &ny, &nx);
strcpy(comments,"dc_ftmm_1D");
JLP_D_WRITEIMAG(ftr, &nx, &ny, &nx, ftb_name, comments);
fprintf(stderr,"\n\n");

JLP_END();
}
