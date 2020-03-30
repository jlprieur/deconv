/***************************************************************
*                                                     
* PROGRAM: dc_ftmm (diane_ftm)                                   
* Same as dc_ftm but takes square modulus in input instead of the PSF.
*                                                       
* PURPOSE: normalize Transfer Function with L1 normalization     
*          and bound it 
*                                                          
* INPUT:  argv[1] = modsq_name                            
*         argv[2] = FTB_name                            
*         argv[3] = Lower thresold for Transfer Function in fourier domain 
*                                                             
* NB:     *.FTB bounded imodulus FT of .PSF (output) 
*                                                             
* From Karim's version: December 1992                                      
*                                                             
* AUTHOR: JLP, SR, JV                                         
*         translated to C by Karim BOUYOUCEF                  
*                                                             
* JLP 
* Version 10/03/08
***************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

int  main(int argc, char *argv[])
{
/*    DECLARATIONS
*/
INT4 nx, ny;
INT_PNTR pntr_ima;
int      status;
char     comments[81], ftb_name[60], modsq_name[61];
float    *modsq, max;
double   alpha, cumul, ww;
register int  i;

 printf(" dc_ftmm, version 10/03/08 \n");
 printf("NB: To determine the correct threshold, first use this program with 0.\n");
/*    TEST OF COMMAND LINE
*/
if (argc == 7) argc = 4;
if (argc != 4)
  {
  fprintf(stderr,"\nUnexpected number of arguments \n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndc_ftmm  modsq_file output_FTB  Lower_thresold\n\n");
  fprintf(stderr,"\nLower_thresold = lower thresold of Transfer Function \n\n");
  return(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
strcpy(modsq_name,argv[1]);
strcpy(ftb_name,argv[2]);
status = sscanf(argv[3],"%lf",&alpha);
if ((status != 1) || (alpha > 1.))
  {
  fprintf(stderr,"\nFATAL ERROR: Lower thresold [%s] incorrect\n\n",argv[2]);
  return(-1);
  }

fprintf(stderr,"\n***************************");
fprintf(stderr,"\nPROGRAM :         dc_ftm");
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nlower thresold for Transfer Function = %f \n",alpha);

/*    INPUT MODSQ 
*/
JLP_INQUIFMT();
JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, modsq_name, comments);
modsq = (float *)pntr_ima;

/* Square root of modsq: */
for (i = 0; i < nx * ny; i++) 
   {
   ww = modsq[i];
   if(ww > 0)
      modsq[i] = sqrt(ww);
   else
      modsq[i] = 0.;
   }

/* Normalize (L1) the Transfer Function
*/
cumul = modsq[(nx/2) + (ny/2) * nx];
for (i = 0; i < nx * ny; i++) modsq[i] /= cumul;

/* Compute the maximum of the image: */
max = modsq[0];
for (i = 0; i < nx * ny; i++) max = MAXI(max, modsq[i]);

printf("Maximum of Transfer Function is %f\n", max);

if(alpha >= max) 
  {
  printf("Fatal error/Threshold is too high: max < alpha\n");
  return(-1);
  }

/* Apply the threshold
*/
for (i = 0; i < nx * ny; i++) 
    if (modsq[i] < alpha) modsq[i] = 0.0;

printf(" JLP99 \n");
/*    STORAGE OF THE RESULT   *.FTB
*/
strcpy(comments,"dc_ftmm");
JLP_WRITEIMAG(modsq, &nx, &ny, &nx, ftb_name, comments);
fprintf(stderr,"\n\n");

return(0);
}
