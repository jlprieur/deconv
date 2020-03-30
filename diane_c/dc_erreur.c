/******************************************************************************
*                                                                    
* PROGRAM: dc_erreur (diane_erreur)
*                                                                      
* PURPOSE: compute all errors necessary for error analysis              
*                                                                        
* INPUT:  argv[1] = generic name                                          
*         argv[2] *.FTB  bounded transfert function 
*         argv[3] = alpha_T
*                                                                          
* INPUT:
*     *.G    regularisation filter  
*     *.FLI  smoothing function    
*     *.SNR  signal to noise ratio in fourier domain
*     *.SIG  sigma_i                               
*     *.KR   stabilization filter                 
*     *.PCR real part of FT of first approximation of the object
*     *.PCI maginary part of FT of first approximation of the object 
*                                                                       
* From Karim's version of March 1993 
*                                                                         
* AUTHOR: JLP, SR, JV                                                      
*         translated to C by Karim BOUYOUCEF                                
*                                                                            
* JLP 
* Version 10/03/2008
******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

int  main(int argc, char *argv[])
{
/*    DECLARATIONS
*/
char     generic_name[61], FTB_name[61], filename[61], comments[81];
INT4     nx, ny, nx1, ny1;
INT_PNTR pntr_ima;
int      err;
float    *G, *FLI, *FTB, *SNR, *SIG, *PCR, *PCI, *KR;
float    eps1, eps2, eps3, x, alphat;
register int i;

/*    TEST OF COMMAND LINE
*/
/* Carefull: 7 parameters always, using JLP "runs" */
if(argc == 7 && *argv[3]) argc = 4;
if(argc != 4) {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndc_erreur  generic_name  transfer_function(*.FTB) alphat\n\n");
  return(-1);
  }

/*    READ COMMAND LINE PARAMETERS
*/
strcpy(generic_name,argv[1]);
strcpy(FTB_name,argv[2]);
err = sscanf(argv[3],"%f",&alphat);
if (err !=1)
  {
  fprintf(stderr,"\nFATAL ERROR: lower thresold for SNR [%s] incorrect\n\n",
          argv[3]);
  return(-1);
  }

/*    DISPLAY informations
*/
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nPROGRAM :       dc_erreur");
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nalphat (lower thresold of reliability of SNR in fourier domain) = %f \n",alphat);

JLP_INQUIFMT();
/*    READ file *.G
*/
sprintf(filename,"%sG",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, filename, comments);
G = (float *)pntr_ima; 
RECENT_FFT(G,G,&nx,&ny,&nx);

/*    READ file *.FLI
*/
sprintf(filename,"%sFLI",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
FLI = (float *)pntr_ima; 
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",filename);
  return(-1);
  }
RECENT_FFT(FLI,FLI,&nx,&ny,&nx);

/*    READ file *.FTB
*/
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, FTB_name, comments);
FTB = (float *)pntr_ima; 
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",FTB_name);
  return(-1);
  }
RECENT_FFT(FTB,FTB,&nx,&ny,&nx);


/*    READ file *.SNR
*/
sprintf(filename,"%sSNR",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
SNR = (float *)pntr_ima; 
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",filename);
  return(-1);
  }
RECENT_FFT(SNR,SNR,&nx,&ny,&nx);

/*    READ file *.SIG
*/
sprintf(filename,"%sSIG",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
SIG = (float *)pntr_ima; 
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",filename);
  return(-1);
  }
RECENT_FFT(SIG,SIG,&nx,&ny,&nx);

/*    READ file *.PCR
*/
sprintf(filename,"%sPCR",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
PCR = (float *)pntr_ima; 
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",filename);
  return(-1);
  }
RECENT_FFT(PCR,PCR,&nx,&ny,&nx);

/*    READ file *.PCI
*/
sprintf(filename,"%sPCI",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
PCI = (float *)pntr_ima; 
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",filename);
  return(-1);
  }
RECENT_FFT(PCI,PCI,&nx,&ny,&nx);

/*    READ file *.KR
*/
sprintf(filename,"%sKR",generic_name);
JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, filename, comments);
KR = (float *)pntr_ima; 
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",filename);
  return(-1);
  }
RECENT_FFT(KR,KR,&nx,&ny,&nx);

/*****************************************
     Calcul des erreurs
  
    eps1 = (eps2 ** 2  +  eps3 ** 2)
    eps2 = eps_i = ||G*St*sigma_i||
    eps3 = eps_o = ||Kr*S*sigma_o||

   sigma_o = | PHITc |
******************************************/
/*    COMPUTE eps_i
*/
eps2 = 0.0;
for (i = 0; i < nx * ny; i++)
  {
  if ((SNR[i] >= alphat) && (FTB[i] != 0.0))
    x = G[i]*(FLI[i]/FTB[i])*SIG[i];
  else
    x = 0.0;
  eps2 += x*x;
  }
eps2 /= (float)(ny*nx*ny*nx);

/*    COMPUTE eps_o
*/
eps3 = 0.0;
for (i = 0; i < nx * ny; i++)
  {
  x = sqrt(PCR[i]*PCR[i]+PCI[i]*PCI[i]);
  x *= FLI[i]*KR[i];
  eps3 += x*x;
  }
eps3 /= (float)(ny*nx);  

/*    COMPUTE eps
*/
eps1 = eps2+eps3;

/*    FREE MEMORY */
free(G);
free(FLI);
free(FTB);
free(SNR);
free(SIG);
free(KR);
free(PCR);
free(PCI);

/*    PRINT results
*/
fprintf(stderr,"\nepsilon = %f",eps1);
fprintf(stderr,"\nepsilon_i = %f",eps2);
fprintf(stderr,"\nepsilon_o = %f",eps3);
fprintf(stderr,"\n\n\n\n");

return(0);
}
