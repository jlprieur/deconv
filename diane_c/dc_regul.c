/******************************************************************************
*                                                                   
* PROGRAM: diane_regul                                               
*                                                                     
* PURPOSE: compute interpolation function *.G for a given support      
*          of the object and estimate if the problem is well conditioned
*                                                 
* INPUT:  argv[1] = generic name                   
*         argv[2] = lower thresold for SNR (alphat) 
*         argv[3] = upper thresold for SNR (alphapt) 
*         argv[4] = PSF radius
*                            
* INPUT:
*     *.SNR  signal to noise ratio in fourier domain
*     *.PHIT first approximation of the object
*     *.HR   synthetic aperture
*
* OUTPUT:
*     *.G    regularisation filter
*     *.KR   stabilization filter
*     *.V    support of the object 
*                               
* From Karim's version of December 1992       
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

#ifndef PI
#define PI 3.14159265358979323846
#endif

int  main(int argc, char *argv[])
{
/*    DECLARATIONS
*/
register    int      i, j;
char     filename[61], generic_name[61], mask_name[61], comments[81];
float    *SNR, *PHIT, *HR, *V, *G, *KR;
float    tau, eta, du, nu, su;
float    X, Y, Z, alphat, alphapt;
int      Ni, Nj, err;
INT_PNTR pntr_ima;
INT4     nx, ny, nx1, ny1;

/*    TEST OF COMMAND LINE
*/
/* Five parameters are needed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
if(argc == 7 && *argv[4]) argc = 5;
if(argc != 5) {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndiane_regul  lower_snr_thresold(alphat)  generic_name(with dot) mask_name upper_snr_thresold(alphapt)\n\n");
  fprintf(stderr,"lower_snr_thresold = alphat");
  fprintf(stderr,"\nupper_snr_thresold = alphapt\n\n");
  return(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
err = sscanf(argv[1],"%f",&alphat);
if (err !=1)
  {
  fprintf(stderr,"\nFATAL ERROR: lower thresold for SNR [%s] incorrect\n\n",argv[2]);
  return(-1);
  }
strcpy(generic_name,argv[2]);
strcpy(mask_name,argv[3]);
err = sscanf(argv[4],"%f",&alphapt);
if ((err !=1) || (alphapt < alphat))
  {
  fprintf(stderr,"\nFATAL ERROR: upper thresold for SNR [%s] incorrect\n\n",argv[3]);
  fprintf(stderr,"\nalpha_t = %f should be smaller than alpha_pt=%f\n\n",alphat,alphapt);
  return(-1);
  }

/*    DISPLAY alphat alphapt
*/
printf("\n***************************");
printf("\nPROGRAM :       dc_regul");
printf("\n***************************");
printf("\nalphat (lower thresold of reliability of SNR in fourier domain) = %f",alphat);
printf("\nalphapt (upper thresold of reliability of SNR in fourier domain) = %f",alphapt);

JLP_INQUIFMT();
/*    READ FILE *.PHIT
*/
sprintf(filename,"%sPHIT",generic_name);
JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
PHIT = (float *)pntr_ima;
RECENT_FFT(PHIT,PHIT,&nx,&ny,&nx);

/*    READ mask file (which will become *.V)
*/
JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,mask_name,comments);
V = (float *)pntr_ima;
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",mask_name);
  return(-1);
  }

/*    READ FILE *.SNR
*/
sprintf(filename,"%sSNR",generic_name);
JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,filename,comments);
SNR = (float *)pntr_ima;
if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",filename);
  return(-1);
  }
RECENT_FFT(SNR,SNR,&nx,&ny,&nx);

/*    READ FILE *.HR
*/
sprintf(filename,"%sHR",generic_name);
JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,filename,comments);
HR = (float *)pntr_ima;
if ((nx != nx1) || (ny != ny1))
  {
  printf("\nFATAL ERROR: file [%s] have incompatible size\n\n",filename);
  return(-1);
  }
RECENT_FFT(HR,HR,&nx,&ny,&nx);


/*    COMPUTE A MESURE OF V
*/
tau = 0.0;
for (i = 0; i < nx * ny; i++) tau += V[i];
tau = (float)sqrt((double)tau);
printf("\nmesure of object support = %f",tau);

/*    MEMORY ALLOCATION
*/
KR = (float*) malloc(ny * nx * sizeof(float));
G = (float*) malloc(ny * nx * sizeof(float));

/*    COMPUTING REGULARIZATION FUNCTION
*/
for (i = 0; i < nx * ny; i++) 
   {
   KR[i] = 0.;
   G[i] = 1.0-HR[i];
   if (G[i] == 0.0)
      {
      if (SNR[i] >= alphapt)
         G[i] = 1.0;
      else if (SNR[i] >= alphat)
         {
         G[i] = (SNR[i]-alphat)/(alphapt-alphat);
         G[i] = sqrt(G[i]);
         }
      }
   else if (SNR[i] < alphat)
      KR[i] = 1.0;
   G[i] = 1.0-G[i]*G[i];
   }

Nj = ny-1;
Ni = nx-1;
du = 1.0/(float)ny;
X = G[0]+G[Nj * nx]+G[Ni + Nj * nx]+G[Ni];
Y = 0.0;
for (i = 1; i < Ni; i++)
    Y += (G[i + 1 * nx]+G[i + Nj * nx]);
for (j = 1; j < Nj; j++)
    Y += (G[1 + j * nx]+G[Ni + j * nx]);
Y *= (du*su/2.0);
Z = 0.0;
for (i=1; i<Ni; i++)
for (j=1; j<Nj; j++)
    Z += G[i + j * nx];
Z *= (du*du/2.0);
nu = (float)sqrt((double)(X+Y+Z));
eta = tau*nu;
printf("\nmesure of interpolation function = %f",nu);
printf("\ninterpolation parameter = %f",eta);

/* restore G */
for (i = 0; i < nx*ny; i++)
    G[i] = (float)sqrt(fabs((double)(1.0-G[i])));

/* interpolation parameter */
if (eta < 3.4)
   {
   printf("\nThe problem is well conditioned.");
   printf("\nwe can start iterative reconstruction.\n");
   }
else if (eta < 8.5)
   {
   printf("\nThe problem has medium stability. Reconstruction error could be important.");
   printf("\nwe can start iterative reconstruction or reduce object support or gain in resolution.\n");
   }
else
   {
   printf("\nThe problem is ill posed. The solution unstable.");
   printf("\nWe must reduce object support or gain in resolution\n\n");
   }

/*    FREE MEMORY */
free(HR);
free(PHIT);
free(SNR);

/*    STORE RESULTS   .G  .KR  .V
*/
RECENT_FFT(G,G,&nx,&ny,&nx);
sprintf(filename,"%sG",generic_name);
sprintf(comments,"diane_regul");
JLP_WRITEIMAG(G,&nx,&ny,&nx,filename,comments);
free(G);

RECENT_FFT(KR,KR,&nx,&ny,&nx);
sprintf(filename,"%sKR",generic_name);
sprintf(comments,"diane_regul");
JLP_WRITEIMAG(KR,&nx,&ny,&nx,filename,comments);
free(KR);

sprintf(filename,"%sV",generic_name);
sprintf(comments,"diane_regul");
JLP_WRITEIMAG(V,&nx,&ny,&nx,filename,comments);
free(V);

printf("\n\n");

return(0);
}
