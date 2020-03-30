/******************************************************************************/
/*                                                                            */
/* PROGRAM: diane_regul                                                       */
/*                                                                            */
/* PURPOSE: compute interpolation function *.G for a given support            */
/*          of the object and estimate if the problem is well conditioned     */
/*                                                                            */
/* INPUT:  argv[1] = generic name                                             */
/*         argv[2] = lower thresold for SNR (alphat)                          */
/*         argv[3] = upper thresold for SNR (alphapt)                         */
/*         argv[4] = PSF radius                                               */
/*                                                                            */
/* NB: *.SNR  signal to noise ratio in fourier domain                 (input) */
/*     *.PHIT first approximation of the object                       (input) */
/*     *.HR   synthetic aperture                                      (input) */
/*     *.V    support of the object                                   (input) */
/*     *.G    regularisation filter                                  (output) */
/*     *.KR   stabilization filter                                   (output) */
/*                                                                            */
/* VERSION: December 1992                                                     */
/*                                                                            */
/* AUTHOR: JLP, SR, JV                                                        */
/*         translated to C by Karim BOUYOUCEF                                 */
/*                                                                            */
/******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>

#define pi 3.14159265358979323846

extern  long int  read_data();
extern  void      write_data();
extern  float **  alloc_square_float();
extern  void      free_square_float();
extern  void      split_image();

void  main(argc,argv)
int   argc;
char  *argv[];
{
/*    DECLARATIONS
*/
auto        int      dim, size[4];
auto        char     type[8], mode[4], nature[8], comments[1024];
auto        char     name[100], file[100];
register    int      i, j;
auto        int      nbr_lin, nbr_col, err;
auto        float    **SNR, **PHIT, **HR, **V, **G, **KR;
auto        float    tau, eta, du, nu, su;
auto        float    X, Y, Z, alphat, alphapt;
auto        int      Ni, Nj;

/*    TEST OF COMMAND LINE
*/
if (argc != 4)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndiane_regul  generic_name  lower_snr_thresold(alphat)  upper_snr_thresold(alphapt)\n\n");
  fprintf(stderr,"lower_snr_thresold = alphat");
  fprintf(stderr,"\nupper_snr_thresold = alphapt\n\n");
  exit(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
strcpy(name,argv[1]);
err = sscanf(argv[2],"%f",&alphat);
if (err !=1)
  {
  fprintf(stderr,"\nFATAL ERROR: lower thresold for SNR [%f] incorrect\n\n",argv[2]);
  exit(-1);
  }
err = sscanf(argv[3],"%f",&alphapt);
if ((err !=1) || (alphapt < alphat))
  {
  fprintf(stderr,"\nFATAL ERROR: upper thresold for SNR [%f] incorrect\n\n",argv[3]);
  exit(-1);
  }

/*    DISPLAY alphat alphapt
*/
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nPROGRAM :       diane_regul");
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nalphat (lower thresold of reliability of SNR in fourier domain) = %f",alphat);
fprintf(stderr,"\nalphapt (upper thresold of reliability of SNR in fourier domain) = %f",alphapt);

/*    READ FILE *.PHIT
*/
dim = 2;
strcpy(type,"float");
strcpy(nature,"real");
sprintf(file,"%s.PHIT",name);
PHIT = read_data(file,dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];
split_image(PHIT,nbr_lin,nbr_col);

/*    READ FILE *.V
*/
sprintf(file,"%s.V",name);
V = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }

/*    READ FILE *.SNR
*/
sprintf(file,"%s.SNR",name);
SNR = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(SNR,nbr_lin,nbr_col);

/*    READ FILE *.HR
*/
sprintf(file,"%s.HR",name);
HR = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(HR,nbr_lin,nbr_col);


/*    COMPUTE A MESURE OF V
*/
tau = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
    tau += V[i][j];
tau = (float)sqrt((double)tau);
fprintf(stderr,"\nmesure of object support = %f",tau);
free_square_float(V,nbr_lin);

/*    MEMORY ALLOCATION
*/
KR = (float**) alloc_square_float(nbr_lin,nbr_col);
G = (float**) alloc_square_float(nbr_lin,nbr_col);

/*    COMPUTING REGULARIZATION FUNCTION
*/
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   {
   G[i][j] = 1.0-HR[i][j];
   if (G[i][j] == 0.0)
      {
      if (SNR[i][j] >= alphapt)
         G[i][j] = 1.0;
      else if (SNR[i][j] >= alphat)
         {
         G[i][j] = (SNR[i][j]-alphat)/(alphapt-alphat);
         G[i][j] = sqrt(G[i][j]);
         }
      }
   else if (SNR[i][j] < alphat)
      KR[i][j] = 1.0;
   G[i][j] = 1.0-G[i][j]*G[i][j];
   }

Ni = nbr_lin-1;
Nj = nbr_col-1;
du = 1.0/(float)nbr_lin;
X = G[0][0]+G[0][Nj]+G[Ni][Nj]+G[Ni][0];
Y = 0.0;
for (i=1; i<Ni; i++)
    Y += (G[i][1]+G[i][Nj]);
for (j=1; j<Nj; j++)
    Y += (G[1][j]+G[Ni][j]);
Y *= (du*su/2.0);
Z = 0.0;
for (i=1; i<Ni; i++)
for (j=1; j<Nj; j++)
    Z += G[i][j];
Z *= (du*du/2.0);
nu = (float)sqrt((double)(X+Y+Z));
eta = tau*nu;
fprintf(stderr,"\nmesure of interpolation function = %f",nu);
fprintf(stderr,"\ninterpolation parameter = %f",eta);

/* restitute G */
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
    G[i][j] = (float)sqrt(fabs((double)(1.0-G[i][j])));

/* interpolation parameter */
if (eta < 3.4)
   {
   fprintf(stderr,"\nThe problem is well conditioned.");
   fprintf(stderr,"\nwe can start iterative reconstruction.\n");
   }
else if (eta < 8.5)
   {
   fprintf(stderr,"\nThe problem has medium stability. Reconstruction error could be important.");
   fprintf(stderr,"\nwe can start iterative reconstruction or reduce object support or gain in resolution.\n");
   }
else
   {
   fprintf(stderr,"\nThe problem is ill posed. The solution unstable.");
   fprintf(stderr,"\nWe must reduce object support or gain in resolution\n\n");
   }

/*    FREE MEMORY
free_square_float(HR,nbr_lin);
free_square_float(PHIT,nbr_lin);
free_square_float(SNR,nbr_lin);

/*    STORE RESULTS   .G  .KR
*/
split_image(G,nbr_lin,nbr_col);
sprintf(file,"%s.G",name);
sprintf(comments,"diane_regul");
write_data(file,G,dim,size,type,mode,nature,comments);
free_square_float(G,nbr_lin);

split_image(KR,nbr_lin,nbr_col);
sprintf(file,"%s.KR",name);
sprintf(comments,"diane_regul");
write_data(file,KR,dim,size,type,mode,nature,comments);
free_square_float(KR,nbr_lin);
fprintf(stderr,"\n\n");

}
