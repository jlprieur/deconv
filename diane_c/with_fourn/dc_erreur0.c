/******************************************************************************/
/*                                                                            */
/* PROGRAM: erreur                                                            */
/*                                                                            */
/* PURPOSE: compute all errors necessary for error analysis                   */
/*                                                                            */
/* INPUT:  argv[1] = generic name                                             */
/*                                                                            */
/* NB: *.G    regularisation filter                                   (input) */
/*     *.PCR  real part of FT first approximation of the object       (input) */
/*     *.PCI  imag. part of FT first approximation of the object      (input) */
/*     *.PHIN final solution                                          (input) */
/*                                                                            */
/* VERSION: March 1993                                                        */
/*                                                                            */
/* AUTHOR: JLP, SR, JV                                                        */
/*         translated to C by Karim BOUYOUCEF                                 */
/*                                                                            */
/******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>

extern  long int  read_data();
extern  void      free_square_float();
extern  void      split_image();
extern  void      fourn();
extern  long int  alloc_square_float();

void  main(argc,argv)
int   argc;
char  *argv[];
{
/*    DECLARATIONS
*/
auto        int      dim, size[4], ndim, isign, nn[4];
auto        char     type[8], mode[4], nature[8], comments[1024];
auto        char     name[100], file[100];
register    int      i, j, k;
auto        int      nbr_lin, nbr_col, err;
auto        float    **G, **PCR, **PCI, **PHIN, **PHIN_re, **PHIN_im, *PHIN_fft;
auto        float    epsilon, mu;
auto        double   delta2, norm_phin, error;

/*    TEST OF COMMAND LINE
*/
if (argc != 4)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\nerreur  generic_name  epsilon  mu_min\n\n");
  exit(-1);
  }

/*    READ COMMAND LINE PARAMETERS
*/
strcpy(name,argv[1]);
err = sscanf(argv[2],"%f",&epsilon);
if (err !=1)
  {
  fprintf(stderr,"\nFATAL ERROR: bad line command\n\n");
  exit(-1);
  }

err = sscanf(argv[3],"%f",&mu);
if (err !=1)
  {
  fprintf(stderr,"\nFATAL ERROR: bad line command\n\n");
  exit(-1);
  }


/*    READ FILE *.G
*/
dim = 2;
strcpy(type,"float");
strcpy(nature,"real");
sprintf(file,"%s.G",name);
G = read_data(file,dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];
split_image(G,nbr_lin,nbr_col);

/*    READ FILE *.PCR
*/
sprintf(file,"%s.PCR",name);
PCR = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(PCR,nbr_lin,nbr_col);

/*    READ FILE *.PCI
*/
sprintf(file,"%s.PCI",name);
PCI = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(PCI,nbr_lin,nbr_col);


/*    READ FILE *.PHIN
*/
sprintf(file,"%s.PHIN",name);
PHIN = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }

/*****************************************
     Calcul de l'erreur
******************************************/

/****   calcul de la TF de PHIN  ***/
PHIN_fft = (float*)malloc(2*nbr_lin*nbr_col*sizeof(float));
PHIN_re = (float **) alloc_square_float(nbr_lin,nbr_col);
PHIN_im = (float **) alloc_square_float(nbr_lin,nbr_col);
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   {
   PHIN_fft[2*(i*nbr_col+j)] = PHIN[i][j];
   PHIN_fft[2*(i*nbr_col+j)+1] = 0.0;
   }
ndim = 2;
isign = -1;
nn[0] = nbr_lin;
nn[1] = nbr_col;
fourn (PHIN_fft-1,nn-1,ndim,isign);

for (k=0; k<nbr_lin*nbr_col; k++)
  {
  i = k/nbr_col;
  j = k-i*nbr_col;
  PHIN_re[i][j] = PHIN_fft[2*k];
  PHIN_im[i][j] = PHIN_fft[2*k+1];
  }

delta2 = 0.0;
norm_phin = 0.0; 
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
    {
    PCR[i][j] -= PHIN_re[i][j];
    PCR[i][j] *= G[i][j];
    PCI[i][j] -= PHIN_im[i][j];
    PCI[i][j] *= G[i][j];
    delta2 += (PCR[i][j] * PCR[i][j] + PCI[i][j] * PCI[i][j]);
    norm_phin += PHIN[i][j] * PHIN[i][j];
    }
delta2 /= (float)(nbr_lin*nbr_col*nbr_lin*nbr_col);
error = sqrt(fabs(epsilon - delta2)/norm_phin);
/*
error *= sqrt(.5*(1.+1./mu));
*/
error *= sqrt(1./mu);

/*    PRINT results
*/
fprintf(stderr,"\nepsilon = %f",epsilon);
fprintf(stderr,"\nmu = %f",mu);
fprintf(stderr,"\ndelta2 = %lf",delta2);
fprintf(stderr,"\nnorm_phin **2 = %lf",norm_phin);
fprintf(stderr,"\nerror = %lf",error);
fprintf(stderr,"\n\n\n\n");
}
