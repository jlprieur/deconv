/******************************************************************************/
/*                                                                            */
/* PROGRAM: diane_erreur                                                      */
/*                                                                            */
/* PURPOSE: compute all errors necessary for error analysis                   */
/*                                                                            */
/* INPUT:  argv[1] = generic name                                             */
/*                                                                            */
/* NB: *.G    regularisation filter                                   (input) */
/*     *.FLI  smoothing function                                      (input) */
/*     *.FTB  bounded transfert function                              (input) */
/*     *.SNR  signal to noise ratio in fourier domain                 (input) */
/*     *.SIG  sigma_i                                                 (input) */
/*     *.KR   stabilization filter                                    (input) */
/*     *.PCR real part of FT of first approximation of the object     (input) */
/*     *.PCI maginary part of FT of first approximation of the object (input) */
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
auto        float    **G, **FLI, **FTB, **SNR, **SIG, **PCR, **PCI, **KR;
auto        float    eps1, eps2, eps3, x, alphat;

/*    TEST OF COMMAND LINE
*/
if (argc != 3)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndiane_erreur  generic_name  alphat\n\n");
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

/*    DISPLAY informations
*/
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nPROGRAM :       diane_erreur");
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nalphat (lower thresold of reliability of SNR in fourier domain) = %f",alphat);

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

/*    READ FILE *.FLI
*/
sprintf(file,"%s.FLI",name);
FLI = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(FLI,nbr_lin,nbr_col);

/*    READ FILE *.FTB
*/
sprintf(file,"%s.FTB",name);
FTB = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(FTB,nbr_lin,nbr_col);


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

/*    READ FILE *.SIG
*/
sprintf(file,"%s.SIG",name);
SIG = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(SIG,nbr_lin,nbr_col);

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
sprintf(file,"%s.PCR",name);
PCI = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(PCI,nbr_lin,nbr_col);


/*    READ FILE *.KR
*/
sprintf(file,"%s.KR",name);
KR = read_data(file,dim,size,type,mode,nature,comments);
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(KR,nbr_lin,nbr_col);

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
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  if ((SNR[i][j] >= alphat) && (FTB[i][j] != 0.0))
    x = G[i][j]*(FLI[i][j]/FTB[i][j])*SIG[i][j];
  else
    x = 0.0;
  eps2 += x*x;
  }
eps2 /= (float)(nbr_lin*nbr_col*nbr_lin*nbr_col);

/*    COMPUTE eps_o
*/
eps3 = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  x = sqrt(PCR[i][j]*PCR[i][j]+PCI[i][j]*PCI[i][j]);
  x *= FLI[i][j]*KR[i][j];
  eps3 += x*x;
  }
eps3 /= (float)(nbr_lin*nbr_col);  

/*    COMPUTE eps
*/
eps1 = eps2+eps3;

/*    FREE MEMORY
free_square_float(G,nbr_lin);
free_square_float(FLI,nbr_lin);
free_square_float(FTB,nbr_lin);
free_square_float(SNR,nbr_lin);
free_square_float(SIG,nbr_lin);
free_square_float(KR,nbr_lin);
free_square_float(PCR,nbr_lin);
free_square_float(PCI,nbr_lin);

/*    PRINT results
*/
fprintf(stderr,"\nepsilon = %f",eps1);
fprintf(stderr,"\nepsilon_i = %f",eps2);
fprintf(stderr,"\nepsilon_o = %f",eps3);
fprintf(stderr,"\n\n\n\n");
}
