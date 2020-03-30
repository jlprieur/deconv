/******************************************************************************/
/*                                                                            */
/* PROGRAM: diane_lissage                                                     */
/*                                                                            */
/* PURPOSE: compute smoothing function and raw approximation                  */
/*          of the object                                                     */
/*                                                                            */
/* INPUT:  argv[1] = generic name                                             */
/*         argv[2] = lower thresold for SNR (alphat)                          */
/*         argv[3] = PSF radius                                               */
/*                                                                            */
/* NB: *.SNR  signal to noise ratio in fourier domain                 (input) */
/*     *.FTB  bounded imodulus FT of .PSF                             (input) */
/*     *.TFR  real part of fourier transform of input image           (input) */
/*     *.TFI  imaginary part of fourier transform of input image      (input) */
/*     *.HR   synthetic aperture                                     (output) */
/*     *.FLI  smoothing filter                                       (output) */
/*     *.PCR  FT real part of first approximation of the object      (output) */
/*     *.PCR  FT imaginary part of first approximation of the object (output) */
/*     *.PHIT first approximation of the object                      (output) */
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
extern  int       fourierC();
extern  int       read_header();

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
auto        int      nbr_lin, nbr_col, err, flag;
auto        float    **SNR, **FTB, **TFR, **TFI, **sectr, **elipr, **stock, **tempr, **tempi;
auto        float    alphat, radius, cumul;
auto        float    somme, error, eta0, epsilon, X2, Y2;
auto        float    ratio, norm, norme, rr, rp, test_real, p_axe, gr_axe;
auto        int      itere, c2, d2, is, it_max;
auto        int      np, irp, r2p, ic;

/*    TEST OF COMMAND LINE
*/
if (argc != 4)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndiane_lissage  generic_name  lower_snr_thresold(alphat)  PSF_radius\n\n");
  fprintf(stderr,"lower_snr_thresold = alphat");
  fprintf(stderr,"\nPSF_radius = radius for synthetic aperture *.HR\n\n");
  exit(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
strcpy(name,argv[1]);
err = sscanf(argv[2],"%f",&alphat);
if (err !=1)
  {
  fprintf(stderr,"\nFATAL ERROR: Lower thresold [%f] incorrect\n\n",argv[2]);
  exit(-1);
  }

/*    READ FILE SIZE
*/
sprintf(file,"%s.FTB",name);
read_header(file,&dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];
if ((nbr_lin < 1) || (nbr_lin > 1024) || (nbr_col < 1) || (nbr_col > 1024))
   {
   fprintf(stderr,"\nFATAL ERROR: invalid dimensions [%d,%d] from %s\n\n",nbr_lin,nbr_col,file);
   exit(-1);
   }

err = sscanf(argv[3],"%f",&radius);
if ((err !=1) || (radius < 1.) || (radius > (float)(nbr_lin/2-1)) || (radius > (float)(nbr_col/2-1)))
  {
  fprintf(stderr,"\nFATAL ERROR: radius for synthetic aperture [%f] incorrect\n\n",argv[3]);
  exit(-1);
  }

/*   DISPLAY  alphat radius
*/
fprintf(stderr,"\n**************************");
fprintf(stderr,"\nPROGRAM :    diane_lissage");
fprintf(stderr,"\n**************************");
fprintf(stderr,"\nalphat = %f",alphat);
fprintf(stderr,"\nradius of Hr (fourier space) = %f",radius);

/*    MEMORY ALLOCATION
*/
tempr = (float**) alloc_square_float(nbr_lin,nbr_col);
sectr = (float**)alloc_square_float(nbr_lin,nbr_col);

/*    COMPUTING SMOOTHING FUNCTION
*/
eta0 = 1.55;
it_max = 50;
epsilon = 0.00001;
np = nbr_lin/2+1;
rr = eta0*nbr_lin/(radius*pi);
fprintf(stderr,"\nradius of smoothing function (direct space) = %f",rr);
rp = radius;
irp = (int)rp+1;
r2p = irp*irp;
is = 0;
for (j=0; j<nbr_col; j++)
   {
   c2 = (j-np)*(j-np);
   for (i=0; i<nbr_lin; i++)
     {
     if ((c2+(i-np)*(i-np)) <= r2p)
       {
       sectr[i][j] = 1.0;
       is++;
       }
     }
   }
split_image(sectr,nbr_lin,nbr_col);

/* ellipsoide of resolution */
p_axe = rr;
gr_axe = p_axe;
ic = 0;
p_axe *= p_axe;
gr_axe *= gr_axe;
elipr = (float**)alloc_square_float(nbr_lin,nbr_col);
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   {
   c2 = (j-np)*(j-np);
   d2 = (i-np)*(i-np);
   test_real = (float)c2/p_axe+(float)d2/gr_axe;
   if (test_real <= 1.0)
      {
      elipr[i][j] = 1.0;
      ic++;
      }
   tempr[i][j] = elipr[i][j];
   }

eta0 = sqrt((float)(ic*is)/(float)nbr_lin);

/* compute smoothing function by power method */
stock = (float **) alloc_square_float(nbr_lin,nbr_col);
tempi = (float**)alloc_square_float(nbr_lin,nbr_col);
itere = 0;
error = epsilon + 1.0;

do
   {
   for (i=0; i<nbr_lin; i++)
   for (j=0; j<nbr_col; j++)
     {
     stock[i][j] = tempr[i][j];
     tempi[i][j] = 0.0;
     }
   flag = 1;
   err = fourierC(tempr,tempi,nbr_lin,nbr_col,flag);
   if (err != 0)
     {
     fprintf(stderr,"\nFATAL ERROR: problem with complex fourier transform\n\n");
     exit(-1);
     }
   for (i=0; i<nbr_lin; i++)
   for (j=0; j<nbr_col; j++)
     if (sectr[i][j] == 0.0)
       {
       tempr[i][j] = 0.0;
       tempi[i][j] = 0.0;
       }
   flag = -1;
   err = fourierC(tempr,tempi,nbr_lin,nbr_col,flag);
   if (err != 0)
     {
     fprintf(stderr,"\nFATAL ERROR: problem with complex fourier transform\n\n");
     exit(-1);
     }
   for (i=0; i<nbr_lin; i++)
   for (j=0; j<nbr_col; j++)
     if (elipr[i][j] == 0.0)
       {
       tempr[i][j] = 0.0;
       tempi[i][j] = 0.0;
       }
   /* normalization and convergence test */
   norme = 0.0;
   error = 0.0;
   for (i=0; i<nbr_lin; i++)
   for (j=0; j<nbr_col; j++)
     {
     X2 = tempr[i][j] * tempr[i][j];
     Y2 = tempi[i][j] * tempi[i][j];
     norme += (X2+Y2);
     tempr[i][j] = (float)sqrt(X2+Y2);
     }
   norme = (float)sqrt((double)norme);
   for (i=0; i<nbr_lin; i++)
   for (j=0; j<nbr_col; j++)
     {
     tempr[i][j] /= norme;
     X2 = (tempr[i][j] - stock[i][j]) * (tempr[i][j] - stock[i][j]);
     error = error + X2;
     }
   error = (float)sqrt((double)error);
   itere++;
   }
while ((itere < it_max) && (error > epsilon));

fprintf(stderr,"\nnumber of iterations = %d",itere);
fprintf(stderr,"\nError = %f",error);
fprintf(stderr,"\nKhi = %f",norme);

/*   FREE MEMORY
*/
free_square_float(elipr,nbr_lin);
free_square_float(stock,nbr_lin);

/* computes smoothing function  */
somme = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   {
   tempi[i][j] = 0.0;
   somme += tempr[i][j];
   }
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   tempr[i][j] /= somme;

flag = 1;
err = fourierC(tempr,tempi,nbr_lin,nbr_col,flag);
if (err != 0)
  {
  fprintf(stderr,"\nFATAL ERROR: problem with complex fourier transform\n\n");
  exit(-1);
  }
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  tempr[i][j] = (float)sqrt((double)(tempr[i][j]*tempr[i][j]+tempi[i][j]*tempi[i][j]));  

/*   FREE MEMORY
*/
free_square_float(tempi,nbr_lin);

/*   STORE RESULT  *.HR
*/
split_image(sectr,nbr_lin,nbr_col);
dim = 2;
size[0] = nbr_lin;
size[1] = nbr_col;
strcpy(type,"float");
strcpy(mode,"bin");
strcpy(nature,"real");
sprintf(file,"%s.HR",name);
sprintf(comments,"diane_lissage");
write_data(file,sectr,dim,size,type,mode,nature,comments);
free_square_float(sectr,nbr_lin);

/*    COMPUTES PHIT
      first approximation of the object
*/
sprintf(file,"%s.SNR",name);
SNR = (float **)read_data(file,dim,size,type,mode,nature,comments);
if (SNR == NULL)
  {
  fprintf(stderr,"FATAL ERROR: when loading file [%s]\n\n",file);
  exit(-1);
  }
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  } 
split_image(SNR,nbr_lin,nbr_col);

sprintf(file,"%s.TFR",name);
TFR = (float **)read_data(file,dim,size,type,mode,nature,comments);
if (TFR == NULL)
  {
  fprintf(stderr,"FATAL ERROR: when loading file [%s]\n\n",file);
  exit(-1);
  }
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(TFR,nbr_lin,nbr_col);

sprintf(file,"%s.TFI",name);
TFI = (float **)read_data(file,dim,size,type,mode,nature,comments);
if (TFI == NULL)
  {
  fprintf(stderr,"FATAL ERROR: when loading file [%s]\n\n",file);
  exit(-1);
  }
if ((size[0] != nbr_lin) || (size[1] != nbr_col))
  {
  fprintf(stderr,"\nFATAL ERROR: file [%s] have incompatible size\n\n",file);
  exit(-1);
  }
split_image(TFI,nbr_lin,nbr_col);

sprintf(file,"%s.FTB",name);
FTB = (float **)read_data(file,dim,size,type,mode,nature,comments);
if (FTB == NULL)
  {
  fprintf(stderr,"FATAL ERROR: when loading file [%s]\n\n",file);
  exit(-1);
  }
split_image(FTB,nbr_lin,nbr_col);


/*    COMPUTE  .PCR  .PCI
*/
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   {
   if ((SNR[i][j] >= alphat) && (FTB[i][j] != 0.0))
      {
      ratio = tempr[i][j]/FTB[i][j];
      TFR[i][j] *= ratio;
      TFI[i][j] *= ratio;
      }
   else
      {
      TFR[i][j] = 0.0;
      TFI[i][j] = 0.0;
      }
    }

/*   STORE RESULTS  .FLI  .PCR  .PCI
*/
split_image(tempr,nbr_lin,nbr_col);
sprintf(file,"%s.FLI",name);
sprintf(comments,"diane_lissage");
write_data(file,tempr,dim,size,type,mode,nature,comments);

split_image(TFR,nbr_lin,nbr_col);
sprintf(file,"%s.PCR",name);
sprintf(comments,"diane_lissage");
write_data(file,TFR,dim,size,type,mode,nature,comments);

split_image(TFI,nbr_lin,nbr_col);
sprintf(file,"%s.PCI",name);
sprintf(comments,"diane_lissage");
write_data(file,TFI,dim,size,type,mode,nature,comments);

/*    FREE MEMORY
*/
free_square_float(tempr,nbr_lin);
free_square_float(SNR,nbr_lin);
free_square_float(FTB,nbr_lin);


/*    COMPUTE PHIT
*/
norm = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   norm += (TFR[i][j]*TFR[i][j]+TFI[i][j]*TFI[i][j]);
norm = (float)sqrt((double)norm/(double)(nbr_lin*nbr_col));
fprintf(stderr,"\nnorm of smoothing function = %f",norm);
flag = -1;
split_image(TFR,nbr_lin,nbr_col);
split_image(TFI,nbr_lin,nbr_col);
err = fourierC(TFR,TFI,nbr_lin,nbr_col,flag);
if (err != 0)
  {
  fprintf(stderr,"\nFATAL ERROR: problem with complex fourier transform\n\n");
  exit(-1);
  }
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   TFR[i][j] = (float)sqrt((double)(TFR[i][j]*TFR[i][j]+TFI[i][j]*TFI[i][j]));


/*   STORE RESULT  .PHIT
*/
sprintf(file,"%s.PHIT",name);
sprintf(comments,"diane_lissage");
write_data(file,TFR,dim,size,type,mode,nature,comments);
free_square_float(TFR,nbr_lin);
free_square_float(TFI,nbr_lin);
fprintf(stderr,"\n\n");

}
