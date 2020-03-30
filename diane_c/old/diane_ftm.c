/***************************************************************/
/*                                                             */
/* PROGRAM: diane_ftm                                          */
/*                                                             */
/* PURPOSE: normalize PSF respecting L1 normalization          */
/*          performs fourier transform of PSF and bound it     */
/*                                                             */
/* INPUT:  argv[1] = generic name                              */
/*         argv[2] = Lower thresold for PSF in fourier domain  */
/*                                                             */
/* NB:   generic_name.PSF point spread fuction (input)         */
/*       generic_name.FTB bounded imodulus FT of .PSF (output) */
/*                                                             */
/* VERSION: December 1992                                      */
/*                                                             */
/* AUTHOR: JLP, SR, JV                                         */
/*         translated to C by Karim BOUYOUCEF                  */
/*                                                             */
/***************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>

extern  long int  read_data();
extern  void      write_data();
extern  float **  alloc_square_float();
extern  void      free_square_float();
extern  void      split_image();
extern  int       fourierC();

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
auto        float    **PSF, **tfi;
auto        float    alpha, cumul;

/*    TEST OF COMMAND LINE
*/
if (argc != 3)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndiane_ftm  generic_name  Lower_thresold\n\n");
  fprintf(stderr,"generic_name = *.PSF (input)  *.FTB (output)");
  fprintf(stderr,"\nLower_thresold = lower thresold of modulus fourier transform of *.PSF\n\n");
  exit(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
strcpy(name,argv[1]);
err = sscanf(argv[2],"%f",&alpha);
if ((err != 1) || (alpha > 1.))
  {
  fprintf(stderr,"\nFATAL ERROR: Lower thresold [%f] incorrect\n\n",argv[2]);
  exit(-1);
  }

fprintf(stderr,"\n***************************");
fprintf(stderr,"\nPROGRAM :         diane_ftm");
fprintf(stderr,"\n***************************");
fprintf(stderr,"\nlower thresold for FT of PSF = %f",alpha);

/*    READS INPUT PSF
*/
dim = 2;
strcpy(type,"float");
strcpy(nature,"real");
sprintf(file,"%s.PSF", argv[1]);
PSF = (float **)read_data(file,dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];


/*    NORMALIZE L1 THE PSF
*/
cumul = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  cumul += PSF[i][j];
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  PSF[i][j] /= cumul;
tfi = alloc_square_float(nbr_lin,nbr_col);
flag = 1;
err = fourierC(PSF,tfi,nbr_lin,nbr_col,flag);
if (err != 0)
  {
  fprintf(stderr,"\nFATAL ERROR: problem with complex fourier transform\n\n");
  exit(-1);
  }

/*    MODULUS OF PSF FOURIER TRANSFORM
*/
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  PSF[i][j] = (float)sqrt((double)(PSF[i][j]*PSF[i][j]+tfi[i][j]*tfi[i][j]));

/*    THRESOLD
*/
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  if (PSF[i][j] < alpha)
    PSF[i][j] = 0.0;

split_image(PSF,nbr_lin,nbr_col);
free_square_float(tfi,nbr_lin);

/*    STORAGE OF THE RESULT   *.FTB
*/
sprintf(file,"%s.FTB",name);
strcpy(comments,"diane_ftm");
write_data(file,PSF,dim,size,type,mode,nature,comments);
free_square_float(PSF,nbr_lin);
fprintf(stderr,"\n\n");

}
