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
/* NB:   generic_name.SNR point spread fuction (input)         */
/*       generic_name.FTR real part of FT of input image(input)*/
/*       generic_name.FTI imag. part of FT of input image(input)*/
/*                                                             */
/* VERSION: March  1993                                        */
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
auto        float    **TFR, **TFI, **SNR;
auto        float    alpha, cumul;

/*    TEST OF COMMAND LINE
*/
if (argc != 2)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndiane_sigmai  generic_name\n\n");
  exit(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
strcpy(name,argv[1]);

/*    READS INPUT SNR
*/
dim = 2;
strcpy(type,"float");
strcpy(nature,"real");
sprintf(file,"%s.SNR", argv[1]);
SNR = (float **)read_data(file,dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];

/*    READS INPUT TFR
*/
dim = 2;
strcpy(type,"float");
strcpy(nature,"real");
sprintf(file,"%s.TFR", argv[1]);
TFR = (float **)read_data(file,dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];

/*    READS INPUT TFI
*/
dim = 2;
strcpy(type,"float");
strcpy(nature,"real");
sprintf(file,"%s.TFI", argv[1]);
TFI = (float **)read_data(file,dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];


for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
    {
    TFR[i][j] = sqrt(TFR[i][j]*TFR[i][j] + TFI[i][j]*TFI[i][j]);
if (SNR[i][j] == 0.0)
   {
   fprintf(stderr,"\nFATAL ERROR: panique a bord on evacue!\n\n");
   exit(-1);
   }
    TFR[i][j] /= SNR[i][j];
    }


/*    STORAGE OF THE RESULT   *.SIG
*/
sprintf(file,"%s.SIG",name);
strcpy(comments,"diane_sigmai");
write_data(file,TFR,dim,size,type,mode,nature,comments);
fprintf(stderr,"\n\n");

}
