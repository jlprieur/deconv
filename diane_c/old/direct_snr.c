/***************************************************************/
/*                                                             */
/* FUNCTION: direct_snr                                        */
/*                                                             */
/* PURPOSE: compute signal to noise ratio in direct space      */
/*          assuming that noise is of additive nature          */
/*          and taking into account a support constraint       */
/*                                                             */
/* INPUT:  argv[1] = input_image		               */
/*         argv[2] = support_constraint                        */
/*                                                             */
/* VERSION: March  1993                                        */
/*                                                             */
/* AUTHOR: Karim BOUYOUCEF                                     */
/*                                                             */
/***************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>

extern    long int    read_data();
extern    void        write_data();
extern    void       free_square_float();

void  main(argc,argv)
int   argc;
char  *argv[];
{
/*
    DECLARATIONS
*/
auto        int    dim, size[4], nbr_lin, nbr_col, nb_o, nb_n;   
auto        float  mean_o, mean_n, var_o, var_n, SNR;
auto        char   type[16], mode[4], nature[8], comments[1024]; 
register    int    i, j;
auto        float  **image, **support;

/*
    TEST OF LINE COMMAND
*/
if (argc != 3)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:");
  fprintf(stderr,"\ndirect_snr  input_image  support_constraint\n\n");
  exit(-1);
  }

/*
    READS INPUT FILES
*/
dim = 2;
strcpy(type,"float");
strcpy(nature,"real");
image = (float **) read_data(argv[1],dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];

support = (float **) read_data(argv[2],dim,size,type,mode,nature,comments);

if ((nbr_lin != size[0]) || (nbr_col != size[1]))
  {
  fprintf(stderr,"\n\aFATAL ERROR: [%s] & [%s] of different sizes\n",argv[1],argv[2]);
  exit(-1);
  }

/*
    PERFORMS OPERATION
*/
nb_o = nb_n = 0;
mean_o = mean_n = 0.0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  if (support[i][j] == 0.0)
     {
     nb_n ++;
     mean_n += image[i][j];
     }
  else 
     {
     nb_o ++;
     mean_o += image[i][j];
     }
  }
mean_n /= (float)nb_n;
mean_o /= (float)nb_o;

nb_o = nb_n = 0;
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  if (support[i][j] == 0.0)
     {
     nb_n ++;
     var_n = (image[i][j]-mean_n)*(image[i][j]-mean_n);
     }
  else 
     {
     nb_o ++;
     var_o = (image[i][j]-mean_o)*(image[i][j]-mean_o);
     }
   }
var_n /= (float)nb_n;
var_o /= (float)nb_o;

SNR = var_o/var_n;
fprintf(stderr,"\nSNR in direct space = %f\n\n",SNR);
free_square_float(image,nbr_lin);
free_square_float(support,nbr_lin);
}
