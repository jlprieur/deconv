/***************************************************************
*                                                    
* FUNCTION: direct_snr                                
*                                                      
* PURPOSE: compute signal to noise ratio in direct space
*          assuming that noise is of additive nature  
*          and taking into account a support constraint
*                             
* INPUT:  argv[1] = input_image
*         argv[2] = support_constraint 
*                      
* From Karim's version of March 1993  
*                        
* AUTHOR: Karim BOUYOUCEF 
*
* JLP
* Version 20/05/1999
***************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

void  main(argc,argv)
int   argc;
char  *argv[];
{
/*
    DECLARATIONS
*/
int    nb_o, nb_n;   
INT4   nx, ny, nx1, ny1;
INT_PNTR pntr_ima;
float  mean_o, mean_n, var_o, var_n, SNR;
char   filename[61], comments[81]; 
register    int    i, j;
float  *image, *support;

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
JLP_BEGIN();
JLP_INQUIFMT();
strcpy(filename,argv[1]);
JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
image = (float *)pntr_ima;

strcpy(filename,argv[2]);
JLP_VM_READIMAG1(&pntr_ima,&nx1,&ny1,filename,comments);
support = (float *)pntr_ima;

if ((nx != nx1) || (ny != ny1))
  {
  fprintf(stderr,"\n\aFATAL ERROR: [%s] & [%s] of different sizes\n",argv[1],argv[2]);
  exit(-1);
  }

/*
    PERFORMS OPERATION
*/
nb_o = nb_n = 0;
mean_o = mean_n = 0.0;
for (i = 0; i < nx*ny; i++)
  {
  if (support[i] == 0.0)
     {
     nb_n ++;
     mean_n += image[i];
     }
  else 
     {
     nb_o ++;
     mean_o += image[i];
     }
  }
mean_n /= (float)nb_n;
mean_o /= (float)nb_o;

nb_o = nb_n = 0;
for (i = 0; i < nx*ny; i++)
  {
  if (support[i] == 0.0)
     {
     nb_n ++;
     var_n = (image[i]-mean_n)*(image[i]-mean_n);
     }
  else 
     {
     nb_o ++;
     var_o = (image[i]-mean_o)*(image[i]-mean_o);
     }
   }
var_n /= (float)nb_n;
var_o /= (float)nb_o;

SNR = var_o/var_n;
fprintf(stderr,"\nSNR in direct space = %f\n\n",SNR);
free(image);
free(support);

JLP_END();
}
