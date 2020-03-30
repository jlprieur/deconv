/***************************************************************
* 
* PROGRAM: diane_snr
* 
* PURPOSE: estimates signal to noise ratio of an image
*          in Fourier domain assuming additive gaussian noise
* 
* INPUT:  argv[1] = input image 
*         argv[2] = signal to noise ratio in direct space
*         argv[3] = reliability (1 good, 2 medium, 3 bad)
*         argv[4] = generic name for output 
*                  
* From Karim's version of December 1992 
*                     
* AUTHOR: JLP, SR, JV  
*         translated to C by Karim BOUYOUCEF 
* 
* JLP
* Version 20/05/1999
* Version 01/02/2011
***************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <jlp_ftoc.h>

int  main(int argc, char *argv[])
{
/*    DECLARATIONS
*/
INT4  nx, ny;
INT_PNTR pntr_ima;
char  name[60], filename[61], comments[81];
register int i, k, l;
int   err, alphat, flag;
float *image, *tfr, *tfi, *sigma_i, *noise_r, *noise_i, *snr;
float snrd, alk, sigma_n, work;

/*    TEST OF COMMAND LINE
*/
/* Carefull: 7 parameters always, using JLP "runs" */
if(argc == 7 && *argv[4]) argc = 5;
if (argc != 5) {
  printf("\nUnexpected number of arguments\n");
  printf("\nUSAGE:\n");
  printf("\ndc_snr  input_image  snr  reliability  generic_name\n\n");
  printf("snr = signal to noise ratio in direct space : mean(signal)/var(signal)");
  printf("\nreliability = reliability of snr   1:good  2:medium  3:bad");
  printf("\ngeneric_name = generic name of results  .SNR  .SIG  .RFT  .IFT\n\n");
  return(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
err = sscanf(argv[2],"%f",&snrd);
if ((err !=1) || (snrd < 0.0) || (snrd > 100000.0))
  {
  printf("\nFATAL ERROR: snr in direct space [%s] incorrect\n\n",argv[2]);
  exit(-1);
  }
err = sscanf(argv[3],"%d",&alphat);
if ((err !=1) || (alphat < 1) || (alphat > 3))
  {
  printf("\nFATAL ERROR: reliability [%s] incorrect\n\n",argv[3]);
  exit(-1);
  }
strcpy(name,argv[4]);

printf("\n************************");
printf("\nPROGRAM :      diane_snr");
printf("\n************************");
printf("\nsnr in direct space = %f",snrd);
printf("\nreliability = %d",alphat);

/*    READS INPUT IMAGE
*/
JLP_INQUIFMT();
strcpy(filename,argv[1]);
JLP_VM_READIMAG1(&pntr_ima,&nx,&ny,filename,comments);
image = (float *)pntr_ima;

/*    FOURIER TRANSFORM OF IMAGE
*/
tfr = (float *) malloc(ny * nx *sizeof(float));
tfi = (float *) malloc(ny * nx *sizeof(float));
for (i = 0; i < nx*ny; i++) tfr[i] = image[i];
/* JLP99: Need recentring here, for some unknown reason... */
err = RECENT_FFT(tfr,tfr,&nx,&ny,&nx);
for (i = 0; i < nx*ny; i++) tfi[i] = 0.;
// Direct FFT if flag=1:
flag = 1;
err = fftw_2D_float(tfr,tfi,(int)nx,(int)ny,(int)flag);
if (err != 0)
  {
  printf("\nFATAL ERROR: problem with complex fourier transform\n\n");
  exit(-1);
  }

/*    MEMORY ALLOCATION
*/
sigma_i = (float *) malloc(ny * nx *sizeof(float));
noise_r = (float *) malloc(ny * nx *sizeof(float));
noise_i = (float *) malloc(ny * nx *sizeof(float));
snr = (float *) malloc(ny * nx *sizeof(float));

for (i = 0; i < nx*ny; i++) sigma_i[i] = 0.0;
alk = (float)alphat/snrd;

/****  random generator initialization  ****/
srand(9);

for (k=1; k<=5; k++)
  {
  for (i = 0; i < nx*ny; i++)
    {
    sigma_n = image[i]*alk;

    /****  gaussian law  ****/
    work = 0.0;
    for (l=0; l<48; l++)
       work += random()/2147483647.0-0.5;
    work *= 0.5;

    noise_r[i] = work*sigma_n;
    noise_i[i] = 0.0;
    image[i] += noise_r[i];
    }
  
  err = fftw_2D_float(noise_r,noise_i,nx,ny,flag);
  if (err != 0)
    {
    printf("\nFATAL ERROR: problem with complex fourier transform\n\n");
    exit(-1);
    }
   
for (i = 0; i < nx*ny; i++)
    sigma_i[i] += sqrt(noise_r[i] * noise_r[i] + noise_i[i] * noise_i[i]);
  }

for (i = 0; i < nx*ny; i++)
  {
  sigma_i[i] *= 2.0/5.0;
  snr[i] = (sqrt(tfr[i]*tfr[i]+tfi[i]*tfi[i])/sigma_i[i]);
  }

err = RECENT_FFT(tfr,tfr,&nx,&ny,&nx);
if (err != 0)
  {
  printf("\nFATAL ERROR: problem when splitting tfr\n\n");
  exit(-1);
  }
err = RECENT_FFT(tfi,tfi,&nx,&ny,&nx);
if (err != 0)
  {
  printf("\nFATAL ERROR: problem when splitting tfi\n\n");
  exit(-1);
  }

err = RECENT_FFT(snr,snr,&nx,&ny,&nx);
if (err != 0)
  {
  printf("\nFATAL ERROR: problem when splitting snr\n\n");
  exit(-1);
  }
err = RECENT_FFT(sigma_i,sigma_i,&nx,&ny,&nx);
if (err != 0)
  {
  printf("\nFATAL ERROR: problem when splitting sigma_i\n\n");
  exit(-1);
  }

/*    STORAGE OF RESULTS
      .SIG  .SNR  .RFT  .IFT
*/
sprintf(filename,"%s.SNR",name);
JLP_WRITEIMAG(snr,&nx,&ny,&nx,filename,comments);
sprintf(filename,"%s.SIG",name);
JLP_WRITEIMAG(sigma_i,&nx,&ny,&nx,filename,comments);
sprintf(filename,"%s.RFT",name);
JLP_WRITEIMAG(tfr,&nx,&ny,&nx,filename,comments);
sprintf(filename,"%s.IFT",name);
JLP_WRITEIMAG(tfi,&nx,&ny,&nx,filename,comments);

/*    FREE MEMORY
*/
free(snr);
free(sigma_i);
free(noise_i);
free(noise_r);
free(tfr);
free(tfi);
free(image);
printf("\n\n");

return(0);
}
