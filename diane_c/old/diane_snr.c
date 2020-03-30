/***************************************************************/
/*                                                             */
/* PROGRAM: diane_snr                                          */
/*                                                             */
/* PURPOSE: estimates signal to noise ratio of an image        */
/*          in Fourier domain assuming additive gaussian noise */
/*                                                             */
/* INPUT:  argv[1] = input image                               */
/*         argv[2] = signal to noise ratio in direct space     */
/*         argv[3] = reliability (1 good, 2 medium, 3 bad)     */
/*         argv[4] = generic name for output                   */
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
extern  int       split_image();
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
register    int      i, j, k, l;
auto        int      nbr_lin, nbr_col, err, alphat, flag;
auto        float    **image, **tfr, **tfi, **sigma_i, **noise_r, **noise_i, **snr;
auto        float    snrd, alk, sigma_n, work;

/*    TEST OF COMMAND LINE
*/
if (argc != 5)
  {
  fprintf(stderr,"\nUnexpected number of arguments\n");
  fprintf(stderr,"\nUSAGE:\n");
  fprintf(stderr,"\ndiane_snr  input_image  snr  reliability  generic_name\n\n");
  fprintf(stderr,"snr = signal to noise ratio in direct space : mean(signal)/var(signal)");
  fprintf(stderr,"\nreliability = reliability of snr   1:good  2:medium  3:bad");
  fprintf(stderr,"\ngeneric_name = generic name of results  .SNR  .SIG  .TFR  .TFI\n\n");
  exit(-1);
  }


/*    READ COMMAND LINE PARAMETERS
*/
err = sscanf(argv[2],"%f",&snrd);
if ((err !=1) || (snrd < 0.0) || (snrd > 100000.0))
  {
  fprintf(stderr,"\nFATAL ERROR: snr in direct space [%f] incorrect\n\n",argv[2]);
  exit(-1);
  }
err = sscanf(argv[3],"%d",&alphat);
if ((err !=1) || (alphat < 1) || (alphat > 3))
  {
  fprintf(stderr,"\nFATAL ERROR: reliability [%d] incorrect\n\n",argv[3]);
  exit(-1);
  }
strcpy(name,argv[4]);

fprintf(stderr,"\n************************");
fprintf(stderr,"\nPROGRAM :      diane_snr");
fprintf(stderr,"\n************************");
fprintf(stderr,"\nsnr in direct space = %f",snrd);
fprintf(stderr,"\nreliability = %d",alphat);

/*    READS INPUT IMAGE
*/
dim = 2;
strcpy(type,"float");
strcpy(nature,"real");
image = (float **)read_data(argv[1],dim,size,type,mode,nature,comments);
nbr_lin = size[0];
nbr_col = size[1];


/*    FOURIER TRANSFORM OF IMAGE
*/
tfr = alloc_square_float(nbr_lin,nbr_col);
tfi = alloc_square_float(nbr_lin,nbr_col);
for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  tfr[i][j] = image[i][j];
flag = 1;
err = fourierC(tfr,tfi,nbr_lin,nbr_col,flag);
if (err != 0)
  {
  fprintf(stderr,"\nFATAL ERROR: problem with complex fourier transform\n\n");
  exit(-1);
  }

/*    MEMORY ALLOCATION
*/
sigma_i = alloc_square_float(nbr_lin,nbr_col);
noise_r = alloc_square_float(nbr_lin,nbr_col);
noise_i = alloc_square_float(nbr_lin,nbr_col);
snr = alloc_square_float(nbr_lin,nbr_col);

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
   sigma_i[i][j] = 0.0;
alk = (float)alphat/snrd;

/****  random generator initialization  ****/
srand(9);

for (k=1; k<=5; k++)
  {
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
    {
    sigma_n = image[i][j]*alk;

    /****  gaussian law  ****/
    work = 0.0;
    for (l=0; l<48; l++)
       work += random()/2147483647.0-0.5;
    work *= 0.5;

    noise_r[i][j] = work*sigma_n;
    noise_i[i][j] = 0.0;
    image[i][j] += noise_r[i][j];
    }
  
  err = fourierC(noise_r,noise_i,nbr_lin,nbr_col,flag);
  if (err != 0)
    {
    fprintf(stderr,"\nFATAL ERROR: problem with complex fourier transform\n\n");
    exit(-1);
    }
   
  for (i=0; i<nbr_lin; i++)
  for (j=0; j<nbr_col; j++)
    sigma_i[i][j] += sqrt(noise_r[i][j] * noise_r[i][j] + noise_i[i][j] * noise_i[i][j]);
  }

for (i=0; i<nbr_lin; i++)
for (j=0; j<nbr_col; j++)
  {
  sigma_i[i][j] *= 2.0/5.0;
  snr[i][j] = (sqrt(tfr[i][j]*tfr[i][j]+tfi[i][j]*tfi[i][j])/sigma_i[i][j]);
  }

err = (int) split_image(tfr,nbr_lin,nbr_col);
if (err != 0)
  {
  fprintf(stderr,"\nFATAL ERROR: problem when splitting tfr\n\n");
  exit(-1);
  }
err = (int)split_image(tfi,nbr_lin,nbr_col);
if (err != 0)
  {
  fprintf(stderr,"\nFATAL ERROR: problem when splitting tfi\n\n");
  exit(-1);
  }
err = (int)split_image(snr,nbr_lin,nbr_col);
if (err != 0)
  {
  fprintf(stderr,"\nFATAL ERROR: problem when splitting snr\n\n");
  exit(-1);
  }
err = (int)split_image(sigma_i,nbr_lin,nbr_col);
if (err != 0)
  {
  fprintf(stderr,"\nFATAL ERROR: problem when splitting sigma_i\n\n");
  exit(-1);
  }

/*    STORAGE OF RESULTS
      .SIG  .SNR  .TFR  .TFI
*/
dim = 2;
size[0] = nbr_lin;
size[1] = nbr_col;
strcpy(type,"float");
strcpy(mode,"bin");
strcpy(nature,"real");
strcpy(comments,"diane_snr");

sprintf(file,"%s.SNR",name);
write_data(file,snr,dim,size,type,mode,nature,comments);
sprintf(file,"%s.SIG",name);
write_data(file,sigma_i,dim,size,type,mode,nature,comments);
sprintf(file,"%s.TFR",name);
write_data(file,tfr,dim,size,type,mode,nature,comments);
sprintf(file,"%s.TFI",name);
write_data(file,tfi,dim,size,type,mode,nature,comments);

/*    FREE MEMORY
*/
free_square_float(snr,nbr_lin);
free_square_float(sigma_i,nbr_lin);
free_square_float(tfr,nbr_lin);
free_square_float(tfi,nbr_lin);
free_square_float(image,nbr_lin);
free_square_float(noise_r,nbr_lin);
free_square_float(noise_i,nbr_lin);
fprintf(stderr,"\n\n");
}
