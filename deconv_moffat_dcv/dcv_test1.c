/****************************************************************** 
* dcv_test1.c 
* From the program in Fortran 90
* To emulate Herve Carfantan's simulations PI1 and PI2 
*
* JLP
* Version 19/03/2001
*******************************************************************/ 
#include <jlp_dcv.h>

#define NN1 256

/* Prototypes of functions contained here: */
int output_fits(double *yy0, double *yy, double *yyb,
                double *hh, int nn, char *extension);
int herve_signal1(double *yy0, int nn);
int herve_psf1(double *hh, int nn);
int herve_signal2(double *yy0, int nn, int nn1);
int herve_psf2(double *hh, int nn, int *hh_nx);
int dcv_conv1D(double *yy0, double *hh, double *yy, int nn, char *option);
int dcv_conv_1D(double *yy0, double *hh, double *yy, int nn);

int main(int argc, char *argv[])
{
int nn1=NN1, display_all, simu1;
int nn, hh_nx;
INT4 nx, ny;
double xx[NN1+1], yy0[NN1+1], yy0_pw[NN1+1], hh[NN1+1], yy[NN1+1]; 
double yyb[NN1+1], hh_pw[NN1+1], w0[NN1+1], w1[NN1+1], yy_noise[NN1+1];
double snr, noise_level, sum;
float ww;
char ans[1], plotdev[32], xlabel[40], ylabel[40], title[40]; 
char extension[20]; 
register int i;

/* snr in dB (20 db = 100.) */
snr=20.;

printf(" Program to generate a signal for simulations\n");
printf(" JLP version 22/03/2001 \n");

/* Reduce the number of arguments when "runs is used and  argc=7 */
if(argc == 7)
 {
 for(i = 6; i > 0; i--)
  {
   printf("*argv[%d] = >%c< (argc=%d)\n", i, *argv[i], argc);
   if(*argv[i] == ' ' || *argv[i] == '\0')argc--;
   else break;
  }
 }

if(argc == 3)
{
sscanf(argv[1],"%d",&simu1);
sscanf(argv[2],"%d",&display_all);
}
else
{
printf(" Usage: dcv_test1 simu_number display_all \n");
printf(" example: dcv_test1 1 0 \n");
printf(" \n\n Interactive input of parameters:\n");
printf(" Menu: \n");
printf(" 1. Simulation PI1 \n");
printf(" 2. Simulation PI2 \n");
printf(" Enter your choice: \n");
scanf("%d",&simu1);
printf(" Display all figures? 1 or 0 (1=yes 0=no) :\n");
scanf("%d",&display_all);
}

if(simu1 == 1)
  {
  printf(" OK: Simulation PI1 \n");
  nn = 100;
  }
else
  {
  simu1 = 2;
  printf(" OK: Simulation PI2 \n");
  nn = 250;
  }

strcpy(xlabel," "); 
strcpy(ylabel," "); 
/* xterm_s: small size 
*  xterm: normal size
*/
strcpy(plotdev,"xterm_s"); 
/* Set to zero all arrays: */
for(i = 0; i <= nn; i++)
 { 
 yy0[i] = 0.; hh[i] = 0.; yy[i] = 0.; yyb[i] = 0.;
 }

/* To prepare FFT's: */
fftw_setup(nn1, 1);

/* X axis: */
for(i = 0; i <= nn1; i++) xx[i] = i;

/* Generate Herve's signal */
if(simu1 == 1)
  herve_signal1(yy0, nn);
else
  herve_signal2(yy0, nn, nn1);

/* Display original signal: */
if(display_all)
  {
  printf("Display original signal:\n");
  strcpy(title, "Original signal"); 
  jlp_display1(xx, yy0, 1, nn, xlabel, ylabel, title, plotdev);
  }

/* Compute power spectrum of original signal: */
for(i = 0; i <= nn1; i++) 
  {yy0_pw[i] = yy0[i]; w0[i] = 0.;}

fftw_double(&yy0_pw[1], &w0[1], nn1, 1, 1);

sum = 0.;
for(i = 1; i <= nn1; i++) sum += yy0[i]; 
printf(" Sum(signal): %.2f\n", sum);

printf(" Central value of Real part of FT: %.2f\n",yy0_pw[1]);
printf(" Ratio: %.2f\n",sum/yy0_pw[1]);
for(i = 1; i <= nn1; i++) 
  yy0_pw[i] = SQUARE(yy0_pw[i]) + SQUARE(w0[i]);

  ny = 1; nx = nn1;
  RECENT_FFT_1D_X(&yy0_pw[1], &yy0_pw[1], &nx, &ny, &nx);

/* Display power spectrum of original signal: */

if(display_all)
  {
  printf("Display power spectrum of original signal:\n");
  strcpy(title, "Power spectrum of original signal"); 
  for(i = 0; i <= nn1; i++) 
     {
     w0[i] = ((double)i - nn1/2)/(double)nn1;
     w1[i] = MY_LOG10(yy0_pw[i]);
     }
  jlp_display1(w0,w1,1,nn1,xlabel,ylabel,title,plotdev);
  }

/* Generate Herve's PSF: */
if(simu1 == 1)
  {
  herve_psf1(hh,nn);
  hh_nx = nn;
  }
else
  herve_psf2(hh,nn,&hh_nx);

/* Display PSF: */
if(display_all)
  {
  printf("Display PSF of filter:\n");
  strcpy(title, "PSF of the filter"); 
  jlp_display1(xx,hh,1,hh_nx,xlabel,ylabel,title,plotdev);
  }

/* Compute transfer function: */
  for(i = 1; i <= nn1; i++) 
    {
    hh_pw[i] = hh[i]; 
    w0[i] = 0.;
    }
/* Arguments: (real,imag,nx,ny,direct=1 or inverse=-1) */
 fftw_double(&hh_pw[1], &w0[1], nn1, 1, 1);

  for(i = 1; i <= nn1; i++) 
     hh_pw[i] = hh_pw[i] * hh_pw[i] + w0[i] * w0[i];
  ny = 1; nx = nn1;
  RECENT_FFT_1D_X(&hh_pw[1], &hh_pw[1], &nx, &ny, &nx);

/* Display Transfer function: */
if(display_all) {
  printf("Frequency response of the filter:\n");
  strcpy(title,"Frequency response of the filter"); 
  for(i = 1; i <= nn1; i++) 
     w0[i] = ((double)i - nn1/2)/(double)nn1;
  for(i = 1; i <= nn1; i++) w1[i] = MY_LOG10(hh_pw[i]);
  jlp_display1(w0,w1,1,nn1,xlabel,ylabel,title,plotdev);
  }

/* Generate filtered signal: */
 dcv_conv1D(yy0,hh,yy,nn1,"full");

/* Display filtered signal: */
if(display_all) {
  printf("Display filtered signal:\n");
  strcpy(title,"Signal filtered by the PSF"); 
  jlp_display1(xx,yy,1,nn,xlabel,ylabel,title,plotdev);
  }

/* Compute power spectrum of filtered signal: */ 
if(display_all) {
  for(i = 1; i <= nn1; i++) {w0[i] = yy[i]; w1[i] = 0.;}
  fftw_double(&w0[1], &w1[1], nn1, 1, 1);
  for(i = 1; i <= nn1; i++) w0[i] = w0[i]*w0[i] + w1[i]*w1[i];
  ny = 1; nx = nn1;
  RECENT_FFT_1D_X(&w0[1], &w0[1], &nx, &ny, &nx);
  printf("Display power spectrum of filtered signal:\n");
  strcpy(title,"Power spectrum of filtered signal"); 
  for(i = 1; i <= nn1; i++) w1[i] = MY_LOG10(w0[i]);
  for(i = 1; i <= nn1; i++) 
     w0[i] = ((double)i - nn1/2)/(double)nn1;
  jlp_display1(w0,w1,1,nn1,xlabel,ylabel,title,plotdev);
  }

/* Generate noisy signal:
* Multiply by 12, since variance of uniform random function is 1/12:
*/
  sum = 0;
  for(i = 1; i <= nn1; i++) sum += yy[i] * yy[i]; 
noise_level = 12. * sum * pow(10.,-snr/10.) / (double)nn;
noise_level = sqrt(noise_level);
printf("Simulated signal with SNR= %.2f dB => noise level: %.3f\n",
        snr, noise_level);

/* Generate random vector (between 0 and 1) */
/* rand is a random generator between 0 and RAND_MAX: */
for(i = 1; i <= nn; i++) yy_noise[i] = (double)rand() / (double)RAND_MAX;

/* Center the noise (since mean=0.5) */
  for(i = 1; i <= nn1; i++) yyb[i] = yy[i] + (yy_noise[i] - 0.5) * noise_level;

/* Display the simulated (noisy and filtered) signal: */
if(display_all) {
  printf("Display simulated signal:\n");
  strcpy(title,"Filtered signal with noise (SNR=20)"); 
  jlp_display1(xx,yyb,1,nn,xlabel,ylabel,title,plotdev);
  }

/* Compute power spectrum of filtered signal: */
if(display_all) {
  for(i = 1; i <= nn1; i++) {w0[i] = yyb[i]; w1[i] = 0.;}
  fftw_double(&w0[1], &w1[1], nn1, 1, 1);
  for(i = 1; i <= nn1; i++) w0[i] = w0[i]*w0[i] + w1[i]*w1[i];
  ny = 1; nx = nn1;
  RECENT_FFT_1D_X(&w0[1], &w0[1], &nx, &ny, &nx);
  printf("Display power spectrum of noisy filtered signal:\n");
  strcpy(title,"Pow. sp. of noisy filtered signal"); 
  for(i = 1; i <= nn1; i++) w1[i] = MY_LOG10(w0[i]);
  for(i = 1; i <= nn1; i++) 
     w0[i] = ((double)i - nn1/2)/(double)nn1;
  jlp_display1(w0,w1,1,nn1,xlabel,ylabel,title,plotdev);
  }

if(simu1 == 1)
  strcpy(extension,"simu1");
else
  strcpy(extension,"simu2");
output_fits(&yy0[1], &yy[1], &yyb[1], &hh[1], nn, extension);

return(0);
}
/************************************************
* To emulate Herve Carfantan's simulation PI1 
************************************************/
int herve_signal1(double *yy0, int nn)
{
register int i;

/* Y axis (Cf. Herve Carfantan's simulation): */
for(i = 0; i <= nn; i++) yy0[i] = 0.;
for(i = 1; i <= 10; i++) yy0[i] = 1.;
for(i = 11; i <= 26; i++) yy0[i] = 1.+0.25*(double)(i-10);
for(i = 27; i <= 46; i++) yy0[i] = 2.;
for(i = 47; i <= 67; i++) yy0[i] = 2.-0.1*(double)(i-47);
for(i = 68; i <= 77; i++) yy0[i] = 4.;

return(0);
}
/************************************************
* To emulate Herve Carfantan's simulation PI1 
* PSF of the filter (rectangle) 
************************************************/
int herve_psf1(double *hh, int nn)
{
register int i;
/* Y axis (Cf. Herve Carfantan's simulation): */ 
for(i = 0; i <= nn; i++) hh[i] = 0.;
for(i = 1; i <= 5; i++) hh[i] = 1.;

return(0);
}
/************************************************
* To emulate Herve Carfantan's simulation PI2 
************************************************/
int herve_signal2(double *yy0, int nn, int nn1)
{
register int i;
FILE *fp;

printf("herve_signal2: simulation #2  nn= %d nn1=%d\n", nn, nn1);
if((fp = fopen("bgg.asc","r")) == NULL)
 {
 printf("herve_signal2/Fatal error opening bgg.asc input file\n");
 exit(-1);
 }
for(i = 0; i <= nn1; i++) yy0[i] = 0.;
/* format(2X,E14.7) */
for(i = 1; i <= nn; i++) fscanf(fp,"  %lf",&yy0[i]);

fclose(fp);
return(0);
}
/************************************************
* To emulate Herve Carfantan's simulation PI2 
* PSF of the filter (rectangle) 
************************************************/
int herve_psf2(double *hh, int nn, int *hh_nx)
{
register int i;
FILE *fp;

printf(" signal #2 \n");
if((fp = fopen("ricker_ri.asc","r")) == NULL)
 {
 printf("herve_psf2/Fatal error opening ricker_ri.asc input file\n");
 exit(-1);
 }

for(i = 0; i <= 256; i++) hh[i] = 0.;
*hh_nx = 32;
for(i = 1; i <= *hh_nx; i++) fscanf(fp,"  %lf",&hh[i]);

fclose(fp);
return(0);
}
/************************************************
* To emulate Herve Carfantan's conv1D
* IN:
*  yy0: signal to filter
*  hh: filter
*  nn: size of the original signal
*  option: "full", "pre ", "post"
* OUT:
*  yy: filtered signal
************************************************/
int dcv_conv1D(double *yy0, double *hh, double *yy, int nn, char *option)
{

switch (option[0])
   {
/* "full" */
   case 'f':
     dcv_conv_1D(yy0, hh, yy, nn);
     break;
/* "pre ", "post", default: */
   default:
    printf(" dcv_conv1D/Fatal: bad option\n");
    exit(-1);
   }   

return(0);
}
/************************************************
* Convolution product
* IN:
*  yy0: signal to filter
*  hh: filter
*  nn: size of the yy0, hh, arrays (here should be equal to nn1=256) 
*  option: "full", "pre ", "post"
* OUT:
*  yy: filtered signal
************************************************/
int dcv_conv_1D(double *yy0, double *hh, double *yy, int nn)
{
register int i;
double *yy0_re, *yy0_im, *hh_re, *hh_im, *ww;

if((yy0_re = (double *)malloc((nn+1)*sizeof(double))) == NULL
 || (yy0_im = (double *)malloc((nn+1)*sizeof(double))) == NULL
 || (hh_re = (double *)malloc((nn+1)*sizeof(double))) == NULL
 || (hh_im = (double *)malloc((nn+1)*sizeof(double))) == NULL
 || (ww = (double *)malloc((nn+1)*sizeof(double))) == NULL)
  {
  printf("dcv_conv_1D/Fatal error allocating memory: nn=%d\n",nn);
  exit(-1);
  }
for(i = 0; i <= nn; i++) yy0_re[i] = yy0[i];
for(i = 0; i <= nn; i++) yy0_im[i] = 0.;
fftw_double(&yy0_re[1], &yy0_im[1], nn, 1, 1);

for(i = 0; i <= nn; i++) hh_re[i] = hh[i];
for(i = 0; i <= nn; i++) hh_im[i] = 0.;
fftw_double(&hh_re[1], &hh_im[1], nn, 1, 1);

for(i = 0; i <= nn; i++) yy[i] = yy0_re[i]*hh_re[i] - yy0_im[i]*hh_im[i];
for(i = 0; i <= nn; i++) ww[i] = yy0_re[i]*hh_im[i] + yy0_im[i]*hh_re[i];
fftw_double(&yy[1], &ww[1], nn, 1, -1);

return(0);
}
/***********************************************
*
************************************************/
int output_fits(double *yy0, double *yy, double *yyb,
                double *hh, int nn, char *extension) 
{
 char filename[60], comments[80]; 
 INT4 nx, ny;

nx=nn; ny=1;
sprintf(filename,"%s_yy0",extension);
strcpy(comments,"Original signal");
JLP_D_WRITEIMAG(yy0, &nx, &ny, &nx, filename, comments);

sprintf(filename,"%s_yy",extension);
strcpy(comments,"Filtered signal");
JLP_D_WRITEIMAG(yy, &nx, &ny, &nx, filename, comments);

sprintf(filename,"%s_yyb",extension);
strcpy(comments,"Simulated signal");
JLP_D_WRITEIMAG(yyb, &nx, &ny, &nx, filename, comments);

sprintf(filename,"%s_psf",extension);
strcpy(comments,"PSF");
JLP_D_WRITEIMAG(hh, &nx, &ny, &nx, filename, comments);

return(0);
}

