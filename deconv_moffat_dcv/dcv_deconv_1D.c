/****************************************************************** 
* Test program 
* To emulate Herve Carfantan's simulation 
*  spdiv => Direct spectral division
*  wiene => Wiener filter 
*  tikho => Tikhonov regularisation
*  maxen => Maximum Entropy method
*
* JLP
* Version 17/04/2002
*******************************************************************/ 
#include "jlp_ftoc.h"
#include "jlp_dcv.h"

extern double *yy0, *yy, *yyd, *w0, *w1;

int main(int argc, char *argv[])
{
double *hh, rms_err, ss, sum, alpha, ftol=1.e-8;
float *float_tab;
int n_simu, nn, nn1, hh_nx, iter, display_all, positive;
INT_PNTR pntr;
INT4 nx, ny;
register int i;
char prefix[20], title[40], filename[60], comments[80];
char option[12];

printf(" Program to deconvolve 1D signals\n");
printf(" JLP version 17/04/2002\n");

/* Reduce the number of arguments when "runs" is used and  argc=7 */
if(argc == 7)
 {
 for(i = 6; i > 0; i--)
  {
   if(*argv[i] == ' ' || *argv[i] == '\0')argc--;
   else break;
  }
 }

if(argc == 3)
{
sscanf(argv[1],"%d",&n_simu);
sscanf(argv[2],"%d",&display_all);
}
else
{
printf(" Usage: dcv_deconv_1D simu_number display_all \n");
printf(" example: dcv_deconv_1D 1 0 \n");
exit(-1);
}

printf(" OK: simulation #%d \n",n_simu);

/* Read it first to determine the size:
*/
if(n_simu == 1)
  strcpy(prefix, "simu1");
else
  strcpy(prefix, "simu2");

sprintf(filename,"%s_yy0",prefix);

JLP_INQUIFMT();

JLP_VM_READIMAG1(&pntr, &nx, &ny, filename, comments);
  if((ny != 1) || (n_simu == 1 && nx != 100) 
     || (n_simu == 2 && nx != 250))
  {
   printf(" Fatal error/input file has wrong size: nx=%d ny=%d \n", nx, ny);
   if(n_simu == 1)
      printf("For simulation #1: nx = 100 and ny = 1!\n");
   else
      printf("For simulation #2: nx = 250 and ny = 1!\n");
   exit(-1);
  }

nn = nx;

switch (n_simu)
  {
   case 1:
/**** nn=100 idim = 128 */
      nn1 = 128;
/*
* Length of the support of hh is hh_nx=5 (for this simulation)
*  hh_nx = 5
*/
      hh_nx = 5; 
      break;
    case 2:
/**** nn=250 idim = 256 */
       nn1 = 256;
/* Length of the support of hh is hh_nx=32 (for this simulation) */
       hh_nx = 32;
       break;
     default:
       printf(" Error: n_simu=%d\n",n_simu);
       exit(-1);
   }

/* Then allocate memory space: */
dcv_1D_init(nn1, nn, hh_nx);

/**************** 
*  Read input files:
*****************/

/*----------------------
* Original signal:
*/
sprintf(filename,"%s_yy0",prefix);
JLP_VM_READIMAG1(&pntr, &nx, &ny, filename, comments);
float_tab = (float *) pntr;
for(i = 0; i <= nn1; i++) yy0[i] = 0.;
TO_DOUBLE(float_tab, &yy0[1], &nx, &ny, &nx);

/*----------------------
* Signal to deconvolve:
*/
sprintf(filename,"%s_yyb",prefix);
for(i = 0; i <= nn1; i++) yy[i] = 0.;
JLP_VM_READIMAG1(&pntr, &nx, &ny, filename, comments);
   if(ny != 1 || nx != nn) {
   printf(" Fatal error: nx=%d (!= 1) ny=%d (!=nn)\n", nx, ny);
   exit(-1);
   }
float_tab = (float *) pntr;
for(i = 0; i <= nn1; i++) yy[i] = 0.;
TO_DOUBLE(float_tab, &yy[1], &nx, &ny, &nx);

/*----------------------
* PSF:
*/
sprintf(filename,"%s_psf",prefix);
if((hh = (double *)malloc((nn1 + 1) * sizeof(double))) == NULL)
  {
  printf("Fatal error allocating memory for hh: nn1=%d\n", nn1);
  exit(-1);
  }
JLP_VM_READIMAG1(&pntr, &nx, &ny, filename, comments);
   if(ny != 1 || nx != nn) {
   printf(" Fatal error: nx=%d (!= 1) ny=%d (!=nn)\n", nx, ny);
   exit(-1);
   }
float_tab = (float *) pntr;
for(i = 0; i <= nn1; i++) hh[i] = 0.;
TO_DOUBLE(float_tab, &hh[1], &nx, &ny, &nx);

/* Compute transfer function: */
dcv_1D_transfert(hh);

/* Set to zero all arrays: */
for(i = 0; i <= nn1; i++) 
  {w0[i] = 0.; w1[i] = 0.;}

/*--------------------------------
* Check if input is OK:
*/
if(display_all) dcv_1D_display_input(nn);

/*--------------------------------
* Deconvolution
*--------------------------------
*/

printf(" ==================== Options: ================\n");
printf(" - spdiv: spectral division \n");
printf(" - wiene: Wiener filter \n");
printf(" Regularisation on x_i only \n");
printf(" - tikho: Tikhonov's regularisation\n");
printf(" - ggaus: Generalized Gauss's regularisation (p=1.1)\n");
printf(" - maxen: regularisation with maximum entropy\n");
printf(" - sqrtr: convex sqrt(s2+x2) regularisation\n");
printf(" Regularisation on x_i+1 - x_i \n");
printf(" - gmark: Gauss-Markov''s regularisation (p=2)\n");
printf(" ======== Enter the option you want: ===========\n");
scanf("%s", option);
printf(" OK: option= %s \n",option);

printf(" Positive constraint? (1 if yes, 0 otherwise)\n");
scanf("%d", &positive);

/*--------------------------------
* Deconvolution by spectral division:
*/
  if(!strncmp(option,"spdiv",5)) {
      printf(" calling dcv_1D_spdiv \n");
      dcv_1D_spdiv();
      printf(" exit from dcv_1D_spdiv \n");
      sum = 0.;
      for(i = 1; i <= nn; i++) sum += SQUARE(yyd[i] - yy0[i]);
      rms_err = sqrt(sum)/ (double)nn;
      printf("rms error=%f\n", rms_err);
      printf("Display signal deconvolved by spectral division:\n");
      sprintf(title, "Spectral division:  rms=%.4f", rms_err);
      dcv_1D_plot2(yyd, yy0, nn, title);
      }
/*--------------------------------
* Deconvolution by Wiener filter:
*/
  else if(!strncmp(option,"wiene",5)) {
      alpha=1.;
      dcv_1D_wiener(alpha);
      sum = 0.;
      for(i = 1; i <= nn; i++) sum += SQUARE(yyd[i] - yy0[i]);
      rms_err = sqrt(sum)/ (double)nn;
      printf("rms error=%f\n", rms_err);
      printf("Display Wiener deconvolved signal:\n");
      sprintf(title, "Wiener: alpha=%.1f:  rms=%.4f", alpha, rms_err);
      if(n_simu == 1)
        dcv_1D_plot2(yyd, yy0, nn, title);
      else
/* Bad visibility for the simulation #2, 
* so only the deconvolved signal is displayed
*/
        dcv_1D_plot1(yyd, nn, title);
      }
/*--------------------------------
* Deconvolution with Tikhonov's regularisation 
* Generalized Gauss with p=2
*/
  else if(!strncmp(option,"tikho",5)) {
      printf(" Enter alpha value: "); scanf("%lf",&alpha);
      dcv_1D_tikhonov(alpha, ftol, &iter, positive);
      sum = 0.;
      for(i = 1; i <= nn; i++) sum += SQUARE(yyd[i] - yy0[i]);
      rms_err = sqrt(sum)/ (double)nn;
      printf("rms error=%f\n", rms_err);
      printf("Display Tikhonov deconvolved signal:\n");
      sprintf(title, "Tikhonov: alpha=%.1f rms=%.4f iter=%d", 
              alpha, rms_err, iter);
      dcv_1D_plot2(yyd, yy0, nn, title);
      }
/*--------------------------------
* Deconvolution with generalized Gauss' regularisation 
* with p=1.1
*/
  else if(!strncmp(option,"ggaus",5)) {
      printf(" Enter alpha value: "); scanf("%lf",&alpha);
      dcv_1D_ggauss(alpha, ftol, &iter, positive);
      sum = 0.;
      for(i = 1; i <= nn; i++) sum += SQUARE(yyd[i] - yy0[i]);
      rms_err = sqrt(sum)/ (double)nn;
      printf("rms error=%f\n", rms_err);
      printf("Display generalized Gauss deconvolved signal:\n");
      sprintf(title, "Gen. Gauss: alpha=%.1f rms=%.4f it=%d", 
              alpha, rms_err, iter);
      dcv_1D_plot2(yyd, yy0, nn, title);
      }
/*--------------------------------
* Deconvolution with convex sqrt(s2+x2) regularisation 
*/
  else if(!strncmp(option,"sqrtr",5)) {
      printf(" Enter alpha value: "); scanf("%lf",&alpha);
      printf(" Enter ss value (0.1,0.001): "); scanf("%lf",&ss);
      dcv_1D_sqrtrf(alpha, ss, ftol, &iter, positive);
      sum = 0.;
      for(i = 1; i <= nn; i++) sum += SQUARE(yyd[i] - yy0[i]);
      rms_err = sqrt(sum)/ (double)nn;
      printf("rms error=%f\n", rms_err);
      printf("Display convex sqrt regul. deconvolved signal:\n");
      sprintf(title, "sqrt: s=%.3f alpha=%.1f rms=%.4f it=%d", 
              ss, alpha, rms_err, iter);
      dcv_1D_plot2(yyd, yy0, nn, title);
      }
/*--------------------------------
* Deconvolution with Maximum Entropy method 
*/
  else if(!strncmp(option,"maxen",5)) {
      printf(" Enter alpha value: "); scanf("%lf",&alpha);
      dcv_1D_mem(alpha, ftol, &iter, positive);
      sum = 0.;
      for(i = 1; i <= nn; i++) sum += SQUARE(yyd[i] - yy0[i]);
      rms_err = sqrt(sum)/ (double)nn;
      printf("rms error=%f\n", rms_err);
      printf("Display MEM deconvolved signal:\n");
      sprintf(title, "MEM: alpha=%.1f rms=%.4f it=%d", 
              alpha, rms_err, iter);
      dcv_1D_plot2(yyd, yy0, nn, title);
      }
/*--------------------------------
* Deconvolution with Gauss Markov regularisation (on x_i+1 - x_i) 
*/
  else if(!strncmp(option,"gmark",5)) {
      printf(" Enter alpha value: "); scanf("%lf",&alpha);
      dcv_1D_gmark(alpha, ftol, &iter, positive);
      sum = 0.;
      for(i = 1; i <= nn; i++) sum += SQUARE(yyd[i] - yy0[i]);
      rms_err = sqrt(sum)/ (double)nn;
      printf("rms error=%f\n", rms_err);
      printf("Display Gauss-Markov deconvolved signal:\n");
      sprintf(title, "Gauss-Markov: alpha=%.1f rms=%.4f it=%d", 
              alpha, rms_err, iter);
      dcv_1D_plot2(yyd, yy0, nn, title);
      }
/*---------------------------------
* Test with Banana function:
*/
  else if(!strncmp(option,"banan",5)) 
      dcv_1D_banana(ftol);
  else
      printf(" Fatal: invalid option\n");

dcv_1D_free();
free(hh);
return(0);
}
