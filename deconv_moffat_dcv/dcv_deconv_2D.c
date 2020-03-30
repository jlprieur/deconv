/****************************************************************** 
* dcv_deconv_2D
* Deconvolution program
*
*  spdiv => Direct spectral division
*  wiene => Wiener filter 
*  tikho => Tikhonov regularisation
*  maxen => Maximum Entropy method
*
* JLP
* Version 17/04/2006
*******************************************************************/ 
#include "jlp_dcv.h"

extern double *yy0, *yy, *yyd;

int dcv_remove_background(double *hh, int nx, int ny, 
                          double *backg, double *noise);
int dcv_2D_get_param(double *yy, int nx, int ny, double *alpha, 
                     double *ss, int *positivity, double *ftol,
                     double noise, char *option);

int main(int argc, char *argv[])
{
double *hh, rms_err, ss, ssum, alpha, ftol; 
double maxval, backg, noise, snr, L_phi, L_y;
float *float_tab;
int simu, iter, positivity;
INT_PNTR pntr;
INT4 nx, ny, nx0, ny0;
register int i;
char prefix[20], filename[60], comments[80];
char option[12];

printf(" Program to deconvolve 2D signals\n");
printf(" JLP version 23/04/2002\n");

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
sscanf(argv[1],"%d",&simu);
sscanf(argv[2],"%s",prefix);
}
else
{
printf(" Usage: dcv_deconv_2D simu prefix \n");
printf(" example: dcv_deconv_2D 1 s5 \n\n");
printf(" Simulation (simu=1): *_psf,*_yy0 are needed (prefix: s3, s5, ...)\n");
printf(" Real data (simu=0): *_psf,*_raw are needed\n");
exit(-1);
}

if(simu)
  printf(" OK: simulation\n");
else
  printf(" OK: real data\n");

JLP_INQUIFMT();

/* Read psf first to determine the size:
*/
sprintf(filename,"%s_psf",prefix);
JLP_VM_READIMAG1(&pntr, &nx, &ny, filename, comments);

/* Then allocate memory space: */
dcv_2D_init(nx, ny);

/**************** 
*  Read input files:
*****************/

/*----------------------
* Original signal:
*/
if(simu) {
   sprintf(filename,"%s_yy0",prefix);
   JLP_VM_READIMAG1(&pntr, &nx0, &ny0, filename, comments);
      if(ny != ny0 || nx != nx0) {
      printf(" Fatal error: incompatible size! (nx=%d,ny=%d) (nx0=%d, ny0=%d)\n", 
              nx, ny, nx0, ny0);
      exit(-1);
      }
   float_tab = (float *) pntr;
   TO_DOUBLE(float_tab, &yy0[1], &nx, &ny, &nx);
   } else {
/*----------------------
* Signal to deconvolve:
*/
   sprintf(filename,"%s_raw",prefix);
   JLP_VM_READIMAG1(&pntr, &nx0, &ny0, filename, comments);
      if(ny != ny0 || nx != nx0) {
      printf(" Fatal error: incompatible size! (nx=%d,ny=%d) (nx0=%d, ny0=%d)\n", 
              nx, ny, nx0, ny0);
      exit(-1);
      }
   float_tab = (float *) pntr;
   TO_DOUBLE(float_tab, &yy[1], &nx, &ny, &nx);
   }

/*----------------------
* PSF:
*/
sprintf(filename,"%s_psf",prefix);
JLP_VM_READIMAG1(&pntr, &nx, &ny, filename, comments);
float_tab = (float *) pntr;
if((hh = (double *)malloc((nx * ny + 1) * sizeof(double))) == NULL)
  {
  printf("fatal error allocating memory space for hh\n");
  exit(-1);
  }
TO_DOUBLE(float_tab, &hh[1], &nx, &ny, &nx);

/* Remove background from the psf: */
 dcv_remove_background(hh, nx, ny, &backg, &noise);

#ifdef DEBUG
      sprintf(filename,"%s_psf_flat",prefix);
      sprintf(comments,"psf with background removed (back=%.4e)", backg);
      JLP_D_WRITEIMAG(hh, &nx, &ny, &nx, filename, comments);
#endif

/* Normalize PSF (to avoid numerical pb...)
*/
maxval = hh[1];
for( i = 1; i <= nx * ny; i++) if(maxval < hh[i]) maxval = hh[i];
if(maxval == 0) { printf("Fatal error: PSF is null!\n"); exit(-1); }
for( i = 1; i <= nx * ny; i++) hh[i] /= maxval;


/* Compute transfer function:
*/
dcv_2D_transfert(hh);

/* Signal to deconvolve: */
if(simu) {
/* Either input it from a file */
    if(0) {
      sprintf(filename,"%s_yyb",prefix);
      JLP_VM_READIMAG1(&pntr, &nx0, &ny0, filename, comments);
        if(ny != ny0 || nx != nx0) {
        printf(" Fatal error: incompatible size! (nx=%d,ny=%d) (nx0=%d, ny0=%d)\n",
                nx, ny, nx0, ny0);
        exit(-1);
        }
      float_tab = (float *) pntr;
      TO_DOUBLE(float_tab, &yy[1], &nx, &ny, &nx);
    } else {
/* Or create it:
* snr is in dB
* snr=20 or 13dB
* snr=100 or 20dB : simu 3: with disk
* snr=1000 or 30dB: simu 4: with disk
* snr=20 or 13dB: simu 5: without disk
*/
      snr = 13.;
      dcv_2D_noisy_signal(snr);
      sprintf(filename,"%s_yyb",prefix);
      sprintf(comments,"%s: convolution of yy0 with psf, SNR=%.1f", prefix, snr);
      JLP_D_WRITEIMAG(yy, &nx, &ny, &nx, filename, comments);
    }
/* End of simulation case: */
}

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

/* Remove background from the input signal: */
 dcv_remove_background(yy, nx, ny, &backg, &noise);
#ifdef DEBUG
      sprintf(filename,"%s_raw_flat",prefix);
      sprintf(comments,"with background removed (back=%.4e)", backg);
      JLP_D_WRITEIMAG(yy, &nx, &ny, &nx, filename, comments);
#endif
/* Obtain parameters for minimization: */
   dcv_2D_get_param(yy, nx, ny, &alpha, &ss, &positivity, &ftol,
                    noise, option); 

L_phi = 0.;
/*--------------------------------
* Deconvolution by spectral division:
*/
  if(!strncmp(option,"spdiv",5)) {
      strcpy(comments," Spectral division");
      dcv_2D_spdiv();
      }
/*--------------------------------
* Deconvolution by Wiener filter:
*/
  else if(!strncmp(option,"wiene",5)) {
      dcv_2D_wiener(alpha);
      sprintf(comments," Wiener: alpha=%.1f",alpha);
      }
/*--------------------------------
* Deconvolution with Tikhonov's regularisation 
* Generalized Gauss with p=2
*/
  else if(!strncmp(option,"tikho",5)) {
      dcv_2D_tikhonov(alpha, positivity, ftol, &iter, &L_phi, &L_y);
      sprintf(comments," Wiener: alpha=%.1f, iter=%d", alpha, iter);
      }
/*--------------------------------
* Deconvolution with generalized Gauss' regularisation 
* with p=1.1
*/
  else if(!strncmp(option,"ggaus",5)) {
      dcv_2D_ggauss(alpha, positivity, ftol, &iter, &L_phi, &L_y);
      sprintf(comments," Gen. Gauss: alpha=%.1f, iter=%d", alpha, iter);
      }
/*--------------------------------
* Deconvolution with convex sqrt(s2+x2) regularisation 
* ss is the level which is a kind of (intensity) threshold 
* between "good signal"  and "medium" or "bad" signal 
* (penalty is less severe than tikhonov, for intensities > ss)
*/
  else if(!strncmp(option,"sqrtr",5)) {
      dcv_2D_sqrtrf(alpha, ss, positivity, ftol, &iter, &L_phi, &L_y);
      sprintf(comments," sqrtr: alpha=%.1f, ss=%f iter=%d", alpha, ss, iter);
      }
/*--------------------------------
* Deconvolution with Maximum Entropy method 
*/
  else if(!strncmp(option,"maxen",5)) {
      dcv_2D_mem(alpha, ftol, &iter, &L_phi, &L_y);
      sprintf(comments," MEM: alpha=%.1f, iter=%d", alpha, iter);
      }
/*--------------------------------
* Deconvolution with Gauss Markov regularisation (on x_i+1 - x_i) 
*/
  else if(!strncmp(option,"gmark",5)) {
      dcv_2D_gmark(alpha, positivity, ftol, &iter, &L_phi, &L_y);
      sprintf(comments," Gauss-Markov: alpha=%.1f, iter=%d", alpha, iter);
      }
  else {
      printf(" Fatal: invalid option\n");
      exit(-1);
      }

if(simu) {
      ssum = 0.;
      for(i = 1; i <= nx * ny; i++) ssum += SQUARE(yyd[i] - yy0[i]);
      rms_err = sqrt(ssum)/ (double)(nx*ny);
      printf("rms error=%f\n", rms_err);
      }

/*----------------- L Curve:
* L_phi = phi(x)
* L_y = || y - H x ||^2
*/
if(L_phi != 0.) 
   printf("L curve: L_phi=%.4e L_y=%.4e\n", L_phi, L_y);

sprintf(filename,"%s_yyd",prefix);
RECENT_FFT_DOUBLE(yyd, yyd, &nx, &ny, &nx);
sprintf(comments,"deconvolution with %s and %.1f", option, alpha);
JLP_D_WRITEIMAG(yyd, &nx, &ny, &nx, filename, comments);

dcv_2D_free();
free(hh);
return(0);
}
/**************************************************************
*
**************************************************************/
int dcv_remove_background(double *hh, int nx, int ny, 
                          double *backg, double *noise)
{
register int i;

/* in "auto_scale1.c": */
auto_sky_d(&hh[1], nx, ny, nx, backg, noise);

for(i = 1; i <= nx*ny; i++) hh[i] -= (*backg);

return(0);
}
/****************************************************************
* Routine to select parameters for deconvolution
*
* INPUT:
*  yy[]: signal array
*  nx, ny: size of yy[] 
*  noise: noise measured in the corners of the yy[] array
*
* OUTPUT:
* alpha: regularization parameter 
* ss: used by sqrtr 
* positivity: 1 if positivity constraint
* ftol: tolerance for minimization
****************************************************************/
int dcv_2D_get_param(double *yy, int nx, int ny, double *alpha, 
                     double *ss, int *positivity, double *ftol,
                     double noise, char *option) 
{
register int i;
double maxval, ssum, alpha0, norm_l2;

/* Look for maximum value:
*/
maxval = yy[1];
for( i = 1; i <= nx * ny; i++) if(maxval < yy[i]) maxval = yy[i];

 ssum = 0.;
 for(i = 1; i <= nx * ny; i++) ssum += SQUARE(yy[i]);
 norm_l2 = ssum / (double)(nx * ny);
 noise = MAXI(noise, 1.e-12);

 printf(" Maximum value: %.4e variance: %.4e, noise: %.4e, SNR: %.4e\n",
          maxval, norm_l2, noise, sqrt(norm_l2)/noise);

/* Input of ss needed for sqrtr: */
/* ss is the level which is a kind of (intensity) threshold 
* between "good signal"  and "medium" or "bad" signal 
* (penalty is less severe than tikhonov, for intensities > ss)
*/
   if(!strncmp(option, "sqrtr", 5)) {
      printf(" sqrt(s^2+x^2)/Enter ss value (0.1 or 0.001 for instance):\n"); 
      scanf("%lf", ss);
      }

 alpha0 = 1.;

 if(!strncmp(option, "tikho", 5)) 
    alpha0 = noise / norm_l2;
 else if(!strncmp(option, "ggaus", 5)) 
    alpha0 = noise / pow(norm_l2, 0.55);
 else if(!strncmp(option, "maxen", 5)) 
    alpha0 = noise / (2. * sqrt(norm_l2) * log(norm_l2));
 else if(!strncmp(option, "sqrtr", 5)) 
    alpha0 = noise / sqrt(SQUARE(*ss) + norm_l2);
 else if(!strncmp(option, "gmark", 5)) 
    alpha0 = noise / norm_l2;
/* For wiener: alpha is the ratio variance of noise / variance of signal 
*    alpha = pow(10.,-snr/(double)10.); for simulations where snr is in dB */
 else if(!strncmp(option, "wiene", 5))
    alpha0 = noise / norm_l2;

/* Now input of the other parameters: */
if(!strncmp(option, "tikho", 5) || !strncmp(option, "ggaus", 5) 
   || !strncmp(option, "sqrtr", 5) || !strncmp(option, "gmark", 5)
   || !strncmp(option, "maxen", 5)) { 
     printf(" tolerance ? [e.g., 1.e-5]\n"); 
     scanf("%lf", ftol);
     if(!strncmp(option, "maxen", 5)) 
       *positivity = 1; 
     else 
       {
       printf(" positivity constraint ? [1,0]\n"); 
       scanf("%d", positivity);
       }
     }

   if(!strncmp(option, "wiene", 5)) {
     printf("NB: Wiener filter is 1. / (power spectrum of transfer fct + alpha)\n");
     }

if(!strncmp(option, "tikho", 5) || !strncmp(option, "ggaus", 5) 
   || !strncmp(option, "sqrtr", 5) || !strncmp(option, "gmark", 5)
   || !strncmp(option, "maxen", 5) || !strncmp(option, "wiene", 5))
   {
    printf(" Value for alpha ? (indicative value: %.4e)\n", alpha0);
    scanf("%lf", alpha);
   }

return(0);
}
