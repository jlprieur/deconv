/*********************************************************************
* clean_deconv.c
* set of routines used for implementing the CLEAN deconvolution
* or independant program used to remove artefacts from bispectrum restored
* images...
*
* October 2008: I multiply both the spectrum of the input image
*               and that of the PSF by the uv_mask
*               in order to make them compatible (since
*               the PSF entered by the user is not the good one...) 
* 
* JLP
* Version 16/10/2008
*********************************************************************/
#include <stdio.h>
#include "jlp_ftoc.h"       /* "jlp_incl/jlp_ftoc.h" */
#include "jlp_fftw.h"       /* "jlp_incl/jlp_fftw.h" */

int clean_setup(double **dirty_map, double **dirty_beam, double **clean_map, 
                double **clean_beam, double **uv_mask, double clean_beam_fwhm, 
                char *in_name, char *transf_name, char *uv_mask_name,
                int *nx, int *ny, char *pt_option);
int clean_save_to_file(double *dirty_map, double *dirty_beam, 
                       double *clean_map, double *clean_beam, 
                       double omega, double noise, long n_iterations,
                       double clean_beam_fwhm, char *in_name, 
                       char *out_prefix, int nx, int ny);
int clean_deconv(double *dirty_map, double *dirty_beam, double *clean_map,
                 double *clean_beam, int nx, int ny, double omega, 
                 double noise, double *max_dirty_map, long *n_iterations);
int create_uv_mask(double **uv_mask, char *uv_mask_name, int nx, int ny);
int dirty_map_from_input_map(double **dirty_map, float *in_map, int nx, int ny);
int dirty_beam_from_transfer(double **dirty_beam, float *transf_fct,
                             double *uv_mask, int nx, int ny);
int dirty_beam_from_psf(double **dirty_beam, float *psf_fct, int nx, int ny);
int gaussian_clean_beam(double **clean_beam, int nx, int ny, double fwhm);
static int max_of_array(double *dirty_map, int nx, int ny, double *max_val,
                         int *i_max, int *j_max);
int filter_with_uv_mask(double *in_array, double *uv_mask, int nx, int ny);

#define MAIN_PROGRAM
#ifdef MAIN_PROGRAM
int main(int argc, char *argv[])
{
double *dirty_map, *dirty_beam, *clean_map, *clean_beam, *uv_mask; 
double omega, noise, max_dirty_map, clean_beam_fwhm;
int ival, status, nx, ny;
long n_iterations;
char in_name[60], out_prefix[40], transf_name[60], uv_mask_name[60]; 
char pt_option[1];

/* Solving the problem with "runs" which has always 7 arguments 
*/
if(argc == 7) {
 if(*argv[6]) argc = 7;
 else if(*argv[5]) argc = 6;
 else if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
 else argc = 1;
}
printf("DEBUG: argc= %d\n", argc);
if(argc != 5 && argc != 6) {
 fprintf(stderr,"********** clean_deconv/Bad syntax *********\nShould be:\n");
 fprintf(stderr,"  clean_deconv omega,noise,clean_beam_fwhm,P/T/S in_image out_prefix in_psf/in_transfer/in_power_spect [uv_mask]\n");
 fprintf(stderr,"Example:\n  clean_deconv 0.3,0.1,5,T im0 dec_im0 transfer0\n");
 fprintf(stderr,"or:  clean_deconv 0.3,0.1,5,P im0 dec_im0 psf0\n");
 fprintf(stderr,"omega = gain of the clean loops\n");
 fprintf(stderr,"noise = final level of the residuals as a fraction of the maximum of the input image\n");
 fprintf(stderr,"clean_beam_fwhm = fwhm of clean beam (in pixels)\n");
 fprintf(stderr,"P=PSF, T=Transfer function, S=square modulus\n");
 fprintf(stderr,"uv_mask: mask to be applied to Fourier domain (not necessary)\n");
 return(-1);
}

ival = sscanf(argv[1], "%lf,%lf,%lf,%c", &omega, &noise, &clean_beam_fwhm, pt_option);
if(ival != 4 || (*pt_option != 'P' && *pt_option != 'p' && *pt_option != 'S'
     && *pt_option != 's' && *pt_option != 'T' && *pt_option != 't')) {
    fprintf(stderr,"omega=%f noise=%f clean_beam_fwhm=%f PSF, Transfer or SquareModulus=%c\n", 
        omega, noise, clean_beam_fwhm, *pt_option);
    fprintf(stderr,"Fatal error: bad format for omega, noise, clean_beam_fwhm, pt_option\n");
    return(-1);
    }
strcpy(in_name, argv[2]);
strcpy(out_prefix, argv[3]);
strcpy(transf_name, argv[4]);
if(argc == 6) { 
  strcpy(uv_mask_name, argv[5]);
  printf("OK: uv_mask is >%s< \n", uv_mask_name);
  } else {
  uv_mask_name[0] = '\0';
  }

printf("OK: from >%s< to >%s< with >%s< \n", in_name, out_prefix, transf_name);
printf("OK: omega=%f noise=%f clean_beam_fwhm=%f PSF, Transfer or SquareModulus=%c\n", 
        omega, noise, clean_beam_fwhm, *pt_option);

JLP_INQUIFMT();

clean_beam = NULL;
clean_map = NULL;
dirty_beam = NULL;
dirty_map = NULL;

status = clean_setup(&dirty_map, &dirty_beam, &clean_map, &clean_beam, 
                     &uv_mask, clean_beam_fwhm, in_name, transf_name, 
                     uv_mask_name, &nx, &ny, pt_option);
if(!status) { 
  status = clean_deconv(dirty_map, dirty_beam, clean_map, clean_beam, nx, ny, 
                        omega, noise, &max_dirty_map, &n_iterations);
  printf("End of deconvolution: %ld iterations, max_dirty_map=%f\n",
          n_iterations, max_dirty_map);

  if(!status) clean_save_to_file(dirty_map, dirty_beam, clean_map, clean_beam, 
                                 omega, noise, n_iterations, clean_beam_fwhm, 
                                 in_name, out_prefix, nx, ny);
}

if(clean_beam != NULL) free(clean_beam);
if(dirty_beam != NULL) free(dirty_beam);
if(clean_map != NULL) free(clean_map);
if(dirty_map != NULL) free(dirty_map);
return(0);
}
#endif
/*************************************************************************
* Setup the arrays necessary for CLEAN deconvolution
*
* INPUT:
*  clean_beam_fwhm: Full Width at Half Maximum of the Gaussian to be used
*                   for the clean beam
*  pt_option: flag set to T if transfer function, S for Square modulus
*              or P if psf
*
* OUTPUT:  
*  dirty_map: dirty map (copy of input image) 
*  clean_map: clean map (filled with zeroes)
*  dirty_beam: dirty beam (centered on nx/2,ny/2)
*  clean_beam: clean beam (centered on nx/2,ny/2)
*  nx = number of lines
*  ny = number of pixels per line
*
*************************************************************************/
int clean_setup(double **dirty_map, double **dirty_beam, double **clean_map, 
                double **clean_beam, double **uv_mask, double clean_beam_fwhm, 
                char *in_name, char *transf_name, char *uv_mask_name,
                int *nx, int *ny, char *pt_option)
{
INT_PNTR pntr; 
INT4 nx1, ny1, nx2, ny2; 
float *in_map, *transf_fct;
char comments[80];
int status;
register int i;

/* Read the input image: */
  printf(" Input image to be cleaned\n");
  status = JLP_VM_READIMAG1(&pntr, &nx1, &ny1, in_name, comments);
  if(status) return(-1);
  in_map = (float *)pntr;
  *nx = nx1;
  *ny = ny1;

/* Create uv mask (from file or set it to unity): */
create_uv_mask(uv_mask, uv_mask_name, *nx, *ny);

printf("pt_option= %c\n", *pt_option);
/* Read the Point Spread Function or the transfer function: */
if(*pt_option == 'T' || *pt_option == 't') 
  printf(" Transfer function (centered in the frame)\n");
else if(*pt_option == 'S' || *pt_option == 's') 
  printf(" Square modulus of the transfer function (centered in the frame)\n");
else
  printf(" Point Spread Function (centered in the frame)\n");

  status = JLP_VM_READIMAG1(&pntr, &nx2, &ny2, transf_name, comments);
  if(status) {free(in_map); return(-1);}
  transf_fct = (float *)pntr;

if(nx1 != nx2 || ny1 != ny2){
  fprintf(stderr, "Fatal error: incompatible size between image and PSF/Transfer function!\n");
  fprintf(stderr, "Image: nx1,ny1 = %d,%d Tranfer function: nx2,ny2=%d,%d\n",
          nx1, ny1, nx2, ny2);
  free(in_map);
  free(transf_fct);
  return(-1);
  } 

if(*pt_option == 'S' || *pt_option == 's') {
  for(i = 0; i < (*nx)*(*ny); i++) {
     if(transf_fct[i] > 0) transf_fct[i] = sqrt(transf_fct[i]);
     else transf_fct[i] = 0.;
     }
 }

/* Compute dirty beam from PSF: */
if(*pt_option == 'P' || *pt_option == 'p') { 
  dirty_beam_from_psf(dirty_beam, transf_fct, *nx, *ny);
/* I multiply both the spectrum of the input image
*             and that of the PSF by the uv_mask
*             in order to make them compatible (since
*             the PSF entered by the user is not the good one...) 
*/
  if(*uv_mask_name) filter_with_uv_mask(*dirty_beam, *uv_mask, *nx, *ny);
} else {
/* Compute dirty beam from transfer function: */
  dirty_beam_from_transfer(dirty_beam, transf_fct, *uv_mask, *nx, *ny);
}

/* Compute clean beam: gaussian of FWHM equal to clean_beam_fwhm: */
gaussian_clean_beam(clean_beam, *nx, *ny, clean_beam_fwhm);

/* Transfer input image to residuals */
dirty_map_from_input_map(dirty_map, in_map, *nx, *ny);
if(*uv_mask_name) filter_with_uv_mask(*dirty_map, *uv_mask, *nx, *ny);

/* Clean map allocated and initialized to zero: */
*clean_map = (double *)malloc((*nx) * (*ny) * sizeof(double));
for(i = 0; i < (*nx) * (*ny); i++) (*clean_map)[i] = 0.;

free(in_map);
free(transf_fct);
return(0);
}
/********************************************************************
* Save results to FITS files
*
********************************************************************/
int clean_save_to_file(double *dirty_map, double *dirty_beam, 
                       double *clean_map, double *clean_beam, 
                       double omega, double noise, long n_iterations,
                       double clean_beam_fwhm, char *in_name, 
                       char *out_prefix, int nx, int ny)
{
char out_name[80], comments[80];
register int i;

/* Clean map is saved to file: */
sprintf(out_name,"%s_cm",out_prefix);
sprintf(comments,"Clean map of %.20s (omega=%5.3f noise=%g it=%ld)", in_name,
        omega, noise, n_iterations);
JLP_D_WRITEIMAG(clean_map, &nx, &ny, &nx, out_name, comments);

/* Clean map + dirty map is saved to file: */
for(i = 0; i < nx * ny; i++) clean_map[i] += dirty_map[i];
sprintf(out_name,"%s_cleaned",out_prefix);
sprintf(comments,"Clean + Dirty map of %.20s (omega=%5.3f noise=%g it=%ld)", 
        in_name, omega, noise, n_iterations);
JLP_D_WRITEIMAG(clean_map, &nx, &ny, &nx, out_name, comments);

/* Dirty map is saved to file: */
sprintf(out_name,"%s_dm",out_prefix);
sprintf(comments,"Dirty map of %.20s (omega=%5.3f noise=%g it=%ld)", 
        in_name, omega, noise, n_iterations);
JLP_D_WRITEIMAG(dirty_map, &nx, &ny, &nx, out_name, comments);

/* Dirty beam is saved to file: */
sprintf(out_name,"%s_db",out_prefix);
sprintf(comments,"Dirty beam of %.30s",in_name);
JLP_D_WRITEIMAG(dirty_beam, &nx, &ny, &nx, out_name, comments);

/* Clean beam is saved to file: */
sprintf(out_name,"%s_cb",out_prefix);
sprintf(comments,"Gaussian Clean beam with fwhm=%.3f pixels", clean_beam_fwhm);
JLP_D_WRITEIMAG(clean_beam, &nx, &ny, &nx, out_name, comments);

return(0);
}
/***********************************************************************
* Compute dirty map from input map
* INPUT:
* in_map: input image to be cleaned
* nx, ny: size of input/output arrays
*
* OUTPUT:
* dirty_map: first value of the residuals 
***********************************************************************/
int dirty_map_from_input_map(double **dirty_map, float *in_map, int nx, int ny)
{
register int i;

*dirty_map = (double *)malloc(nx * ny * sizeof(double));

for(i = 0; i < nx * ny; i++) (*dirty_map)[i] = in_map[i];

return(0);
}
/***********************************************************************
* Compute clean beam: gaussian of FWHM equal to clean_beam_fwhm: 
*
* INPUT:
* nx, ny: size of input/output arrays
* fwhm: Full Width at Half Maximum of the Gaussian
*
* OUTPUT:
* clean_beam: dirty beam (centered on nx/2,ny/2, 
*             with a dynamical range of 1/EPSILON)
***********************************************************************/
int gaussian_clean_beam(double **clean_beam, int nx, int ny, double fwhm)
{
int nx2, ny2;
double EPSILON = 1.e-3, ww, sigma2;
register int i, j;

*clean_beam = (double *)malloc(nx * ny * sizeof(double));

/* FWHM of a Gaussian of sigma=1 is 2*x with: 
* exp(- x^2 / sigma^2) = 0.5:
* Hence x^2 = 0.693 and 2*x = 1.665
*/
sigma2 = SQUARE(fwhm / 1.665);
nx2 = nx/2;
ny2 = ny/2;
for(j = 0; j < ny; j++) {
  for(i = 0; i < nx; i++) {
   ww = SQUARE(i - nx2) + SQUARE(j - ny2);
   ww = exp(- ww / sigma2);
/* Limit the dynamical range to 1/EPSILON */
   if(ww > EPSILON) (*clean_beam)[i + j * nx] = ww; 
   else (*clean_beam)[i + j * nx] = 0.;
  }
}
return(0);
}
/***********************************************************************
* Compute dirty beam from transfer function
*
* I multiply both the spectrum of the input image
* and the transfer function by the uv_mask
* in order to make them compatible (since
* the transfer function entered by the user is not the good one...) 
*
* INPUT:
* transf_fct: transfer function (centered on nx/2,ny/2)
* uv_mask: mask to be applied on the Fourier domain
* nx, ny: size of input/output arrays
*
* OUTPUT:
* dirty_beam: dirty beam (centered on nx/2,ny/2)
***********************************************************************/
int dirty_beam_from_transfer(double **dirty_beam, float *transf_fct,
                             double *uv_mask, int nx, int ny)
{
double *re, *im, ww;
register int i;

*dirty_beam = (double *)malloc(nx * ny * sizeof(double));
for(i = 0; i < nx * ny; i++) (*dirty_beam)[i] = 0.;

re = (double *)malloc(nx * ny * sizeof(double));
im = (double *)malloc(nx * ny * sizeof(double));

/* Multiplication of the transfer function with the uv_mask: */
for(i = 0; i < nx * ny; i++) re[i] = transf_fct[i] * uv_mask[i];
for(i = 0; i < nx * ny; i++) im[i] = 0.;

RECENT_FFT_DOUBLE(re,re,&nx,&ny,&nx);

/*
int fftw_2D_double(double *re, double *im, int nx, int ny, int direct);
*/
fftw_2D_double(re, im, nx, ny, 1);

/* Normalisation to one: */
ww = re[0];
if(ww == 0.) {
   fprintf(stderr, "dirty_beam_from_transfer/Error re[0) is null \n");
   return(-1);
   }
for(i = 0; i < nx * ny; i++) re[i] /= ww;

RECENT_FFT_DOUBLE(re, re, &nx, &ny, &nx);

/* Assuming that transfer function is real and symmetric
* its Fourier transform is real */
for(i = 0; i < nx * ny; i++) (*dirty_beam)[i] = re[i];

free(re);
free(im);
return(0);
}
/***********************************************************************
* Compute dirty beam from the Point Spread Function
*
* INPUT:
* psf_fct: PSFunction (centered on nx/2,ny/2)
* uv_mask: mask to be applied on the Fourier domain
* nx, ny: size of input/output arrays
*
* OUTPUT:
* dirty_beam: dirty beam (centered on nx/2,ny/2)
***********************************************************************/
int dirty_beam_from_psf(double **dirty_beam, float *psf_fct, int nx, int ny)
{
double max_val;
int i_max, j_max;
register int i;

*dirty_beam = (double *)malloc(nx * ny * sizeof(double));
for(i = 0; i < nx * ny; i++) (*dirty_beam)[i] = psf_fct[i];

/* Check if maximum is positive: */
max_of_array(*dirty_beam, nx, ny, &max_val, &i_max, &j_max);

/* Normalisation to unity: */
if(max_val <= 0.) {
   fprintf(stderr, "dirty_beam_from_PSF/Error, central value is not positive (max=%f)\n", max_val);
   exit(-1);
   }
for(i = 0; i < nx * ny; i++) (*dirty_beam)[i] /= max_val;

return(0);
}
/*********************************************************************
* clean_deconv
* Clean algorithm used for deconvolution.
* Derived from a routine of Eric ANTERRIEU, Lib_an11.c Version July 1990
*
* INPUT:  
*  dirty_map[nx*ny] = dirty map 
*  dirty_beam[nx*ny] = dirty beam
*  clean_beam[nx*ny] = clean beam
*  nx = number of lines
*  ny = number of pixels per line
*  omega = loop gain (percentage of the peak to be removed at each iteration)
*  noise = noise (percentage of the maximum value of the dirty map)
*
* OUTPUT: 
*  clean_map[nx*ny] = clean map
*  dirty_map[nx*ny] = residual map
*  max_dirty_map : maximum of the dirty map (residuals)
*  n_iterations = number of iterations
*
***************************************************************/
int clean_deconv(double *dirty_map, double *dirty_beam, double *clean_map,
                 double *clean_beam, int nx, int ny, double omega, 
                 double noise, double *max_dirty_map, long *n_iterations)
{
int    nx2, ny2, dm_xmax, dm_ymax, db_xmax, db_ymax, cb_xmax, cb_ymax;
double  max_val, threshold, old_max_val;
register  int  i, j, ixcm, iycm, ixdm, iydm;

nx2 = nx/2;  ny2 = ny/2;

for(i = 0; i < nx * ny; i++) clean_map [i] = 0.0;

/* Look for the maximum of the dirty beam: */
max_of_array(dirty_beam, nx, ny, &max_val, &db_xmax, &db_ymax);
if(max_val > 0) {
 for(i = 0; i < nx * ny; i++) dirty_beam[i] /= max_val;
 } else {
 fprintf(stderr, "clean_deconv/fatal error: max val of dirty beam = %f\n",
         max_val);
 }

/* Look for the maximum of the clean beam: */
max_of_array(clean_beam, nx, ny, &max_val, &cb_xmax, &cb_ymax);
if(max_val > 0) {
 for(i = 0; i < nx * ny; i++) clean_beam[i] /= max_val;
 } else {
 fprintf(stderr, "clean_deconv/fatal error: max val of clean beam = %f\n",
         max_val);
 }

/* Look for the maximum of the map: */
max_of_array(dirty_map, nx, ny, &max_val, &dm_xmax, &dm_ymax);

if(max_val <=0) {
  fprintf(stderr,
          "clean_deconv/error: maximum of residuals is negative: max=%f\n",
          max_val);
  *max_dirty_map = max_val;
  return(-1);
  }

/******************** MAIN LOOP ****************************************/
*n_iterations = 0;
threshold = noise * max_val;
printf("clean_deconv: omega=%f threshold=%f\n", omega, threshold);
printf("Dirty beam: ixc,iyc=%d,%d Clean beam: ixc, iyc=%d,%d\n", 
        db_xmax, db_ymax, cb_xmax, cb_ymax);
while (max_val >= threshold) {
  (*n_iterations)++;
  ixdm = dm_xmax - db_xmax;
  iydm = dm_ymax - db_ymax;
  ixcm = dm_xmax - cb_xmax;
  iycm = dm_ymax - cb_ymax;
/* DEBUG:
printf("clean_deconv: max_val=%f omega=%f threshold=%f\n", max_val, omega, threshold);
*/
/* Removes max_val * omega * clean_beam  centered at (i_max, j_max)
* from the dirty map (residuals) and update the clean map: */
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++) {
      if ( ((i-ixdm) >= 0) && ((j-iydm) >= 0)
        && ((i-ixdm) < nx) && ((j-iydm) < ny)
        && ((i-ixcm) >= 0) && ((j-iycm) >= 0)
        && ((i-ixcm) < nx) && ((j-iycm) < ny) ) {
        clean_map[i + j * nx] += (max_val * omega 
                                  * clean_beam[(i-ixcm) + (j - iycm) * nx]);
        dirty_map[i + j * nx] -= (max_val * omega
                                  * dirty_beam[(i-ixdm) + (j - iydm) * nx]);
        }
      }
    }
  old_max_val = max_val;
  max_of_array(dirty_map, nx, ny, &max_val, &dm_xmax, &dm_ymax);
/* Protection against bad dirty maps...
  if(max_val > old_max_val) {
    fprintf(stderr,"Error: residuals are not decreasing! (max_val=%f, old_max=%f iter=%ld)\n", 
            max_val, old_max_val, *n_iterations);
    break;
  }
*/
} /* End of while loop */

*max_dirty_map = max_val;

return(0);
}
/**********************************************************************
* Look for the maximum value in a 2D array
*
* INPUT:
*  dirty_map: input array
*  nx, ny: size of input array
*
* OUTPUT:
*  max_val: maximum value found in the input array
*  i_max, j_max: x, y location of the maximum 
**********************************************************************/
static int max_of_array(double *array, int nx, int ny, double *max_val,
                         int *i_max, int *j_max)
{
register int i, j;

*max_val = array[0];
*i_max = 0;
*j_max = 0;
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++) {
      if (array[i + j * nx] > (*max_val)) {
        *max_val = array[i + j * nx];
        *i_max = i;
        *j_max = j;
        }
      }
    }
return(0);
}
/***********************************************************************
* Create uv mask from file if available
* otherwise create a mask set to unity
*
***********************************************************************/
int create_uv_mask(double **uv_mask, char *uv_mask_name, int nx, int ny)
{
INT_PNTR pntr;
INT4 nx1, ny1;
int i_max, j_max, status;
double max_val;
float *uv_mask_flt;
char comments[80];
register int i;

/* uv_mask allocated and initialized to unity: */
*uv_mask = (double *)malloc(nx * ny * sizeof(double));
for(i = 0; i < nx * ny; i++) (*uv_mask)[i] = 1.;

/* Read uv mask if any: */
if(*uv_mask_name) {
  status = JLP_VM_READIMAG1(&pntr, &nx1, &ny1, uv_mask_name, comments);
  if(status) return(-1);
  uv_mask_flt = (float *)pntr;

  if(nx1 != nx || ny1 != ny){
    fprintf(stderr, "create_uv_mask/Error: incompatible size between input image and uv_mask!\n");
    fprintf(stderr, "Image: nx,ny = %d,%d uv_mask: nx1,ny1=%d,%d\n",
            nx, ny, nx1, ny1);
    free(uv_mask_flt);
    return(-1);
    }
/* Transfer to uv_mask array: */
  for(i = 0; i < nx * ny; i++) (*uv_mask)[i] = uv_mask_flt[i];

/* Free memory: */
  free(uv_mask_flt);

/* Check if maximum is one: */
   max_of_array(*uv_mask, nx, ny, &max_val, &i_max, &j_max);

/* If maximum is not one, exit from here: 
*/
  if(max_val != 1.) {
   fprintf(stderr, "create_uv_mask/Fatal error: maximum of mask is %f (!= 1.0)\n",
            max_val);
   exit(-1);
   }
}

return(0);
}
/****************************************************************************
*
* INPUT:
* in_out_array: input array to be filtered
* uv_mask: mask in Fourier domain
*
* OUTPUT:
* in_out_array: array filtered in Fourier domain with uv_mask
****************************************************************************/
int filter_with_uv_mask(double *in_out_array, double *uv_mask, int nx, int ny)
{
double *re, *im;
register int i;

re = (double *)malloc(nx * ny * sizeof(double));
im = (double *)malloc(nx * ny * sizeof(double));

/* Copy input image to real part: */
  for(i = 0; i < nx * ny; i++) re[i] = in_out_array[i];
  for(i = 0; i < nx * ny; i++) im[i] = 0.;

/* Shift image to avoid artefacts due to non-continuities at the edges */
RECENT_FFT_DOUBLE(re,re,&nx,&ny,&nx);

/*
int fftw_2D_double(double *re, double *im, int nx, int ny, int direct);
*/
fftw_2D_double(re, im, nx, ny, 1);

RECENT_FFT_DOUBLE(re,re,&nx,&ny,&nx);
RECENT_FFT_DOUBLE(im,im,&nx,&ny,&nx);

/* Multiplication of the transfer function with the uv_mask: */
  for(i = 0; i < nx * ny; i++) re[i] *= uv_mask[i];
  for(i = 0; i < nx * ny; i++) im[i] *= uv_mask[i];

RECENT_FFT_DOUBLE(re,re,&nx,&ny,&nx);
RECENT_FFT_DOUBLE(im,im,&nx,&ny,&nx);

fftw_2D_double(re, im, nx, ny, -1);

RECENT_FFT_DOUBLE(re,re,&nx,&ny,&nx);

/* Copy real part to output image: */
  for(i = 0; i < nx * ny; i++) in_out_array[i] = re[i];


free(re);
free(im);
return(0);
}
