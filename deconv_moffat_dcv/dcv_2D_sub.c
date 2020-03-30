/****************************************************************** 
* Routines used by dcv_deconv_2D.c
*
*  spdiv => Direct spectral division
*  wiene => Wiener filter 
*  tikho => Tikhonov regularisation
*  maxen => Maximum Entropy method
*
* JLP
* Version 22/04/2002
*******************************************************************/
#include <jlp_dcv.h>
#include "jlp_fftw.h"

/* Public: */
double *yy, *yy0, *yyd;

/* Private: */
/*******************
* ss2 is a parameter needed for convex sqrt regularization
*/
static double ss2;
static int nx, ny, positivity;
static double *w0, *w1, *hh_re, *hh_im;
static double *hht_re, *hht_im;
static double *w_re, *w_im;
static double *yy_re, *yy_im;
static double alpha;
static char option[12];
/* gscale: used for calibrating constants (e.g. for Maximum Entropy) */
static double gscale=1.e-1;
/* Complex array (for fast FFT's) */
static FFTW_COMPLEX *cp0;

/* Public:  (declared in jlp_dcv.h)
int dcv_2D_gmark(double alpha_0, int positivity_0, double ftol, int *iter, 
                 double *L_phi, double *L_y);
int dcv_2D_tikhonov(double alpha_0, int positivity_0, double ftol, int *iter,
                    double *L_phi, double *L_y);
int dcv_2D_ggauss(double alpha_0, int positivity_0, double ftol, int *iter,
                  double *L_phi, double *L_y);
int dcv_2D_sqrtrf(double alpha_0, double ss, int positivity_0, double ftol, 
                  int *iter, double *L_phi, double *L_y);
int dcv_2D_mem(double alpha_0, double ftol, int *iter,
               double *L_phi, double *L_y);
int dcv_2D_spdiv();
int dcv_2D_wiener(double alpha_0);
int dcv_2D_init(int nx_0, int ny_0);
void dcv_2D_free();
int dcv_2D_noisy_signal(double snr);
int dcv_2D_transfert(double *hh);
*/


/* Private functions defined here: */
static double func_2D(double *xx);
static void dfunc_2D(double *xx, double *dx);
static int dcv_2D_check_grad(double *x1, double (*func)(double []), 
                             void (*dfunc)(double[], double[]));
static int dcv_2D_conv_hh(double *xx, double *hh_re, double *hh_im, double *ww);
static double phi_2D_gmark(double *xx);
static void dphi_2D_gmark(double *xx, double *dx);
static double phi_2D_mem(double *xx);
static void dphi_2D_mem(double *xx, double *dx);
static double phi_2D_tikho(double *xx);
static void dphi_2D_tikho(double *xx, double *dx);
static double phi_2D_ggauss(double *xx);
static void dphi_2D_ggauss(double *xx, double *dx);
static double phi_2D_sqrtrf(double *xx);
static void dphi_2D_sqrtrf(double *xx, double *dx);

/**************************************************************
* Function E to be minimized for deconvolution
*
* This subroutine returns   E = |Y-Ax|^2 + alpha  phi
**************************************************************/
static double func_2D(double *xx)
{
double sum2_y, phi;
register int i;

/* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_2D_conv_hh(xx, hh_re, hh_im, w0);

/* |Y-Ax|^2 */
sum2_y = 0.;
for(i = 1; i <= nx * ny; i++) sum2_y += SQUARE(w0[i] - yy[i]);

if(!strncmp(option,"tikho",5)) 
  phi = phi_2D_tikho(xx); 
else if(!strncmp(option,"maxen",5))
  phi = phi_2D_mem(xx); 
else if(!strncmp(option,"ggaus",5))
  phi = phi_2D_ggauss(xx); 
else if(!strncmp(option,"gmark",5))
  phi = phi_2D_gmark(xx); 
else if(!strncmp(option,"sqrtr",5))
  phi = phi_2D_sqrtrf(xx); 
else 
  {
  printf("func_2D/Fatal error: invalid option\n");
  exit(-1);
  }

/* Penalty for negative points: */
if(positivity) {
  for(i = 1; i <= nx * ny; i++) if(xx[i] < 0.) phi += 1000. * gscale;
  }

return(sum2_y + alpha * phi);
}
/**************************************************************
* Gradient of potential function for minimisation:
* grad_Er : value of the gradient of the criterium in X
*                                       (                 \---
* This subroutine returns  grad(E) =grad( |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
*                                       (                 /___
*                                                      (s,r) in C
*
* d/dx (HX - Y)^T (HX - Y) = 2 H^T (HX - Y)
**************************************************************/
static void dfunc_2D(double *xx, double *dx)
{
register int i;

/* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_2D_conv_hh(xx, hh_re, hh_im, w0);

for(i = 1; i <= nx * ny; i++) w0[i] -= yy[i];

/* w0 = conv(w0,hht):
*/
dcv_2D_conv_hh(w0, hht_re, hht_im, w0);

/* Now compute gradient dx: */
if(!strncmp(option,"tikho",5))
  dphi_2D_tikho(xx, dx);
else if(!strncmp(option,"maxen",5))
  dphi_2D_mem(xx, dx);
else if(!strncmp(option,"ggaus",5))
  dphi_2D_ggauss(xx, dx);
else if(!strncmp(option,"gmark",5))
  dphi_2D_gmark(xx, dx);
else if(!strncmp(option,"sqrtr",5))
  dphi_2D_sqrtrf(xx, dx);
else
  {
  printf("dfunc_2D/Fatal error: invalid option\n");
  exit(-1);
  }

if(positivity)
   {
   for(i = 1; i <= nx * ny; i++) 
       if(xx[i] < 0.) dx[i] = -1000. * gscale; 
   }
for(i = 1; i <= nx * ny; i++) dx[i] = 2. * w0[i] + alpha * dx[i]; 

}
/*************************************************************
* Deconvolution with Gauss-Markov regularisation 
* Criterium to minimize is:
* || y - H x ||^2 + alpha phi
*
* with:
*        _
*        \                     p                      p
* phi =  / | x        - x     |  + | x       - x     |
*        -    (i+1,j)    (i,j)        (i,j+1)   (i,j)
*        i
* (Smooth function, which reduces the variation between two successive pixels)
*
* INPUT:
* ftol: tolerance
**************************************************************/
int dcv_2D_gmark(double alpha_0, int positivity_0, double ftol, int *iter,
                 double *L_phi, double *L_y)
{
double fret;
int nn, iter_max = 2000, lbfgs = 0;
register int i;

/* Initialize static parameters: */
alpha = alpha_0;
positivity = positivity_0;
strcpy(option, "gmark");

/* Check if gradient is OK */
dcv_2D_check_grad(yyd, phi_2D_gmark, dphi_2D_gmark);

/* Starting point: */
for(i = 1; i <= nx * ny; i++) yyd[i] = 0.;

  nn = nx * ny;
  jlp_minimize(yyd, nn, nn, ftol, iter, iter_max, &fret, func_2D, dfunc_2D, 
               0, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* L-curve: */
*L_phi = phi_2D_gmark(yyd);
*L_y = func_2D(yyd) - alpha * (*L_phi);

return(0);
}
/**************************************************************
* Regularisation function to be minimized for Gauss-Markov
* Bauman & Sauer' function with p=2:
*        _
*        \                     p                      p
* phi =  / | x        - x     |  + | x       - x     |
*        -    (i+1,j)    (i,j)        (i,j+1)   (i,j)
*        i
*
* (Smooth function, which reduces the variation between two successive pixels)
*
**************************************************************/
static double phi_2D_gmark(double *xx)
{
double ssum, dx, dy;
register int i, j;

/* Central part: */
for(j = 1; j <= ny-1; j++) { 
   for(i = 1; i <= nx-1; i++) { 
     dx = xx[i+1 + (j-1) * nx] - xx[i + (j-1) * nx]; 
     dy = xx[i + j * nx] - xx[i + (j-1) * nx]; 
     ssum += SQUARE(dx) + SQUARE(dy);
     }
   } 

/* Wrap the image for the edges:
* Last line
*/
   for(i = 1; i <= nx-1; i++) { 
     dx = xx[i+1 + (ny-1) * nx] - xx[i + (ny-1) * nx]; 
     dy = xx[i] - xx[i + (ny-1) * nx]; 
     ssum += SQUARE(dx) + SQUARE(dy);
     }

/* Last column: */
for(j = 1; j <= ny-1; j++) { 
     dx = xx[1 + (j-1) * nx] - xx[nx + (j-1) * nx]; 
     dy = xx[nx + j * nx] - xx[nx + (j-1) * nx]; 
     ssum += SQUARE(dx) + SQUARE(dy);
   } 

/* Last pixel: */
     dx = xx[1 + (ny-1) * nx] - xx[nx + (ny-1) * nx]; 
     dy = xx[nx] - xx[nx + (ny-1) * nx]; 
     ssum += SQUARE(dx) + SQUARE(dy);

return(ssum);
}
/**************************************************************
* Gradient of potential function for Gauss-Markov
* Bauman & Sauer' function with p=2:
*        _
*        \                     p                      p
* phi =  / | x        - x     |  + | x       - x     |
*        -    (i+1,j)    (i,j)        (i,j+1)   (i,j)
*        i
* (Smooth function, which reduces the variation between two successive pixels)
*
**************************************************************/
static void dphi_2D_gmark(double *xx, double *dx)
{
register int i, j;

/* Gradient is linear, so
* First go along the lines
*/

for(j = 1; j <= ny; j++)  
   for(i = 2; i <= nx-1; i++)  
       dx[i + (j-1) * nx] = - 2. * ( xx[i+1 + (j-1) * nx] - xx[i + (j-1) * nx] ) 
                            + 2. * ( xx[i + (j-1) * nx] - xx[i-1 + (j-1) * nx] );

/* First and last pixels of each line (wrapping up the lines):
*/
for(j = 1; j <= ny; j++)  
   {
    dx[1 + (j-1) * nx] = - 2. * ( xx[2 + (j-1) * nx] - xx[1 + (j-1) * nx] ) 
                         + 2. * ( xx[1 + (j-1) * nx] - xx[nx + (j-1) * nx] );
    dx[nx + (j-1) * nx] = - 2. * ( xx[1 + (j-1) * nx] - xx[nx + (j-1) * nx] ) 
                         + 2. * ( xx[nx + (j-1) * nx] - xx[nx-1 + (j-1) * nx] );
   }


/* Then go along the columns: 
*/
for(j = 2; j <= ny-1; j++)  
   for(i = 1; i <= nx; i++)  
       dx[i + (j-1) * nx] += - 2. * ( xx[i + j * nx] - xx[i + (j-1) * nx] ) 
                            + 2. * ( xx[i + (j-1) * nx] - xx[i + (j-2) * nx] );

/* First and last pixels of each column (wrapping up the columns):
*/
for(i = 1; i <= nx; i++)  
   {
    dx[i] += - 2. * ( xx[i + nx] - xx[i] ) 
             + 2. * ( xx[i] - xx[i + (ny-1) * nx] );
    dx[i + (ny-1) * nx] += - 2. * ( xx[i] - xx[i + (ny-1) * nx] ) 
                         + 2. * ( xx[i + (ny-1) * nx] - xx[i + (ny-2) * nx] );
   }

}
/**************************************************************
* Deconvolution with MEM 
* Criterium to minimize is:
* || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^2
**************************************************************/
int dcv_2D_mem(double alpha_0, double ftol, int *iter,
               double *L_phi, double *L_y)
{
double fret;
int nn, iter_max = 2000, lbfgs = 0;
register int i;

/* Initialize static parameters: */
alpha = alpha_0;
strcpy(option, "maxen");
positivity = 1;

/* Check if gradient is OK */
dcv_2D_check_grad(yyd, phi_2D_mem, dphi_2D_mem);

/* Starting point: */
for(i = 1; i <= nx * ny; i++) yyd[i] = 0.;

  nn = nx * ny;
  jlp_minimize(yyd, nn, nn, ftol, iter, iter_max, &fret, func_2D, dfunc_2D, 
               0, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* L-curve: */
*L_phi = phi_2D_mem(yyd);
*L_y = func_2D(yyd) - alpha * (*L_phi);

return(0);
}
/**************************************************************
* Entropy potential function:
*        _
*        \     
* phi =  /   x log(x )
*        -    i     i
*        i
*
* Penalty of 1000. if negative...
**************************************************************/
static double phi_2D_mem(double *xx)
{
register int i;
double ssum;

ssum = 0.;
for(i = 1; i <= nx * ny; i++) 
  {
  if(xx[i] > 0.)
    ssum += xx[i] * log(xx[i]);
  else
    ssum += 1000. * gscale; 
  }

return(ssum);
}
/**************************************************************
* Gradient for MEM 
* grad_Er : value of the gradient of the criterium in X
*                                     (                 \---
* This function returns  grad(E) =grad( |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
*                                     (                 /___
*                                                      (s,r) in C
*
* d/dx (HX - Y)^T (HX - Y) = 2 H^T (HX - Y)
**************************************************************/
static void dphi_2D_mem(double *xx, double *dx)
{
register int i;

/*
* Gradient de la fonction potentiel entropique
*        _
*        \     
* phi =  /   x log(x )
*        -    i     i
* => dphi = sum(1 + log(x)):
*/
for(i = 1; i <= nx * ny; i++) 
   {
   if(xx[i] > 0.) 
     dx[i] = 1. + log(xx[i]);
   else
     dx[i] = -1000. * gscale;
   }

}
/**************************************************************
* Deconvolution by Tikhonov's regularisation 
* Criterium to minimize is:
* || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^2
* (Carfantan: grf)
**************************************************************/
int dcv_2D_tikhonov(double alpha_0, int positivity_0, double ftol, int *iter,
                    double *L_phi, double *L_y)
{
double fret;
int nn, iter_max = 2000, lbfgs = 0;
register int i;

/* Initialize static parameters: */
alpha = alpha_0;
positivity = positivity_0;
strcpy(option, "tikho");

/* Check if gradient is OK */
dcv_2D_check_grad(yyd, phi_2D_tikho, dphi_2D_tikho);

/* Starting point: */
for(i = 1; i <= nx * ny; i++) yyd[i] = 0.;

  nn = nx * ny;
  jlp_minimize(yyd, nn, nn, ftol, iter, iter_max, &fret, func_2D, dfunc_2D, 
               0, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* L-curve: */
*L_phi = phi_2D_tikho(yyd);
*L_y = func_2D(yyd) - alpha * (*L_phi);

return(0);
}
/**************************************************************
* Regularisation function to be minimized for Tikhonov
* Bauman & Sauer's potential function, with p=2
*        _
*        \        p
* phi =  /  | x  |
*        -     i
*        i
*  phi = sum(abs(real(x(:))).^p);
*  dphi = p*sign(x).*(abs(real(x)).^(p-1));
**************************************************************/
static double phi_2D_tikho(double *xx)
{
register int i;
double ssum;

ssum = 0.;
for(i = 1; i <= nx * ny; i++) ssum += SQUARE(xx[i]); 

return(ssum); 
}
/**************************************************************
* Gradient of potential function for Tikhonov's regularisation:
* Bauman & Sauer's potential function, with p=2
*        _
*        \        p
* phi =  /  | x  |
*        -     i
*        i
*  dphi = p*sign(x).*(abs(real(x)).^(p-1));
* here: dphi = 2 x
**************************************************************/
static void dphi_2D_tikho(double *xx, double *dx)
{
register int i;

for(i = 1; i <= nx * ny; i++) dx[i] = 2. * xx[i];

}
/**************************************************************
* Deconvolution by generalized Gauss' regularisation 
* Criterium to minimize is:
* || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^p
* with p=1.1
* (Carfantan: ggrf)
**************************************************************/
int dcv_2D_ggauss(double alpha_0, int positivity_0, double ftol, int *iter,
                  double *L_phi, double *L_y)
{
double fret;
int nn, iter_max = 2000, lbfgs = 0;
register int i;

/* Initialize static parameters: */
alpha = alpha_0;
positivity = positivity_0;
strcpy(option, "ggaus");

/* Check if gradient is OK */
dcv_2D_check_grad(yyd, phi_2D_ggauss, dphi_2D_ggauss);

/* Starting point: */
for(i = 1; i <= nx * ny; i++) yyd[i] = 0.;

  nn = nx * ny;
  jlp_minimize(yyd, nn, nn, ftol, iter, iter_max, &fret, func_2D, dfunc_2D, 
               0, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* L-curve: */
*L_phi = phi_2D_ggauss(yyd);
*L_y = func_2D(yyd) - alpha * (*L_phi);

return(0);
}
/**************************************************************
* Regularisation function to be minimized for generalized Gauss
* Bauman & Sauer's potential function, with p=1.1
*        _
*        \        p
* phi =  /  | x  |
*        -     i
*        i
*  phi = sum(abs(real(x(:))).^p);
*  dphi = p*sign(x).*(abs(real(x)).^(p-1));
*
**************************************************************/
static double phi_2D_ggauss(double *xx)
{
register int i;
double ssum;

ssum = 0.;
for(i = 1; i <= nx * ny; i++) ssum += pow(ABS(xx[i]),1.1); 

return(ssum);
}
/**************************************************************
* Gradient for generalized Gauss' regularisation:
* Bauman & Sauer's potential function, with p=1.1
*        _
*        \        p
* phi =  /  | x  |
*        -     i
*        i
*  dphi = p*sign(x).*(abs(real(x)).^(p-1));
* here: dphi = 1.1 x**0.1
**************************************************************/
static void dphi_2D_ggauss(double *xx, double *dx)
{
register int i;

for(i = 1; i <= nx * ny; i++) 
  {
  if(xx[i] > 0.)
     dx[i] = 1.1 * pow(xx[i],0.1);
  else
     dx[i] = -1.1 *  pow(-xx[i],0.1);
  }

}
/**************************************************************
* Deconvolution by generalized Gauss' regularisation 
* Criterium to minimize is:
*             || y - H x ||^2 + alpha phi 
* with:
*        _      __________
*        \     /  2     2 |
* phi =  / \  / ss  +  x
*        -  \/           i
*        i
*
* ss is the level which is a kind of (intensity) threshold 
* between "good signal"  and "medium" or "bad" signal
* (penalty is less severe than tikhonov, for intensities > ss)
**************************************************************/
int dcv_2D_sqrtrf(double alpha_0, double ss, int positivity_0, 
                  double ftol, int *iter, double *L_phi, double *L_y)
{
register int i;
int nn, iter_max = 2000, lbfgs = 0;
double fret;

/* Initialize static parameters: */
alpha = alpha_0;
positivity = positivity_0;
strcpy(option, "sqrtr");

/*
* Set ss2 (private parameter) to the value chosen when calling the routine:
*/
ss2 = ss*ss;

/* Check if gradient is OK */
dcv_2D_check_grad(yyd, phi_2D_sqrtrf, dphi_2D_sqrtrf);

/* Starting point: */
for(i = 1; i <= nx * ny; i++) yyd[i] = 0.;

  nn = nx * ny;
  jlp_minimize(yyd, nn, nn, ftol, iter, iter_max, &fret, func_2D, dfunc_2D, 
               0, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* L-curve: */
*L_phi = phi_2D_sqrtrf(yyd);
*L_y = func_2D(yyd) - alpha * (*L_phi);

return(0);
}
/**************************************************************
* Regularisation function to be minimized for sqrt(s^2+x^2)
*        _      __________
*        \     /  2     2 |
* phi =  / \  / ss  +  x
*        -  \/           i
*        i
* (Carfantan: sqrtrf, with ss=0.01 or ss=0.001)
**************************************************************/
static double phi_2D_sqrtrf(double *xx)
{
double ssum;
register int i;

ssum = 0.;
for(i = 1; i <= nx * ny; i++) ssum += sqrt(ss2 + SQUARE(xx[i]));

return(ssum);
}
/**************************************************************
* Gradient for convex sqrt(s^2+x^2) regularisation:
* with:
*        _      __________
*        \     /  2     2 |
* phi =  / \  / ss  +  x
*        -  \/           i
*        i
* (Carfantan: sqrtrf, with ss=0.01 or ss=0.001)
**************************************************************/
static void dphi_2D_sqrtrf(double *xx, double *dx)
{
register int i;

for(i = 1; i <= nx * ny; i++) 
   dx[i] = xx[i] / sqrt(ss2 + SQUARE(xx[i]));

}
/**************************************************************
* Deconvolution by Wiener filter
**************************************************************/
int dcv_2D_wiener(double alpha_0)
{
register int i;

/* Initialize static parameters: */
alpha = alpha_0;
strcpy(option, "wiene");
positivity = 0;

/* Compute FFT of signal: */ 
   for(i = 1; i <= nx * ny; i++) {
     c_re(cp0[i]) = yy[i];
     c_im(cp0[i]) = 0.;
     }
   fftw_fast(&cp0[1], nx, ny, 1);
   for(i = 1; i <= nx * ny; i++) {
     yy_re[i] = c_re(cp0[i]);
     yy_im[i] = c_im(cp0[i]);
     }

for(i = 1; i <= nx*ny; i++) w0[i] = SQUARE(hh_re[i]) + SQUARE(hh_im[i]);

/****************************************
* 1/(a+ib) = (a-ib)/(a2+b2)
* Wiener filter = ((a2+b2) / ((a2+b2) + alpha)) / (a+ib)
* Hence:        = (a-ib) / ((a2+b2) + alpha))
****************************************/
for(i = 1; i <= nx*ny; i++) 
  if((w0[i] + alpha) != 0.){
  w_re[i] = hh_re[i] / (w0[i] + alpha);
  w_im[i] = - hh_im[i] / (w0[i] + alpha);
  }

/* Deconvolution: */
   for(i = 1; i <= nx*ny; i++) {
     c_re(cp0[i]) = w_re[i] * yy_re[i] - w_im[i] * yy_im[i];
     c_im(cp0[i]) = w_re[i] * yy_im[i] + w_im[i] * yy_re[i];
     }
   fftw_fast(&cp0[1], nx, ny, -1);
   for(i = 1; i <= nx * ny; i++) {
     yyd[i] = c_re(cp0[i]);
     }

return(0);
} 
/**************************************************************
* Deconvolution by spectral division 
**************************************************************/
int dcv_2D_spdiv()
{
register int i;

/* Compute FFT of signal: */ 
   for(i = 1; i <= nx * ny; i++) {
     c_re(cp0[i]) = yy[i];
     c_im(cp0[i]) = 0.;
     }
   fftw_fast(&cp0[1], nx, ny, 1);
   for(i = 1; i <= nx * ny; i++) {
     yy_re[i] = c_re(cp0[i]);
     yy_im[i] = c_im(cp0[i]);
     }

for(i = 1; i <= nx*ny; i++)
w0[i] = hh_re[i] * hh_re[i] + hh_im[i] * hh_im[i];

/* 1/(a+ib) = (a-ib)/(a2+b2)
*/
for(i = 1; i <= nx*ny; i++) {w_re[i] = 0.; w_im[i] = 0.;}
for(i = 1; i <= nx*ny; i++) 
  if(w0[i] != 0.)
   {
   w_re[i] = hh_re[i] / w0[i];
   w_im[i] = - hh_im[i] / w0[i];
   }

/* Deconvolution: */
   for(i = 1; i <= nx*ny; i++) {
     c_re(cp0[i]) = w_re[i] * yy_re[i] - w_im[i] * yy_im[i];
     c_im(cp0[i]) = w_re[i] * yy_im[i] + w_im[i] * yy_re[i];
     }
   fftw_fast(&cp0[1], nx, ny, -1);
   for(i = 1; i <= nx * ny; i++) {
     yyd[i] = c_re(cp0[i]);
     }

return(0); 
}
/**************************************************************
* Convolution H*X :  w0 = conv(xx,hh)
*
* IN:
* xx: function to be convolved by hh
* hh_re, hh_im: TF of hh
*
* OUT:
* ww = conv(hh,xx)  WARNING: often ww=w0 in the calling program...
*
**************************************************************/
static int dcv_2D_conv_hh(double *xx, double *hh_re, double *hh_im, double *ww)
{
register int i;
/* Convolution = product in Fourier domain: */

/* Compute FFT of xx: */ 
   for(i = 1; i <= nx * ny; i++) {
     c_re(cp0[i]) = xx[i];
     c_im(cp0[i]) = 0.;
     }
   fftw_fast(&cp0[1], nx, ny, 1);
   for(i = 1; i <= nx * ny; i++) {
     w_re[i] = c_re(cp0[i]);
     w_im[i] = c_im(cp0[i]);
     }

/* Product in Fourier domain: */
   for(i = 1; i <= nx*ny; i++) {
     c_re(cp0[i]) = w_re[i] * hh_re[i] - w_im[i] * hh_im[i];
     c_im(cp0[i]) = w_re[i] * hh_im[i] + w_im[i] * hh_re[i];
     }
   fftw_fast(&cp0[1], nx, ny, -1);
   for(i = 1; i <= nx * ny; i++) {
     ww[i] = c_re(cp0[i]);
     }

return(0);
}
/**************************************************************
* Check the validity of the gradient
* and compare f(x+dx)-f(x)/dx with grad_f(x)
*
* x1,x2,dx: work space of dimension nx * ny
**************************************************************/
static int dcv_2D_check_grad(double *x1, double (*func)(double []), 
                             void (*dfunc)(double[], double[]))
{
double *x2, *dx;
register int i, j;
int nw, idim;
double eps=1.e-4, tolerance=1.e-4;
double f_x2,f_x1,error;

idim = nx * ny + 1;
if( (x2 = (double *)malloc(idim * sizeof(double))) == NULL
    || (dx = (double *)malloc(idim * sizeof(double))) == NULL)
  {
   printf("dcv_2D_check_grad/Fatal error allocating memory: idim=%d\n",
           idim);
   exit(-1);
  }

printf(" dcv_2D_check_grad/Start checking gradient \n");

/* Generate random vector (between 0 and 1) */
for(i = 1; i <= nx * ny; i++) x1[i] = (double)rand() / (double)RAND_MAX;

f_x1 = (*func)(x1); 

/* Loop on a subsample of components (too long otherwise): */
nw = (nx * ny)/2;
for(i = nw; i < nw + 10; i++)
   {
   for(j = 1; j <= nx * ny; j++) x2[j] = x1[j];
/* Small variation of component #i: */
   x2[i] = x1[i] + eps;
   f_x2 = (*func)(x2); 
   (*dfunc)(x1,dx);
   error = (f_x2 - f_x1)/eps - dx[i];
#ifdef DEBUG
   printf(" x1[%d]=%.3f f_x1=%e f_x2=%e (f_x2 - f_x1)/eps=%e dx[i]=%e error=%e\n",
            i, x1[i], f_x1, f_x2, (f_x2 - f_x1)/eps, dx[i], error);
#endif
/* Relative error: (error with "abs", OK with "ABS") */
   error /= (ABS(dx[i]) + 1.e-12);
    if(error > tolerance) {
      printf("dcv_2D_check_grad/Error! \n");
      printf("component #i=%d:  relative error =%.4e\n", i, error);
      }
   }

printf("dcv_2D_check_grad/gradient is OK. \n");

free(x2);
free(dx);
return(0);
}
/***********************************************************************
* To initialize all static (private) arrays
*
***********************************************************************/
int dcv_2D_init(int nx_0, int ny_0)
{
int idim;

/* Transfer to static variables: */
nx = nx_0;
ny = ny_0;

/* Allocation of memory space: */
idim = 1 + nx_0 * ny_0;
if( (yy = (double *)malloc(idim * sizeof(double))) == NULL
 || (yy0 = (double *)malloc(idim * sizeof(double))) == NULL
 || (yyd = (double *)malloc(idim * sizeof(double))) == NULL
 || (yy_re = (double *)malloc(idim * sizeof(double))) == NULL
 || (yy_im = (double *)malloc(idim * sizeof(double))) == NULL
 || (hh_re = (double *)malloc(idim * sizeof(double))) == NULL
 || (hh_im = (double *)malloc(idim * sizeof(double))) == NULL
 || (hht_re = (double *)malloc(idim * sizeof(double))) == NULL
 || (hht_im = (double *)malloc(idim * sizeof(double))) == NULL
 || (w_re = (double *)malloc(idim * sizeof(double))) == NULL
 || (w_im = (double *)malloc(idim * sizeof(double))) == NULL
 || (w0 = (double *)malloc(idim * sizeof(double))) == NULL
 || (w1 = (double *)malloc(idim * sizeof(double))) == NULL
 || (cp0 = (FFTW_COMPLEX *) malloc(idim * sizeof(FFTW_COMPLEX))) == NULL)
 {
  printf("dcv_2D_init/Fatal error allocating memory space: idim=%d\n",idim);
  exit(-1);
 }

/* Initialization for fast FFT's: */
fftw_setup(nx, ny);

return(0);
}
/************************************************************
* To free memory
************************************************************/
void dcv_2D_free()
{
free(yy);
free(yy0);
free(yyd);
free(yy_re);
free(yy_im);
free(hh_re);
free(hh_im);
free(hht_re);
free(hht_im);
free(w_re);
free(w_im);
free(w0);
free(w1);
free(cp0);

FFTW_SHUTDOWN();

return;
}
/****************************************************************
* Generates a noisy signal with a given SNR
* snr: in dB
* INPUT:
*   yy0[nx*ny]
* OUTPUT:
*   yy[nx*ny]
****************************************************************/
int dcv_2D_noisy_signal(double snr)
{
double ssum, noise_level;
register int i;
 
/* dcv_conv_hh returns: yy = conv(hh,yy0)
*/
  dcv_2D_conv_hh(yy0, hh_re, hh_im, yy);
  RECENT_FFT_DOUBLE(yy, yy, &nx, &ny, &nx);

/* Generate noisy signal with SNR in dB:
* Multiply by 12, since variance of uniform random function is 1/12.
*/
  ssum = 0.;
  for(i = 1; i <= nx * ny; i++) ssum += SQUARE(yy[i]);
 
  noise_level = 12. * ssum * pow(10.,-snr/(double)10.) / (double)(nx * ny);
  noise_level = sqrt(noise_level);
  printf("Simulated signal with SNR=%.2f dB => noise level:%.3f\n",
         snr, noise_level);

/* Generate random vector (between 0 and 1) */
  for(i = 1; i <= nx * ny; i++) w0[i] = (double)rand() / (double)RAND_MAX;

/* Center noise and normalize it: */
  for(i = 1; i <= nx * ny; i++) yy[i] += (w0[i] - 0.5) * noise_level;

return(0);
}
/**************************************************************
* Compute transfer function: 
*
* INPUT:
* hh
* OUTPUT:
* hh_re, hh_im
* hht_re, hht_im
**************************************************************/
int dcv_2D_transfert(double *hh)
{
register int i, j;

   for(i = 1; i <= nx * ny; i++) {
     c_re(cp0[i]) = hh[i];
     c_im(cp0[i]) = 0.;
     }
   fftw_fast(&cp0[1], nx, ny, 1);
   for(i = 1; i <= nx * ny; i++) {
     hh_re[i] = c_re(cp0[i]);
     hh_im[i] = c_im(cp0[i]);
     }

/* Compute FFT function of reversed PSF (i.e., PSF(-x)):
*/
   for(j = 1; j <= ny; j++) {
     for(i = 1; i <= nx; i++) {
       c_re(cp0[i + (j-1) * nx]) = hh[(nx - (i-1)) + (ny - j) * nx];
       c_im(cp0[i + (j-1) * nx]) = 0.;
     }
   }
   fftw_fast(&cp0[1], nx, ny, 1);
   for(i = 1; i <= nx * ny; i++) {
     hht_re[i] = c_re(cp0[i]);
     hht_im[i] = c_im(cp0[i]);
     }

return(0);
}
