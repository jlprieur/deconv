/****************************************************************** 
* Module used by dcv_deconv_1D.f90
*
*  spdiv => Direct spectral division
*  wiene => Wiener filter 
*  tikho => Tikhonov regularisation
*  maxen => Maximum Entropy method
*
* JLP
* Version 23/03/2001
*******************************************************************/
#include <jlp_dcv.h>
#include "jlp_fftw.h"

/* Public: */
double *yy, *yy0, *yyd, *w0, *w1;

/* Private: */
/*******************
* For dcv_1D_mod:
* ss2 is a parameter needed for convex sqrt regularization
*/
static double ss2;
static int nn1, nn, hh_nx;
static double *hh_re, *hh_im, *hht_re, *hht_im;
static double *xx1, *xx2, *w_re, *w_im;
static double *yy_re, *yy_im;
static double alpha;
/* Complex array (for fast FFT's) */
static FFTW_COMPLEX *cp0;

/* Public:  (declared in jlp_dcv.h)
int dcv_1D_gmark(double alpha_0, double ftol, int *iter, int positive);
int dcv_1D_tikhonov(double alpha_0, double ftol, int *iter, int positive);
int dcv_1D_ggauss(double alpha_0, double ftol, int *iter, int positive);
int dcv_1D_sqrtrf(double alpha_0, double ss, double ftol, int *iter, int positive);
int dcv_1D_mem(double alpha_0, double ftol, int *iter, int positive);
int dcv_1D_spdiv();
int dcv_1D_wiener(double alpha_0);
int dcv_1D_banana(double ftol);
int dcv_1D_init(int nn1_0, int nn_0, int hh_nx_0);
void dcv_1D_free();
int dcv_1D_plot1(double *yy1, int nn, char *title);
int dcv_1D_plot1_log(double *yy1, int nn1, char *title);
int dcv_1D_plot2(double *yy1, double *yy2, int nn, char *title);
int dcv_1D_display_input(int nn);
int dcv_1D_transfert(double *hh);
*/

/* Private functions defined here: */
static int dcv_1D_check_grad(double *x1, double (*func)(double []), 
                          void (*dfunc)(double[], double[]));
static int dcv_1D_conv_hh(double *xx, double *hh_re, double *hh_im, double *ww);
static double func_1D_gmark(double *xx);
static void dfunc_1D_gmark(double *xx, double *dx);
static double func_1D_mem(double *xx);
static void dfunc_1D_mem(double *xx, double *dx);
static double func_1D_tikho(double *xx);
static void dfunc_1D_tikho(double *xx, double *dx);
static double func_1D_ggauss(double *xx);
static void dfunc_1D_ggauss(double *xx, double *dx);
static double func_1D_sqrtrf(double *xx);
static void dfunc_1D_sqrtrf(double *xx, double *dx);
static double func_1D_banana(double *x);
static void dfunc_1D_banana(double *xx, double *dx);

/*************************************************************
* To display a curve (linear) 
**************************************************************/
int dcv_1D_plot1(double *yy1, int nn, char *title)
{
register int i;
char xlabel[40], ylabel[40], plotdev[40];

/* X axis: */
 for(i = 1; i <= nn; i++) xx1[i] = (double)i;

strcpy(xlabel," ");
strcpy(ylabel," ");
strcpy(plotdev,"xterm_s");
jlp_display1(xx1, yy1, 1, nn, xlabel, ylabel, title, plotdev);

return(0);
}
/*************************************************************
* To display a curve (semi-log)
**************************************************************/
int dcv_1D_plot1_log(double *yy1, int nn1, char *title)
{
register int i;

/* X axis: */
 for(i = 1; i <= nn1; i++) xx1[i] = (double)(i - nn1/2)/ (double)nn1;
/* Y axis: */
 for(i = 1; i <= nn1; i++) xx2[i] = MY_LOG10(yy1[i]);

 jlp_display1(xx1, xx2, 1, nn1, " ", " ", title, "xterm_s");

return(0);
}
/*************************************************************
* To display two curves (linear) 
**************************************************************/
int dcv_1D_plot2(double *yy1, double *yy2, int nn, char *title)
{
register int i;

/* X axis: */
 for(i = 1; i <= nn; i++) xx1[i] = (double)i;

 jlp_display2(xx1, yy1, 1, nn, xx1, yy2, 1, nn, " ", " ", title, "xterm_s",
              "L0", "L1");

return(0);
}
/*************************************************************
* Deconvolution with Gauss-Markov regularisation 
* Criterium to minimize is:
* || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^2
*
* INPUT:
* ftol: tolerance
**************************************************************/
int dcv_1D_gmark(double alpha_0, double ftol, int *iter,
                 int positive)
{
double fret;
int iter_max = 2000, lbfgs = 0;
register int i;

alpha = alpha_0;

/* Check if gradient is OK */
dcv_1D_check_grad(yyd, func_1D_gmark, dfunc_1D_gmark);

/* Starting point: */
for(i = 1; i <= nn1; i++) yyd[i] = 0.;

  jlp_minimize(yyd, nn, nn1, ftol, iter, iter_max, &fret, 
               func_1D_gmark, dfunc_1D_gmark, positive, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* Pb with last point sometimes, so: */
yyd[nn] = yyd[nn-1]; 

return(0);
}
/**************************************************************
* Function E to be minimized for Gauss-Markov 
* E : value of the criterium in X
*                     \---             2
* E = |Y-Ax|^2+ alpha  >   | x   - x  | 
*                     /___    i+1   i
*                      i
*
**************************************************************/
static double func_1D_gmark(double *xx)
{
double sum2_y, sum2_dx;
register int i;

/* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

sum2_y = 0.;
for(i = 1; i <= nn; i++) sum2_y += SQUARE(w0[i] - yy[i]); 

/* Bauman & Sauer' function with p=2:
*        _
*        \             p
* phi =  / | x   - x  |
*        -    i+1   i
*        i
* (Smooth function, which reduces the variation between two successive pixels)
*/
sum2_dx = 0.;
for(i = 2; i <= nn; i++) sum2_dx += SQUARE(xx[i] - xx[i-1]); 

return(sum2_y + alpha * sum2_dx);

}
/**************************************************************
* Gradient for Gauss-Markov 
* grad_Er : value of the gradient of the criterium in X
*                                     (                 \---
* This function returns  grad(E) =grad( |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
*                                     (                 /___
*                                                      (s,r) in C
*
* d/dx (HX - Y)^T (HX - Y) = 2 H^T (HX - Y)
**************************************************************/
static void dfunc_1D_gmark(double *xx, double *dx)
{
register int i;

/* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

for(i = 1; i <= nn; i++) dx[i] = w0[i] - yy[i];

dcv_1D_conv_hh(dx, hht_re, hht_im, w0);

/* Shift hh_nx pixels to origin, 
* since convolution by hht has shifted the input signal by that amount: */
for(i = 1; i <= nn; i++) w0[i] = w0[hh_nx + i - 1];
for(i = nn + 1; i <= nn1; i++) w0[i] = 0.;

/* Gradient of Bauman & Sauer' function with p=2:
*        _
*        \             2
* phi =  / | x   - x  |
*        -    i+1   i
*        i
* (Smooth function, which reduces the variation between two successive pixels)
*/

dx[1] = 2. * w0[1] - 2. * alpha * (xx[2]-xx[1]);

for(i = 2; i <= nn - 1; i++) 
 dx[i] = 2. * w0[i] - 2. * alpha * (xx[i+1]-xx[i]) 
         + 2. * alpha * (xx[i]-xx[i-1]);

dx[nn] = 2. * w0[nn] + 2. * alpha*(xx[nn]-xx[nn-1]);

}
/**************************************************************
* Deconvolution with MEM 
* Criterium to minimize is:
* || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^2
**************************************************************/
int dcv_1D_mem(double alpha_0, double ftol, int *iter,
               int positive)
{
double fret;
int iter_max = 2000, lbfgs = 0;
register int i;

alpha = alpha_0;

/* Check if gradient is OK */
dcv_1D_check_grad(yyd, func_1D_mem, dfunc_1D_mem);

/* Starting point: */
for(i = 1; i <= nn1; i++) yyd[i] = 0.;

  jlp_minimize(yyd, nn, nn1, ftol, iter, iter_max, &fret, 
               func_1D_mem, dfunc_1D_mem, positive, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* Pb with last point sometimes, so: */
yyd[nn] = yyd[nn-1]; 

return(0);
}
/**************************************************************
* Function E to be minimized for MEM 
* E : value of the criterium in X
*                                             \---
* This function returns   E = |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
*                                             /___
*                                             (s,r) in C
*
**************************************************************/
static double func_1D_mem(double *xx)
{
register int i;
double sum2_y, sum_x;

/*
* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

sum2_y = 0.;
for(i = 1; i <= nn; i++) sum2_y += SQUARE(w0[i] - yy[i]); 

/*
* Entropy potential function:
*        _
*        \     
* phi =  /   x log(x )
*        -    i     i
*        i
*
* Penalty of 1000. if negative...
*/
sum_x = 0.;
for(i = 1; i <= nn; i++) 
  {
  if(xx[i] > 0.)
    sum_x += xx[i] * log(xx[i]);
  else
    sum_x += 1000.; 
  }

return(sum2_y + alpha * sum_x);
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
static void dfunc_1D_mem(double *xx, double *dx)
{
register int i;

/*
* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

for(i = 1; i <= nn; i++) dx[i] = w0[i] - yy[i]; 

dcv_1D_conv_hh(dx, hht_re, hht_im, w0);

/* Shift hh_nx pixels to origin, 
* since convolution by hht has shifted the input signal by that amount: */
for(i = 1; i <= nn; i++) w0[i] = w0[hh_nx + i - 1];
for(i = nn + 1; i <= nn1; i++) w0[i] = 0.;

/*
* Gradient de la fonction potentiel entropique
*        _
*        \     
* phi =  /   x log(x )
*        -    i     i
* => dphi = (1 + log(x)):
*/
for(i = 1; i <= nn; i++) 
   {
   if(xx[i] > 0.) 
     w1[i] = 1. + log(xx[i]);
   else
     w1[i] = -1000.;
   }

for(i = 1; i <= nn; i++) dx[i] = 2. * w0[i] + alpha * w1[i];

}
/**************************************************************
* Deconvolution by Tikhonov's regularisation 
* Criterium to minimize is:
* || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^2
* (Carfantan: grf)
**************************************************************/
int dcv_1D_tikhonov(double alpha_0, double ftol, int *iter,
                    int positive)
{
double fret;
int iter_max = 2000, lbfgs = 0;
register int i;

alpha = alpha_0;

/* Check if gradient is OK */
dcv_1D_check_grad(yyd, func_1D_tikho, dfunc_1D_tikho);

/* Starting point: */
for(i = 1; i <= nn; i++) yyd[i] = 0.;

  jlp_minimize(yyd, nn, nn1, ftol, iter, iter_max, &fret, 
               func_1D_tikho, dfunc_1D_tikho, positive, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* Pb with last point sometimes, so: */
yyd[nn] = yyd[nn-1]; 

return(0);
}
/**************************************************************
* Function E to be minimized for Tikhonov's regularisation:
* E : value of the criterium in X
*                                            \---
* This funtion returns   E = |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
*                                            /___
*                                           (s,r) in C
*
**************************************************************/
static double func_1D_tikho(double *xx)
{
register int i;
double sum2_y, sum_x;

/*
* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

sum2_y = 0.;
for(i = 1; i <= nn; i++) sum2_y += SQUARE(w0[i] - yy[i]); 

/*
* Bauman & Sauer's potential function, with p=2
*        _
*        \        p
* phi =  /  | x  |
*        -     i
*        i
*  phi = sum(abs(real(x(:))).^p);
*  dphi = p*sign(x).*(abs(real(x)).^(p-1));
*/
sum_x = 0.;
for(i = 1; i <= nn; i++) sum_x += SQUARE(xx[i]); 

return(sum2_y + alpha * sum_x); 
}
/**************************************************************
* Gradient for Tikhonov's regularisation:
* grad_Er : value of the gradient of the criterium in X
*                                    (                 \---
* This funtion returns  grad(E) =grad( |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
*                                    (                 /___
*                                                      (s,r) in C
*
* d/dx (HX - Y)^T (HX - Y) = 2 H^T (HX - Y)
**************************************************************/
static void dfunc_1D_tikho(double *xx, double *dx)
{
register int i;

/*
* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

for(i = 1; i <= nn; i++) dx[i] = w0[i] - yy[i]; 

dcv_1D_conv_hh(dx, hht_re, hht_im, w0);

/* Shift hh_nx pixels to origin, 
* since convolution by hht has shifted the input signal by that amount: */
for(i = 1; i <= nn; i++) w0[i] = w0[hh_nx + i - 1];
for(i = nn + 1; i <= nn1; i++) w0[i] = 0.;

/*
* Bauman & Sauer's potential function, with p=2
*        _
*        \        p
* phi =  /  | x  |
*        -     i
*        i
*  dphi = p*sign(x).*(abs(real(x)).^(p-1));
* here: dphi = 2 x
*/
for(i = 1; i <= nn; i++) dx[i] = 2. * w0[i] + 2. * alpha * xx[i];

}
/**************************************************************
* Deconvolution by generalized Gauss' regularisation 
* Criterium to minimize is:
* || y - H x ||^2 + alpha Sigma_{k=1}^{k=N} x_k^p
* with p=1.1
* (Carfantan: ggrf)
**************************************************************/
int dcv_1D_ggauss(double alpha_0, double ftol, int *iter,
                  int positive)
{
double fret;
int iter_max = 2000, lbfgs = 0;
register int i;

alpha = alpha_0;

/* Check if gradient is OK */
dcv_1D_check_grad(yyd, func_1D_ggauss, dfunc_1D_ggauss);

/* Starting point: */
for(i = 1; i <= nn; i++) yyd[i] = 0.;

  jlp_minimize(yyd, nn, nn1, ftol, iter, iter_max, &fret, 
               func_1D_ggauss, dfunc_1D_ggauss, positive, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* Pb with last point sometimes, so: */
yyd[nn] = yyd[nn-1]; 

return(0);
}
/**************************************************************
* Function E to be minimized for generalized Gauss' regularisation:
* E : value of the criterium in X
*                                             \---
* This function returns   E = |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
*                                             /___
*                                            (s,r) in C
*
**************************************************************/
static double func_1D_ggauss(double *xx)
{
register int i;
double sum2_y, sum_x;

/*
* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

sum2_y = 0.;
for(i = 1; i <= nn; i++) sum2_y += SQUARE(w0[i] - yy[i]); 

/*
* Bauman & Sauer's potential function, with p=1.1
*        _
*        \        p
* phi =  /  | x  |
*        -     i
*        i
*  phi = sum(abs(real(x(:))).^p);
*  dphi = p*sign(x).*(abs(real(x)).^(p-1));
*/
sum_x = 0.;
for(i = 1; i <= nn; i++) sum_x += pow(ABS(xx[i]),1.1); 

return(sum2_y + alpha * sum_x);
}
/**************************************************************
* Gradient for generalized Gauss' regularisation:
* grad_Er : value of the gradient of the criterium in X
*                                     (                 \---
* This function returns  grad(E) =grad( |Y-Ax|^2+ alpha  >    phi(Xs-Xr)
*                                     (                 /___
*                                                     (s,r) in C
*
* d/dx (HX - Y)^T (HX - Y) = 2 H^T (HX - Y)
**************************************************************/
static void dfunc_1D_ggauss(double *xx, double *dx)
{
register int i;

/*
* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

for(i = 1; i <= nn; i++) dx[i] = w0[i] - yy[i]; 

dcv_1D_conv_hh(dx, hht_re, hht_im, w0);

/* Shift hh_nx pixels to origin, 
* since convolution by hht has shifted the input signal by that amount: */
for(i = 1; i <= nn; i++) w0[i] = w0[hh_nx + i - 1];
for(i = nn + 1; i <= nn1; i++) w0[i] = 0.;

/*
* Bauman & Sauer's potential function, with p=1.1
*        _
*        \        p
* phi =  /  | x  |
*        -     i
*        i
*  dphi = p*sign(x).*(abs(real(x)).^(p-1));
* here: dphi = 1.1 x**0.1
*/
for(i = 1; i <= nn; i++) 
  {
  if(xx[i] > 0.)
     w1[i] = pow(xx[i],0.1);
  else
     w1[i] = - pow(-xx[i],0.1);
  }

for(i = 1; i <= nn; i++) dx[i] = 2. * w0[i] + 1.1 * alpha * w1[i];

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
int dcv_1D_sqrtrf(double alpha_0, double ss, double ftol, int *iter,
                  int positive)
{
register int i;
int iter_max = 2000, lbfgs = 0;
double fret;

alpha = alpha_0;

/*
* Set ss2 (private parameter) to the value chosen when calling the routine:
*/
ss2 = ss*ss;

/* Check if gradient is OK */
dcv_1D_check_grad(yyd, func_1D_sqrtrf, dfunc_1D_sqrtrf);

/* Starting point: */
for(i = 1; i <= nn; i++) yyd[i] = 0.;

  jlp_minimize(yyd, nn, nn1, ftol, iter, iter_max, &fret, 
               func_1D_sqrtrf, dfunc_1D_sqrtrf, positive, lbfgs);
  printf(" Number of iterations: %d \n",*iter);
  printf(" Value of the minimum: %.5g\n",fret);

/* Pb with last point sometimes, so: */
yyd[nn] = yyd[nn-1]; 

  printf(" alpha=%f ss=%f \n",alpha_0, sqrt(ss2));
return(0);
}
/**************************************************************
* Function E to be minimized for convex sqrt(s^2+x^2) regularisation:
* E : value of the criterium in X
*                                               
* This function returns   E = |Y-Ax|^2+ alpha  phi
*                                              
* with:
*        _      __________
*        \     /  2     2 |
* phi =  / \  / ss  +  x
*        -  \/           i
*        i
* (Carfantan: sqrtrf, with ss=0.01 or ss=0.001)
**************************************************************/
static double func_1D_sqrtrf(double *xx)
{
double sum2_y, sum_x;
register int i;

/*
* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

sum2_y = 0.;
for(i = 1; i <= nn; i++) sum2_y += SQUARE(w0[i] - yy[i]);

sum_x = 0.;
for(i = 1; i <= nn; i++) sum_x += sqrt(ss2 + SQUARE(xx[i]));

return(sum2_y + alpha * sum_x);
}
/**************************************************************
* Gradient for convex sqrt(s^2+x^2) regularisation:
* grad_Er : value of the gradient of the criterium in X
*
* This function returns  grad(E) =grad( |Y-Ax|^2+ alpha phi )
*                                                      
* with:
*        _      __________
*        \     /  2     2 |
* phi =  / \  / ss  +  x
*        -  \/           i
*        i
* (Carfantan: sqrtrf, with ss=0.01 or ss=0.001)
*
* d/dx (HX - Y)^T (HX - Y) = 2 H^T (HX - Y)
**************************************************************/
static void dfunc_1D_sqrtrf(double *xx, double *dx)
{
register int i;

/*
* Convolution = product in Fourier domain:
* w0 = conv(xx,hh):
*/
dcv_1D_conv_hh(xx, hh_re, hh_im, w0);

for(i = 1; i <= nn; i++) dx[i] = w0[i] - yy[i]; 

dcv_1D_conv_hh(dx, hht_re, hht_im, w0);

/* Shift hh_nx pixels to origin, 
* since convolution by hht has shifted the input signal by that amount: */
for(i = 1; i <= nn; i++) w0[i] = w0[hh_nx + i - 1];
for(i = nn + 1; i <= nn1; i++) w0[i] = 0.;

for(i = 1; i <= nn; i++) 
   dx[i] = 2. * w0[i] + alpha * xx[i] / sqrt(ss2 + xx[i]*xx[i]);

}
/**************************************************************
* Banana function : start at X=[-1.9;2]. minimum at X=[1;1] : f(X)=0;
* z = 100*(x(2)^2-x(1))^2 + (1-x(1))^2
*  dx = [  -200*(x(2)^2 - x(1))-2*(1-x(1)) ] 
        [   100*(4 * x(2)^3 - 4 x(1))*x(2) ]
*
**************************************************************/
int dcv_1D_banana(double ftol)
{
int n=2, iter;
int positive = 0;
int iter_max = 2000, lbfgs = 0;
double p[3], fret;

printf("Test with Banana function\n");

/* Starting point: */
  p[1]= -1.9;
  p[2]= 2.;

  jlp_minimize(p, n, n, ftol, &iter, iter_max, &fret, 
               func_1D_banana, dfunc_1D_banana, positive, lbfgs);
  printf(" Location of the minimum: %.3f %.3f\n", p[1], p[2]);
  printf(" Number of iterations: %d\n",iter);
  printf(" Value of the minimum: %.5e\n",fret);

return(0); 
}
/**********************************************************
* Banana function:
* z = 100*(x(2)^2-x(1))^2 + (1-x(1))^2
**********************************************************/
static double func_1D_banana(double *x)
{
 return(SQUARE(1. - x[1]) + 100. * SQUARE(SQUARE(x[2]) - x[1]));
}
/**********************************************************
* Banana function:
*  dx = [  -200*(x(2)^2 - x(1))-2*(1-x(1)) ] 
        [   100*(4 * x(2)^3 - 4 x(1))*x(2) ]
**********************************************************/
static void dfunc_1D_banana(double *x, double *dx)
{
 dx[1] = -200. * (SQUARE(x[2]) - x[1]) - 2. * (1. - x[1]);
 dx[2] = 400. * (x[2]*x[2]*x[2] - x[1] * x[2]);
}
/**************************************************************
* Deconvolution by Wiener filter
**************************************************************/
int dcv_1D_wiener(double alpha_0)
{
register int i;

alpha = alpha_0;

/* Compute FFT of signal: */
   for(i = 1; i <= nn1; i++) {
     c_re(cp0[i]) = yy[i];
     c_im(cp0[i]) = 0.;
     }
   fftw_fast(&cp0[1], nn1, 1, 1);
   for(i = 1; i <= nn1; i++) {
     yy_re[i] = c_re(cp0[i]);
     yy_im[i] = c_im(cp0[i]);
     }

/* Compute square modulus of transfer function: */
  for(i = 1; i <= nn1; i++) 
     w0[i] = SQUARE(hh_re[i]) + SQUARE(hh_im[i]);

/****************************************
* 1/(a+ib) = (a-ib)/(a2+b2)
* Wiener filter = ((a2+b2) / ((a2+b2) + alpha)) / (a+ib)
* Hence:        = (a-ib) / ((a2+b2) + alpha))
****************************************/
for(i = 1; i <= nn1; i++) {w_re[i] = 0.; w_im[i] = 0.;}
for(i = 1; i <= nn1; i++) 
  if((w0[i] + alpha) != 0.){
  w_re[i] = hh_re[i] / (w0[i] + alpha);
  w_im[i] = - hh_im[i] / (w0[i] + alpha);
  }

/* Plot Wiener filter: */
if(0){
  printf(" Display Wiener filter\n");
  for(i = 1; i <= nn1; i++) w1[i] = SQUARE(w_re[i]) + SQUARE(w_im[i]);
  dcv_1D_plot1_log(w1, nn1, "Wiener filter");
  }

/* Deconvolution: */
  for(i = 1; i <= nn1; i++) {
    c_re(cp0[i]) = w_re[i] * yy_re[i] - w_im[i] * yy_im[i];
    c_im(cp0[i]) = w_re[i] * yy_im[i] + w_im[i] * yy_re[i];
    }
   fftw_fast(&cp0[1], nn1, 1, -1);
   for(i = 1; i <= nn; i++) yyd[i] = c_re(cp0[i]);
   for(i = nn+1; i <= nn1; i++) yyd[i] = 0.; 

return(0);
} 
/**************************************************************
* Deconvolution by spectral division 
**************************************************************/
int dcv_1D_spdiv()
{
register int i;

/* Compute FFT of signal: */
   for(i = 1; i <= nn1; i++) {
     c_re(cp0[i]) = yy[i];
     c_im(cp0[i]) = 0.;
     }
   fftw_fast(&cp0[1], nn1, 1, 1);
   for(i = 1; i <= nn1; i++) {
     yy_re[i] = c_re(cp0[i]);
     yy_im[i] = c_im(cp0[i]);
     }

/* Compute square modulus of transfer function: */
  for(i = 1; i <= nn1; i++) 
     w0[i] = SQUARE(hh_re[i]) + SQUARE(hh_im[i]);

/* 1/(a+ib) = (a-ib)/(a2+b2)
*/
for(i = 1; i <= nn1; i++) {w_re[i] = 0.; w_im[i] = 0.;}
for(i = 1; i <= nn1; i++) 
  if(w0[i] != 0.)
   {
   w_re[i] = hh_re[i] / w0[i];
   w_im[i] = - hh_im[i] / w0[i];
   }

/* Plot inverse filter: */
if(1){
  printf(" Display inverse filter\n");
  for(i = 1; i <= nn1; i++) w1[i] = SQUARE(w_re[i]) + SQUARE(w_im[i]);
  dcv_1D_plot1_log(w1, nn1, "Inverse filter");
  }

/* Deconvolution: */
  for(i = 1; i <= nn1; i++) {
    c_re(cp0[i]) = w_re[i] * yy_re[i] - w_im[i] * yy_im[i];
    c_im(cp0[i]) = w_re[i] * yy_im[i] + w_im[i] * yy_re[i];
    }
   fftw_fast(&cp0[1], nn1, 1, -1);
   for(i = 1; i <= nn; i++) yyd[i] = c_re(cp0[i]);
   for(i = nn+1; i <= nn1; i++) yyd[i] = 0.;

return(0); 
}
/**************************************************************
* Display input files (to check if OK)
**************************************************************/
int dcv_1D_display_input(int nn)
{
INT4 nx, ny;
int display_all = 1;
register int i;
char title[40];

/* Display original signal: */
if(display_all){
  printf("Display original signal:\n");
  strcpy(title,"Original signal"); 
  dcv_1D_plot1(yy0, nn, title);
  }

/* Display noisy signal to deconvolve: */
if(display_all){
  printf("Display noisy signal to deconvolve:\n");
  strcpy(title,"Noisy signal to deconvolve"); 
  dcv_1D_plot1(yy, nn, title);
  }

/* Compute power spectrum of original signal: */
  for(i = 1; i <= nn1; i++) {
    c_re(cp0[i]) = yy0[i];
    c_im(cp0[i]) = 0.;
    }
  fftw_fast(&cp0[1], nn1, 1, 1);
  for(i = 1; i <= nn1; i++) {
    w0[i] = SQUARE(c_re(cp0[i])) + SQUARE(c_im(cp0[i]));
    }
  nx = nn1; ny = 1;
  RECENT_FFT_1D_X(&w0[1], &w0[1], &nx, &ny, &nx);

 printf("nn1 = %d\n", nn1);
/* Display power spectrum of original signal: */
if(display_all){
  printf("Display power spectrum of original signal: \n");
  strcpy(title,"Power spectrum of original signal"); 
  dcv_1D_plot1_log(w0, nn1, title);
  }

/* Compute square modulus of transfer function: */
  for(i = 1; i <= nn1; i++) 
     w0[i] = SQUARE(hh_re[i]) + SQUARE(hh_im[i]);
  nx = nn1; ny = 1;
  RECENT_FFT_1D_X(&w0[1], &w0[1], &nx, &ny, &nx);

/* Display Transfer function: */
if(display_all){
  printf("Frequency response of the filter:\n");
  strcpy(title,"Frequency response of the filter"); 
  dcv_1D_plot1_log(w0, nn1, title);
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
static int dcv_1D_conv_hh(double *xx, double *hh_re, double *hh_im, double *ww)
{
register int i;

/* Just in case: */
for(i = nn+1; i <= nn1; i++) xx[i] = 0.;

/* Compute FFT of xx: */ 
  for(i = 1; i <= nn1; i++) {
    c_re(cp0[i]) = xx[i];
    c_im(cp0[i]) = 0.;
    }
  fftw_fast(&cp0[1], nn1, 1, 1);

/* Product in Fourier domain: */
for(i = 1; i <= nn1; i++)
  {
  w_re[i] = c_re(cp0[i]) * hh_re[i] - c_im(cp0[i]) * hh_im[i];
  w_im[i] = c_re(cp0[i]) * hh_im[i] + c_im(cp0[i]) * hh_re[i];
  }

/* Back to direct space: */
  for(i = 1; i <= nn1; i++) {
    c_re(cp0[i]) = w_re[i];
    c_im(cp0[i]) = w_im[i];
    }
  fftw_fast(&cp0[1], nn1, 1, -1);
  for(i = 1; i <= nn1; i++) ww[i] = c_re(cp0[i]);

return(0);
}
/**************************************************************
* Check the validity of the gradient
* and compare f(x+dx)-f(x)/dx with grad_f(x)
*
* x1: work space of dimension nn
**************************************************************/
static int dcv_1D_check_grad(double *x1, double (*func)(double []), 
                          void (*dfunc)(double[], double[]))
{
register int i, j;
double eps=1.e-4, tolerance=1.e-4;
double f_xx1,f_x1,error;

/* Generate random vector (between 0 and 1) */
for(i = 1; i <= nn1; i++) x1[i] = 0.;
for(i = 1; i <= nn; i++) x1[i] = (double)rand() / (double)RAND_MAX;
for(i = 1; i <= nn1; i++) xx1[i] = 0.;

f_x1 = (*func)(x1); 

/* Loop on all components: */
for(i = 1; i <= nn; i++)
   {
   for(j = 1; j <= nn; j++) xx1[j] = x1[j];
/* Small variation of component #i: */
   xx1[i] = x1[i] + eps;
   f_xx1 = (*func)(xx1); 
/* Gradient stored in xx2: */
   (*dfunc)(xx1,xx2);
   error = (f_xx1 - f_x1)/eps - xx2[i];
   error = error / (ABS(f_xx1 - f_x1)/eps + 1.e-12);
#ifdef DEBUG
   if(i < 10)
   printf(" f_x1=%e f_xx1=%e (f_xx1-f_x1)/eps = %e dx[%d] = %e error=%e\n",
            f_x1, f_xx1, (f_xx1 - f_x1)/eps, i, xx2[i], error);
#endif
    if(error > tolerance) {
      printf("dcv_1D_check_grad/Error! \n");
      printf("component #i=%d:  relative error =%.4e\n", i, error);
      }
   }

printf("dcv_1D_check_grad/End: gradient has been checked.\n");
return(0);
}
/***********************************************************************
* To initialize all static (private) arrays
*
* xx1, xx2: arrays used by dcv_1D_plot1, dcv_1D_plot_log, dcv_1D_plot2
*           and check_grad
***********************************************************************/
int dcv_1D_init(int nn1_0, int nn_0, int hh_nx_0)
{
int idim;

/* Transfer to static variables: */
nn1 = nn1_0;
nn = nn_0;
hh_nx = hh_nx_0;

/* Allocation of memory space: */
idim = nn1+1;
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
 || (xx1 = (double *)malloc(idim * sizeof(double))) == NULL
 || (xx2 = (double *)malloc(idim * sizeof(double))) == NULL
 || (w0 = (double *)malloc(idim * sizeof(double))) == NULL
 || (w1 = (double *)malloc(idim * sizeof(double))) == NULL
 || (cp0 = (FFTW_COMPLEX *) malloc(idim * sizeof(FFTW_COMPLEX))) == NULL)
 {
  printf("dcv_1D_init/Fatal error allocating memory space: idim=%d\n",idim);
  exit(-1);
 }

/* To prepare fast FFT's: */
printf(" fftw_setup with (%d,1)\n", nn1);
fftw_setup(nn1,1);

return(0);
}
/************************************************************
* To free memory
************************************************************/
void dcv_1D_free()
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
free(xx1);
free(xx2);
free(w0);
free(w1);
FFTW_SHUTDOWN();
return;
}
/*************************************************************
* Compute transfer function: 
**************************************************************/
int dcv_1D_transfert(double *hh)
{
register int i;

   for(i = 1; i <= nn1; i++) {
     c_re(cp0[i]) = hh[i];
     c_im(cp0[i]) = 0.;
     }
   fftw_fast(&cp0[1], nn1, 1, 1);
   for(i = 1; i <= nn1; i++) {
     hh_re[i] = c_re(cp0[i]);
     hh_im[i] = c_im(cp0[i]);
     }

/* Compute FFT function of reversed PSF (i.e., PSF(-x)):
* Length of the support of hh is hh_nx */
     for(i = 1; i <= nn1; i++) {
       c_re(cp0[i]) = 0.;
       c_im(cp0[i]) = 0.;
     }
     for(i = 1; i <= hh_nx; i++) 
       c_re(cp0[i]) = hh[hh_nx - i + 1];
   fftw_fast(&cp0[1], nn1, 1, 1);
   for(i = 1; i <= nn1; i++) {
     hht_re[i] = c_re(cp0[i]);
     hht_im[i] = c_im(cp0[i]);
     }

return(0);
}
